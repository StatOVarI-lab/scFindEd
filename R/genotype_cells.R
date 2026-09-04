#' Calling Cell Genotypes on Each Amplicon
#'
#' This function is used to assign a 'most-probable haploid genotype' for every
#' cell barcode in each amplicon of an experiment. It applies a Dirichlet
#' multinomial model and determines the variant in each cell that is the
#' most-probable to be the 'haploid genotype'.
#'
#' The function also calculates credible intervals for the posterior means
#' calculated, based on the provided `credible_interval` threshold. The lower
#' bound may be used to set a minimum probability for accepting an assigned
#' most-probable genotype.
#'
#' A numeric `alpha_0` can be provided directly and used across every amplicon.
#' However, if left blank (or `alpha_0 = "empirical`), then a value for `alpha_0`
#' will be estimated for each amplicon using an Empirical Bayes optimiser.
#'
#' @param alpha_0 Determines the alpha_0 scale parameter used to generate
#' pseudocounts in the Dirichlet multinomial model. Defaults to estimating alpha_0
#' through Empirical Bayes optimisation (also `alpha_0 = "empirical"`), but also
#' accepts a numeric value that will then be applied to all amplicons.
#' @param credible_interval Float that sets the threshold at which lower and
#' upper credible intervals are calculated. Must be in the range `[0, 1]`.
#'
#' @inheritParams global_params source_dir sample_name panel_ID output_dir verbose
#'
#' @return For each amplicon, exports the haploid genotypes file, containing the
#' most probable variant as the haploid genotype in each cell barcode, along with
#' reported posterior means and credible intervals.
#'
#' @export
genotype_cells <- function(
  source_dir,
  sample_name,
  panel_ID,
  alpha_0,
  credible_interval = 0.95,
  output_dir = NULL,
  verbose = TRUE
) {
  # Import SNP summary
  design_summary_file_name <- paste0(panel_ID, "_design_summary_joined.csv")
  design_summary_file <- import_design_summary_info(source_dir,
                                                    design_summary_file_name,
                                                    verbose)

  # Extract unique amplicon IDs
  amplicon_list <- unique(design_summary_file[["AmpID"]])
  n_amplicons <- length(amplicon_list)

  # Calculate credible interval tail size
  if (!(credible_interval > 0 & credible_interval < 1)) {
    cli::cli_abort("credible_interval must be a numeric value in the range [0, 1].")
  }
  CI_tail <- (1 - credible_interval) / 2

  for (i in 1:n_amplicons) {
    AmpID <- amplicon_list[i]

    if (verbose) cli::cli_h2(AmpID)

    # Import true mutations
    true_mutations_file_name <- paste0(sample_name, "_", AmpID,
                                       "_true_mutations.tsv")
    variant_counts <- import_true_mutations(source_dir,
                                            true_mutations_file_name,
                                            verbose) |>
      dplyr::rename(variant = merged_HGVSc)

    # Estimate global_prior
    global_prior <- variant_counts |>
      dplyr::group_by(variant) |>
      dplyr::summarise(variant_count = sum(count), .groups = "drop") |>
      dplyr::mutate(q_global = variant_count / sum(variant_count)) |>
      dplyr::select(variant, q_global)

    # Determine alpha_0 (i.e. parse provided argument value)
    if (is.null(alpha_0) || identical(alpha_0, "empirical")) {
      alpha_0_estimate <- empirical_bayes_alpha_0(variant_counts, global_prior, verbose)
    } else if (is.numeric(alpha_0)) {
      alpha_0_estimate <- alpha_0
    } else {
      cli::cli_abort('alpha_0 must be either "empirical" or numeric.')
    }

    # -- Haploidy --

    if (verbose) cli::cli_process_start(msg = "Calculating most probable haploid genotypes...",
                                        msg_done = "Most probable haploid genotypes calculated.")

    haploid_genotypes <- variant_counts |>
      # Join the global prior values
      dplyr::inner_join(global_prior, by = "variant") |>
      dplyr::group_by(dplyr::across(dplyr::any_of(c("cell_barcode",
                                                    "identity_barcode",
                                                    "cell_count")))) |>
      dplyr::mutate(
        # Scale prior by alpha_0_estimate to obtain pseudocounts
        alpha_local = alpha_0_estimate * q_global,

        # Beta posterior shape parameters for marginal distribution of theta_i
        alpha_posterior = alpha_local + count,
        beta_posterior = (alpha_0_estimate + cell_count) - alpha_posterior,

        # Conjugate Dirichlet-Multinomial posterior expectation
        posterior_mean = alpha_posterior / (alpha_0_estimate + cell_count),

        # Credible interval bounds derived from Beta quantile function
        CI_lower = qbeta(CI_tail, shape1 = alpha_posterior, shape2 = beta_posterior),
        CI_upper = qbeta(1 - CI_tail, shape1 = alpha_posterior, shape2 = beta_posterior),
        CI_confidence = 1 - (CI_upper - CI_lower)
      ) |>

      # Select only the maximum probability variant for the final genotype call
      dplyr::slice_max(order_by = posterior_mean, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::rename(
        variant_genotype = variant,
        variant_count = count
      ) |>
      dplyr::select(dplyr::any_of(c("sample_id", "cell_barcode", "amp_id",
                                    "identity_barcode", "variant_genotype",
                                    "variant_count", "cell_count", "posterior_mean",
                                    "CI_confidence", "CI_lower", "CI_upper")))

    if (verbose) cli::cli_process_done()

    # -- Optional file saving (only if output_dir is specified) --

    if (!is.null(output_dir)) {
      haploid_genotypes_file_name <- paste0(sample_name, "_", AmpID,
                                            "_haploid_genotypes.tsv")

      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      readr::write_tsv(haploid_genotypes, file = file.path(
        output_dir, haploid_genotypes_file_name))

      cli::cli_alert_success(
        "{AmpID} haploid genotypes saved to: {normalizePath(output_dir)}.")
    }
  }
}


#' @keywords internal
empirical_bayes_alpha_0 <- function(
  variant_counts,
  global_prior,
  verbose = TRUE,
  lower = 1e-2,
  upper = 1e4
) {
  if (verbose) cli::cli_process_start(msg = "Estimating alpha_0 with Empirical Bayes...")

  # Attach q_global to each row once — no dense matrix needed
  variant_q <- global_prior[["q_global"]][match(variant_counts[["variant"]], global_prior[["variant"]])]
  counts    <- variant_counts[["count"]]

  # Per-cell totals (N_c) and cell count, computed on the long data directly
  N_c     <- rowsum(counts, variant_counts[["cell_barcode"]])[, 1]
  n_cells <- length(N_c)

  negative_log_likelihood <- function(alpha_0) {
    # Term 1: Cell-level coverage penalty
    term1 <- n_cells * lgamma(alpha_0) - sum(lgamma(N_c + alpha_0))

    # Term 2: Variant-level gamma sums (only over non-zero counts)
    alpha_row <- alpha_0 * variant_q
    term2 <- sum(lgamma(counts + alpha_row) - lgamma(alpha_row))

    -(term1 + term2)
  }

  # 1D bounded optimisation using Brent's method
  optimisation <- optimize(negative_log_likelihood, interval = c(lower, upper))
  estimated_alpha_0 <- optimisation[["minimum"]]

  if (verbose) cli::cli_process_done(msg_done = "Empirical Bayes alpha_0 estimate: {round(estimated_alpha_0, 3)}.")

  if (estimated_alpha_0 == lower) {
    cli::cli_alert_warning("The empirical Bayes estimate for alpha_0 returns the lower optimisation bound.")
  } else if (estimated_alpha_0 == upper) {
    cli::cli_alert_warning("The empirical Bayes estimate for alpha_0 returns the upper optimisation bound.")
  }

  return(estimated_alpha_0)
}
