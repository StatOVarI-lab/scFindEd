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
#' The per-cell genotypes are then collated across amplicons to determine the
#' most probable mutations in each cell, as well as the identifying the presence
#' of inserted edits and assigning confident WTs.
#'
#' @param alpha_0 Determines the alpha_0 scale parameter used to generate
#' pseudocounts in the Dirichlet multinomial model. Defaults to estimating alpha_0
#' through Empirical Bayes optimisation (also `alpha_0 = "empirical"`), but also
#' accepts a numeric value that will then be applied to all amplicons.
#' @param credible_interval Float that sets the threshold at which lower and
#' upper credible intervals are calculated. Must be in the range `[0, 1]`.
#'
#' @inheritParams global_params source_dir sample_name panel_ID ploidy lower_CI_assignment_threshold output_dir plots_output_dir verbose
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
  ploidy = "haploid",
  credible_interval = 0.95,
  lower_CI_assignment_threshold = NULL,
  output_dir = NULL,
  plots_output_dir = NULL,
  verbose = TRUE
) {
  # -- Ensure all provided arguments are valid --

  validate_genotype_cells_args(source_dir,
                               sample_name,
                               panel_ID,
                               alpha_0,
                               ploidy,
                               credible_interval,
                               lower_CI_assignment_threshold,
                               output_dir,
                               plots_output_dir,
                               verbose)

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

  if (ploidy == "haploid") {
    CI_tail <- (1 - credible_interval) / 2
  } else if (ploidy == "diploid") {
    CI_tail <- 1 - credible_interval
  }

  assigned_genotypes <- list()

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

    # -- Assign Cells Genotypes --

    if (verbose) cli::cli_process_start(msg = "Calculating most probable {ploidy} genotypes...",
                                        msg_done = "Most probable {ploidy} genotypes calculated.")

    genotypes <- variant_counts |>
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
        CI_upper = if (ploidy == "diploid") qbeta(1 - CI_tail, shape1 = alpha_posterior, shape2 = beta_posterior) else NULL,
        CI_confidence = if (ploidy == "diploid") 1 - (CI_upper - CI_lower) else 1 - (1 - CI_lower)
      ) |>
      dplyr::ungroup() |>
      dplyr::rename(variant_genotype = variant, variant_count = count) |>
      dplyr::select(dplyr::any_of(c("sample_id", "cell_barcode", "amp_id",
                                    "identity_barcode", "variant_genotype",
                                    "variant_count", "cell_count", "posterior_mean",
                                    "CI_confidence", "CI_lower", "CI_upper")))

    if (verbose) cli::cli_process_done()

    assigned_genotypes[[AmpID]] <- genotypes

    # -- Optional file saving (only if output_dir is specified) --

    if (!is.null(output_dir)) {
      genotypes_file_name <- paste0(sample_name, "_", AmpID, "_", ploidy,
                                    "_genotypes.tsv")

      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      readr::write_tsv(genotypes, file = file.path(
        output_dir, genotypes_file_name))

      cli::cli_alert_success(
        "{AmpID} {ploidy} genotypes saved to: {normalizePath(output_dir)}.")
    }

    # -- Per Amplicon Posterior Mean Distribution (only if plots_output_dir is specified) --
    if (!is.null(plots_output_dir)) {
      AmpID_posterior_mean_plot <- posterior_mean_distribution_plot(AmpID,
                                                                    genotypes)

      AmpID_posterior_mean_plot_file_name <- paste0(sample_name, "_", AmpID,
                                                    "_posterior_mean_distribution.png")

      # Optional file saving (only if plots_output_dir is provided)
      if (!dir.exists(plots_output_dir)) dir.create(plots_output_dir, recursive = TRUE)
      ggplot2::ggsave(file.path(plots_output_dir, AmpID_posterior_mean_plot_file_name),
                      plot = AmpID_posterior_mean_plot,
                      width = 8, height = 5, dpi = 300)

      cli::cli_alert_success("{AmpID} posterior mean distribution plot exported to: {normalizePath(plots_output_dir)}")
    }
  }

  # -- Haploid Genotype Summary --

  haploid_summary <- haploid_genotype_assignment(design_summary_file,
                                                 assigned_genotypes,
                                                 lower_CI_assignment_threshold,
                                                 verbose)

  # -- Optional file saving (only if output_dir is specified) --

  if (!is.null(output_dir)) {
    haploid_summary_file_name <- paste0(sample_name, "_", ploidy, "_genotypes_summary.csv")

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    readr::write_csv(haploid_summary, file = file.path(
      output_dir, haploid_summary_file_name))

    cli::cli_alert_success(
      "{stringr::str_to_title(ploidy)} genotypes summary saved to: {normalizePath(output_dir)}.")
  }

  # -- Distribution of Cell Counts per Identity Barcode (only if plots_output_dir is specified) --
  if (!is.null(plots_output_dir)) {
    cell_distribution_plot <- cells_per_identity_barcode_plot(haploid_summary)

    cell_distribution_plot_file_name <- paste0(sample_name, "_cells_per_identity_barcode.png")

    # Optional file saving (only if plots_output_dir is provided)
    if (!dir.exists(plots_output_dir)) dir.create(plots_output_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(plots_output_dir, cell_distribution_plot_file_name),
                    plot = cell_distribution_plot,
                    width = 8, height = 5, dpi = 300)

    cli::cli_alert_success("{AmpID} posterior mean distribution plot exported to: {normalizePath(plots_output_dir)}")
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


#' @keywords internal
haploid_genotype_assignment <- function(
  design_summary_file,
  assigned_genotypes,
  lower_CI_assignment_threshold,
  verbose
) {
  if (verbose) {
    cli::cli_h2("Assigning Genotypes Across Amplicons")
    cli::cli_process_start(msg = "Collating genotypes...",
                           msg_done = "Genotypes collated.")
  }

  # Extract a list of edit genotypes
  edit_list <- design_summary_to_genotypes(design_summary_file)

  # Assemble a data frame of all cell_barcodes and their identity_barcode
  per_cell_genotypes <- dplyr::bind_rows(assigned_genotypes) |>
    dplyr::distinct(cell_barcode, identity_barcode)
  if (anyDuplicated(per_cell_genotypes[["cell_barcode"]]) > 0) {
    cli::cli_alert_warning("Some cell_barcodes have multiple identity_barcode assignments.")
  }

  # -- Pull across edit genotype assignments on every amplicon --

  WT_calls   <- NULL
  mutation_calls <- NULL

  for (ampID in names(assigned_genotypes)) {
    # Take only the most probable variant for cell genotype calling
    assigned_genotypes[[ampID]] <- assigned_genotypes[[ampID]] |>
      dplyr::group_by(cell_barcode) |>
      dplyr::slice_max(order_by = posterior_mean, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()

    amp_mutation_calls <- assigned_genotypes[[ampID]] |>
      # If variant_genotype is one of the edits in amplicon, and the lower
      # credible interval is >= the acceptable assignment threshold
      dplyr::mutate(not_WT = variant_genotype != "WT") |>
      dplyr::transmute(
        cell_barcode,
        identity_barcode,
        genotype = dplyr::if_else(not_WT, paste0(ampID, ":", variant_genotype), NA),
        CI_lower = dplyr::if_else(not_WT, CI_lower, NA)
      )
    colnames(amp_mutation_calls)[3:4] <- paste0(ampID, c("_genotype", "_CI_lower"))

    mutation_calls <- if(is.null(mutation_calls)) {
      amp_mutation_calls
    } else {
      dplyr::full_join(mutation_calls, amp_mutation_calls,
                       by = c("cell_barcode", "identity_barcode"))
    }

    # Checks if the cell is wildtype on this amplicon
    amp_WT_calls <- assigned_genotypes[[ampID]] |>
      dplyr::transmute(cell_barcode, WT = variant_genotype == "WT")
    colnames(amp_WT_calls)[2] <- ampID

    WT_calls <- if (is.null(WT_calls)) {
      amp_WT_calls
    } else {
      dplyr::full_join(WT_calls, amp_WT_calls, by = "cell_barcode")
    }
  }

  # -- Summarise confident WTs across amplicons --

  # Confident WT: no non-WT calls, and WT on at least 50% of amplicons
  WT_calls <- WT_calls |>
    dplyr::mutate(confident_WT =
                    rowSums(dplyr::pick(-cell_barcode) == FALSE, na.rm = TRUE) == 0 &
                    rowSums(dplyr::pick(-cell_barcode), na.rm = TRUE) /
                    ncol(dplyr::pick(-cell_barcode)) >= 0.5
    )

  # -- Summarise mutations and edits across amplicons --

  genotype_cols <- grep("_genotype$", colnames(mutation_calls), value = TRUE)
  CI_cols       <-  sub("_genotype$", "_CI_lower", genotype_cols)

  genotype_matrix <- as.matrix(mutation_calls[genotype_cols])
  CI_matrix       <- as.matrix(mutation_calls[CI_cols])

  # TRUE where the call on that amplicon is an inserted edit
  is_edit <- matrix(genotype_matrix %in% edit_list,
                    nrow = nrow(genotype_matrix),
                    ncol = ncol(genotype_matrix))

  # Index of the amplicon with the highest-CI *edit* in each cell
  top_edit <- max.col(replace(CI_matrix, !is_edit, -Inf), ties.method = "first")
  top_idx  <- cbind(seq_len(nrow(mutation_calls)), top_edit)

  # Pipe-separated inserted edits only
  edit_only <- genotype_matrix
  edit_only[!is_edit] <- NA
  all_edits <- as.data.frame(edit_only) |>
    tidyr::unite("all_edits", dplyr::everything(), sep = " | ", na.rm = TRUE) |>
    dplyr::pull(all_edits)

  # Feed back in to the output summary table
  mutation_calls <- mutation_calls |>
    # Attach in WT_calls tracker
    dplyr::left_join(dplyr::select(WT_calls, cell_barcode, confident_WT),
                     by = "cell_barcode") |>
    # all_mutations: pipe-separated mutations, "" when there are none
    tidyr::unite("all_mutations", dplyr::all_of(genotype_cols),
                 sep = " | ", na.rm = TRUE, remove = FALSE) |>
    dplyr::mutate(
      all_mutations = dplyr::na_if(all_mutations, ""),
      all_edits     = dplyr::na_if(all_edits, ""),
      n_edits       = rowSums(is_edit),
      most_probable_edit     = dplyr::if_else(n_edits == 0, NA, genotype_matrix[top_idx]),
      most_probable_CI_lower = dplyr::if_else(n_edits == 0, NA, CI_matrix[top_idx]),
      confident_edit = most_probable_CI_lower >= lower_CI_assignment_threshold
    )

  # Remove per-amplicon columns
  mutation_calls <- mutation_calls |>
    dplyr::select(-dplyr::all_of(as.vector(rbind(genotype_cols, CI_cols))))

  # each_amp_cols <- as.vector(rbind(genotype_cols, CI_cols))
  # mutation_calls <- mutation_calls |>
  #   dplyr::relocate(dplyr::all_of(each_amp_cols), .after = dplyr::last_col())

  if (verbose) {
    cli::cli_process_done()

    total_cells        <- nrow(mutation_calls)
    confident_WT_cells <- sum(mutation_calls[["confident_WT"]])
    mutant_cells       <- sum(!is.na(mutation_calls[["all_mutations"]]))
    single_edit_cells  <- sum(mutation_calls[["n_edits"]] == 1)
    multi_edit_cells   <- sum(mutation_calls[["n_edits"]] > 1)

    cli::cli_alert_info("Confident Wildtype Cells: {confident_WT_cells} / {total_cells}")
    cli::cli_alert_info("Mutant Cells: {mutant_cells} / {total_cells}")
    cli::cli_alert_info("Cells with One edit: {single_edit_cells} / {total_cells}")
    cli::cli_alert_info("Cells with Multiple edits: {multi_edit_cells} / {total_cells}")
  }

  return(mutation_calls)
}



#' @keywords internal
posterior_mean_distribution_plot <- function(
  AmpID,
  genotypes
) {
  plot_colour <- scales::brewer_pal(palette = "Set2")(1)

  # Plot the distribution of the posterior mean across AmpID cell_barcodes
  plot <- ggplot2::ggplot(genotypes, ggplot2::aes(x = posterior_mean)) +
    ggplot2::geom_histogram(binwidth = 0.025, boundary = 1,
                            fill = plot_colour, colour = "black") +
    ggplot2::labs(
      title = paste0(AmpID, ": Posterior Mean Distribution"),
      x = "Posterior Mean",
      y = "Frequency (Number of Cell Barcodes)"
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  return(plot)
}


#' @keywords internal
cells_per_identity_barcode_plot <- function(summary) {
  # Count the number of unique cell_barcodes per identity_barcode
  identity_barcode_cells <- summary |>
    dplyr::group_by(identity_barcode) |>
    dplyr::summarise(cells = length(unique(cell_barcode)))

  plot_colour <- scales::brewer_pal(palette = "Set2")(1)

  # Plot the distribution of the number of cell_barcodes per identity_barcode
  plot <- ggplot2::ggplot(identity_barcode_cells, ggplot2::aes(x = cells)) +
    ggplot2::geom_histogram(binwidth = 1, boundary = 0,
                            fill = plot_colour, colour = "black") +
    ggplot2::labs(
      title = "Number of Cells per Identity Barcode",
      x = "Number of Unique Cell Barcodes",
      y = "Frequency (Number of Identity Barcode Combinatios)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  return(plot)
}


