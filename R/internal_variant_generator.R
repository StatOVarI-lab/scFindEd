#' @keywords internal
#' Generates all SNP and 1bp indel variants for a given WT/HDR input and
#' amplicon sequence template.
generate_all_variants <- function(
  WT_input,
  HDR_input,
  amplicon_input,
  all_reads = NULL,
  flank_len = 8,
  editing_window_extension = 0,
  verbose
) {
  BASES <- c("A", "C", "G", "T")
  
  if (verbose) cli::cli_process_start(msg = "Generating all variants...",
                                      msg_done = "All variants generated.")
  
  # Convert inputs to uppercase for consistency
  WT_template <- toupper(WT_input)
  HDR_template <- toupper(HDR_input)
  amplicon_template <- toupper(amplicon_input)
  
  # -- Optional: Apply editing_window_extension --
  
  if (editing_window_extension > 0) {
    # Locate where WT_template sits inside amplicon_template
    match_pos <- regexpr(WT_template, amplicon_template)[1]
    if (match_pos == -1) {
      cli::cli_abort("WT_input not found in amplicon_input.")
    }
    
    # Define and extract editing_window_extension flanks
    upstream_extend_start <- max(1, match_pos - editing_window_extension)
    upstream_extend_end <- match_pos - 1
    downstream_extend_start <- match_pos + nchar(WT_template)
    downstream_extend_end <- min(
      nchar(amplicon_template), 
      downstream_extend_start + editing_window_extension - 1)
    
    upstream_extend <- substr(amplicon_template, upstream_extend_start, upstream_extend_end)
    downstream_extend <- substr(amplicon_template, downstream_extend_start, downstream_extend_end)
    
    # Redefine WT_ and HDR_templates to include editing_window_extension flanks
    WT_template <- paste0(upstream_extend, WT_template, downstream_extend)
    HDR_template <- paste0(upstream_extend, HDR_template, downstream_extend)
  }
  
  # -- Extract Position of (extended) WT_template
  
 # Locate where WT_template sits inside amplicon_template
  match_pos <- regexpr(WT_template, amplicon_template)[1]
  if (match_pos == -1) {
    cli::cli_abort("WT_input not found in amplicon_input.")
  } 
  
  # Define and extract editing_window flanks
  upstream_start <- max(1, match_pos - flank_len)
  upstream_end <- match_pos - 1
  downstream_start <- match_pos + nchar(WT_template)
  downstream_end <- 
    min(nchar(amplicon_template), downstream_start + flank_len - 1)
  
  upstream_flank <- substr(amplicon_template, upstream_start, upstream_end)
  downstream_flank <- substr(amplicon_template, downstream_start, downstream_end)
  
  # -- Identify SNP and Protected Region --
  
  amplicon_length <- nchar(amplicon_template)
  
  SNP_pos <- which(strsplit(WT_template, "")[[1]] != strsplit(HDR_template, "")[[1]])
  SNP_window <- lapply(SNP_pos, function(p) {
    seq(from = max(1, p-1),
        to   = min(amplicon_length, p+1))
  })
  protected_pos <- sort(unique(unlist(SNP_window)))
  
  # -- Generate WT & HDR Variants --
  
  df_WT_variants <- generate_variants(WT_template, "WT", protected_pos, verbose)
  df_HDR_variants <- generate_variants(HDR_template, "HDR", protected_pos, verbose)
  
  # Combine WT & HDR variant data frames, and generate complete sequences
  df_variants <- rbind(df_WT_variants, df_HDR_variants)
  df_variants$full_sequence <- paste0(upstream_flank, df_variants$sequence, downstream_flank)
  df_variants <- df_variants[!duplicated(df_variants$full_sequence),]
  
  # -- Add Sequences from all_reads if Provided --
  
  if (!is.null(all_reads)) {
    inner_seqs <- sapply(all_reads$merged_sequence,
                         extract_between_flanks,
                         upstream_flank,
                         downstream_flank)
    
    # Remove NA instances, and duplicates
    inner_seqs <- unique(inner_seqs[!is.na(inner_seqs)])
    
    # Add unique new sequences to df_variants
    unique_inner_seqs <- setdiff(inner_seqs, df_variants$sequence)
    if (length(unique_inner_seqs) > 0) {
      df_random <- data.frame(
        sequence = unique_inner_seqs,
        genotype_group = "random",
        mutation = unique_inner_seqs,
        full_sequence = paste0(upstream_flank, unique_inner_seqs, downstream_flank)
      )
      
      df_variants <- rbind(df_variants, df_random)
    }
  }
  
  # -- Add genotype Column --
  df_variants$genotype <- ifelse(
    is.na(df_variants$genotype_group),
    df_variants$full_sequence,
    paste(df_variants$genotype_group, df_variants$mutation, sep = "_")
  )
  
  if (verbose) {
    cli::cli_process_done()
    cli::cli_alert_info("Total variants: {nrow(df_variants)}")
  }
  
  return(df_variants)
}


#' @keywords internal
#' Generates all variants of a given editing window, assuming a genotype and a
#' set of protected positions (a window around the SNP) are provided.
generate_variants <- function(
  template, 
  genotype_prefix, 
  protected_pos,
  verbose
) {
  BASES <- c("A", "C", "G", "T")
  
  template_length <- nchar(template)
  
  # -- Create reference data frame --
  
  df_ref <- data.frame(
    sequence = template,
    genotype_group = genotype_prefix,
    mutation = "reference",
    stringsAsFactors = FALSE
  )
  
  # -- Generate SNP Possibilities
  
  SNP_list <- c()
  SNP_mutations <- c()
  
  for (i in 1:template_length) {
    if (!(i %in% protected_pos)) {
      ref_base <- substr(template, i, i) 
      for (base in setdiff(BASES, ref_base)) {
        SNP_seq <- paste0(substr(template, 1, i-1), base, substr(template, i+1, template_length))
        SNP_list <- c(SNP_list, SNP_seq)
        SNP_mutations <- c(SNP_mutations, paste0(ref_base, "->", base, " at ", i))
      }
    }
  }
  
  df_SNP <- data.frame(
    sequence = SNP_list,
    genotype_group = paste0(genotype_prefix, "_SNP"),
    mutation = SNP_mutations,
    stringsAsFactors = FALSE
  )
  
  # -- Generate 1bp Indel Possibilities
  
  del_list <- c()
  del_mutations <- c()
  
  for (i in 1:template_length) {
    if (!(i %in% protected_pos)) {
      del_seq <- paste0(substr(template, 1, i-1), substr(template, i+1, template_length))
      del_list <- c(del_list, del_seq)
      del_mutations <- c(del_mutations, paste0("del at ", i))
    }
  }
  
  ins_list <- c()
  ins_mutations <- c()
  
  for (i in 1:template_length) {
    if (i < min(protected_pos) || i >= max(protected_pos)) {
      for (b in BASES) {
        ins_seq <- paste0(substr(template, 1, i), b, substr(template, i+1, template_length))
        ins_list <- c(ins_list, ins_seq)
        ins_mutations <- c(ins_mutations, paste0("ins ", b, " at ", i+1))
      }
    }
  }
  
  df_indel <- data.frame(
    sequence = c(del_list, ins_list),
    genotype_group = paste0(genotype_prefix, "_indel"),
    mutation = c(del_mutations, ins_mutations),
    stringsAsFactors = FALSE
  )
  
  # -- Combine Variant Data Frames and Add Flanks
  
  df_variants <- rbind(df_ref, df_SNP, df_indel)
  
  return(df_variants)
}

#' @keywords internal
#' Wrapper function for assigning all given read sequences a genotype by matching
#' with a table of generated variants.
map_genotypes <- function(
  sequences_df,
  variants_df,
  sequence_col = "merged_sequence",
  n_cores = 1,
  verbose
) {
  BASES <- c("A", "C", "G", "T")
  
  if (verbose) cli::cli_progress_bar("Assigning genotypes to reads", total = 20)
  
  # Assign the number of threads
  setDTthreads(n_cores)
  
  # Convert to data frames to data.tables
  sequences_dt <- as.data.table(sequences_df)
  variants_dt <- as.data.table(variants_df)
  
  # Ensure sequences are upper case
  set(sequences_dt, j = sequence_col, value = toupper(sequences_dt[[sequence_col]]))
  set(variants_dt, j = "full_sequence", value = toupper(variants_dt$full_sequence))
  
  # Create an empty genotype column
  set(sequences_dt, j = "genotype", value = NA_character_)
  
  n_variants <- nrow(variants_dt)
  progress_bar_interval <- ceiling(n_variants / 20) # 5% increments
  
  for (i in 1:n_variants) {
    target_sequence <- variants_dt$full_sequence[i]
    genotype_value <- variants_dt$genotype[i]
    
    # Don't check empty target_sequences
    if (is.na(target_sequence) || target_sequence == "") next
    
    # Parallelised matching and assignment
    sequences_dt[stringi::stri_detect_fixed(get(sequence_col), 
                                            target_sequence) & is.na(genotype),
                 genotype := genotype_value]
    
    # Update progress_bar
    if (i %% progress_bar_interval == 0) cli::cli_progress_update()
  }
  
  cli::cli_progress_done()

  if (verbose) cli::cli_alert_success("Genotypes assigned to reads.")
  
  return(as.data.frame(sequences_dt))
}


#' @keywords internal
#' Helper function for extracting the sequence of bases that lie between the 
#' upstream and downstream flanks.
extract_between_flanks <- function(
  sequence,
  upstream_flank,
  downstream_flank
) {
  # Ensure sequence is upper case
  sequence <- toupper(sequence)
  
  upstream_hits <- gregexpr(upstream_flank, sequence, fixed = TRUE)[[1]]
  downstream_hits <- gregexpr(downstream_flank, sequence, fixed = TRUE)[[1]]
  
  if (upstream_hits[1] == -1 || downstream_hits[1] == -1) {
    return(NA)
  }
  
  if (length(upstream_hits) == 1 && length(downstream_hits) == 1 &&
      upstream_hits[1] > 0 && downstream_hits[1] > upstream_hits[1]) {
    start_pos <- upstream_hits[1] + nchar(upstream_flank)
    end_pos <- downstream_hits[1] - 1
    return(substr(sequence, start_pos, end_pos))
  } else {
    return(NA)
  }
}

