#' Global Parameter Documentation
#'
#' @param source_dir String containing the relative path to the directory that
#' source data should be imported from.
#' @param sample_name String containing the exact file name (including the file
#' type) for the experiment sample.
#' @param panel_ID String containing the ID number for the design summary panel.
#' @param identity_AmpID String or vector containing the amplicon ID(s) (AmpID)
#' for the amplicon(s) that contains the inserted identity barcode(s).
#' @param min_reads Quality control threshold: minimum number of reads.
#' @param barcoded Boolean that tells the function whether the experiment
#' utilises an identity barcode system.
#' @param ploidy Character string that defines whether genotyping shoud be
#' consider as haploid or diploid. Note that only haploidy is fully implemented
#' currently.
#' @param lower_CI_assignment_threshold Float that sets the minimum lower credible
#' interval bound for assigning a confident haploid genotype.
#' @param output_dir String containing the relative path to the directory that
#' the output data should be exported to.
#' @param plots_output_dir String containing the relative path to the directory
#' that the output quality control plots should be exported to.
#' @param verbose Boolean that determines whether detailed outputs should be
#' printed.
#'
#' @name global_params
#' @keywords internal
NULL
