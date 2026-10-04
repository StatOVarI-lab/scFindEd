[![Stable Version - v0.1.0](https://img.shields.io/badge/stable_version-v.0.1.0-blue?logo=github)](https://github.com/StatOVarI-lab/scFindEd)
[![Dev Version - v0.1.0](https://img.shields.io/badge/dev_version-v.0.1.0-blue?logo=github)](https://github.com/StatOVarI-lab/scFindEd/tree/dev)

[![Documentation](https://img.shields.io/badge/documentation-white?logo=googledocs)](https://github.com/StatOVarI-lab/scFindEd/blob/main/vignettes/primary-vignette.Rmd)

## Overview 

*A toolkit for identity barcode extraction and genotyping of single-cell-based base editing experiments.* 

`scFindEd` provides a suite of functions that form a pipeline for barcoding and genotyping paired-end amplicon sequencing data. It is the complementary R package for `scFindEd-nf` - the upstream processing pipeline that aligns the raw sequencing reads and transforms the data ready for `scFindEd` analysis. There are three primary purposes that `scFindEd` serves, each forming a dedicated branch of the package. 

1. Barcode Extraction - based upon experiment summary data and pre-processed identity amplicon sequencing reads, each reads identity barcodes are extracted and used to assign each droplet barcode its identity barcode as a 'barcode system'. 
2. Genotyping - again based upon experiment summary data and all remaining pre-processed sequencing reads, each read is analysed for consistent variations from the reference. and ultimately each droplet barcode is assigned a genotype for each possible amplicon. Genotyping has been made agnostic to Barcode Extraction, meaning `scFindEd` can be used to genotype single-cell paired-end sequencing experiments that do not use a lineage tracing barcode system. 
3. Synthetic Experiment Generation - in order to validate the performance of `scFindEd`, a synthetic data generation pipeline has been implemented. The synthetic experiment generated is tailor-made to be run through `scFindEd-nf` pre-processing, and generates a complete ground-truth table for validating both identity barcode assignment and genotyping.

Users of `scFindEd` may also opt for several quality-control plots to be produced throughout the analysis pipeline. Examples of these may be found in `/plots/`.

Key experimental methodologies that `scFindEd` has been developed as an accompaniment for include both scSNV-seq $^{[1]}$ and Saturation-seq $^{[2]}$.

## Installation 

The package remains in active development and is not yet public. It relies on the `devtools` and `usethis` packages for installation and setup.  
### For developers working locally: 
``` 
install.packages("devtools")
library(devtools)

devtools::install()
library(scFindEd)

# To learn more about an exported function
help("function_name")
```
### Dependencies
Please see below for a complete list of `scFindEd` package dependencies. This information can also be found in `DESCRIPTION`.

```
CRAN-based:
- cli
- dplyr
- data.table
- ggplot2
- patchwork
- purrr
- readr
- rlang
- scales
- stringdist
- stringi
- stringr
- tibble
- tidyr
- truncnorm

Bioconductor-based:
- Biostrings
- cigarillo
- IRanges
- S4Vectors
```

## scFindEd Workflow

A vignette that walks through the application of each subsequent function in the Barcode Extraction $\rightarrow$ Genotyping pipeline can be found in `/vignettes/primary-vignette.Rmd`.

```
scFindEd
1. Barcode Extraction                   2. Genotyping  
├── extract_identity_barcode_info()     ├── summarise_variants()
├── evaluate_optimal_barcode_flanks()   ├── identify_true_mutations()
├── process_identity_barcodes()         ├── plot_variant_summary()
├── identify_whitelisted_barocdes()     └── genotype_cells()
└── assign_identity_barcodes()

2. Synthetic Data Generation
└── generate_synthetic_experiment()
```

## Data Structure

In order for `scFindEd` to access both input data files as well as the intermediary files that `scFindEd` produces, we advise creating a directory such as the below, and instructing all functions to use this directory with both `source_dir` and `output_dir`.

```
data/ 
  ├── fragment_mutations/ 
  │     ├── {sample_name}_{AmpID}_fragment_mutations.tsv.gz
  │     ├── {sample_name}_{AmpID}_fragment_mutations.tsv.gz
  │     ├── ...
  │     ├── {sample_name}_{AmpID}_fragment_mutations.tsv.gz
  │     └── {sample_name}_{AmpID}_fragment_mutations.tsv.gz
  │
  ├── {sample_name}_{identity_AmpID}_reads.tsv.gz 
  │
  ├── {panel_ID}_design_summary_joined.csv
  │
  ├── ref_amplicons_sorted.fa
  │
  ├── {BC1_whitelist_file_name}
  ├── {BC2_whitelist_file_name}
  ├── {BC3_whitelist_file_name}
  └── ... # etc. for however many barcodes are in your identity barcode system
```

All files except for `BC1_whitelist_file_name`, `BC2_whitelist_file_name` etc. should be produced by the Nextflow, and therefore should have standardised data structure. `BC1_whitelist_file_name`, `BC2_whitelist_file_name` etc. should be `.csv` files with the following structures:

```
{BC1_whitelist_file_name} e.g. "75K_barcode.csv"

identity_barcode_1A
identity_barcode_1B
identity_barcode_1C
identity_barcode_1D
...


{BC2_whitelist_file_name} e.g. "BC2_whitelist_barcodes.csv"

identity_barcode_2A
identity_barcode_2B
identity_barcode_2C
identity_barcode_2D
...
```

## Development Notes 

This package is currently under active development; the licence is pending and will be added soon. 

When adding new code: 
- Run `devtools::document()` to update documentation and NAMESPACE. 
- Run `devtools::load_all()` to test interactively. 
- Run `usethis::use_package("package_name")` to add a new package dependency.
- Run `usethis::use_test("function_name")` to add new tests. 
- Run `devtools::check()` to run all checks before committing.

Please feel free to suggest changes or point out bugs using GitHub's [Issues](https://github.com/StatOVarI-lab/scFindEd/issues) system, and we will aim to address them where possible.

## References
[^1] Cooper SE, Coelho MA, Strauss ME, Gontarczyk AM, Wu Q, Garnett MJ, Marioni JC, Bassett AR. scSNV-seq: high-throughput phenotyping of single nucleotide variants by coupled single-cell genotyping and transcriptomics. Genome Biol. 2024 Jan 15;25(1):20. doi: 10.1186/s13059-024-03169-y.
[^2] Strauss ME, Waters AJ, Robertson H, Brendler-Spaeth T, Gontarczyk A, Gupta P, Kataria S, Gitterman D, Ntereke T, Wells L, Billington J, Bassett A, Cooper SE, Adams DJ. Saturation-seq integrates single-cell saturation genome editing and RNA-seq to quantify _NFE2L2_ (NRF2) variant effects. bioRxiv. 2026 Jul. doi: 10.64898/2026.06.30.735631.


