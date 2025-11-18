# scFindEd

## Overview 

*A toolkit for analysing barcodes.* 

`scFindEd` provides a set of functions for processing and analysing barcode amplicon sequencing data. It currently supports tasks such as reading FASTQ files and will later include additional barcode analysis and visualisation tools.  

## Installation 

The package is in active development and not yet public. It relies on the `devtools` and `usethis` packages for installation and setup.  

For developers working locally: 
 
``` 
devtools::install()
library(scFindEd)
```

## Functionalitites 

`process_BC.R` 

Overview: Briefly describe what this function does (e.g. processes barcode FASTQ files and extracts relevant sequences).

How to use: Add a short example or command

Output: Describe the expected result

Issues: Mention any current limitations or known issues.

Future implementations: List planned extensions or improvements.

`other_function.R` 

 

## Example Data 

The package includes a lightweight dataset (`example_barcodes`) for demonstration and testing. 

Larger datasets used for internal testing (e.g. FASTQ files and barcode reference tables) are not included in the GitHub repository due to size but can be accessed locally by the developers. 

You can add a short note here explaining the expected input data structure e.g. file format, column names, etc.

## Development Notes 

- The package is currently under active development. 

- For reproducibility, dependencies are listed in `DESCRIPTION`; we may adopt `renv` later. 

- The licence is pending and will be added soon. 

- Local testing uses the path/to/directories directories, which are not tracked on GitHub. 

- When adding new code: 

      - run `devtools::document()` to update documentation and NAMESPACE. 

      - run `devtools::load_all()` to test interactively. 

      - use `usethis::use_test("function_name")` to add new tests. 

      - use `devtools::check()` to run all checks before commiting. 

 

[Add other notes, thoughts, or future plans here.] 



