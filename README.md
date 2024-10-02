# CCAFE
Case Control Allele Frequency (AF) Reconstruction R Package

This repository contains the source code for the CaseControlAF R package which can be used to reconstruct the allele frequency (AF) for cases and controls separately given commonly available summary statistics. 

The package contains two functions:

1) CaseControl_AF
2) CaseControl_SE

See full documentation, vignettes, and examples here: (https://wolffha.github.io/CCAFE/)

## Download the package

To install this package using BioConductor (Not yet available):

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CCAFE")
```

To download this package using *devtools* in R:

```R
require(devtools)
devtools::install_github("https://github.com/wolffha/CCAFE/")
```

## CaseControl_AF

Use this function when you have the following statistics (for each variant)

* Number of cases
* Number of controls
* Odds Ratio (OR) or beta coefficient
* **AF** (allele frequency) for the total sample (cases and controls combined)

### Usage
**data**: a dataframe with a row for each variant and columns for OR and total AF

**N_case**: an integer for the number of case samples

**N_control**: an integer for the number of control samples

**OR_colname**: a string containing the exact column name in 'data' with the OR

**AF_total_colname**: a string containing the exact column name in 'data' with the total AF

Returns a dataframe with two columns: AF_case and AF_control. The number of rows is equal to the number of variants.

## CaseControl_SE
Use this function when you have the following statistics (for each variant)

* Number of cases
* Number of controls
* Odds Ratio (OR) or beta coefficient
* **SE** of the log(OR) for each variant

*Code adapted from ReACt GroupFreq function available here: (https://github.com/Paschou-Lab/ReAct/blob/main/GrpPRS_src/CountConstruct.c)*

### Usage
**data**: a dataframe where each row is a variant and columns for the OR, SE, chromosome, and position

**N_case**: an integer for the number of case samples

**N_control**: an integer for the number of control samples

**OR_colname**: a string containing the exact column name in *data* with the odds ratios

**SE_colname**: a string containing the exact column name in *data* with the standard errors

**position_colname**: a string containing the exact column name in *data* with the positions of the variants

**chromosome_colname**: a string containing the exact column name in *data* with the chromosome of the variants. 
Note, sex chromosomes can be either characters ('X', 'x', 'Y', 'y') or numeric where X=23 and Y=24

**sex_chromosomes**: boolean, TRUE if variants from sex chromosome(s) are included in the dataset

**do_correction**: boolean, TRUE if data is provided to correct the estimates using proxy MAFs

**remove_sex_chromosomes**: boolean, TRUE if variants on sex chromosomes should be removed. This is only necessary if *sex_chromosomes* == TRUE and the number of XX/XY individuals per case and control sample is NOT known

CaseControl_SE has the following optional inputs: 

If *sex_chromosomes* == TRUE and *remove_sex_chromosomes* == FALSE, then the following inputs are required:

**N_XX_case**: the number of XX chromosome case individuals

**N_XX_control**: the number of XX chromosome control individuals

**N_XY_case**: the number of XY chromosome case individuals

**N_XY_control**: the number of XY chromosome control individuals

If *do_correction* == TRUE, then data must be provided that includes harmonized data with proxy MAFs

**correction_data**: a dataframe with the following EXACT column names: CHR, POS, proxy_MAF, containing data for variants harmonized between the observed and proxy datasets

Returns the *data* dataframe with three additional columns with names: MAF_case, MAF_control and MAF_total containing the estimated minor allele frequency in the cases, controls, and total sample. The number of rows is equal to the number of variants. If proxyMAFs_colname is not NA, will include three additional columns containing the adjusted estimated MAFs (MAF_case_adj, MAF_control_adj, MAF_total_adj)

**NOTE:** This method assumes we are estimating the minor allele frequency (MAF)

### Examples and documentation

See full documentation, vignettes, and examples here: (https://wolffha.github.io/CaseControlAF/)

