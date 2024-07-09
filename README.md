# CaseControlAF
Case Control Allele Frequency (AF) Reconstruction R Package

This repository contains the source code for the CaseControlAF R package which can be used to reconstruct the allele frequency (AF) for cases and controls separately given commonly available summary statistics. 

The package contains two functions:

1) CaseControl_AF
2) CaseControl_SE

## Download the package

To install this package using BioConductor:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CaseControlAF")
```

To download this package using *devtools* in R:

```R
require(devtools)
devtools::install_github("https://github.com/wolffha/CaseControlAF")
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
**data**: a dataframe with a row for each variant and columns for OR, SE, and optionally proxy MAFs

**N_case**: an integer for the number of case samples

**N_control**: an integer for the number of control samples

**OR_colname**: a string containing the exact column name in 'data' with the OR

**SE_colname**: a string containing the exact column name in 'data' with the SE

**proxyMAFs_colname**:  a string containing the exact column name in 'data' with the proxy MAFs to be used to correct the bias in the estimates. Default is NA - will only produce adjusted MAFs if not NA

Returns a dataframe with three columns with names: MAF_case, MAF_control and MAF_total containing the estimated minor allele frequency in the cases, controls, and total sample. The number of rows is equal to the number of variants. If proxyMAFs_colname is not NA, will include three additional columns containing the adjusted estimated MAFs (MAF_case_adj, MAF_control_adj, MAF_total_adj)

**NOTE:** This method assumes we are estimating the minor allele frequency (MAF)

### Examples

See examples here: (https://wolffha.github.io/CaseControlAF/)

