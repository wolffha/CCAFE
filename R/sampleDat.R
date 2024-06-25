#' PanUKBB and gnomAD Diabetes Data
#'
#' This is a subset of 500 variants on chromosome 1 from the Pan-UKBB diabetes GWAS with the whole sample (pop), case, and control minor allele frequency (MAF) for those classified in Pan-UKBB as European (EUR). This variants (which are mapped to GRCh37) have been harmonized with gnomAD non-Finnish European (NFE) MAFs.
#'
#' @format ## `sampleDat`
#' A data frame with 500 rows and 11 columns:
#' \describe{
#'   \item{CHR}{chromosome number}
#'   \item{POS}{base-pair position of variant (GRCh37 coordinates)}
#'   \item{REF}{Reference allele}
#'   \item{ALT}{Alternate allele}
#'   \item{true_maf_case}{MAF in EUR cases in Pan-UKBB Diabetes}
#'   \item{true_maf_control}{MAF in EUR controls in Pan-UKBB Diabetes}
#'   \item{true_maf_pop}{MAF in whole EUR sample in Pan-UKBB Diabetes}
#'   \item{beta}{beta from EUR GWAS in Pan-UKBB Diabetes}
#'   \item{SE}{SE of beta from EUR GWAS in Pan-UKBB Diabetes}
#'   \item{OR}{OR from EUR GWAS in Pan-UKBB Diabetes}
#'   \item{gnomad_maf}{MAF in gnomAD NFE}
#' }
#' @source <https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files>
#' @source <https://gnomad.broadinstitute.org/downloads>
"sampleDat"
