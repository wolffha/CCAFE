#' @title CCAFE_convertVCF
#' @description
#' Formats information from a VCF object for use in CCAFE methods as follows:
#' From the rowRanges object: seqnames (chromosome), ranges (position),
#; REF, ALT, RSID
#' From the geno object: ES (effect size of ALT), SE, AF (allele frequency
#' of ALT)
#'
#' @param vcf a Variant Call Format (VCF) file read in using VariantAnnotation
#' BioConductor package
#' @return a dataframe object with columns Position, RSID, Chromosome, REF,
#' ALT, beta, SE, AF, OR
#' @export
#' @author Hayley Wolff (Stoneman), \email{hayley.wolff@cuanschutz.edu}
#' @examples
#' library(VariantAnnotation)
#' library(CCAFE)
#'
#' # load the data
#' data("vcf_sample")
#'
#' # run the method
#' df_sample <- CCAFE_convertVCF(vcf_sample)
#' print(head(df_sample))
#'
#' # can then use in CCAFE methods
#' # since we have total AF, will use CaseControl_AF
#' df_sample <- CaseControl_AF(data = df_sample,
#'                           N_case = 48286,
#'                           N_control = 250671,
#'                           OR_colname = "OR",
#'                           AF_total_colname = "AF")
#' head(df_sample)

CCAFE_convertVCF <- function(vcf) {
  # make sure object is VariantAnnotation VCF
  if(!inherits(vcf, "CollapsedVCF") && !inherits(vcf, "ExpandedVCF")) {
    stop("The input object is not a valid VCF (CollapsedVCF or ExandedVCF)
    from the VariantAnnotation package")
  } else {
    message("Valid VCF object, converting to dataframe for CCAFE...")
  }

  # first we will get the info from GRanges object (position, RSID)
  meta <- as.data.frame(IRanges::ranges(vcf))
  meta <- meta[,c(1, 4)]
  colnames(meta) <- c("Position", "RSID")
  meta$Chromosome <- as.vector(GenomicRanges::seqnames(
    SummarizedExperiment::rowRanges(vcf)))

  # now we can also get the meta data (REF, ALT) from the GRanges object
  meta <- cbind(meta, S4Vectors::mcols(vcf)[,c(2,3)])
  rownames(meta) <- seq(1, nrow(meta))

  # now we will get the info from the geno object
  geno_dat <- data.frame(
    beta = unlist(VariantAnnotation::geno(vcf)$ES),
    SE = unlist(VariantAnnotation::geno(vcf)$SE),
    AF = unlist(VariantAnnotation::geno(vcf)$AF)
  )

  df_data <- cbind(meta, geno_dat)

  # use OR (not beta) in CCAFE so we will create that column
  df_data$OR <- exp(df_data$beta)

  df_data <- as.data.frame(df_data)

  return(df_data)
}
