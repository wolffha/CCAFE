
#' @title CaseControl_AF
#' @description
#' This is a function to derive the case and control AFs from GWAS summary statistics when
#' the user has access to the whole sample AF, the sample sizes, and the OR (or beta).
#' If user has SE instead of sample AF use [CaseControlAF::CaseControl_SE()]
#'
#' @param data dataframe with each row being a variant and columns for AF_total and OR
#' @param N_case the number of cases in the sample
#' @param N_control the number of controls in the sample
#' @param OR_colname a string containing the exact column name in 'data' with the OR
#' @param AF_total_colname a string containing the exact column name in 'data' with the whole sample AF
#'
#' @return returns a dataframe with two columns (AF_case, AF_control) and rows
#' equal to the number of variants
#'
#' @author Hayley Wolff (Stoneman), \email{hayley.wolff@cuanschutz.edu}
#'
#' @references https://github.com/wolffha/CaseControlAF
#'
#' @seealso \url{https://github.com/wolffha/CaseControlAF} for further documentation
#'
#' @examples
#' library(CaseControlAF)
#'
#' data("sampleDat")
#' sampleDat <- as.data.frame(sampleDat)
#'
#' nCase_sample = 16550
#' nControl_sample = 403923
#'
#' # get the estimated case and control AFs
#' af_method_results <- CaseControl_AF(data = sampleDat,
#'                                     N_case = nCase_sample,
#'                                     N_control = nControl_sample,
#'                                     OR_colname = "OR",
#'                                     AF_total_colname = "true_maf_pop")
#'
#' head(af_method_results)
#'
#' @export
CaseControl_AF <- function(data, N_case = 0, N_control = 0, OR_colname = "OR", AF_total_colname = "AF"){

  # do input checking

  # check valid input for case/control sample size
  if(N_case <= 0) {
    stop("ERROR: 'N_case' needs to be a number > 0")
  }

  if(N_control <= 0) {
    stop("ERROR: 'N_control' needs to be a number > 0")
  }

  # check valid input data type
  if(!is.data.frame(data)) {
    stop("ERROR: 'data' must be a dataframe")
  }

  # check that OR and AF_total columns exist
  if(!OR_colname %in% colnames(data)) {
    stop("ERROR: 'OR_colname' does not exist in 'data'")
  }

  if(!AF_total_colname %in% colnames(data)) {
    stop("ERROR: 'AF_total_colname' does not exist in 'data'")
  }

  OR <- data[,c(OR_colname)]
  AF_total <- data[,AF_total_colname]

  # check for valid input for OR and AF_total
  if(typeof(OR) != "double") {
    stop("ERROR: 'OR' values must all be numbers (hint: check for NAs)")
  }

  if(typeof(AF_total) != "double") {
    stop("ERROR: 'AF_total' values must all be numbers (hint: check for NAs)")
  }

  if(any(AF_total < 0)) {
    stop("ERROR: 'AF_total' cannot contain negative AFs")
  }

  if(any(AF_total > 1)) {
    stop("ERROR: 'AF_total' cannot contain values > 1")
  }

  # solve for real roots of a quadratic using the quadratic formula
  quad_roots<-function(a,b,c){
    c(((-b-sqrt(b^2-4*a*c))/(2*a)),((-b+sqrt(b^2-4*a*c))/(2*a)))
    # You could add a return statement to this function and then return a vector of the roots
  }

  # calculate total sample size
  N_total <- N_control+N_case

  # set a, b, c of quadratic equation derived in manuscript
  a <- (N_control/N_case)*(OR-1)
  b <- (OR-((N_total/N_case)*AF_total*OR))+((N_control/N_case)+(N_total*AF_total/N_case))
  c <- -(N_total/N_case)*AF_total

  AF_control <- rep(0, length(a)) # vector to store the control AFs

  for(i in seq_len(length(a))) {
    # find roots of quadratic eq and choose root between 0 and 1 as AF_control
    AF_control_opts <-  quad_roots(a[i], b[i], c[i])
    if(AF_control_opts[1]>1 | AF_control_opts[1]<0){
      AF_control[i] <- AF_control_opts[2]
    }else{
      AF_control[i] <- AF_control_opts[1]
    }
  }

  # calculate AF_case with known relationship shown in manuscript
  AF_case <- (N_total/N_case)*AF_total - (N_control/N_case)*AF_control

  data$AF_case <- AF_case
  data$AF_control <- AF_control

  # Output shows case AF first, then control AF
  return(data)
}
