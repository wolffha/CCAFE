
#' @title CaseControl_AF
#' @description
#' This is a function to derive the case and control AFs from GWAS summary
#' statistics when the user has access to the whole sample AF, the sample sizes,
#' and the OR (or beta).
#' If user has SE instead of sample AF use [CCAFE::CaseControl_SE()]
#'
#' @param data dataframe with each row being a variant and columns for
#' AF_total and OR
#' @param N_case the number of cases in the sample
#' @param N_control the number of controls in the sample
#' @param OR_colname a string containing the exact column name in 'data'
#' with the OR
#' @param AF_total_colname a string containing the exact column name in 'data'
#' with the whole sample AF
#'
#' @return returns a dataframe with two columns (AF_case, AF_control) and rows
#' equal to the number of variants
#'
#' @author Hayley Wolff (Stoneman), \email{hayley.wolff@cuanschutz.edu}
#'
#' @references https://github.com/wolffha/CCAFE
#'
#' @seealso \url{https://github.com/wolffha/CCAFE} for further documentation
#'
#' @examples
#' library(CCAFE)
#'
# data("sampleDat")
# sampleDat <- as.data.frame(sampleDat)
# 
# nCase_sample = 16550
# nControl_sample = 403923
# 
# sampleDat$N_case <- sample(1000:2000, nrow(sampleDat), replace = TRUE)
# sampleDat$N_control <- sample(1000:2000, nrow(sampleDat), replace = TRUE)
# 
# 
# # get the estimated case and control AFs
# af_method_results <- CaseControl_AF(data = sampleDat,
#                                     N_case = nCase_sample,
#                                     N_control = nControl_sample,
#                                     OR_colname = "OR",
#                                     AF_total_colname = "true_maf_pop",
#                                     N_variant_colname = TRUE,
#                                     N_case_colname = "N_case",
#                                     N_control_colname = "N_control")
# 
# 
# head(af_method_results)
#
#' @export
#' 


CaseControl_AF <- function(data,
                           N_case = 0,
                           N_control = 0,
                           OR_colname = "OR",
                           AF_total_colname = "AF",
                           N_variant_colname = FALSE,
                           N_case_colname = NULL,
                           N_control_colname = NULL){

  data <- as.data.frame(data)
  # do input checking
  n_var <- nrow(data)
  
  if (N_variant_colname) {
    
    if (is.null(N_case_colname) || is.null(N_control_colname)) {
      stop("When setting the variant column name value to TRUE, 'N_case_colname' and 'N_control_colname' must be provided")
    }
    
    if (!N_case_colname %in% colnames(data)) {
      stop(paste0("'", N_case_colname, "' not found in data"))
    }
    
    if (!N_control_colname %in% colnames(data)) {
      stop(paste0("'", N_control_colname, "' not found in data"))
    }
    
    N_case <- data[[N_case_colname]]
    N_control <- data[[N_control_colname]]
  }
  
  # N_case
  if(!is.numeric(N_case)) {
    stop("'N_case' must be numeric")
  }
  if(length(N_case) == 1) {
    N_case <- rep(N_case, n_var)
  }
  if(length(N_case) != n_var) {
    stop("'N_case' must be length 1 or nrow(data)")
  }
  if(any(N_case <= 0, na.rm = TRUE)) {
    stop("'N_case' must contain values > 0")
  }
  
  # N_control
  if(!is.numeric(N_control)) {
    stop("'N_control' must be numeric")
  }
  if(length(N_control) == 1) {
    N_control <- rep(N_control, n_var)
  }
  if(length(N_control) != n_var) {
    stop("'N_control' must be length 1 or nrow(data)")
  }
  if(any(N_control <= 0, na.rm = TRUE)) {
    stop("'N_control' must contain values > 0")
  }

  # check valid input data type
  if(!is.data.frame(data)) {
    stop("'data' must be a dataframe")
  }

  # check that OR and AF_total columns exist
  if(!OR_colname %in% colnames(data)) {
    stop("'OR_colname' does not exist in 'data'")
  }

  if(!AF_total_colname %in% colnames(data)) {
    stop("'AF_total_colname' does not exist in 'data'")
  }

  OR <- data[,c(OR_colname)]
  AF_total <- data[,AF_total_colname]

  # check for valid input for OR and AF_total
  if(typeof(OR) != "double") {
    stop("'OR' values must all be numbers (hint: check for NAs)")
  }

  if(typeof(AF_total) != "double") {
    stop("'AF_total' values must all be numbers (hint: check for NAs)")
  }

  if(any(AF_total < 0)) {
    stop("'AF_total' cannot contain negative AFs")
  }

  if(any(AF_total > 1)) {
    stop("'AF_total' cannot contain values > 1")
  }

    #calculate total sample size
  N_total <- N_control+N_case

  #set a, b, c of quadratic equation derived in manuscript
  a <- (N_control/N_case)*(OR-1)
  b <- (OR-((N_total/N_case)*AF_total*OR))+((N_control/N_case)+
                                              (N_total*AF_total/N_case))
  c <- -(N_total/N_case)*AF_total

  # determine AF_control by applying quadratic formula, selecting root [0,1]
  # as AF_control
  AF_control <- vapply(seq_along(a), function(i) {
    AF_control_opts <- quad_roots(a[i], b[i], c[i])
    if (AF_control_opts[1] > 1 | AF_control_opts[1] < 0) {
      return(AF_control_opts[2])
    } else {
      return(AF_control_opts[1])
    }
  }, c(0.0))

  #calculate AF_case with known relationship shown in manuscript
  AF_case <- (N_total/N_case)*AF_total - (N_control/N_case)*AF_control

  data$AF_case <- AF_case
  data$AF_control <- AF_control

  #Output shows case AF first, then control AF
  return(data)
}

# solve for real roots of a quadratic using the quadratic formula
quad_roots<-function(a,b,c){
  c(((-b-sqrt(b^2-4*a*c))/(2*a)),((-b+sqrt(b^2-4*a*c))/(2*a)))
}

