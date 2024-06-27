#' @title CaseControl_SE
#' @description This is a function to derive the case and control AFs from GWAS summary statistics when
#' the user has access to the whole sample AF, the sample sizes, and the OR (or beta).
#' If user has SE instead of sample AF use [CaseControlAF::CaseControl_SE()]
#' This code uses the GroupFreq function adapted from C from <https://github.com/Paschou-Lab/ReAct/blob/main/GrpPRS_src/CountConstruct.c>
#'
#' @param N_case an integer of the number of Case individuals
#' @param N_control an integer of the number of Control individuals
#' @param OR a vector of odds ratios
#' @param SE a vector of standard errors for the OR calculation *make sure the indices match between the two vectors
#' @param proxyMAFs a vector of minor allele frequencies (MAFs) for the true whole sample MAF to be used to correct the bias in the estimates. Default is NA - will only produce adjusted MAFs if not NA
#'
#' @note
#'  Make sure the vectors are in the same order for your variants
#'
#' @return a dataframe with 3 columns: MAF_case, MAF_control, MAF_pop for the estimated MAFs for each variant. If proxyMAFs is not NA, then will output 3 additional columns (MAF_case_adj, MAF_control_adj, MAF_pop_adj) with the adjusted estimates.
#'
#' @author Hayley Wolff (Stoneman), \email{hayley.wolff@cuanschutz.edu}
#'
#' @references https://github.com/wolffha/CaseControlAF
#'
#' @seealso \url{https://github.com/wolffha/CaseControlAF} for further documentation
#'
#' @examples
#' data("sampleDat")
#'
#' nCase_sample = 16550
#' nControl_sample = 403923
#'
#' # get the estimated case and control MAFs
#' se_method_results <- CaseControl_SE(N_case = nCase_sample,
#'                                     N_control = nControl_sample,
#'                                     OR = sampleDat$OR,
#'                                     SE = sampleDat$SE)
#'
#' head(se_method_results)
#'
#' @import genpwr
#' @import dplyr
#' @export
CaseControl_SE <- function(N_case, N_control, OR, SE, proxyMAFs = NA) {
  # genpwr package is used to solve for the roots of the quadratic equation
  #require(genpwr)

  # do input checking
  if(length(OR) != length(SE)) {
    stop("ERROR: 'OR' and 'SE' vectors need to be the same length")
  }
  if(typeof(OR) != "double") {
    stop("ERROR: 'OR' vector needs to be a vector of numbers (hint: check for NAs)")
  }
  if(typeof(SE) != "double") {
    stop("ERROR: 'SE' vector needs to be a vector of numbers (hint: check for NAs)")
  }
  if(N_case <= 0) {
    stop("ERROR: 'N_case' needs to be a number > 0")
  }
  if(N_control <= 0) {
    stop("ERROR: 'N_control' needs to be a number > 0")
  }
  if(length(proxyMAFs) != 1) {
    if(length(proxyMAFs) != length(OR)) {
      stop("ERROR: 'proxyMAFs' needs to be the same length as the vector of ORs. If you don't have proxy true MAFs then set proxyMAFs = NA")
    }
  }
  if(any(SE < 0)) {
    stop("ERROR: 'SE' vector cannot contain negative values")
  }

  # this function uses w, x, y, z as the derivation in the ReACt paper does - see their supplement
  w = SE^2
  x = 2*N_case
  y = 2*N_control
  z = OR

  AC_control = rep(0, length(OR))
  AC_case = rep(0, length(OR))

  for(i in 1:length(OR)) {
    # need to make sure the discriminate will be positive or it won't solve
  #inflate w(se) by 1.001, for at max 49 times (usually can be done within 5 iterations.
  #maximum inflation is 0.050195, ~5%), if disc still < 0 then give up
    for(j in 0:99) {
      w[i] = w[i] * (1.001^j)
      disc = (((2*z[i]*y*(1-z[i])) - (w[i]*x*y*z[i]))^2) - 4*(w[i]*x*z[i] + (1-z[i])^2)*(y*z[i]*(x + y*z[i]))
      if (!is.na(disc) & disc >= 0) {
        break
      }
    }
    if (is.na(disc) | disc < 0) { # if discriminate is <0 cannot solve and return 0s
      AC_control[i] = 0
      AC_case[i] = 0
    } else { # now actually solve the quadratic for allele counts (AC)
        # solve for the a, b, and c of the quadratic equation
        # this quadratic is solving for the allele count (AC) of the controls
        # overall their derivation relies on AC (rather than AF) and then calculates AF
        a = (w[i]*x*z[i]) + (1-z[i])^2
        b = 2*y*z[i]*(1-z[i]) - w[i]*x*y*z[i]
        c = y*z[i]*(x + y*z[i])

         #find roots of quadratic equation
        AF_control_opts =  genpwr::quad_roots(a, b, c)
        # in order to select which root, we need to use each option (d1 and d2) to calculate a, b, c of the 2x2 table of allele counts
        d1 = AF_control_opts[1]
        c1 = y - d1
        b1 = (x*d1)/(y*z[i] - z[i]*d1 + d1)
        a1 = x - b1

        d2 = AF_control_opts[2]
        c2 = y - d2
        b2 = (x*d2)/(y*z[i] - z[i]*d2 + d2)
        a2 = x - b2

        vec1 = c(a1, b1, c1, d1) # vector of a,b,c,d using root 1
        vec2 = c(a2, b2, c2, d2) # vector of a,b,c,d using root 2

        # if both roots allow for all values to be positive, choose the larger
        if(!any(vec1 < 0) & !(any(vec2 < 0))) {
          if(d1 > d2) { # if d1 is the larger root, then the AC_control = c1 and AC_case = a1
            AC_control[i] = c1
            AC_case[i] = a1
          } else {
            AC_control[i] = c2
            AC_case[i] = a2
          }
        } else if(!any(vec1 < 0)) { # if d1 allows all values of 2x2 to be positive but NOT d2, use c1 and a1
          AC_control[i] = c1
          AC_case[i] = a1
        } else { # if d2 allows all values to be positive but NOT d1, use c2 and a2
          AC_control[i] = c2
          AC_case[i] = a2
        }
      }
    }

  # calculate and return the MAF
  MAF_case = AC_case/(2*N_case)
  MAF_control = AC_control/(2*N_control)
  MAF_pop = (MAF_case*N_case + MAF_control*N_control)/(N_case + N_control)

  # correct the MAFs using the proxy MAFs
  if(length(proxyMAFs) == 1){
    return(data.frame(MAF_case = MAF_case,
                      MAF_control = MAF_control,
                      MAF_pop = MAF_pop))
  } else {
    get_adjusted <- function(bins = "large", estimated, true) {
      if(bins == "large") {
        binDat <- data.frame(bins = c("[0.0, 0.1)", "[0.1, 0.2)", "[0.2, 0.3)",
                                      "[0.3, 0.4)", "[0.4, 0.5]"),
                             min = c(0, 0.1, 0.2, 0.3, 0.4),
                             max = c(0.1, 0.2, 0.3, 0.4, 0.5),
                             x = rep(0, 5),
                             x2 = rep(0,5),
                             intercept = rep(0, 5))
      } else {
        binDat <- data.frame(bins = c("[0.0, 0.05)", "[0.05, 0.1)", "[0.1, 0.15)","[0.15, 0.2)",
                                      "[0.2, 0.25)", "[0.25, 0.3)", "[0.3, 0.35)", "[0.35, 0.4)",
                                      "[0.4, 0.45)", "[0.45, 0.5]"),
                             min = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45),
                             max = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5),
                             x = rep(0, 10),
                             x2 = rep(0, 10),
                             intercept = rep(0, 10))
      }

      # store the models and adjusted MAFs
      dat <- data.frame(estimated = estimated, true = true)
      adjusted <- rep(0, length(estimated))

      for(i in 1:nrow(binDat)) {
        # filter data for model fitting to those within MAF bin
        subdat <- dplyr::filter(dat, true >= binDat[i,]$min &
                                  true < binDat[i,]$max)
        mod <- stats::lm(data = subdat, formula = estimated ~ stats::poly(true, 2, raw = T))

        # get studentized residuals use to filter outliers and refit model
        stud_res <- mod$residuals * stats::sd(mod$residuals)
        tokeep <- which(abs(stud_res) < 3)

        mod2 <- stats::lm(data = subdat[tokeep,], formula = estimated ~ stats::poly(true,2, raw = T))

        # save refit model coefficients
        binDat[i,]$intercept <- mod2$coefficients[1]
        binDat[i,]$x <- mod2$coefficients[2]
        binDat[i,]$x2 <- mod2$coefficients[3]

        # find indices of variants within this MAF bin
        keeps <- which(dat$true >= binDat[i,]$min & dat$true < binDat[i,]$max)
        pred = dat$true*mod2$coefficients[2] + dat$true^2*mod2$coefficients[3] + mod2$coefficients[1]
        bias = dat$true - pred
        # calculate and save adjusted MAF
        adjusted[keeps] <- bias[keeps]
      }
      return(adjusted)
    }
    bias = get_adjusted(estimated = MAF_pop, true = proxyMAFs)

    MAF_case_adj = MAF_case + bias
    MAF_control_adj = MAF_control + bias
    MAF_pop_adj = MAF_pop + bias

    # round any negatives to 0
    MAF_case_adj[MAF_case_adj < 0] = 0
    MAF_control_adj[MAF_control_adj < 0] = 0
    MAF_pop_adj[MAF_pop_adj < 0] = 0

    return(data.frame(MAF_case = MAF_case,
                      MAF_control = MAF_control,
                      MAF_pop = MAF_pop,
                      MAF_case_adj = MAF_case_adj,
                      MAF_control_adj = MAF_control_adj,
                      MAF_pop_adj = MAF_pop_adj))
  }
}
