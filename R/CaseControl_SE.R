#' @title CaseControl_SE
#' @description This is a function to derive the case, control, and total MAFs from GWAS summary statistics when
#' the user has access to the sample sizes, and the OR (or beta), and SE for the log(OR) for each variant.
#' If user has total AF instead of SE use [CCAFE::CaseControl_AF()]
#' This code uses the GroupFreq function adapted from C from <https://github.com/Paschou-Lab/ReAct/blob/main/GrpPRS_src/CountConstruct.c>
#'
#' @param data dataframe where each row is a variant and columns contain the OR, SE, chromosome and positions
#' @param N_case an integer of the number of Case individuals
#' @param N_control an integer of the number of Control individuals
#' @param OR_colname a string containing the exact column name in 'data' with the OR
#' @param SE_colname a string containing the exact column name in 'data' with the SE
#' @param chromosome_colname a string containing the exact column name in 'data' with the chromosomes, default "chr"
#' @param position_colname a string containing the exact column name in 'data' with the position, default "pos"
#' @param sex_chromosomes boolean, TRUE if variants from sex chromosomes are included in the dataset. Sex chromosomes can be numeric (23, 24) or character (X, Y). If numeric, assumes X=23 and Y=24.
#' @param N_XX_case the number of XX chromosome case individuals (REQUIRED if sex_chromosomes == TRUE)
#' @param N_XX_control the number of XX chromosome control individuals (REQUIRED if sex_chromosomes == TRUE)
#' @param N_XY_case the number of XY chromosome case individuals (REQUIRED if sex_chromosomes == TRUE)
#' @param N_XY_control the number of XY chromosome control individuals (REQUIRED if sex_chromosomes == TRUE)
#' @param do_correction boolean, TRUE if data is provided to perform correction
#' @param correction_data a dataframe with the following exact columns: CHR, POS, proxy_MAF with data that is harmonized between the proxy true datasets and the observed dataset
#' @param remove_sex_chromosomes boolean, TRUE if should keep autosomes only. This is needed when the number of biological sex males/females per case and control group is not known.
#'
#' @return returns data as a dataframe with three additional columns: MAF_case, MAF_control, MAF_total for the estimated MAFs for each variant. If do_correction = TRUE, then will output 3 additional columns (MAF_case_adj, MAF_control_adj, MAF_total_adj) with the adjusted estimates.
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
#' data("sampleDat")
#' sampleDat <- as.data.frame(sampleDat)
#'
#' nCase_sample = 16550
#' nControl_sample = 403923
#'
#' # get the estimated case and control MAFs
#' se_method_results <- CaseControl_SE(data = sampleDat,
#'                                     N_case = nCase_sample,
#'                                     N_control = nControl_sample,
#'                                     OR_colname = "OR",
#'                                     SE_colname = "SE",
#'                                     chromosome_colname = "CHR",
#'                                     position_colname = "POS")
#'
#' head(se_method_results)
#'
#' @import dplyr
#' @export
CaseControl_SE <- function(data, N_case = 0, N_control = 0, OR_colname = "OR", SE_colname = "SE",
                           chromosome_colname = "chr", sex_chromosomes = FALSE,
                           position_colname = "pos",
                           N_XX_case = NA,
                           N_XX_control = NA,
                           N_XY_case = NA,
                           N_XY_control = NA,
                           do_correction = FALSE,
                           correction_data = NA,
                           remove_sex_chromosomes = TRUE) {
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

  # check that OR and SE columns exist
  if(!OR_colname %in% colnames(data)) {
    stop("ERROR: 'OR_colname' does not exist in 'data'")
  }

  if(!SE_colname %in% colnames(data)) {
    stop("ERROR: 'SE_colname' does not exist in 'data'")
  }

  if(!chromosome_colname %in% colnames(data)) {
    stop("ERROR: 'chromosome_colname does not exist in 'data'")
  }

  if(!position_colname %in% colnames(data)) {
    stop("ERROR: 'position_colname does not exist in 'data'")
  }

  OR <- data[,OR_colname]
  SE <- data[,SE_colname]

  # check for valid input for OR and SE
  if(typeof(OR) != "double") {
    stop("ERROR: 'OR' values must all be numbers (hint: check for NAs)")
  }

  if(typeof(SE) != "double") {
    stop("ERROR: 'SE' values must all be numbers (hint: check for NAs)")
  }

  if(any(SE < 0)) {
    stop("ERROR: 'SE' cannot contain negative values")
  }

  sexchr_type <- NA

  # check sex chromosome info
  if(typeof(data[,c(chromosome_colname)]) == "character") { # if column is character
    # check if can be coerced to all numeric
    if(any(suppressWarnings(is.na(as.numeric(data[,c(chromosome_colname)]))))) {
      sexchr_type <- "chr"
      # there were characters in chromosome column so make sure sex chromosomes is true
      if(!sex_chromosomes & !remove_sex_chromosomes) {
        stop("ERROR: Non-numeric chromosomes (sex chromosomes) present but 'sex_chromosomes' = FALSE")
      }
    } else { # if there aren't any non-numeric
      # make sure there are none > 22 which could be numeric representations of sex chromosomes
      if(any(as.numeric(data[,c(chromosome_colname)]) > 22)) {
        if(!sex_chromosomes & !remove_sex_chromosomes) {
          sexchr_type <- "num"
          stop("ERROR: Chromosomes > 22 (sex chromosomes) present but 'sex_chromosomes' = FALSE")
        }
      }
    }
  } else {
    # if chromosome column is all numeric check for those > 22 which are likely sex chromosomes
    if(any(data[,c(chromosome_colname)]) > 22) {
      sexchr_type <- "num"
      if(!sex_chromosomes) {
        stop("ERROR: Chromosomes > 22 (sex chromosomes) present but 'sex_chromosomes' = FALSE")
      }
    }
  }

  # make sure N_XX and N_XY are included if sex_chromosomes is true
  if(sex_chromosomes & !remove_sex_chromosomes) {
    if(is.na(N_XX_case) | is.na(N_XX_control) | is.na(N_XY_case) | is.na(N_XY_control)) {
      stop("ERROR: 'sex_chromosomes' = TRUE but N_XX_case, N_XX_control, N_XY_case and N_XY_control are not included")
    }
  }
  # input checking for correction
  if(do_correction) {
    if(!is.data.frame(correction_data)) {
      stop("ERROR: 'correction_data' must be a dataframe")
    }

    if(!"CHR" %in% colnames(correction_data)) {
      stop("ERROR: column name 'CHR' must be in 'correction_data'")
    }

    if(!"POS" %in% colnames(correction_data)) {
      stop("ERROR: column name 'POS' must be in 'correction_data'")
    }

    if(!"proxy_MAF" %in% colnames(correction_data)) {
      stop("ERROR: column name 'proxy_MAF' must be in 'correction_data'")
    }
  }

  # solve for real roots of a quadratic using the quadratic formula
  quad_roots<-function(a,b,c){
    c(((-b-sqrt(b^2-4*a*c))/(2*a)),((-b+sqrt(b^2-4*a*c))/(2*a)))
  }

  # this function contains code to calculate the MAF using code adapted from ReACt. For explanations of w, x, y, z, see original ReACt paper (Supplement 5.2)
  solve_maf <- function(w, x, y, z, N_case, N_control) {
    AC_control <- rep(0, length(OR))
    AC_case <- rep(0, length(OR))

    for(i in seq_len(length(OR))) {
      # need to make sure the discriminate will be positive or it won't solve
      # inflate w(se) by 1.001, for at max 49 times (usually can be done within 5 iterations.
      #maximum inflation is 0.050195, ~5%), if disc still < 0 then give up
      for(j in seq(from = 0, to = 99, by = 1)) {
        w[i] <- w[i] * (1.001^j)
        disc <- (((2*z[i]*y*(1-z[i])) - (w[i]*x*y*z[i]))^2) - 4*(w[i]*x*z[i] + (1-z[i])^2)*(y*z[i]*(x + y*z[i]))
        if (!is.na(disc) & disc >= 0) {
          break
        }
      }
      if (is.na(disc) | disc < 0) { # if discriminate is <0 cannot solve and return 0s
        AC_control[i] <- 0
        AC_case[i] <- 0
      } else { # now actually solve the quadratic for allele counts (AC)
        # solve for the a, b, and c of the quadratic equation
        # this quadratic is solving for the allele count (AC) of the controls
        # overall their derivation relies on AC (rather than AF) and then calculates AF
        a <- (w[i]*x*z[i]) + (1-z[i])^2
        b <- 2*y*z[i]*(1-z[i]) - w[i]*x*y*z[i]
        c <- y*z[i]*(x + y*z[i])

        #find roots of quadratic equation
        AF_control_opts <- suppressWarnings(quad_roots(a, b, c))
        # in order to select which root, we need to use each option (d1 and d2) to calculate a, b, c of the 2x2 table of allele counts
        d1 <- AF_control_opts[1]
        c1 <- y - d1
        b1 <- (x*d1)/(y*z[i] - z[i]*d1 + d1)
        a1 <- x - b1

        d2 <- AF_control_opts[2]
        c2 <- y - d2
        b2 <- (x*d2)/(y*z[i] - z[i]*d2 + d2)
        a2 <- x - b2

        vec1 <- c(a1, b1, c1, d1) # vector of a,b,c,d using root 1
        vec2 <- c(a2, b2, c2, d2) # vector of a,b,c,d using root 2

        # if both roots allow for all values to be positive, choose the larger
        if(!any(vec1 < 0) & !(any(vec2 < 0))) {
          if(d1 > d2) { # if d1 is the larger root, then the AC_control = c1 and AC_case = a1
            AC_control[i] <- c1
            AC_case[i] <- a1
          } else {
            AC_control[i] <- c2
            AC_case[i] <- a2
          }
        } else if(!any(vec1 < 0)) { # if d1 allows all values of 2x2 to be positive but NOT d2, use c1 and a1
          AC_control[i] <- c1
          AC_case[i] <- a1
        } else { # if d2 allows all values to be positive but NOT d1, use c2 and a2
          AC_control[i] <- c2
          AC_case[i] <- a2
        }
      }
    }

    # calculate and return the MAF
    MAF_case <- AC_case/(x)
    MAF_control <- AC_control/(y)
    MAF_pop <- (AC_case + AC_control)/(x + y)

    return(data.frame(MAF_case, MAF_control, MAF_pop))
  }

  if(!sex_chromosomes) {
    # this function uses w, x, y, z as the derivation in the ReACt paper does - see their supplement
    w <- SE^2
    x <- 2*N_case
    y <- 2*N_control
    z <- OR

    res <- solve_maf(w, x, y, z, N_case, N_control)
    data$MAF_case <- res$MAF_case
    data$MAF_control <- res$MAF_control
    data$MAF_total <- res$MAF_pop

  } else {
    if(remove_sex_chromosomes) {
      keep_chrs <- seq(1, 22, 1)
      data$temp_chr <- data[,c(chromosome_colname)]
      data <- data[data$temp_chr %in% keep_chrs, ]

      SE <- data[,c(SE_colname)]
      OR <- data[,c(OR_colname)]

      w <- SE^2
      x <- 2*N_case
      y <- 2*N_control
      z <- OR

      res <- solve_maf(w, x, y, z, N_case, N_control)
      data$MAF_case <- res$MAF_case
      data$MAF_control <- res$MAF_control
      data$MAF_total <- res$MAF_pop
    } else {
      data$temp_chr <- data[,c(chromosome_colname)]

      N_XX_total <- N_XX_case + N_XX_control
      N_XY_total <- N_XY_case + N_XY_control
      N_total <- N_case + N_control

      if(sexchr_type == "num") {
        data_nosex <- data[data$temp_chr <= 22, ]
        data_x <- data[data$temp_chr == 23, ]
        data_y <- data[data$temp_chr == 24, ]
      } else {
        data_nosex <- data[!data$temp_chr %in% c("X", "x", "Y", "y"), ]
        data_x <- data[data$temp_chr %in% c("X", "x"), ]
        data_y <- data[data$temp_chr %in% c("Y", "y"), ]

      }
      # do autosome data
      SE <- data_nosex[,c(SE_colname)]
      OR <- data_nosex[,c(OR_colname)]

      w <- SE^2
      x <- 2*N_case
      y <- 2*N_control
      z <- OR

      res <- solve_maf(w, x, y, z, N_case, N_control)
      data_nosex$MAF_case <- res$MAF_case
      data_nosex$MAF_control <- res$MAF_control
      data_nosex$MAF_total <- res$MAF_pop

      # do sex chromosomes
      # start with x
      if(nrow(data_x) > 0) {
        SE <- data_x[,c(SE_colname)]
        OR <- data_x[,c(OR_colname)]

        N_X_control <- (2*N_XX_control + N_XY_control)
        N_X_case <- (2*N_XX_case + N_XY_case)
        N_X_total <- (2*N_XX_total + N_XY_total)

        w <- SE^2
        x <- N_X_case
        y <- N_X_control
        z <- OR
        res <- solve_maf(w, x, y, z, N_X_case, N_X_control)
        data_x$MAF_case <- res$MAF_case
        data_x$MAF_control <- res$MAF_control
        data_x$MAF_total <- res$MAF_pop
      }
      # do y chromosome
      if(nrow(data_y > 0)) {
        SE <- data_y[,c(SE_colname)]
        OR <- data_y[,c(OR_colname)]

        N_Y_control <- (N_XY_control)
        N_Y_case <- (N_XY_case)
        N_Y_total <- (N_XY_total)

        w <- SE^2
        x <- N_Y_case
        y <- N_Y_control
        z <- OR

        res <- solve_maf(w, x, y, z, N_Y_case, N_Y_control)
        data_y$MAF_case <- res$MAF_case
        data_y$MAF_control <- res$MAF_control
        data_y$MAF_total <- res$MAF_pop
      }
      data <- rbind(data_nosex, data_x, data_y)
    }
  }


  # do the correction if do_correction
  if(!do_correction) {
    return(data)
  } else {
    get_model <- function(bins = "large", estimated, true) {
      if(bins == "large") {
        binDat <- data.frame(bins = c("[0.0, 0.1)", "[0.1, 0.2)", "[0.2, 0.3)",
                                      "[0.3, 0.4)", "[0.4, 0.5]"),
                             min = c(0, 0.1, 0.2, 0.3, 0.4),
                             max = c(0.1, 0.2, 0.3, 0.4, 0.5),
                             x = rep(0, 5),
                             x2 = rep(0,5),
                             intercept = rep(0, 5),
                             ymin = rep(0, 5),
                             ymax = rep(0, 5))
      } else {
        binDat <- data.frame(bins = c("[0.0, 0.05)", "[0.05, 0.1)", "[0.1, 0.15)","[0.15, 0.2)",
                                      "[0.2, 0.25)", "[0.25, 0.3)", "[0.3, 0.35)", "[0.35, 0.4)",
                                      "[0.4, 0.45)", "[0.45, 0.5]"),
                             min = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45),
                             max = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5),
                             x = rep(0, 10),
                             x2 = rep(0, 10),
                             intercept = rep(0, 10),
                             ymin = rep(0, 10),
                             ymax = rep(0, 10))
      }

      # store the models and adjusted MAFs
      dat <- data.frame(estimated = estimated, true = true)
      adjusted <- rep(0, length(estimated))

      for(i in seq_len(nrow(binDat))) {
        # filter data for model fitting to those within MAF bin
        subdat <- dplyr::filter(dat, true >= binDat[i,]$min &
                                  true < binDat[i,]$max)
        mod <- stats::lm(data = subdat, formula = estimated ~ stats::poly(true, 2, raw = TRUE))

        # get studentized residuals use to filter outliers and refit model
        stud_res <- mod$residuals * stats::sd(mod$residuals)
        tokeep <- which(abs(stud_res) < 3)

        if(length(tokeep) > 0) {
          mod2 <- stats::lm(data = subdat[tokeep,], formula = estimated ~ stats::poly(true,2, raw = TRUE))

          # save refit model coefficients
          binDat[i,]$intercept <- mod2$coefficients[1]
          binDat[i,]$x <- mod2$coefficients[2]
          binDat[i,]$x2 <- mod2$coefficients[3]
        } else {
          # save refit model coefficients
          binDat[i,]$intercept <- mod$coefficients[1]
          binDat[i,]$x <- mod$coefficients[2]
          binDat[i,]$x2 <- mod$coefficients[3]
        }

        # calculate max y'
        binDat[i,]$ymax <- binDat[i,]$max^2*binDat[i,]$x2 + binDat[i,]$max*binDat[i,]$x + binDat[i,]$intercept

        # assign y' min
        if(i == 1) {
          binDat[i,]$ymin <- 0
        } else {
          binDat[i,]$ymin <- binDat[i-1, ]$ymax
        }
      }
      return(binDat)
    }
    get_bias <- function(binDat, estimated){
      # first create dataframe with estimated, and the appropriate intercept, x, x2 values given the
      # value of estimated
      dat <- data.frame(estimated = estimated)

      breaks <- c(0, binDat$ymax, 1)
      dat$bin <- cut(dat$estimated, breaks = breaks, include.lowest = TRUE,
                     right = FALSE, labels = FALSE)
      dat[dat$bin == 6, ]$bin <- 5

      dat$intercept <- 0
      dat$x2 <- 0
      dat$x <- 0

      for(i in seq(1, 5, 1)) {
        dat[dat$bin == i, ]$intercept <- binDat[i, ]$intercept
        dat[dat$bin == i, ]$x <- binDat[i,]$x
        dat[dat$bin == i, ]$x2 <- binDat[i, ]$x2
      }
      # for each variant calculate the value of x for the regression (gnomAD value the corresponds to
      # that estimated value)
      get_x <- function(dat) {
        # we have only AF' which, in terms of fitting the regression = y'
        # need to calculate x given y'
        # y' = beta_0 + beta_1*x + beta_2*x^2 (solve the quadratic)

        # data columns will be in the following order:
        # estimated, bin, intercept, x2, x
        x_opts <- suppressWarnings(quad_roots(dat[4], dat[5], (dat[3] - dat[1])))
        if(all(is.na(x_opts))) {
          x <- 0
        } else {
          if(!any(x_opts >= 0 & x_opts <= 0.5)) {
            x <- 0
          } else {
            x <- x_opts[which(x_opts >= 0 & x_opts <= .5)]
            if(length(x) > 1) { # if both are between 0 and .5 then select the one closest to the estimated value
              x <- x[which.min(abs(x - dat[1]))]
            }
          }
        }
        return(as.numeric(x))
      }

      dat$gnomad <- apply(dat, 1, function(x) unlist(get_x(x)))

      # calculate the bias
      bias <- dat$gnomad - dat$estimated

      # anything from 0 to 0.05 don't correct
      #bias[which(estimated < 0.05)] <- 0
      # if 0 output as gnomAD estimate don't correct
      bias[which(dat$gnomad == 0)] <- 0

      return(bias)
    }

    merged <- suppressWarnings(left_join(correction_data, data, by = c("CHR" = chromosome_colname, "POS" = position_colname)))

    mod <- get_model(estimated = merged$MAF_total, true = merged$proxy_MAF)
    bias <- get_bias(mod, estimated = data$MAF_total)

    data$MAF_case_adj <- data$MAF_case + bias
    data$MAF_control_adj <- data$MAF_control + bias
    data$MAF_total_adj <- data$MAF_total + bias
    return(data)
  }
}
