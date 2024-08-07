data("sampleDat")

sampleDat <- as.data.frame(sampleDat)

expected_mafcase_SE <- c(1.764905e-01, 3.263835e-01, 1.392353e-01, 2.513875e-06, 1.946826e-05, 1.250050e-02)
expected_mafcontrol_SE <- c(1.764405e-01, 3.266896e-01, 1.387297e-01, 2.794421e-06, 6.821993e-05, 1.269041e-02)
expected_mafpop_SE <- c(1.764425e-01, 3.266776e-01, 1.387496e-01, 2.783379e-06, 6.630105e-05, 1.268294e-02)

expected_mafcase_adj <- c(0.20818449, 0.45723980, 0.16577908, 0.00000000, 0.00000000, 0.00963508)
expected_mafcontrol_adj <- c(0.208134479, 0.457545911, 0.165273492, 0.000000000, 0.000000000, 0.009824987)
expected_mafpop_adj <- c(0.208136448, 0.457533862, 0.165293392, 0.000000000, 0.000000000, 0.009817512)

nCase_sample = 16550
nControl_sample = 403923

res_SE <- CaseControl_SE(data = sampleDat,
                         N_case = nCase_sample,
                         N_control = nControl_sample,
                         OR_colname = "OR",
                         SE_colname = "SE",
                         chromosome_colname = "CHR",
                         position_colname = "POS",
                         do_correction = F,
                         sex_chromosomes = F)

corr_dat <- data.frame(CHR = sampleDat$CHR,
                       POS = sampleDat$POS,
                       proxy_MAF = sampleDat$gnomad_maf)

res_SE_adj <- CaseControl_SE(data = sampleDat,
                             N_case = nCase_sample,
                             N_control = nControl_sample,
                             OR_colname = "OR",
                             SE_colname = "SE",
                             chromosome_colname = "CHR",
                             position_colname = "POS",
                             do_correction = T,
                             sex_chromosomes = F,
                             correction_data = corr_dat)

test_that("Number variants retained CaseControl_SE", {
  expect_equal(nrow(res_SE), 500)
})

test_that("Number variants retained CaseControl_SE with correction", {
  expect_equal(nrow(res_SE_adj), 500)
})

test_that("Correct number of columns CaseControl_SE", {
  expect_equal(ncol(res_SE), (ncol(sampleDat) + 3))
})

test_that("Correct number of columns CaseControl_SE with correction", {
  expect_equal(ncol(res_SE_adj), (ncol(sampleDat) + 6))
})

test_that("Output correct cases CaseControl_SE", {
  expect_equal(round(head(res_SE$MAF_case), 3), round(expected_mafcase_SE, 3))
})

test_that("Output correct controls CaseControl_SE", {
  expect_equal(round(head(res_SE$MAF_control), 3), round(expected_mafcontrol_SE, 3))
})

test_that("Output correct pop CaseControl_SE", {
  expect_equal(round(head(res_SE$MAF_total), 3), round(expected_mafpop_SE, 3))
})

testDat <- data.frame(OR = c(1, 1, 1, 1),
                      SE = c(.1, .2, .3, .4),
                      CHR = c(1, 1, 1, 1),
                      POS = c(1, 2, 3, 4))
test_that("Get no error from proper input" , {
  expect_no_error(CaseControl_SE(data = testDat, N_case = 10, N_control = 10, OR_colname = "OR",
                                 SE_colname = "SE",
                                 chromosome_colname = "CHR",
                                 position_colname = "POS",
                                 do_correction = F,
                                 sex_chromosomes = F))
})

testDat <- data.frame(OR = c(1, NA, 1, 1),
                      SE = c(.1, .2, .3, .4),
                      CHR = c(1, 1, 1, 1),
                      POS = c(1, 2, 3, 4))
test_that("Get error from NA in OR" , {
  expect_error(CaseControl_AF(data = testDat, N_case = 10, N_control = 10, OR_colname = "OR", SE_colname = "SE",
                              chromosome_colname = "CHR",
                              position_colname = "POS",
                              do_correction = F,
                              sex_chromosomes = F))
})

testDat <- data.frame(OR = c(1, 1, 1, 1),
                      SE = c(.1, .2, NA, .4),
                      CHR = c(1, 1, 1, 1),
                      POS = c(1, 2, 3, 4))
test_that("Get error from NA in SE" , {
  expect_error(CaseControl_AF(data = testDat, N_case = 10, N_control = 10, OR_colname = "OR", SE_colname = "SE",
                              chromosome_colname = "CHR",
                              position_colname = "POS",
                              do_correction = F,
                              sex_chromosomes = F))
})

testDat <- data.frame(OR = c(1, 1, 1, 1),
                      SE = c(.1, .2, .3, .4),
                      CHR = c(1, 1, 1, 1),
                      POS = c(1, 2, 3, 4))
test_that("Get error negative case sample size" , {
  expect_error(CaseControl_AF(data = testDat, N_case = -10, N_control = 10, OR_colname = "OR", SE_colname = "SE",
                              chromosome_colname = "CHR",
                              position_colname = "POS",
                              do_correction = F,
                              sex_chromosomes = F))
})

test_that("Get error negative control sample size" , {
  expect_error(CaseControl_AF(data = testDat, N_case = 10, N_control = -10, OR_colname = "OR",
                              SE_colname = "SE",
                              chromosome_colname = "CHR",
                              position_colname = "POS",
                              do_correction = F,
                              sex_chromosomes = F))
})

testDat <- data.frame(OR = c(1, 1, 1, 1),
                      SE = c(-.1, .2, .3, .4),
                      CHR = c(1, 1, 1, 1),
                      POS = c(1, 2, 3, 4))
test_that("Get error negative SE" , {
  expect_error(CaseControl_SE(data = testDat, N_case = 10, N_control = 10, OR_colname = "OR",
                              SE_colname = "SE",
                              chromosome_colname = "CHR",
                              position_colname = "POS",
                              do_correction = F,
                              sex_chromosomes = F))
})

testDat <- data.frame(OR = c(1, 1, 1, 1),
                      SE = c(.1, .2, .3, .4),
                      CHR = c(1, 1, 1, 1),
                      POS = c(1, 2, 3, 4))
testDat <- as.list(testDat)
test_that("Get error if input not dataframe", {
  expect_error(CaseControl_SE(data = testDat, N_case = 10, N_control = 10, OR_colname = "OR",
                              SE_colname = "SE",
                              chromosome_colname = "CHR",
                              position_colname = "POS",
                              do_correction = F,
                              sex_chromosomes = F))
})

test_that("Get error if odds ratio column name is wrong", {
  expect_error(CaseControl_SE(data = testDat, N_case = 10, N_control = 10, OR_colname = "Odds.Ratio",
                              SE_colname = "SE",
                              chromosome_colname = "CHR",
                              position_colname = "POS",
                              do_correction = F,
                              sex_chromosomes = F))
})

test_that("Get error if standard error columnn name is wrong", {
  expect_error(CaseControl_SE(data = testDat, N_case = 10, N_control = 10, OR_colname = "OR",
                              SE_colname = "err",
                              chromosome_colname = "CHR",
                              position_colname = "POS",
                              do_correction = F,
                              sex_chromosomes = F))
})


