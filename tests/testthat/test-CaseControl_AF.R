data("sampleDat")

expected_afcase_AF <- c(2.041813e-01, 4.468049e-01, 1.647775e-01, 0.000000e+00, 1.077692e-05, 1.381395e-02)
expected_afcontrol_AF <- c(2.041254e-01, 4.471490e-01, 1.641968e-01, 0.000000e+00, 3.776488e-05, 1.402353e-02)

nCase_sample = 16550
nControl_sample = 403923


res_AF <- CaseControl_AF(N_case = nCase_sample,
                         N_control = nControl_sample,
                         OR = sampleDat$OR,
                         AF_population = sampleDat$true_maf_pop)

test_that("Number variants retained CaseControl_AF", {
  expect_equal(nrow(res_AF), 500)
})

test_that("Correct number of columns CaseControl_AF", {
  expect_equal(ncol(res_AF), 2)
})

test_that("Output correct cases CaseControl_AF", {
  expect_equal(round(head(res_AF$AF_case),3), round(expected_afcase_AF, 3))
})

test_that("Output correct controls CaseControl_AF", {
  expect_equal(round(head(res_AF$AF_control), 3), round(expected_afcontrol_AF, 3))
})

test_that("Get no error from proper input" , {
  expect_no_error(CaseControl_AF(N_case = 10, N_control = 10, OR = c(1, 1, 1, 1), AF_population = c(.1, .2, .3, .4)))
})

test_that("Get error from input lengths mismatch" , {
  expect_error(CaseControl_AF(N_case = 10, N_control = 10, OR = c(1, 1, 1, 1), AF_population = c(.1, .2, .3, .4, .5)))
})

test_that("Get error from NA in OR" , {
  expect_error(CaseControl_AF(N_case = 10, N_control = 10, OR = c(1, NA, 1, 1), AF_population = c(.1, .2, .3, .4)))
})

test_that("Get error from NA in AF_population" , {
  expect_error(CaseControl_AF(N_case = 10, N_control = 10, OR = c(1, 1, 1, 1), AF_population = c(.1, NA, .3, .4)))
})

test_that("Get error negative case sample size" , {
  expect_error(CaseControl_AF(N_case = -10, N_control = 10, OR = c(1, 1, 1, 1), AF_population = c(.1, .2, .3, .4)))
})

test_that("Get error negative control sample size" , {
  expect_error(CaseControl_AF(N_case = 10, N_control = -10, OR = c(1, 1, 1, 1), AF_population = c(.1, .2, .3, .4)))
})

test_that("Get error negative AF" , {
  expect_error(CaseControl_AF(N_case = 10, N_control = 10, OR = c(1, 1, 1, 1), AF_population = c(-.1, .2, .3, .4)))
})

test_that("Get error AF > 1" , {
  expect_error(CaseControl_AF(N_case = 10, N_control = 10, OR = c(1, 1, 1, 1), AF_population = c(1.1, .2, .3, .4)))
})
