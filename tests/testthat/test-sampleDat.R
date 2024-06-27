data("sampleDat")

test_that("Correct number of variants", {
  expect_equal(dim(sampleDat)[1], 500)
})

test_that("Correct number of columns", {
  expect_equal(dim(sampleDat)[2], 11)
})

