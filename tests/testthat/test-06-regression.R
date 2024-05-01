context("regression")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("naive bart matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "bart", method.trt = "none", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L)
  expect_equal(fitted(fit, "cate"), -0.237484791465709)
})

test_that("bart on p.score matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "bart", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  expect_equal(fitted(fit, "cate"), 0.0213594201614189)
})

test_that("bart w/p.weighting matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  expect_equal(fitted(fit, "pate"), 0.0211982499681833)
})

test_that("bart w/TMLE matches old", {
  # Because this is not documented, to enable this test execute from R
  #   Sys.setenv(NOT_CRAN = "true")
  # or from shell
  #   export NOT_CRAN=true
  skip_on_cran()
  skip_if_not_installed("tmle")
  
  set.seed(22)
  
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "tmle", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)

  tmle_version <- packageVersion("tmle")
  if (tmle_version >= "2.0.1") {
    expect_equal(fitted(fit, "pate"), 0.445429512755897)
  } else if (tmle_version >= "1.5.0") {
    expect_equal(fitted(fit, "pate"), 0.30048319956979)
  } else {
    expect_equal(fitted(fit, "pate"), 0.293195268298812)
  }
})

