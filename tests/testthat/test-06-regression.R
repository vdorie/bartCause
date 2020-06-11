context("regression")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("naive bart matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "bart", method.trt = "none", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L)
  expect_equal(fitted(fit, "cate"), 1.70530256582424e-16)
})

test_that("bart on p.score matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "bart", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  expect_equal(fitted(fit, "cate"), 1.82520665248376e-16)
})

test_that("bart w/p.weighting matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  expect_equal(fitted(fit, "pate"), 1.8371826654854e-16)
})

test_that("bart w/TMLE matches old", {
  skip_if_not_installed("tmle")
  
  set.seed(22)
  
  fit <- bartc(y, z, x, testData,
               method.rsp = "tmle", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  
  if (packageVersion("tmle") >= "1.5.0")
    expect_equal(fitted(fit, "pate"), 0.471054207787337)
  else
    expect_equal(fitted(fit, "pate"), 0.460325271945704)
})

