context("regression")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("naive bart matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "bart", method.trt = "none", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L)
  expect_equal(fitted(fit, "est"), -0.226641057943278)
})

test_that("bart on p.score matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "bart", method.trt = "bart.xval", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  expect_equal(fitted(fit, "est"), 0.216224889131285)
})

test_that("bart w/p.weighting matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "p.weight", method.trt = "bart.xval", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  expect_equal(fitted(fit, "est"), 0.195646445418943)
})

test_that("bart w/TMLE matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "tmle", method.trt = "bart.xval", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  expect_equal(fitted(fit, "est"), 0.21118043658438)
})

