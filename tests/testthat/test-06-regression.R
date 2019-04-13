context("regression")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("naive bart matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "bart", method.trt = "none", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L)
  expect_equal(fitted(fit, "cate"), -0.226641057943278)
})

test_that("bart on p.score matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "bart", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  expect_equal(fitted(fit, "cate"), -0.292368476106707)
})

test_that("bart w/p.weighting matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  expect_equal(fitted(fit, "pate"), -0.245729972409719)
})

test_that("bart w/TMLE matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "tmle", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  expect_equal(fitted(fit, "pate"), 0.775833381080805)
})

