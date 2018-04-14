context("regression")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("naive bart matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "bart", method.trt = "none", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L)
  expect_equal(fitted(fit, "est"), 0.272157631329396)
})

test_that("bart on p.score matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "bart", method.trt = "bart.xval", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L)
  expect_equal(fitted(fit, "est"), 1.04316036256787)
})

test_that("bart w/pweighting matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "pweight", method.trt = "bart.xval", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L)
  expect_equal(fitted(fit, "est"), 1.11933625142046)
})

test_that("bart w/TMLE matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, testData,
               method.rsp = "tmle", method.trt = "bart.xval", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L)
  expect_equal(fitted(fit, "est"), 0.32052378345552)
})

