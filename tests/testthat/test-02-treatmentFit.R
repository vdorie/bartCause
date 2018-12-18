context("bartc treatment fits")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("glm fit matches manual call", {
  res <- bartCause:::getGLMTreatmentFit(z, x, testData)
  expect_equal(res$p.score, fitted(stats::glm(z ~ x, stats::binomial, testData)))
})

test_that("glm fit passes arguments to glm", {
  res <- bartCause:::getGLMTreatmentFit(z, x, testData, start = c(0, 0, 0, 0))
  expect_equal(res$p.score, fitted(stats::glm(z ~ x, stats::binomial, testData, start = c(0, 0, 0, 0))))
})

test_that("bart fit matches manual call", {
  set.seed(22)
  res <- bartCause:::getBartTreatmentFit(z, x, testData, n.chains = 1L, n.threads = 1L, n.burn = 50, n.samples = 150, n.trees = 75L)
  set.seed(22)
  expect_equal(res$p.score, apply(pnorm(dbarts::bart2(z ~ x, testData, n.chains = 1L, n.threads = 1L, n.burn = 50, n.samples = 150, n.trees = 75L, verbose = FALSE)$yhat.train), 2L, mean))
})

test_that("xbart fit matches manual call", {
  set.seed(22)
  res <- bartCause:::getBartXValTreatmentFit(z, x, testData, n.chains = 1L, n.threads = 1L, n.burn = 25, n.samples = 75, n.trees = 25L, n.reps = 10L)
  set.seed(22)
  k <- c(0.5, 1, 2, 4, 8)
  xVal <- dbarts::xbart(z ~ x, testData, k = k, n.threads = 1L, n.burn = 25, n.samples = 75, n.trees = 25L, n.reps = 10L, verbose = FALSE)
  k <- k[which.min(apply(xVal, 2L, mean))]
  
  expect_equal(res$p.score, apply(pnorm(dbarts::bart2(z ~ x, testData, k = k, n.chains = 1L, n.threads = 1L, n.burn = 25, n.samples = 75, n.trees = 25L, verbose = FALSE)$yhat.train), 2L, mean))
  rm(k, xVal)
})

test_that("xbart fit passes on extra args", {
  res <- bartCause:::getBartXValTreatmentFit(z, x, testData, n.chains = 1L, n.threads = 1L, n.burn = 10L,
                                             n.samples = list(15, 8), n.trees = c(25L, 30L), n.reps = 10L)
  expect_equal(ncol(res$samples), 8L)
})
