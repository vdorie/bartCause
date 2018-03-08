context("bartc treatment fits")

source(system.file("common", "linearData.R", package = "BartCause"))

test_that("glm fit matches manual call", {
  res <- BartCause:::getGLMTreatmentFit(z, x, testData)
  expect_equal(res$p.score, fitted(stats::glm(z ~ x, stats::binomial, testData)))
})

test_that("glm fit passes arguments to glm", {
  res <- BartCause:::getGLMTreatmentFit(z, x, testData, start = c(0, 0, 0, 0))
  expect_equal(res$p.score, fitted(stats::glm(z ~ x, stats::binomial, testData, start = c(0, 0, 0, 0))))
})

test_that("bart fit matches manual call", {
  set.seed(22)
  res <- BartCause:::getBartTreatmentFit(z, x, testData, n.chains = 1L, n.threads = 1L, n.burn = 50, n.samples = 150, n.trees = 75L)
  set.seed(22)
  expect_equal(res$p.score, apply(pnorm(dbarts::bart2(z ~ x, testData, n.chains = 1L, n.threads = 1L, n.burn = 50, n.samples = 150, n.trees = 75L, verbose = FALSE)$yhat.train), 2L, mean))
})

test_that("xbart fit matches manual call", {
  set.seed(22)
  res <- BartCause:::getBartXValTreatmentFit(z, x, testData, n.chains = 1L, n.threads = 1L, n.burn = 50, n.samples = 150, n.trees = 75L, n.reps = 10L)
  set.seed(22)
  k <- c(0.5, 1, 2, 4, 8)
  xVal <- dbarts::xbart(z ~ x, testData, k = k, n.threads = 1L, n.burn = 50, n.samples = 150, n.trees = 75L, n.reps = 10L, verbose = FALSE)
  k <- k[which.min(apply(xVal, 2L, mean))]
  
  expect_equal(res$p.score, apply(pnorm(dbarts::bart2(z ~ x, testData, k = k, n.chains = 1L, n.threads = 1L, n.burn = 50, n.samples = 150, n.trees = 75L, verbose = FALSE)$yhat.train), 2L, mean))
  rm(k, xVal)
})

