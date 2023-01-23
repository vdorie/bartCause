context("bartc treatment fits")

source(system.file("common", "groupedData.R", package = "bartCause"))

test_that("glm fit matches manual call", {
  res <- bartCause:::getGLMTreatmentFit(y, z, x, data = testData)
  expect_equal(res$p.score, fitted(stats::glm(z ~ x, testData, family = stats::binomial)))
})

test_that("glm with fixef fit matches manual call", {
  res <- bartCause:::getGLMTreatmentFit(y, z, x, data = testData, group.by = g, use.ranef = FALSE)
  expect_equal(res$p.score, fitted(stats::glm(z ~ x + g, testData, family = stats::binomial)))
})

test_that("glmer fit matches manual call", {
  skip_if_not_installed("lme4")
  res <- bartCause:::getGLMTreatmentFit(y, z, x, data = testData, group.by = g, use.ranef = TRUE)
  expect_equal(res$p.score, fitted(lme4::glmer(z ~ x + (1 | g), testData, family = stats::binomial)))
})

test_that("glm fit passes arguments to glm", {
  res <- bartCause:::getGLMTreatmentFit(y, z, x, data = testData, start = c(0, 0, 0, 0))
  expect_equal(res$p.score, fitted(stats::glm(z ~ x, testData, family = stats::binomial, start = c(0, 0, 0, 0))))
})

test_that("bart fit matches manual call", {
  set.seed(22)
  res <- bartCause:::getBartTreatmentFit(y, z, x, data = testData, n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L)
  set.seed(22)
  expect_equal(res$p.score, apply(pnorm(dbarts::bart2(z ~ x, testData, n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L, verbose = FALSE)$yhat.train), 2L, mean))
})

test_that("bart fit with fixef matches manual call", {
  set.seed(22)
  res <- bartCause:::getBartTreatmentFit(y, z, x, data = testData, n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L, group.by = g, use.ranef = FALSE)
  set.seed(22)
  expect_equal(res$p.score, apply(pnorm(dbarts::bart2(z ~ x + g, testData, n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L, verbose = FALSE)$yhat.train), 2L, mean))
})

test_that("rbart_vi fit matches manual call", {
  set.seed(22)
  res <- bartCause:::getBartTreatmentFit(y, z, x, data = testData, n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L, group.by = g, use.ranef = TRUE)
  set.seed(22)
  rbartFit <- dbarts::rbart_vi(z ~ x, testData, n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L, verbose = FALSE, group.by = g)
  rbartPred <- unname(apply(pnorm(rbartFit$yhat.train + rbartFit$ranef[,as.factor(testData$g)]), 2L, mean))
  expect_equal(res$p.score, rbartPred)
})

test_that("bart fit adds extra defaults", {
  set.seed(22)
  res <- bartCause:::getBartTreatmentFit(y, z, x, data = testData, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L, keepTrees = TRUE, combineChains = FALSE)
  expect_equal(dim(res$samples), c(10L, 13L, nrow(testData$x)))
  expect_true(!is.null(res$fit$fit))
  expect_true(res$fit$fit$control@keepTrees)
})

# commenting this out until more control over how long the crossvalidation runs is baked in
if (FALSE) test_that("xbart fit matches manual call", {
  set.seed(22)
  res <- bartCause:::getBartTreatmentFit(y, z, x, data = testData, n.chains = 1L, n.threads = 1L, n.burn = 25, n.samples = 75, n.trees = 25L, crossvalidate = TRUE)
  set.seed(22)
  k <- c(0.5, 1, 2, 4, 8)
  xVal <- dbarts::xbart(z ~ x, testData, k = k, n.threads = 1L, n.burn = 25, n.samples = 75, n.trees = 25L, n.reps = 10L, verbose = FALSE)
  k <- k[which.min(apply(xVal, 2L, mean))]
  
  expect_equal(res$p.score, apply(pnorm(dbarts::bart2(z ~ x, testData, k = k, n.chains = 1L, n.threads = 1L, n.burn = 25, n.samples = 75, n.trees = 25L, verbose = FALSE)$yhat.train), 2L, mean))
  rm(k, xVal)
})

test_that("glm fit fails for non-binary treatment with literals", {
  testData$z <- testData$z + 1
  expect_error(bartCause:::getGLMTreatmentFit(y, z, x, data = testData))
})

test_that("glm fit fails for non-binary treatment with expressions", {
  expect_error(bartCause:::getGLMTreatmentFit(y, z + 1, x, data = testData))
})

test_that("bart fit fails for non-binary treatment with literals", {
  testData$z <- testData$z + 1
  expect_error(bartCause:::getBartTreatmentFit(y, z, x, data = testData, n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L))
})

test_that("glm fit fails for non-binary treatment with data.frame", {
  data <- data.frame(z = testData$z + 1, y = testData$y, x = testData$x)
  z <- testData$z + 1
  x <- testData$x
  y <- testData$y

  expect_error(bartCause:::getGLMTreatmentFit(response = y, treatment = z, confounders = x.1 + x.2 + x.3, data = data))
})

test_that("glm fit fails for non-binary treatment with data.frame", {
  data <- data.frame(z = testData$z + 1, y = testData$y, x = testData$x)
  z <- testData$z + 1
  x <- testData$x
  y <- testData$y

  expect_error(bartCause:::getBartTreatmentFit(y, z, x, data = data, n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L))
})

