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

test_that("varying intercept fit matches manual call", {
  skip_if_not_installed("stan4bart")
  df <- with(testData, data.frame(z = z, V1 = x[,1L], V2 = x[,2L], V3 = x[,3L], g = g))
  
  set.seed(22)
  res <- bartCause:::getBartTreatmentFit(y, z, x, data = testData, chains = 1L, iter = 16L, warmup = 8L,
                                         bart_args = list(n.trees = 7L), group.by = g, use.ranef = TRUE)
  expect_true(inherits(res$fit, "stan4bartFit"))
  expect_equal(res$fit$call$formula, str2lang("z ~ bart(V1 + V2 + V3) + (1 | g)"))
  
  set.seed(22)
  s4bFit <- stan4bart::stan4bart(z ~ bart(V1 + V2 + V3) + (1 | g), df, chains = 1L, iter = 16L, warmup = 8L,
                                 verbose = -1L, bart_args = list(n.trees = 7L))
  expect_equal(res$p.score, apply(dbarts::extract(s4bFit, combine_chains = FALSE), 1L, mean))
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

test_that("glm and bart treatment fits require treatment and confounders", {
  expect_error(bartCause:::getGLMTreatmentFit(y, confounders = x, data = testData), "'treatment' variable must be specified")
  expect_error(bartCause:::getGLMTreatmentFit(y, z, data = testData), "'confounders' variable must be specified")
  expect_error(bartCause:::getBartTreatmentFit(y, confounders = x, data = testData), "'treatment' variable must be specified")
  expect_error(bartCause:::getBartTreatmentFit(y, z, data = testData), "'confounders' variable must be specified")
})

test_that("glm fit works with '.' as confounders and a data.frame", {
  df <- with(testData, data.frame(y = y, z = z, x = x))
  res <- bartCause:::getGLMTreatmentFit(y, z, ., data = df)
  manual <- fitted(stats::glm(z ~ x.1 + x.2 + x.3, df, family = stats::binomial))
  expect_equal(unname(res$p.score), unname(manual))
})

test_that("bart treatment fit composes group.by with parametric and rejects crossvalidate", {
  skip_if_not_installed("stan4bart")
  ## group.by no longer collides with parametric: both land in one formula
  res <- bartCause:::getTreatmentDataCall(stan4bart::stan4bart, z, x, data = testData,
                                          parametric = w, group.by = g, use.ranef = TRUE, use.lmer = FALSE)
  expect_equal(res$call, str2lang("stan4bart::stan4bart(z ~ w + bart(x) + (1 | g), treatment = z, data = testData)"))
  
  expect_error(
    bartCause:::getBartTreatmentFit(y, z, x, data = testData, group.by = g, crossvalidate = TRUE,
                                    n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L),
    "crossvalidation not yet supported"
  )
})

