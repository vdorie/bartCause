context("counterfactual of a row with a missing response")

set.seed(25)
mrN <- 60L
mrData <- data.frame(x = rnorm(mrN), g = factor(rep_len(c("a", "b", "c"), mrN)))
mrData$z <- rbinom(mrN, 1L, plogis(mrData$x))
mrData$y <- mrData$x + mrData$z * (1 + 0.5 * mrData$x) + rnorm(mrN)
mrData$yb <- as.integer(mrData$y > median(mrData$y))
mrData$w <- runif(mrN, 0.5, 1.5)
mrMissing <- c(1L, 5L, 22L, 37L, mrN)

mrFit <- function(response = "y", estimand = "ate", missing = mrMissing, data = mrData, method.trt = "none", ...) {
  data[[response]][missing] <- NA
  extra <- list(...)
  call <- c(list(quote(bartCause::bartc), as.name(response), quote(z), quote(x), data = quote(data),
                 estimand = estimand, method.trt = method.trt, p.scoreAsCovariate = FALSE, verbose = FALSE,
                 seed = 3L, n.samples = 12L, n.burn = 8L, n.trees = 8L, n.chains = 2L, n.threads = 1L,
                 keepTrees = TRUE), extra)
  eval(as.call(call))
}

## the expected outcomes at the rows of 'data', with the treatment set to
## 'trt', predicted by the fit's own trees at the rows' predictors as the user
## wrote them (not through the fit's test matrix)
mrHand <- function(fit, data, trt, byGroup = FALSE, combine = TRUE) {
  newdata <- data.frame(x = data$x, z = trt)
  if (byGroup) {
    newdata$g <- data$g
  }
  ## the glm propensity score the p.weight and tmle fits take as a covariate
  if ("ps" %in% colnames(fit$data.rsp@x)) newdata$ps <- unname(fitted(fit$fit.trt))
  newdata <- newdata[, colnames(fit$data.rsp@x), drop = FALSE]
  predict(fit$fit.rsp, newdata, type = "ev", combineChains = combine)
}

expectOwnCounterfactual <- function(fit, data = mrData, byGroup = FALSE) {
  mu.0 <- extract(fit, "mu.0", sample = "all")
  mu.1 <- extract(fit, "mu.1", sample = "all")
  expect_equal(unname(mu.0), unname(mrHand(fit, data, 0, byGroup)), tolerance = 1e-8)
  expect_equal(unname(mu.1), unname(mrHand(fit, data, 1, byGroup)), tolerance = 1e-8)
  ## a row's individual effect is not zero by construction
  icate <- extract(fit, "icate", sample = "all")
  expect_equal(unname(icate), unname(mu.1 - mu.0), tolerance = 1e-8)
  rows <- which(fit$missingRows)
  expect_true(all(apply(icate[, rows, drop = FALSE], 2L, sd) > 0))
}

## the response method's two arrays, observed treatment and counterfactual, by chain
expectOwnArrays <- function(fit, data = mrData) {
  trt <- fit$trt
  expect_equal(unname(fit$mu.hat.obs), unname(mrHand(fit, data, trt, combine = FALSE)), tolerance = 1e-8)
  expect_equal(unname(fit$mu.hat.cf), unname(mrHand(fit, data, 1 - trt, combine = FALSE)), tolerance = 1e-8)
  rows <- which(fit$missingRows)
  expect_true(all(apply(fit$mu.hat.cf[, , rows, drop = FALSE] != fit$mu.hat.obs[, , rows, drop = FALSE], 3L, any)))
}

test_that("a missing-response row is predicted at both treatments at its own predictors, for each estimand and outcome", {
  for (response in c("y", "yb")) for (estimand in c("ate", "att", "atc")) {
    fit <- mrFit(response, estimand)
    expect_equal(sum(fit$missingRows), length(mrMissing))
    expectOwnCounterfactual(fit)
  }
})

test_that("a missing-response row at the first and last positions and in a subset is predicted at both treatments at its own predictors", {
  fit <- mrFit(missing = c(1L, 2L, 3L))
  expectOwnCounterfactual(fit)
  fit <- mrFit(missing = mrN - 0:2)
  expectOwnCounterfactual(fit)
  sub <- seq_len(mrN)[-c(2L, 9L, 30L)]
  fit <- mrFit(subset = sub)
  expectOwnCounterfactual(fit, mrData[sub, ])
  fit <- mrFit(missing = c(1L, 2L, 3L), subset = sub)
  expectOwnCounterfactual(fit, mrData[sub, ])
})

test_that("a missing-response row is predicted at both treatments at its own predictors with group.by as a fixed effect", {
  fit <- mrFit(group.by = quote(g), use.ranef = FALSE, group.effects = TRUE)
  expect_true("g" %in% colnames(fit$data.rsp@x))
  expectOwnCounterfactual(fit, byGroup = TRUE)
})

test_that("a missing-response row is predicted at both treatments at its own predictors with a common support rule and weights", {
  for (rule in c("sd", "chisq")) {
    fit <- mrFit(estimand = "att", commonSup.rule = rule, commonSup.cut = 0.5)
    expectOwnCounterfactual(fit)
    ## the rule cuts some rows, and the estimate averages the hand effects of the rest
    rows <- fit$commonSup.sub & fit$trt > 0
    expect_true(any(rows) && !all(fit$commonSup.sub))
    hand <- mrHand(fit, mrData, 1) - mrHand(fit, mrData, 0)
    expect_equal(summary(fit, target = "cate")$estimates$estimate, mean(rowMeans(hand[, rows, drop = FALSE])),
                 tolerance = 1e-8)
  }
  fit <- mrFit(weights = quote(w))
  expectOwnCounterfactual(fit)
  fit <- mrFit(weights = quote(w), method.rsp = "p.weight", method.trt = "glm", estimand = "att")
  expectOwnArrays(fit)
})

test_that("a missing-response row is predicted at both treatments at its own predictors under tmle", {
  skip_on_cran()
  skip_if_not_installed("tmle")
  fit <- mrFit(method.rsp = "tmle", method.trt = "glm", maxIter = 3L)
  expectOwnArrays(fit)
})
