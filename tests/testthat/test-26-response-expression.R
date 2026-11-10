context("a response written as an expression")

set.seed(26)
n <- 60L
exprData <- data.frame(x = rnorm(n))
exprData$z <- rbinom(n, 1L, plogis(exprData$x))
exprData$y <- exp(exprData$x / 2 + exprData$z / 2 + rnorm(n, sd = 0.3))
exprData$ly <- log(exprData$y)

## the data is named in the call, as a user's is, and read in the caller's frame
exprFit <- function(response, data) {
  call <- substitute(bartc(response, z, x, data = data, verbose = FALSE, seed = 3, n.samples = 10L, n.burn = 5L,
                           n.trees = 5L, n.chains = 1L, n.threads = 1L),
                     list(response = substitute(response), data = substitute(data)))
  suppressMessages(eval(call, parent.frame()))
}

test_that("a response that is an expression of the columns is fit as that column is", {
  fit <- exprFit(log(y), exprData)
  expect_is(fit, "bartcFit")
  expect_equal(summary(fit)$estimates, summary(exprFit(ly, exprData))$estimates)
})

test_that("an expression response with missing values has its missing rows predicted", {
  missingData <- exprData
  missingData$y[c(3L, 9L)] <- NA
  missingData$ly <- log(missingData$y)
  fit <- exprFit(log(y), missingData)
  expect_equal(which(fit$missingRows), c(3L, 9L))
  expect_equal(summary(fit)$estimates, summary(exprFit(ly, missingData))$estimates)
  expect_false(anyNA(extract(fit, "icate")))
})

test_that("a group.by of the wrong length says rows, as the weights do", {
  expect_error(bartc(y, z, x, data = exprData, group.by = factor(rep(1:2, 20L)), use.ranef = FALSE, verbose = FALSE),
               "'group.by' must have one value for each row \\(60\\), but its length is 40")
})
