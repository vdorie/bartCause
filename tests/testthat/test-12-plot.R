context("plot methods")

## Regression coverage for the dbarts 1.0-0 finalization: dbarts stores @x as a
## dbartsMixedMatrix, which the common-support plot must coerce before crossprod /
## %*% / data.frame. This path is exercised only here (not by the other tests),
## and it was caught originally by the package example, not the suite.

set.seed(22)
n <- 100L
x <- matrix(rnorm(3 * n), n, 3)
z <- rbinom(n, 1, pnorm(0.5 * x[, 1]))
y <- as.numeric(2 * z + x[, 1] - 0.5 * x[, 2] + rnorm(n))

fit <- bartc(y, z, x, n.samples = 40L, n.burn = 20L, n.chains = 2L,
             n.threads = 1L, verbose = FALSE, commonSup.rule = "sd")

test_that("plot methods run against a dbartsMixedMatrix design matrix", {
  pf <- tempfile(fileext = ".pdf")
  grDevices::pdf(pf)
  on.exit({ grDevices::dev.off(); unlink(pf) }, add = TRUE)

  expect_error(plot_sigma(fit), NA)
  # PCA path exercises crossprod() and %*% on the coerced design matrix
  expect_error(plot_support(fit, xvar = "pca.1", yvar = "pca.2"), NA)
  # tree path exercises data.frame() on the coerced design matrix
  skip_if_not_installed("rpart")
  expect_error(plot_support(fit, xvar = "tree.1", yvar = "css", legend.x = NULL), NA)
})
