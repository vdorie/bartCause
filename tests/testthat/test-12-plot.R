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

test_that("plot_sigma overlays one trace per chain (dbarts 1.0-0 flattens sigma/first.sigma)", {
  ## Regression test: dbarts 1.0-0 stores fit.rsp$sigma/$first.sigma as flat,
  ## sample-major vectors when n.chains > 1 (no per-chain matrix), so the
  ## un-reshaped is.null(dim(.)) check used to always take the single-trace
  ## branch: the burn-in marker landed at n.chains*n.burn and all chains were
  ## concatenated into one line instead of overlaid.
  n.chains  <- fit$n.chains
  n.burn    <- length(fit$fit.rsp$first.sigma) / n.chains
  n.samples <- length(fit$fit.rsp$sigma) / n.chains
  expect_equal(n.chains, 2L)

  captured <- list(lines = list())
  local_mocked_bindings(
    plot   = function(...) { captured$plot <<- list(...) },
    abline = function(...) { captured$abline <<- list(...) },
    lines  = function(...) { captured$lines[[length(captured$lines) + 1L]] <<- list(...) },
    .package = "bartCause"
  )
  plot_sigma(fit)

  expect_equal(captured$plot$xlim, c(1L, n.samples))
  expect_equal(captured$abline$v, n.burn)
  expect_equal(length(captured$lines), n.chains)

  sigma.mat <- matrix(fit$fit.rsp$sigma, nrow = n.chains)
  warmup.mat <- matrix(fit$fit.rsp$first.sigma, nrow = n.chains)
  full <- cbind(warmup.mat, sigma.mat)
  for (i in seq_len(n.chains))
    expect_equal(captured$lines[[i]][[2L]], full[i, ])
})

test_that("plot_est traces the conditional average treatment effect per chain", {
  captured <- list(lines = list())
  local_mocked_bindings(
    plot  = function(...) { captured$plot <<- list(...) },
    lines = function(...) { captured$lines[[length(captured$lines) + 1L]] <<- list(...) },
    .package = "bartCause"
  )
  plot_est(fit)

  cate <- extract(fit, "cate", combineChains = FALSE)
  expect_equal(captured$plot$xlim, c(1L, ncol(cate)))
  expect_equal(length(captured$lines), fit$n.chains)
  for (i in seq_len(fit$n.chains))
    expect_equal(captured$lines[[i]][[2L]], cate[i, ])
})

test_that("plot_indiv histograms the requested fitted quantity", {
  captured <- list()
  local_mocked_bindings(
    hist = function(x, ...) { captured$x <<- x },
    .package = "bartCause"
  )
  plot_indiv(fit, type = "mu.obs")
  expect_equal(captured$x, fitted(fit, type = "mu.obs"))

  plot_indiv(fit, type = "icate")
  expect_equal(captured$x, fitted(fit, type = "icate"))

  expect_error(plot_indiv(fit, type = "not-a-type"), "type must be in")
})

test_that("plot_support validates its arguments", {
  noSupportFit <- bartc(y, z, x, n.samples = 5L, n.burn = 3L, n.chains = 1L,
                        n.threads = 1L, verbose = FALSE, commonSup.rule = "none")
  expect_error(plot_support(noSupportFit), "requires support rule other than 'none'")
  expect_error(plot_support(fit, xvar = "not-a-real-variable"), "unrecognized variable")
  expect_error(plot_support(fit, sample = "not-a-sample-arg"), "sample must be in")
})
