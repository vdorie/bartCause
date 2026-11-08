context("long confounder formulas and plots with missing responses")

test_that("a confounder formula too long for one line is fit by every route", {
  set.seed(1)
  n <- 80L
  p <- 25L
  d <- as.data.frame(matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("xvariable", seq_len(p)))))
  d$z <- rbinom(n, 1, 0.5)
  d$y <- rnorm(n)
  form <- stats::reformulate(paste0("xvariable", seq_len(p)))
  expect_gt(length(deparse(form[[2L]])), 1L)
  
  ## the formula is spliced into the call as a symbol expression, data by name
  fitForm <- function(confounders, data, ...)
    eval(bquote(bartc(y, z, .(confounders), data = data, n.samples = 5L, n.burn = 5L, n.chains = 2L,
                      n.threads = 1L, verbose = FALSE, seed = 1L, ...)))
  
  expect_error(fitForm(form[[2L]], d), NA)
  expect_error(fitForm(form[[2L]], d, method.trt = "bart"), NA)
  expect_error(fitForm(form[[2L]], d, method.trt = "glm"), NA)
  
  ## the same data with a dot, which the grepl in the treatment fit looks for
  expect_error(fitForm(quote(.), d), NA)
})

test_that("plot_support on a fit with missing responses plots every row", {
  set.seed(1)
  n <- 120L
  d <- data.frame(y = rnorm(n), z = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
  missing <- c(3L, 9L, 40L)
  d$y[missing] <- NA
  fit <- bartc(y, z, x1 + x2, data = d, commonSup.rule = "sd", n.samples = 10L, n.burn = 10L,
               n.chains = 2L, n.threads = 1L, verbose = FALSE, seed = 1L)
  expect_equal(nrow(fit$data.rsp@x), n - length(missing))
  
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(plot_support(fit, xvar = "x1", yvar = "y.cf"), NA)
  expect_error(plot_support(fit, xvar = "pca.1", yvar = "pca.2"), NA)
  expect_error(plot_support(fit, xvar = "x1", yvar = "y", sample = "all"), NA)
  
  ## the design plotted is the data's, all rows, treatment not flipped
  plotted <- NULL
  trace <- function(x, y, ...) { plotted <<- list(x = x, y = y); invisible() }
  testthat::local_mocked_bindings(plot = trace, .package = "bartCause")
  plot_support(fit, xvar = "x1", yvar = "x2", sample = "all", legend.x = NULL)
  expect_equal(plotted$x, d$x1)
  expect_equal(plotted$y, d$x2)
  
  ## the treatment is the data's, and the response has NA at the missing rows
  plot_support(fit, xvar = "x1", yvar = "z", sample = "all", legend.x = NULL)
  expect_equal(plotted$y, d$z)
  plot_support(fit, xvar = "x1", yvar = "y", sample = "all", legend.x = NULL)
  expect_equal(plotted$y, d$y)
  expect_true(all(is.na(plotted$y[missing])))
})
