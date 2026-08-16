context("subset propensity score alignment")

## Regression tests for: data = <data.frame>, subset = <...>, and a
## p.score-producing method.trt. Before the fix, getResponseDataCall assigned
## the subset-length score into the FULL-length data.frame column, which
## errors when the lengths don't divide and silently recycles when they do.

test_that("a p.score-producing method.trt aligns to subset when lengths divide (bug1, silent)", {
  n <- 120L
  set.seed(101)
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z = rbinom(n, 1, 0.5), y = rnorm(n))
  sub <- seq(2L, n, by = 2L)             # 60 rows, divides n
  scores <- seq(0.01, 0.60, by = 0.01)   # subset-length, distinguishable from any recycling

  fit <- bartc(y, z, x1 + x2, data = df, subset = sub, method.trt = scores, method.rsp = "bart",
               verbose = FALSE, n.burn = 3L, n.samples = 5L, n.trees = 5L, n.chains = 1L, n.threads = 1L)

  ## the response design must train on the scores at the subset rows, not a
  ## recycled-then-subset column, and must agree with what the fit reports
  expect_equal(as.numeric(fit$data.rsp@x[, "ps"]), scores)
  expect_equal(as.numeric(fit$data.rsp@x[, "ps"]), fit$p.score)
})

test_that("a p.score-producing method.trt no longer errors when subset length doesn't divide nrow (bug1, loud)", {
  n <- 100L
  set.seed(102)
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z = rbinom(n, 1, 0.5), y = rnorm(n))
  sub <- seq_len(73L) # does not divide n = 100

  set.seed(103)
  fit <- bartc(y, z, x1 + x2, data = df, subset = sub, method.trt = "bart", method.rsp = "bart",
               verbose = FALSE, n.burn = 3L, n.samples = 5L, n.trees = 5L, n.chains = 1L, n.threads = 1L)
  expect_is(fit, "bartcFit")
  expect_equal(length(fit$trt), 73L)
  expect_equal(as.numeric(fit$data.rsp@x[, "ps"]), fit$p.score)

  set.seed(103)
  fit <- bartc(y, z, x1 + x2, data = df, subset = sub, method.trt = "bart", method.rsp = "bcf",
               verbose = FALSE, n.burn = 3L, n.samples = 5L, n.trees = 5L, n.chains = 1L, n.threads = 1L)
  expect_is(fit, "bartcFit")
  expect_equal(length(fit$trt), 73L)
})

test_that("p.score = <vector> gives a message naming p.scoreAsCovariate, not a length error (bug2)", {
  df <- data.frame(x1 = rnorm(10), z = rbinom(10, 1, 0.5), y = rnorm(10))
  expect_error(
    bartc(y, z, x1, data = df, p.score = runif(10), method.trt = "none", verbose = FALSE),
    "p.scoreAsCovariate must be a single TRUE or FALSE; to supply propensity scores directly, use method.trt = <vector>",
    fixed = TRUE
  )
})
