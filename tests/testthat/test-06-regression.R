context("regression")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("naive bart matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "bart", method.trt = "none", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L)
  # snapshot refreshed for dbarts fbff1989 (default proposal mixture moved to
  # birth_death 0.6, swap 0, change 0.4); at this seed the 5-tree/5-sample fit
  # draws no split on the treatment, so the estimate is zero to rounding -
  # the pin is a draw-stream tripwire, not an estimate
  expect_equal(fitted(fit, "cate"), 1.71473946153355e-16)
})

test_that("bart on p.score matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "bart", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  # snapshot refreshed for dbarts fbff1989 (default proposal mixture moved);
  # noise-dominated 5-tree/5-sample fit (true tau ~ 0.27)
  expect_equal(fitted(fit, "cate"), 0.0385169225627653)
})

test_that("bart w/p.weighting matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  # snapshot refreshed for dbarts fbff1989 (default proposal mixture moved);
  # noise-dominated fit (true tau ~ 0.27)
  expect_equal(fitted(fit, "pate"), 0.0605258119600773)
})

test_that("bart w/TMLE matches old", {
  # Because this is not documented, to enable this test execute from R
  #   Sys.setenv(NOT_CRAN = "true")
  # or from shell
  #   export NOT_CRAN=true
  skip_on_cran()
  skip_if_not_installed("tmle")
  
  set.seed(22)
  
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "tmle", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)

  tmle_version <- packageVersion("tmle")
  if (tmle_version >= "2.1") {
    # snapshot refreshed for dbarts 8e1e674c (prior-tree-init draw shift);
    # finite estimate for this noise-dominated fit (true tau ~ 0.27). Only the
    # installed tmle branch (2.1.x) was regenerated; older-tmle branches below
    # are historical.
    expect_equal(fitted(fit, "pate"), -0.12495456896179358)
  } else if (tmle_version >= "2.0.1") {
    expect_equal(fitted(fit, "pate"), 0.445429512755897)
  } else if (tmle_version >= "1.5.0") {
    expect_equal(fitted(fit, "pate"), 0.30048319956979)
  } else {
    expect_equal(fitted(fit, "pate"), 0.293195268298812)
  }
})

