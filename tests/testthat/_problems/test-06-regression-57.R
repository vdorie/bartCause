# Extracted from test-06-regression.R:57

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "bartCause", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("regression")
source(system.file("common", "linearData.R", package = "bartCause"))

# test -------------------------------------------------------------------------
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
