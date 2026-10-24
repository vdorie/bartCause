context("regression")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("naive bart matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "bart", method.trt = "none", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L)
  # snapshot refreshed for dbarts f946d728 (the default rule's one cut per gap
  # is weighted by the gap's width, which moves the prior of the treatment's
  # single cut); noise-dominated 5-tree/5-sample fit (true tau ~ 0.27) - the
  # pin is a draw-stream tripwire, not an estimate
  expect_equal(fitted(fit, "cate"), -0.0347119569743732)
})

test_that("bart on p.score matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "bart", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  # snapshot refreshed for dbarts d1de9ad1 (a probit forest is rescaled once a
  # sweep, which moves the propensity score's draws); at this seed the
  # 5-tree/5-sample fit leaves the estimate at zero to rounding - the pin is a
  # draw-stream tripwire, not an estimate
  expect_equal(fitted(fit, "cate"), -3.67234020970386e-16)
})

test_that("bart w/p.weighting matches old", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L, n.trees = 5L, n.reps = 5L)
  # snapshot refreshed for dbarts d1de9ad1 (a probit forest is rescaled once a
  # sweep, which moves the propensity score's draws); at this seed the
  # 5-tree/5-sample fit leaves the estimate at zero to rounding - the pin is a
  # draw-stream tripwire, not an estimate
  expect_equal(fitted(fit, "pate"), -8.33846449515389e-17)
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
    # tmle's default SuperLearner library fits its dbarts learner through
    # dbarts's BayesTree-style door (bartBT), so this value tracks that
    # door's tree-move mixture as well as bart's; moved since dbarts 093dd035
    # by d1de9ad1 (a probit forest is rescaled once a sweep) and f946d728 (the
    # default rule's one cut per gap is weighted by its width), and recorded
    # after both. Skipped on CRAN. Noise-dominated fit (true tau ~ 0.27). Only
    # the installed tmle branch (2.1.x) was regenerated; older-tmle branches
    # below are historical.
    expect_equal(fitted(fit, "pate"), 0.721932327807755)
  } else if (tmle_version >= "2.0.1") {
    expect_equal(fitted(fit, "pate"), 0.445429512755897)
  } else if (tmle_version >= "1.5.0") {
    expect_equal(fitted(fit, "pate"), 0.30048319956979)
  } else {
    expect_equal(fitted(fit, "pate"), 0.293195268298812)
  }
})

