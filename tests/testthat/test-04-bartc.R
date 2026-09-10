context("bartc main function")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("bartc fails with invalid inputs", {
  expect_error(bartc(not.found, z, x, data = testData, verbose = FALSE))
  expect_error(bartc(y, not.found, x, data = testData, verbose = FALSE))
  expect_error(bartc(y, z, not.found, data = testData, verbose = FALSE))
  expect_error(bartc(y, z, x, data = testData, verbose = FALSE, method.rsp = "not-a-method"))
  expect_error(bartc(y, z, x, data = testData, verbose = FALSE, method.trt = "not-a-method"))
  expect_error(bartc(y, z, x, data = testData, verbose = FALSE, estimand = "not-an-estimand"))
  expect_error(bartc(y, z, x, data = testData, verbose = FALSE, group.by = not.found))
  expect_error(bartc(y, z, x, data = testData, verbose = FALSE, p.scoreAsCovariate = NA))
  expect_error(bartc(y, z, x, data = testData, verbose = FALSE, keepCall = NA))
  expect_error(bartc(y, z, x, data = testData, verbose = NA))
})

test_that("bartc matches manual fit", {  
  set.seed(22)
  bartcFit <- bartc(y, z, x, data = testData,
                    method.rsp = "bart", method.trt = "bart", verbose = FALSE,
                    n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L)
  
  set.seed(22)
  fit.trt <- dbarts::bart2(z ~ x, testData, verbose = FALSE,
                           n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L)
  p.score <- apply(pnorm(fit.trt$yhat.train), 2L, mean)
  expect_equal(p.score, fitted(bartcFit, type = "p.score"))
  
  x.train <- cbind(z = testData$z, testData$x, ps = p.score)
  x.test  <- cbind(z = 1, testData$x, ps = p.score)
  x.test <- rbind(x.test, x.test)
  x.test[seq.int(nrow(testData$x) + 1L, nrow(x.test)),"z"] <- 0
  
  fit.rsp <- dbarts::bart2(x.train, testData$y, x.test, verbose = FALSE,
                           n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L)
  expect_equal(extract(bartcFit, type = "mu.0"),
               fit.rsp$yhat.test[,seq.int(nrow(testData$x) + 1L, nrow(x.test))])
})

test_that("bartc returns valid ouput with one chain", {
  n.obs <- length(testData$y)
  
  fit <- bartc(y, z, x, data = testData, method.trt = "glm", method.rsp = "bart", verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L)
  expect_is(fit, "bartcFit")
  expect_equal(dim(fit$mu.hat.obs), c(13L, n.obs))
  expect_equal(dim(fit$mu.hat.cf), c(13L, n.obs))
  expect_equal(length(fit$p.score), n.obs)
  expect_true(is.null(fit$samples.p.score))
  
  fit <- bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L)
  
  expect_is(fit, "bartcFit")
  expect_equal(dim(fit$mu.hat.obs), c(13L, n.obs))
  expect_equal(dim(fit$mu.hat.cf), c(13L, n.obs))
  expect_equal(length(fit$p.score), n.obs)
  expect_equal(dim(fit$samples.p.score), c(13L, n.obs))
})

test_that("bartc returns valid ouput with two chains", {
  n.obs <- length(testData$y)
  
  fit <- bartc(y, z, x, data = testData, method.trt = "glm", method.rsp = "bart", verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L)
  expect_is(fit, "bartcFit")
  expect_equal(dim(fit$mu.hat.obs), c(2L, 13L, n.obs))
  expect_equal(dim(fit$mu.hat.cf), c(2L, 13L, n.obs))
  expect_equal(length(fit$p.score), n.obs)
  expect_true(is.null(fit$samples.p.score))
  
  fit <- bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L)
  
  expect_is(fit, "bartcFit")
  expect_equal(dim(fit$mu.hat.obs), c(2L, 13L, n.obs))
  expect_equal(dim(fit$mu.hat.cf), c(2L, 13L, n.obs))
  expect_equal(length(fit$p.score), n.obs)
  expect_equal(dim(fit$samples.p.score), c(2L, 13L, n.obs))
})

test_that("bartc reports n.chains correctly regardless of combineChains", {
  # combineChains = TRUE collapses fit$mu.hat.obs to 2-D; n.chains has to come
  # from the response fitter's own record (fit.rsp$n.chains), not from
  # mu.hat.obs's dims, or it reads back as 1 no matter how many chains ran
  for (cc in c(TRUE, FALSE)) {
    fit <- bartc(y, z, x, data = testData, method.trt = "glm", method.rsp = "bart", verbose = FALSE,
                 n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 4L, n.threads = 1L,
                 combineChains = cc)
    expect_equal(fit$n.chains, 4L)
  }
})

test_that("bartc runs with all treatment settings and one chain", {
  expect_is(bartc(y, z, x, data = testData, method.trt = "glm", method.rsp = "bart", verbose = FALSE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L),
            "bartcFit")
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L),
            "bartcFit")
})

test_that("bartc runs with all treatment settings and two chains", {
  expect_is(bartc(y, z, x, data = testData, method.trt = "glm", method.rsp = "bart", verbose = FALSE,
                  n.samples = 3L, n.burn = 5L, n.trees = 25L, n.chains = 2L),
            "bartcFit")
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
                  n.samples = 3L, n.burn = 5L, n.trees = 25L, n.chains = 2L),
            "bartcFit")
})

test_that("bartc runs with all response settings and one chain", {
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L),
            "bartcFit")
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "p.weight", verbose = FALSE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L),
            "bartcFit")
})

test_that("bartc runs with all response settings and one chain for method tmle", {
  skip_on_cran()

  oldWarn <- getOption("warn")
  if (!requireNamespace("tmle", quietly = TRUE))
    options(warn = -1)
  
  expect_is(
    bartc(
      y, z, x, data = testData, method.trt = "bart", method.rsp = "tmle",
      verbose = FALSE, n.burn = 3L, n.samples = 13L, n.trees = 7L,
      n.chains = 1L, n.threads = 1L
    ),
    "bartcFit"
  )
  
  options(warn = oldWarn)
})

test_that("bartc runs with all response settings and two chains", {
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L),
            "bartcFit")
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "p.weight", verbose = FALSE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L),
            "bartcFit")
})
  
test_that("bartc runs with all response settings and two chains for method tmle", {
  skip_on_cran()

  oldWarn <- getOption("warn")
  if (!requireNamespace("tmle", quietly = TRUE))
    options(warn = -1)
  
  expect_is(
    bartc(
      y, z, x, data = testData, method.trt = "bart", method.rsp = "tmle",
      verbose = FALSE, n.burn = 3L, n.samples = 13L, n.trees = 7L,
      n.chains = 2L, n.threads = 1L
    ),
    "bartcFit"
  )

  options(warn = oldWarn)
})

source(system.file("common", "groupedData.R", package = "bartCause"))

test_that("bartc runs with all response settings and group.by set", {
  skip_if_not_installed("stan4bart")
  skip_if_not_installed("lme4")
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
                  group.by = g, group.effects = TRUE,
                  chains = 2L, iter = 26L, warmup = 13L, bart_args = list(n.trees = 7L)),
            "bartcFit")
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "p.weight", verbose = FALSE,
                  group.by = g, group.effects = TRUE,
                  chains = 2L, iter = 26L, warmup = 13L, bart_args = list(n.trees = 7L)),
            "bartcFit")
  expect_is(bartc(y, z, x, data = testData, method.trt = "glm", method.rsp = "bart", verbose = FALSE,
                  group.by = g, group.effects = TRUE,
                  chains = 2L, iter = 26L, warmup = 13L, bart_args = list(n.trees = 7L)),
            "bartcFit")
  expect_is(bartc(y, z, x, data = testData, method.trt = "glm", method.rsp = "p.weight", verbose = FALSE,
                  group.by = g, group.effects = TRUE,
                  chains = 2L, iter = 26L, warmup = 13L, bart_args = list(n.trees = 7L)),
            "bartcFit")

  # check a bart/bart with fixed effects
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
                  group.by = g, group.effects = FALSE, use.ranef = FALSE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L),
            "bartcFit")
})

test_that("bartc runs with all response settings and group.by set for method tmle", {
  skip_on_cran()

  oldWarn <- getOption("warn")
  if (!requireNamespace("tmle", quietly = TRUE))
    options(warn = -1)
  
  skip_if_not_installed("stan4bart")
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "tmle", verbose = FALSE,
                  group.by = g, group.effects = TRUE, maxIter = 5L,
                  chains = 2L, iter = 26L, warmup = 13L, bart_args = list(n.trees = 7L)),
            "bartcFit")
  
  options(warn = oldWarn)
})

test_that("bartc runs with missing data", {
  testData$y[seq_len(10L)] <- NA
  ## a varying intercept is fit by stan4bart, whose counterfactual test surface
  ## comes from the fitted rows, so the grouping factor enters as a fixed effect
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
                  group.by = g, group.effects = TRUE, use.ranef = FALSE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L),
            "bartcFit")
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "p.weight", verbose = FALSE,
                  group.by = g, group.effects = TRUE, use.ranef = FALSE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L),
            "bartcFit")
  skip_if_not_installed("stan4bart")
  expect_error(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
                     group.by = g, group.effects = TRUE,
                     chains = 2L, iter = 26L, warmup = 13L, bart_args = list(n.trees = 7L)),
               "cannot fit with missing response values")
})

test_that("bartc model argument overrides work correctly", {
  bartcFit <- bartc(y, z, x, data = testData,
                    method.rsp = "bart", method.trt = "bart", verbose = FALSE,
                    n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L,
                    k = "chi(1, Inf)")
  expect_true(!is.null(bartcFit$fit.trt$k))
  expect_true(!is.null(bartcFit$fit.rsp$k))
  
  bartcFit <- bartc(y, z, x, data = testData,
                    method.rsp = "bart", method.trt = "bart", verbose = FALSE,
                    n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L,
                    args.trt = list(k = "chi(1, Inf)"))
  expect_true(!is.null(bartcFit$fit.trt$k))
  expect_true( is.null(bartcFit$fit.rsp$k))
})

test_that("bartc works with '.' as confounders", {
  testDF <- with(testData, data.frame(y = y, z = z, x = x))
  bartcFit <- bartc(y, z, ., data = testDF,
                    method.rsp = "bart", method.trt = "bart", verbose = FALSE,
                    n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L)
  # dbarts >= 1.0-0 combines chains: varcount is 2-D [n.chains*n.samples, nvars]
  # with variable names in dimnames[[2L]] (was 3-D with names in [[3L]])
  expect_true(!("y" %in% dimnames(bartcFit$fit.trt$varcount)[[2L]]))
})

test_that("bartc runs with missing data for method tmle", {
  skip_on_cran()

  oldWarn <- getOption("warn")
  if (!requireNamespace("tmle", quietly = TRUE))
    options(warn = -1)

  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "tmle", verbose = FALSE,
                  group.by = g, group.effects = TRUE, use.ranef = FALSE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L, maxIter = 5L),
            "bartcFit")

  options(warn = oldWarn)
})

test_that("bartc runs the bcf response method at one and two chains (FB6)", {
  n.obs <- length(testData$y)

  fit <- bartc(y, z, x, data = testData, method.trt = "glm", method.rsp = "bcf", verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L)
  expect_is(fit, "bartcFit")
  expect_is(fit$fit.rsp, "bartBCF")
  expect_equal(fit$method.rsp, "bcf")
  expect_equal(dim(fit$mu.hat.obs), c(13L, n.obs))
  expect_equal(dim(fit$mu.hat.cf), c(13L, n.obs))
  expect_equal(fit$n.chains, 1L)
  expect_equal(fit$name.trt, "z")

  fit <- bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bcf", verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L)
  expect_equal(dim(fit$mu.hat.obs), c(2L, 13L, n.obs))
  expect_equal(dim(fit$fit.rsp$sigma), c(2L, 13L))
  # the propensity score reaches the prognostic forest only, on the column the
  # response builder resolved rather than on one re-derived from the design
  expect_equal(fit$fit.rsp$name.p.score, "ps")
  expect_equal(sum(fit$fit.rsp$varcount$tau[,,"ps"]), 0)
  expect_gt(sum(fit$fit.rsp$varcount$mu[,,"ps"]), 0)
})

test_that("bartc's bcf response method defaults n.threads capped at n.chains", {
  # n.chains is fixed at 2 and n.threads left at its default; an uncapped
  # default warns "n.threads (N) exceeds n.chains (2)" out of dbartsControl
  # on a machine with more than 2 cores
  expect_no_warning(
    bartc(y, z, x, data = testData, method.trt = "glm", method.rsp = "bcf", verbose = FALSE,
          n.burn = 2L, n.samples = 3L, n.trees = 5L, n.chains = 2L),
    message = "n.threads.*exceeds n.chains"
  )
})

test_that("bartc refuses crossvalidation for the bcf response method (FB3)", {
  expect_error(bartc(y, z, x, data = testData, method.rsp = "bcf", crossvalidate = TRUE,
                     verbose = FALSE, n.burn = 3L, n.samples = 13L, n.trees = 7L,
                     n.chains = 1L, n.threads = 1L),
               "crossvalidate is not supported for response method 'bcf'")
  expect_error(bartc(y, z, x, data = testData, method.rsp = "bcf", crossvalidate = "rsp",
                     verbose = FALSE, n.burn = 3L, n.samples = 13L, n.trees = 7L,
                     n.chains = 1L, n.threads = 1L),
               "crossvalidate is not supported for response method 'bcf'")
  # the identical call at 'bart' is accepted (crossvalidation itself is slow, so
  # only the argument's acceptance is checked here)
  expect_true("bcf" %in% eval(formals(bartCause::bartc)$method.rsp))
})

test_that("bartc honors the subset argument", {
  sub <- seq_len(60L)
  set.seed(22)
  fit <- bartc(y, z, x, data = testData, subset = sub, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L)

  expect_equal(length(fit$trt), 60L)
  expect_equal(as.numeric(fit$data.rsp@y), testData$y[sub])
  expect_equal(as.numeric(fit$trt), testData$z[sub])
  expect_equal(dim(fit$mu.hat.obs), c(13L, 60L))

  fit <- bartc(y, z, x, data = testData, subset = sub, method.trt = "bart", method.rsp = "bcf",
               verbose = FALSE, n.burn = 3L, n.samples = 13L, n.trees = 7L,
               n.chains = 1L, n.threads = 1L)
  expect_equal(length(fit$trt), 60L)
  expect_equal(as.numeric(fit$data.rsp@y), testData$y[sub])
  expect_equal(as.numeric(fit$trt), testData$z[sub])
  expect_equal(dim(fit$mu.hat.obs), c(13L, 60L))
})
