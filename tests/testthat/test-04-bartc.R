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
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
                  group.by = g, group.effects = TRUE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L),
            "bartcFit")
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "p.weight", verbose = FALSE,
                  group.by = g, group.effects = TRUE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L),
            "bartcFit")
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
                  group.by = g, group.effects = TRUE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L),
            "bartcFit")
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "p.weight", verbose = FALSE,
                  group.by = g, group.effects = TRUE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L),
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
  
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "tmle", verbose = FALSE,
                  group.by = g, group.effects = TRUE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L, maxIter = 5L),
            "bartcFit")
  
  options(warn = oldWarn)
})

test_that("bartc runs with missing data", {
  testData$y[seq_len(10L)] <- NA
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
                  group.by = g, group.effects = TRUE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L),
            "bartcFit")
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "p.weight", verbose = FALSE,
                  group.by = g, group.effects = TRUE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L),
            "bartcFit")
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
  expect_true(!("y" %in% dimnames(bartcFit$fit.trt$varcount)[[3L]]))
})

test_that("bartc runs with missing data for method tmle", {
  skip_on_cran()

  oldWarn <- getOption("warn")
  if (!requireNamespace("tmle", quietly = TRUE))
    options(warn = -1)
  
  expect_is(bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "tmle", verbose = FALSE,
                  group.by = g, group.effects = TRUE,
                  n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L, maxIter = 5L),
            "bartcFit")
  
  options(warn = oldWarn)
})
