context("bartc main function")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("bartc fails with invalid inputs", {
  expect_error(bartc(not.found, z, x, testData, verbose = FALSE))
  expect_error(bartc(y, not.found, x, testData, verbose = FALSE))
  expect_error(bartc(y, z, not.found, testData, verbose = FALSE))
  expect_error(bartc(y, z, x, testData, verbose = FALSE, method.rsp = "not-a-method"))
  expect_error(bartc(y, z, x, testData, verbose = FALSE, method.trt = "not-a-method"))
  expect_error(bartc(y, z, x, testData, verbose = FALSE, estimand = "not-an-estimand"))
  expect_error(bartc(y, z, x, testData, verbose = FALSE, group.by = not.found))
  expect_error(bartc(y, z, x, testData, verbose = FALSE, p.scoreAsCovariate = NA))
  expect_error(bartc(y, z, x, testData, verbose = FALSE, keepCall = NA))
  expect_error(bartc(y, z, x, testData, verbose = NA))
})

test_that("bartc matches manual fit", {  
  set.seed(22)
  bartcFit <- bartc(y, z, x, testData,
                    method.rsp = "bart", method.trt = "bart", verbose = FALSE,
                    n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L)
  
  set.seed(22)
  fit.trt <- bart2(z ~ x, testData, verbose = FALSE,
                   n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L)
  p.score <- apply(pnorm(fit.trt$yhat.train), 2L, mean)
  expect_equal(p.score, fitted(bartcFit, value = "p.score"))
  
  x.train <- cbind(testData$x, z = testData$z, ps = p.score)
  x.test  <- cbind(testData$x, z = 1, ps = p.score)
  x.test <- rbind(x.test, x.test)
  x.test[seq.int(nrow(testData$x) + 1L, nrow(x.test)),"z"] <- 0
  
  fit.rsp <- bart2(x.train, testData$y, x.test, verbose = FALSE,
                   n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L)
  expect_equal(extract(bartcFit, value = "y0"),
               t(fit.rsp$yhat.test[,seq.int(nrow(testData$x) + 1L, nrow(x.test))]))
})

test_that("bartc returns valid ouput with one chain", {
  n.obs <- length(testData$y)
  
  fit <- bartc(y, z, x, testData, method.trt = "glm", method.rsp = "bart",
               n.samples = 200L, n.burn = 100L, n.chains = 1L, verbose = FALSE)
  expect_is(fit, "bartcFit")
  expect_equal(length(fit$samples.est), 200L)
  expect_equal(dim(fit$yhat.obs), c(n.obs, 200L))
  expect_equal(dim(fit$yhat.cf), c(n.obs, 200L))
  expect_equal(length(fit$p.score), n.obs)
  expect_true(is.null(fit$samples.p.score))
  
  fit <- bartc(y, z, x, testData, method.trt = "bart", method.rsp = "bart",
               n.samples = 200L, n.burn = 100L, n.chains = 1L, verbose = FALSE)
  
  expect_is(fit, "bartcFit")
  expect_equal(length(fit$samples.est), 200L)
  expect_equal(dim(fit$yhat.obs), c(n.obs, 200L))
  expect_equal(dim(fit$yhat.cf), c(n.obs, 200L))
  expect_equal(length(fit$p.score), n.obs)
  expect_equal(dim(fit$samples.p.score), c(n.obs, 200L))
})

test_that("bartc returns valid ouput with two chains", {
  n.obs <- length(testData$y)
  
  fit <- bartc(y, z, x, testData, method.trt = "glm", method.rsp = "bart",
               n.samples = 100L, n.burn = 50L, n.chains = 2L, verbose = FALSE)
  expect_is(fit, "bartcFit")
  expect_equal(length(fit$samples.est), 200L)
  expect_equal(dim(fit$samples.est), c(2L, 100L))
  expect_equal(dim(fit$yhat.obs), c(n.obs, 2L, 100L))
  expect_equal(dim(fit$yhat.cf), c(n.obs, 2L, 100L))
  expect_equal(length(fit$p.score), n.obs)
  expect_true(is.null(fit$samples.p.score))
  
  fit <- bartc(y, z, x, testData, method.trt = "bart", method.rsp = "bart",
               n.samples = 100L, n.burn = 50L, n.chains = 2L, verbose = FALSE)
  
  expect_is(fit, "bartcFit")
  expect_equal(length(fit$samples.est), 200L)
  expect_equal(dim(fit$yhat.obs), c(n.obs, 2L, 100L))
  expect_equal(dim(fit$yhat.cf), c(n.obs, 2L, 100L))
  expect_equal(length(fit$p.score), n.obs)
  expect_equal(dim(fit$samples.p.score), c(n.obs, 2L, 100L))
})

test_that("bartc runs with all treatment settings and one chain", {
  expect_is(bartc(y, z, x, testData, method.trt = "glm", method.rsp = "bart",
                  n.samples = 20L, n.burn = 10L, n.trees = 25L, n.chains = 1L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "bart",
                  n.samples = 20L, n.burn = 10L, n.trees = 25L, n.chains = 1L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart.xval", method.rsp = "bart",
                  n.samples = 20L, n.burn = 10L, n.trees = 25L, n.chains = 1L, n.threads = 2L, verbose = FALSE),
            "bartcFit")
})

test_that("bartc runs with all treatment settings and two chains", {
  expect_is(bartc(y, z, x, testData, method.trt = "glm", method.rsp = "bart",
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "bart",
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart.xval", method.rsp = "bart",
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, n.threads = 2L, verbose = FALSE),
            "bartcFit")
})

test_that("bartc runs with all response settings and one chain", {
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "bart",
                  n.samples = 20L, n.burn = 10L, n.trees = 25L, n.chains = 1L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "p.weight",
                  n.samples = 20L, n.burn = 10L, n.trees = 25L, n.chains = 1L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "tmle",
                  n.samples = 20L, n.burn = 10L, n.trees = 25L, n.chains = 1L, verbose = FALSE),
            "bartcFit")
})

test_that("bartc runs with all response settings and two chains", {
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "bart",
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "p.weight",
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "tmle",
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "bartcFit")
})

test_that("bartc runs with all response settings and group.by set", {
  testData$g <- sample(3L, nrow(testData$x), replace = TRUE)
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "bart", group.by = g,
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "p.weight", group.by = g,
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "tmle", group.by = g,
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE, maxIter = 10),
            "bartcFit")
  
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "bart", group.by = g, use.rbart = TRUE,
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "p.weight", group.by = g, use.rbart = TRUE,
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "tmle", group.by = g, use.rbart = TRUE,
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE, maxIter = 10),
            "bartcFit")
})

test_that("bartc runs with missing data", {
  testData$g <- sample(3L, nrow(testData$x), replace = TRUE)
  testData$y[seq_len(10L)] <- NA
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "bart", group.by = g,
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "p.weight", group.by = g,
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "bartcFit")
  expect_is(bartc(y, z, x, testData, method.trt = "bart", method.rsp = "tmle", group.by = g,
                  n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "bartcFit")
})


