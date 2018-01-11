context("cibart main function")

source(system.file("common", "linearData.R", package = "cibart"))

test_that("cibart fails with invalid inputs", {
  expect_error(cibart(not.found, z, x, testData, verbose = FALSE))
  expect_error(cibart(y, not.found, x, testData, verbose = FALSE))
  expect_error(cibart(y, z, not.found, testData, verbose = FALSE))
  expect_error(cibart(y, z, x, testData, verbose = FALSE, method.rsp = "not-a-method"))
  expect_error(cibart(y, z, x, testData, verbose = FALSE, method.trt = "not-a-method"))
  expect_error(cibart(y, z, x, testData, verbose = FALSE, estimand = "not-an-estimand"))
  expect_error(cibart(y, z, x, testData, verbose = FALSE, group.by = not.found))
  expect_error(cibart(y, z, x, testData, verbose = FALSE, propensityScoreAsCovariate = NA))
  expect_error(cibart(y, z, x, testData, verbose = FALSE, keepCall = NA))
  expect_error(cibart(y, z, x, testData, verbose = NA))
})

test_that("cibart returns valid ouput with one chain", {
  n.obs <- length(testData$y)
  
  fit <- cibart(y, z, x, testData, method.trt = "glm", method.rsp = "bart",
                n.samples = 200L, n.burn = 100L, n.chains = 1L, verbose = FALSE)
  expect_is(fit, "cibartFit")
  expect_equal(length(fit$samples.est), 200L)
  expect_equal(dim(fit$samples.indiv.diff), c(n.obs, 200L))
  expect_equal(length(fit$p.score), n.obs)
  expect_true(is.null(fit$samples.p.score))
  
  fit <- cibart(y, z, x, testData, method.trt = "bart", method.rsp = "bart",
                n.samples = 200L, n.burn = 100L, n.chains = 1L, verbose = FALSE)
  
  expect_is(fit, "cibartFit")
  expect_equal(length(fit$samples.est), 200L)
  expect_equal(dim(fit$samples.indiv.diff), c(n.obs, 200L))
  expect_equal(length(fit$p.score), n.obs)
  expect_equal(dim(fit$samples.p.score), c(n.obs, 200L))
})

test_that("cibart returns valid ouput with two chains", {
  n.obs <- length(testData$y)
  
  fit <- cibart(y, z, x, testData, method.trt = "glm", method.rsp = "bart",
                n.samples = 100L, n.burn = 50L, n.chains = 2L, verbose = FALSE)
  expect_is(fit, "cibartFit")
  expect_equal(length(fit$samples.est), 200L)
  expect_equal(dim(fit$samples.est), c(2L, 100L))
  expect_equal(dim(fit$samples.indiv.diff), c(n.obs, 2L, 100L))
  expect_equal(length(fit$p.score), n.obs)
  expect_true(is.null(fit$samples.p.score))
  
  fit <- cibart(y, z, x, testData, method.trt = "bart", method.rsp = "bart",
                n.samples = 100L, n.burn = 50L, n.chains = 2L, verbose = FALSE)
  
  expect_is(fit, "cibartFit")
  expect_equal(length(fit$samples.est), 200L)
  expect_equal(dim(fit$samples.indiv.diff), c(n.obs, 2L, 100L))
  expect_equal(length(fit$p.score), n.obs)
  expect_equal(dim(fit$samples.p.score), c(n.obs, 2L, 100L))
})

test_that("cibart runs with all treatment settings and one chain", {
  expect_is(cibart(y, z, x, testData, method.trt = "glm", method.rsp = "bart",
                   n.samples = 20L, n.burn = 10L, n.trees = 25L, n.chains = 1L, verbose = FALSE),
            "cibartFit")
  expect_is(cibart(y, z, x, testData, method.trt = "bart", method.rsp = "bart",
                   n.samples = 20L, n.burn = 10L, n.trees = 25L, n.chains = 1L, verbose = FALSE),
            "cibartFit")
  expect_is(cibart(y, z, x, testData, method.trt = "bart.xval", method.rsp = "bart",
                   n.samples = 20L, n.burn = 10L, n.trees = 25L, n.chains = 1L, n.threads = 2L, verbose = FALSE),
            "cibartFit")
})

test_that("cibart runs with all treatment settings and two chains", {
  expect_is(cibart(y, z, x, testData, method.trt = "glm", method.rsp = "bart",
                   n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "cibartFit")
  expect_is(cibart(y, z, x, testData, method.trt = "bart", method.rsp = "bart",
                   n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "cibartFit")
  expect_is(cibart(y, z, x, testData, method.trt = "bart.xval", method.rsp = "bart",
                   n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, n.threads = 2L, verbose = FALSE),
            "cibartFit")
})

test_that("cibart runs with all response settings and one chain", {
  expect_is(cibart(y, z, x, testData, method.trt = "bart", method.rsp = "bart",
                   n.samples = 20L, n.burn = 10L, n.trees = 25L, n.chains = 1L, verbose = FALSE),
            "cibartFit")
  expect_is(cibart(y, z, x, testData, method.trt = "bart", method.rsp = "pweight",
                   n.samples = 20L, n.burn = 10L, n.trees = 25L, n.chains = 1L, verbose = FALSE),
            "cibartFit")
  expect_is(cibart(y, z, x, testData, method.trt = "bart", method.rsp = "tmle",
                   n.samples = 20L, n.burn = 10L, n.trees = 25L, n.chains = 1L, verbose = FALSE),
            "cibartFit")
})

test_that("cibart runs with all response settings and two chains", {
  expect_is(cibart(y, z, x, testData, method.trt = "bart", method.rsp = "bart",
                   n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "cibartFit")
  expect_is(cibart(y, z, x, testData, method.trt = "bart", method.rsp = "pweight",
                   n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "cibartFit")
  expect_is(cibart(y, z, x, testData, method.trt = "bart", method.rsp = "tmle",
                   n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "cibartFit")
})

test_that("cibart runs with all response settings and group.by set", {
  testData$g <- sample(3L, nrow(testData$x), replace = TRUE)
  expect_is(cibart(y, z, x, testData, method.trt = "bart", method.rsp = "bart", group.by = g,
                   n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "cibartFit")
  expect_is(cibart(y, z, x, testData, method.trt = "bart", method.rsp = "pweight", group.by = g,
                   n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE),
            "cibartFit")
  expect_is(cibart(y, z, x, testData, method.trt = "bart", method.rsp = "tmle", group.by = g,
                   n.samples = 10L, n.burn = 5L, n.trees = 25L, n.chains = 2L, verbose = FALSE, maxIter = 10),
            "cibartFit")
})

