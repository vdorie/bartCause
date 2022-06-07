context("common support diagnostics")

source(system.file("common", "linearData.R", package = "bartCause"))
testData$g <- sample(3L, nrow(testData$x), replace = TRUE)

test_that("sd common support diagnostic works", {
  expect_is(bartc(y, z, x, data = testData,
                  method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
                  n.burn = 0L, n.samples = 3L, n.trees = 7L, n.chains = 1L, n.threads = 1L,
                  commonSup.rule = "sd", maxIter = 2L), "bartcFit")
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.burn = 0L, n.samples = 3L, n.trees = 7L, n.chains = 1L, n.threads = 1L,
               n.thin = 1L,
               group.by = g,
               commonSup.rule = "sd", maxIter = 2L)
  expect_is(fit, "bartcFit")
  expect_is(summary(fit), "bartcFit.summary")
  
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.burn = 0L, n.samples = 3L, n.trees = 7L, n.chains = 1L, n.threads = 1L,
               n.thin = 1L,
               group.by = g, group.effects = TRUE,
               commonSup.rule = "sd", maxIter = 2L)
  expect_is(fit, "bartcFit")
  expect_is(summary(fit), "bartcFit.summary")
})

test_that("chisq common support diagnostic works", {
  expect_is(bartc(y, z, x, data = testData,
                  method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
                  n.burn = 0L, n.samples = 3L, n.trees = 7L, n.chains = 1L, n.threads = 1L,
                  commonSup.rule = "chisq", maxIter = 2L), "bartcFit")
  expect_is(bartc(y, z, x, data = testData,
                  method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
                  n.burn = 0L, n.samples = 3L, n.trees = 7L, n.chains = 1L, n.threads = 1L,
                  n.thin = 1L,
                  group.by = g,
                  commonSup.rule = "chisq", maxIter = 2L), "bartcFit")
})

