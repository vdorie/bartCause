context("common support diagnostics")

source(system.file("common", "linearData.R", package = "BartCause"))
testData$g <- sample(3L, nrow(testData$x), replace = TRUE)

test_that("sd common support diagnostic works", {
  expect_is(bartc(y, z, x, testData,
                  method.rsp = "tmle", method.trt = "bart", estimand = "att", verbose = FALSE,
                  n.samples = 5L, n.burn = 5L, n.chains = 1L, n.thread = 1L,
                  commonSup.rule = "sd", maxIter = 5L), "bartcFit")
  fit <- bartc(y, z, x, testData,
               method.rsp = "tmle", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.samples = 5L, n.burn = 5L, n.chains = 1L, n.thread = 1L, group.by = g,
               commonSup.rule = "sd", maxIter = 5L)
  expect_is(fit, "bartcFit")
  expect_is(summary(fit), "bartcFit.summary")
})

test_that("chisq common support diagnostic works", {
  expect_is(bartc(y, z, x, testData,
                  method.rsp = "tmle", method.trt = "bart", estimand = "att", verbose = FALSE,
                  n.samples = 5L, n.burn = 5L, n.chains = 1L, n.thread = 1L,
                  commonSup.rule = "chisq", maxIter = 5L), "bartcFit")
  expect_is(bartc(y, z, x, testData,
                  method.rsp = "tmle", method.trt = "bart", estimand = "att", verbose = FALSE,
                  n.samples = 5L, n.burn = 5L, n.chains = 1L, n.thread = 1L, group.by = g,
                  commonSup.rule = "chisq", maxIter = 5L), "bartcFit")
})

