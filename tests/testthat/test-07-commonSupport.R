context("common support diagnostics")

source(system.file("common", "linearData.R", package = "cibart"))
testData$g <- sample(3L, nrow(testData$x), replace = TRUE)

test_that("sd common support diagnostic works", {
  expect_is(cibart(y, z, x, testData,
                method.rsp = "tmle", method.trt = "bart", estimand = "att", verbose = FALSE,
                n.samples = 5L, n.burn = 5L, n.chains = 1L, n.thread = 1L,
                commonSup.rule = "sd", maxIter = 5L), "cibartFit")
  fit <- cibart(y, z, x, testData,
                method.rsp = "tmle", method.trt = "bart", estimand = "att", verbose = FALSE,
                n.samples = 5L, n.burn = 5L, n.chains = 1L, n.thread = 1L, group.by = g,
                commonSup.rule = "sd", maxIter = 5L)
  expect_is(fit, "cibartFit")
  expect_is(summary(fit), "cibartFit.summary")
})

test_that("chisq common support diagnostic works", {
  expect_is(cibart(y, z, x, testData,
                method.rsp = "tmle", method.trt = "bart", estimand = "att", verbose = FALSE,
                n.samples = 5L, n.burn = 5L, n.chains = 1L, n.thread = 1L,
                commonSup.rule = "chisq", maxIter = 5L), "cibartFit")
  expect_is(cibart(y, z, x, testData,
                method.rsp = "tmle", method.trt = "bart", estimand = "att", verbose = FALSE,
                n.samples = 5L, n.burn = 5L, n.chains = 1L, n.thread = 1L, group.by = g,
                commonSup.rule = "chisq", maxIter = 5L), "cibartFit")
})

