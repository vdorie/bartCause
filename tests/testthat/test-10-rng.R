context("rng")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("bartc with fixed seed is reproducible", {
  fit1 <- bartc(y, z, x, data = testData,
                method.rsp = "bart", method.trt = "bart", verbose = FALSE,
                n.samples = 5L, n.burn = 0L, n.trees = 7L, n.chains = 2L, n.threads = 1L,
                seed = 12345L)
  
  fit2 <- bartc(y, z, x, data = testData,
                method.rsp = "bart", method.trt = "bart", verbose = FALSE,
                n.samples = 5L, n.burn = 0L, n.trees = 7L, n.chains = 2L, n.threads = 1L,
                seed = 12345L)
  
  expect_equal(fit1$mu.hat.obs, fit2$mu.hat.obs)
  expect_equal(fit1$p.score, fit2$p.score)
  
  fit3 <- bartc(y, z, x, data = testData,
                method.rsp = "bart", method.trt = "bart", verbose = FALSE,
                n.samples = 5L, n.burn = 0L, n.trees = 7L, n.chains = 2L, n.threads = 2L,
                seed = 12345L)
  
  fit4 <- bartc(y, z, x, data = testData,
                method.rsp = "bart", method.trt = "bart", verbose = FALSE,
                n.samples = 5L, n.burn = 0L, n.trees = 7L, n.chains = 2L, n.threads = 2L,
                seed = 12345L)
  
  expect_equal(fit3$mu.hat.obs, fit4$mu.hat.obs)
  expect_equal(fit1$p.score, fit2$p.score)
  
  expect_true(any(fit1$mu.hat.obs != fit3$mu.hat.obs))
  expect_equal(fit1$p.score, fit2$p.score)
})

