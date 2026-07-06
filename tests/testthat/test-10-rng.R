context("rng")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("bartc with fixed seed is reproducible", {
  # As of dbarts 1.0-0, a seeded fit no longer depends on the thread count:
  # each chain runs its own generator seeded deterministically from `seed`, so
  # the same seed reproduces bit-for-bit whether run with one thread or two, and
  # different seeds give different draws.
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
  expect_equal(fit3$p.score, fit4$p.score)

  # same seed is now invariant to the thread count
  expect_equal(fit1$mu.hat.obs, fit3$mu.hat.obs)
  expect_equal(fit1$p.score, fit3$p.score)

  # a different seed gives different draws
  fit5 <- bartc(y, z, x, data = testData,
                method.rsp = "bart", method.trt = "bart", verbose = FALSE,
                n.samples = 5L, n.burn = 0L, n.trees = 7L, n.chains = 2L, n.threads = 1L,
                seed = 67890L)

  expect_true(any(fit1$mu.hat.obs != fit5$mu.hat.obs))
  expect_true(any(fit1$p.score != fit5$p.score))
})
