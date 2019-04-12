context("common support diagnostics")

source(system.file("common", "linearData.R", package = "bartCause"))

n.train <- 80L
x <- testData$x[seq_len(n.train),]
y <- testData$y[seq_len(n.train)]
z <- testData$z[seq_len(n.train)]

x.new <- testData$x[seq.int(n.train + 1L, nrow(testData$x)),]
n.test <- nrow(x.new)


test_that("predict gives sane results", {
  n.samples <- 15L
  n.chains  <- 2L
  fit <- bartc(y, z, x, method.trt = "glm", method.rsp = "bart",
               n.samples = n.samples, n.burn = 10L, n.chains = n.chains,
               n.threads = 1L, n.trees = 10L, keepTrees = TRUE,
               verbose = FALSE)
  
  # check predict for single row
  expect_equal(length(predict(fit, x.new[1,], value = "mu.0")), n.samples * n.chains)
  
  p.score <- predict(fit, x.new, value = "p.score")
  mu.1    <- predict(fit, x.new, value = "mu.1", combineChains = FALSE)
  mu.0    <- predict(fit, x.new, value = "mu.0", combineChains = TRUE)
  icate   <- predict(fit, x.new, value = "icate", combineChains = TRUE)
  
  expect_true(is.null(dim(p.score)))
  expect_equal(dim(mu.1), c(n.test, n.chains, n.samples))
  expect_equal(dim(mu.0), c(n.test, n.chains * n.samples))
  expect_equal(as.vector(icate), as.vector(mu.1) - as.vector(mu.0))
})

test_that("predict results matches training data", {
  n.samples <- 15L
  n.chains  <- 2L
  fit <- bartc(y, z, x, method.trt = "bart", method.rsp = "bart",
               n.samples = n.samples, n.burn = 10L, n.chains = n.chains,
               n.threads = 1L, n.trees = 10L, keepTrees = TRUE,
               args.trt = list(k = 1.5), verbose = FALSE)
  
  p.score <- extract(fit, value = "p.score")
  mu.1    <- extract(fit, value = "mu.1")
  mu.0    <- extract(fit, value = "mu.0")
  icate   <- extract(fit, value = "icate")
  
  p.score.new <- predict(fit, x, value = "p.score")
  mu.1.new    <- predict(fit, x, value = "mu.1")
  mu.0.new    <- predict(fit, x, value = "mu.0")
  icate.new   <- predict(fit, x, value = "icate")
  
  expect_equal(p.score, p.score.new)
  expect_equal(mu.0, mu.0.new)
  expect_equal(mu.1, mu.1.new)
  expect_equal(icate, icate.new)
})

rm(n.train, n.test, x, y, z, x.new)

# testData$g <- sample(3L, nrow(testData$x), replace = TRUE)

