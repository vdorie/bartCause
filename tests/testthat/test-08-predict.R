context("predict")

source(system.file("common", "linearData.R", package = "bartCause"), local = TRUE)

n.train <- 80L
x <- testData$x[seq_len(n.train),]
y <- testData$y[seq_len(n.train)]
z <- testData$z[seq_len(n.train)]

x.new <- testData$x[seq.int(n.train + 1L, nrow(testData$x)),]
n.test <- nrow(x.new)

test_that("predict gives sane results", {
  n.samples <- 7L
  n.chains  <- 2L
  fit <- bartc(y, z, x, method.trt = "glm", method.rsp = "bart",
               n.chains = n.chains, n.threads = 1L, n.burn = 0L, n.samples = n.samples, n.trees = 13L,
               keepTrees = TRUE,
               verbose = FALSE)
  
  # check predict for single row
  expect_equal(length(predict(fit, x.new[1,], type = "mu.0")), n.samples * n.chains)
  
  p.score <- predict(fit, x.new, type = "p.score")
  mu.1    <- predict(fit, x.new, type = "mu.1", combineChains = FALSE)
  mu.0    <- predict(fit, x.new, type = "mu.0", combineChains = TRUE)
  icate   <- predict(fit, x.new, type = "icate", combineChains = TRUE)
  
  expect_true(is.null(dim(p.score)))
  expect_equal(dim(mu.1), c(n.chains, n.samples, n.test))
  expect_equal(dim(mu.0), c(n.chains * n.samples, n.test))
  expect_equal(as.vector(icate), as.vector(matrix(aperm(mu.1, c(2L, 1L, 3L)), n.samples * n.chains)) - as.vector(mu.0))
})

test_that("predict results matches training data", {
  n.samples <- 7L
  n.chains  <- 2L
  fit <- bartc(y, z, x, method.trt = "bart", method.rsp = "bart",
               n.chains = n.chains, n.threads = 1L, n.burn = 0L, n.samples = n.samples, n.trees = 13L,
               keepTrees = TRUE,
               args.trt = list(k = 1.5), verbose = FALSE)
  
  p.score <- extract(fit, type = "p.score")
  mu.1    <- extract(fit, type = "mu.1")
  mu.0    <- extract(fit, type = "mu.0")
  icate   <- extract(fit, type = "icate")
  
  p.score.new <- predict(fit, x, type = "p.score")
  mu.1.new    <- predict(fit, x, type = "mu.1")
  mu.0.new    <- predict(fit, x, type = "mu.0")
  icate.new   <- predict(fit, x, type = "icate")
  
  expect_equal(p.score, p.score.new)
  expect_equal(mu.0, mu.0.new)
  expect_equal(mu.1, mu.1.new)
  expect_equal(icate, icate.new)
})

set.seed(22)
g <- sample(3L, nrow(x), replace = TRUE)

n.samples <- 7L
n.chains  <- 2L

test_that("predict works with grouped data, glm trt model", {
 
  fit <- bartc(y, z, x, method.trt = "glm", method.rsp = "bart", group.by = g,
               n.chains = n.chains, n.threads = 1L, n.burn = 0L, n.samples = n.samples, n.trees = 13L,
               keepTrees = TRUE, use.ranef = FALSE,
               args.trt = list(k = 1.5), verbose = FALSE)
  
  p.score <- fitted(fit, type = "p.score")
  mu.1    <- extract(fit, type = "mu.1")
  mu.0    <- extract(fit, type = "mu.0")
  icate   <- extract(fit, type = "icate")
  
  p.score.new <- predict(fit, x, group.by = g, type = "p.score")
  mu.1.new    <- predict(fit, x, group.by = g, type = "mu.1")
  mu.0.new    <- predict(fit, x, group.by = g, type = "mu.0")
  icate.new   <- predict(fit, x, group.by = g, type = "icate")
  
  expect_equal(p.score, p.score.new)
  expect_equal(mu.0, mu.0.new)
  expect_equal(mu.1, mu.1.new)
  expect_equal(icate, icate.new)
})

test_that("predict works with grouped data, glmer trt model", {
  skip_if_not_installed("lme4")

  suppressWarnings(
    fit <- bartc(y, z, x, method.trt = "glm", method.rsp = "bart", group.by = g,
                 n.chains = n.chains, n.threads = 1L, n.burn = 0L, n.samples = n.samples, n.trees = 13L,
                 keepTrees = TRUE, use.ranef = FALSE,
                 args.trt = list(k = 1.5), verbose = FALSE)
  )
  
  p.score <- fitted(fit, type = "p.score")
  mu.1    <- extract(fit, type = "mu.1")
  mu.0    <- extract(fit, type = "mu.0")
  icate   <- extract(fit, type = "icate")
  
  p.score.new <- predict(fit, x, group.by = g, type = "p.score")
  mu.1.new    <- predict(fit, x, group.by = g, type = "mu.1")
  mu.0.new    <- predict(fit, x, group.by = g, type = "mu.0")
  icate.new   <- predict(fit, x, group.by = g, type = "icate")
  
  expect_equal(p.score, p.score.new)
  expect_equal(mu.0, mu.0.new)
  expect_equal(mu.1, mu.1.new)
  expect_equal(icate, icate.new)
})

test_that("predict works with grouped data, bart trt model", {
  fit <- bartc(y, z, x, method.trt = "bart", method.rsp = "bart", group.by = g,
               n.chains = n.chains, n.threads = 1L, n.burn = 0L, n.samples = n.samples, n.trees = 13L,
               keepTrees = TRUE,
               args.trt = list(k = 1.5), verbose = FALSE)
  
  p.score <- extract(fit, type = "p.score")
  mu.1    <- extract(fit, type = "mu.1")
  mu.0    <- extract(fit, type = "mu.0")
  icate   <- extract(fit, type = "icate")
  
  p.score.new <- predict(fit, x, group.by = g, type = "p.score")
  mu.1.new    <- predict(fit, x, group.by = g, type = "mu.1")
  mu.0.new    <- predict(fit, x, group.by = g, type = "mu.0")
  icate.new   <- predict(fit, x, group.by = g, type = "icate")
  
  expect_equal(p.score, p.score.new)
  expect_equal(mu.0, mu.0.new)
  expect_equal(mu.1, mu.1.new)
  expect_equal(icate, icate.new)
  
  fit <- bartc(y, z, x, method.trt = "bart", method.rsp = "bart", group.by = g,
               n.chains = n.chains, n.threads = 1L, n.burn = 0L, n.samples = n.samples, n.trees = 13L,
               keepTrees = TRUE, use.ranef = FALSE,
               args.trt = list(k = 1.5), verbose = FALSE)
  
  p.score <- extract(fit, type = "p.score")
  mu.1    <- extract(fit, type = "mu.1")
  mu.0    <- extract(fit, type = "mu.0")
  icate   <- extract(fit, type = "icate")
  
  p.score.new <- predict(fit, x, group.by = g, type = "p.score")
  mu.1.new    <- predict(fit, x, group.by = g, type = "mu.1")
  mu.0.new    <- predict(fit, x, group.by = g, type = "mu.0")
  icate.new   <- predict(fit, x, group.by = g, type = "icate")
  
  expect_equal(p.score, p.score.new)
  expect_equal(mu.0, mu.0.new)
  expect_equal(mu.1, mu.1.new)
  expect_equal(icate, icate.new)
})

rm(testData, n.train, x, y, z, g, n.samples, n.chains, x.new, n.test)

