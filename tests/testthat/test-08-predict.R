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
  expect_equal(length(predict(fit, x.new[1,], value = "y0")), n.samples * n.chains)
  
  p.score <- predict(fit, x.new, value = "p.score")
  y1      <- predict(fit, x.new, value = "y1", combineChains = FALSE)
  y0      <- predict(fit, x.new, value = "y0", combineChains = TRUE)
  ite     <- predict(fit, x.new, value = "indiv.diff", combineChains = TRUE)
  
  expect_true(is.null(dim(p.score)))
  expect_equal(dim(y1), c(n.test, n.chains, n.samples))
  expect_equal(dim(y0), c(n.test, n.chains * n.samples))
  expect_equal(as.vector(ite), as.vector(y1) - as.vector(y0))
})

test_that("predict results matches training data", {
  n.samples <- 15L
  n.chains  <- 2L
  fit <- bartc(y, z, x, method.trt = "bart", method.rsp = "bart",
               n.samples = n.samples, n.burn = 10L, n.chains = n.chains,
               n.threads = 1L, n.trees = 10L, keepTrees = TRUE,
               args.trt = list(k = 1.5), verbose = FALSE)
  
  p.score <- extract(fit, value = "p.score")
  y1      <- extract(fit, value = "y1")
  y0      <- extract(fit, value = "y0")
  ite     <- extract(fit, value = "indiv.diff")
  
  p.score.new <- predict(fit, x, value = "p.score")
  y1.new      <- predict(fit, x, value = "y1")
  y0.new      <- predict(fit, x, value = "y0")
  ite.new     <- predict(fit, x, value = "indiv.diff")
  
  expect_equal(p.score, p.score.new)
  expect_equal(y0, y0.new)
  expect_equal(y1, y1.new)
  expect_equal(ite, ite.new)
})

rm(n.train, n.test, x, y, z, x.new)

# testData$g <- sample(3L, nrow(testData$x), replace = TRUE)

