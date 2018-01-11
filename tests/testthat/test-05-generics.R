context("generic functions")

source(system.file("common", "linearData.R", package = "cibart"))

set.seed(22)
testData$g <- sample(3L, nrow(testData$x), replace = TRUE)
groups <- unique(testData$g)
set.seed(22)
fit <- cibart(y, z, x, testData, method.trt = "glm", group.by = g, n.samples = 50L, n.burn = 25L, n.thread = 1L, verbose = FALSE)

p.score <- fitted(glm(z ~ x, family = stats::binomial, data = testData))

set.seed(22)
x.train <- cbind(testData$x, z = testData$z, p.score)
x.test  <- x.train; x.test[,"z"] <- 1 - x.test[,"z"]
bartFit <- dbarts::bart2(x.train, testData$y, x.test, n.samples = 50L, n.burn = 25L, n.thread = 1L, verbose = FALSE)

samples.y  <- aperm(bartFit$yhat.train, c(3L, 1L, 2L))
samples.cf <- aperm(bartFit$yhat.test, c(3L, 1L, 2L))
samples.y0 <- samples.y * (1 - testData$z) + samples.cf * testData$z
samples.y1 <- samples.y * testData$z + samples.cf * (1 - testData$z)

test_that("fitted matches manual fit", {
  est <- fitted(fit, "est")
  y1 <- fitted(fit, "y1")
  y0 <- fitted(fit, "y0")
  indiv.diff <- fitted(fit, "indiv.diff")
  y <- fitted(fit, "y")
  
  expect_equal(y0, apply(samples.y0, 1L, mean))
  expect_equal(y1, apply(samples.y1, 1L, mean)) 
  
  expect_equal(indiv.diff, apply(samples.y1 - samples.y0, 1L, mean))
  expect_equal(y, apply(samples.y, 1L, mean))
  expect_equal(fitted(fit, "p.score"), p.score)
  
  expect_equal(length(est), length(groups))
  for (group in groups)
    expect_equal(est[[as.character(group)]], mean(indiv.diff[testData$g == group]))
})

test_that("extract matches manual fit", {
  ## first that combine chains works
  y0 <- extract(fit, "y0")
  expect_equal(y0, matrix(samples.y0, dim(samples.y0)[1L], dim(samples.y0)[2L] * dim(samples.y0)[3L]))
  
  y0 <- extract(fit, "y0", combineChains = FALSE)
  est <- extract(fit, "est", combineChains = FALSE)
  y1 <- extract(fit, "y1", combineChains = FALSE)
  indiv.diff <- extract(fit, "indiv.diff", combineChains = FALSE)
  y <- extract(fit, "y", combineChains = FALSE)
  
  expect_equal(y0, samples.y0) 
  expect_equal(y1, samples.y1)
  
  expect_equal(indiv.diff, samples.y1 - samples.y0)
  expect_equal(y, samples.y)
  
  expect_equal(length(est), length(groups))
  for (group in groups)
    expect_equal(est[[as.character(group)]], apply(indiv.diff[testData$g == group,,], c(2L, 3L), mean))
})

test_that("summary object contain correct information", {
  sum <- summary(fit)
  
  expect_true(sum$call == parse(text = 'cibart(response = y, treatment = z, confounders = x, data = testData, method.trt = "glm", group.by = g, n.samples = 50L, n.burn = 25L, n.thread = 1L, verbose = FALSE)')[[1L]])
  expect_equal(sum$method.rsp, "bart")
  expect_equal(sum$method.trt, "glm")
  expect_equal(sum$ci.style, eval(formals(cibart:::summary.cibartFit)$ci.style)[1L])
  expect_equal(sum$ci.level, eval(formals(cibart:::summary.cibartFit)$ci.level))
  expect_equal(sum$numObservations, length(testData$y))
  expect_equal(sum$numSamples, 50L * 4L)
  expect_equal(sum$n.chains, 4L)
  expect_equal(sum$estimates$estimate, unname(fitted(fit, "est")))
})

test_that("summary works with confidence interval styles", {
  expect_is(summary(fit, ci.style = "norm"),  "cibartFit.summary")
  expect_is(summary(fit, ci.style = "quant"), "cibartFit.summary")
  expect_is(summary(fit, ci.style = "hpd"),   "cibartFit.summary")
})

