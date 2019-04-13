context("generic functions")

source(system.file("common", "linearData.R", package = "bartCause"))

set.seed(22)
testData$g <- sample(3L, nrow(testData$x), replace = TRUE)
groups <- unique(testData$g)
set.seed(22)
fit <- bartc(y, z, x, testData, method.trt = "glm", group.by = g, n.samples = 50L,
             n.burn = 25L, n.chains = 4L, n.threads = 1L, verbose = FALSE)

p.score <- fitted(glm(z ~ x, family = stats::binomial, data = testData))

set.seed(22)
x.train <- cbind(testData$x, z = testData$z, p.score)
x.test  <- x.train; x.test[,"z"] <- 1 - x.test[,"z"]
bartFit <- dbarts::bart2(x.train, testData$y, x.test, n.samples = 50L, n.burn = 25L,
                         n.chains = 4L, n.threads = 1L, verbose = FALSE)

samples.obs <- aperm(bartFit$yhat.train, c(3L, 1L, 2L))
samples.cf  <- aperm(bartFit$yhat.test, c(3L, 1L, 2L))
samples.mu.0 <- samples.obs * (1 - testData$z) + samples.cf * testData$z
samples.mu.1 <- samples.obs * testData$z +       samples.cf * (1 - testData$z)

test_that("fitted matches manual fit", {
  cate <- fitted(fit, "cate")
  mu.1 <- fitted(fit, "mu.1")
  mu.0 <- fitted(fit, "mu.0")
  icate <- fitted(fit, "icate")
  mu.obs <- fitted(fit, "mu")
  
  expect_equal(mu.0, apply(samples.mu.0, 1L, mean))
  expect_equal(mu.1, apply(samples.mu.1, 1L, mean)) 
  
  expect_equal(icate, apply(samples.mu.1 - samples.mu.0, 1L, mean))
  expect_equal(mu.obs, apply(samples.obs, 1L, mean))
  expect_equal(fitted(fit, "p.score"), p.score)
  
  expect_equal(length(cate), length(groups))
  for (group in groups)
    expect_equal(cate[[as.character(group)]], mean(icate[testData$g == group]))
})

test_that("extract matches manual fit", {
  ## first that combine chains works
  mu.0 <- extract(fit, "mu.0")
  expect_equal(mu.0, matrix(samples.mu.0, dim(samples.mu.0)[1L], dim(samples.mu.0)[2L] * dim(samples.mu.0)[3L]))
  
  mu.0   <- extract(fit, "mu.0",  combineChains = FALSE)
  cate   <- extract(fit, "cate",  combineChains = FALSE)
  mu.1   <- extract(fit, "mu.1",  combineChains = FALSE)
  icate  <- extract(fit, "icate", combineChains = FALSE)
  mu.obs <- extract(fit, "mu",    combineChains = FALSE)
  
  expect_equal(mu.0, samples.mu.0) 
  expect_equal(mu.1, samples.mu.1)
  
  expect_equal(icate, samples.mu.1 - samples.mu.0)
  expect_equal(mu.obs, samples.obs)
  
  expect_equal(length(cate), length(groups))
  for (group in groups)
    expect_equal(cate[[as.character(group)]], apply(icate[testData$g == group,,], c(2L, 3L), mean))
})

test_that("generics work for p.weights", {
  oldWarnLevel <- options()$warn
  options(warn = -1L)
  pfit <- bartc(y, z, x, testData, method.trt = "bart", method.rsp = "p.weight",
                estimand = "att", group.by = g, n.samples = 50L,
                n.burn = 25L, n.threads = 1L, verbose = FALSE)
  pfit.sum <- summary(pfit)
  
  p.weights  <- extract(pfit, "p.weights", sample = "all")

  g.sel <- lapply(unique(testData$g), function(level) which(testData$g == level))
  
  #boundValues <- asNamespace("bartCause")$boundValues
  boundValues <- bartCause:::boundValues
  
  ## match internal bounding
  yBounds <- c(.005, .995)
  p.scoreBounds <- c(0.025, 0.975)
  
  for (i in seq_along(unique(testData$g))) {
    m <- min(testData$y[g.sel[[i]]]); M <- max(testData$y[g.sel[[i]]])
    
    yhat.0 <- boundValues((boundValues(extract(pfit, "mu.0", sample = "all")[g.sel[[i]],], c(m, M)) - m) / (M - m), yBounds)
    yhat.1 <- boundValues((boundValues(extract(pfit, "mu.1", sample = "all")[g.sel[[i]],], c(m, M)) - m) / (M - m), yBounds)
    
    icate <- yhat.1 - yhat.0

    expect_equal(pfit.sum$est$estimate[i], mean(apply((icate * p.weights[g.sel[[i]],]), 2L, mean) * (M - m)))
  }
  options(warn = oldWarnLevel)
  
  expect_equal(apply(p.weights, 1L, mean), fitted(pfit, "p.weights", sample = "all"))
})

test_that("summary object contain correct information", {
  sum <- summary(fit)
  
  testCall <- testCall <- parse(text = 'bartc(response = y, treatment = z, confounders = x, data = testData, method.trt = "glm", group.by = g, verbose = FALSE, n.samples = 50L, n.burn = 25L, n.chains = 4L, n.threads = 1L)')[[1L]]
  
  expect_true(length(testCall) == length(sum$call) && sum$call == testCall)
  expect_equal(sum$method.rsp, "bart")
  expect_equal(sum$method.trt, "glm")
  expect_equal(sum$ci.info$ci.style, eval(formals(bartCause:::summary.bartcFit)$ci.style)[1L])
  expect_equal(sum$ci.info$ci.level, eval(formals(bartCause:::summary.bartcFit)$ci.level))
  expect_equal(sum$numObservations, length(testData$y))
  expect_equal(sum$numSamples, 50L * 4L)
  expect_equal(sum$n.chains, 4L)
  expect_equal(sum$estimates$estimate, unname(fitted(fit, "cate")))
})

test_that("summary works with different styles", {
  expect_is(summary(fit, ci.style = "norm"),  "bartcFit.summary")
  expect_is(summary(fit, ci.style = "quant"), "bartcFit.summary")
  expect_is(summary(fit, ci.style = "hpd"),   "bartcFit.summary")
  expect_is(summary(fit, pate.style = "ppd"), "bartcFit.summary")
  
  fit <- bartc(y, z, x, testData, method.trt = "glm", method.rsp = "tmle", group.by = g, n.samples = 50L,
             n.burn = 25L, n.chains = 4L, n.threads = 1L, verbose = FALSE)
  expect_is(summary(fit), "bartcFit.summary")
  expect_is(summary(fit, pate.style = "ppd"), "bartcFit.summary")
})

