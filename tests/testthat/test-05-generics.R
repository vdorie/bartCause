context("generic functions")

source(system.file("common", "groupedData.R", package = "bartCause"))

set.seed(22)
fit <- bartc(y, z, x, testData, method.trt = "glm", n.samples = 50L,
             group.by = g, use.ranef = FALSE, group.effects = TRUE,
             n.burn = 25L, n.chains = 4L, n.threads = 1L, verbose = FALSE)

p.score <- fitted(glm(z ~ x + g, family = stats::binomial, data = testData))

set.seed(22)
x.train <- cbind(testData$x, z = testData$z, p.score, testData$g)
x.test  <- x.train; x.test[,"z"] <- 1 - x.test[,"z"]
bartFit <- dbarts::bart2(x.train, testData$y, x.test, n.samples = 50L, n.burn = 25L,
                         n.chains = 4L, n.threads = 1L, verbose = FALSE)

obsCfToTrtCtl <- function(obs, cf, trt) {
  if (length(dim(obs)) > 2L) {
    aperm(aperm(obs, c(3L, 1L, 2L)) * trt + aperm(cf, c(3L, 1L, 2L)) * (1 - trt), c(2L, 3L, 1L))
  } else {
    t(t(obs) * trt + t(cf) * (1 - trt))
  }
}

samples.obs <- bartFit$yhat.train
samples.cf  <- bartFit$yhat.test
samples.mu.0 <- obsCfToTrtCtl(samples.obs, samples.cf, 1 - testData$z)
samples.mu.1 <- obsCfToTrtCtl(samples.obs, samples.cf, testData$z)

test_that("fitted matches manual fit", {
  cate <- fitted(fit, "cate")
  mu.1 <- fitted(fit, "mu.1")
  mu.0 <- fitted(fit, "mu.0")
  icate <- fitted(fit, "icate")
  mu.obs <- fitted(fit, "mu.obs")
  
  expect_equal(mu.0, apply(samples.mu.0, 3L, mean))
  expect_equal(mu.1, apply(samples.mu.1, 3L, mean)) 
  
  expect_equal(icate, apply(samples.mu.1 - samples.mu.0, 3L, mean))
  expect_equal(mu.obs, apply(samples.obs, 3L, mean))
  expect_equal(fitted(fit, "p.score"), p.score)
  
  groups <- levels(as.factor(testData$g))
  expect_equal(length(cate), length(groups))
  for (group in groups)
    expect_equal(cate[[as.character(group)]], mean(icate[testData$g == group]))
})

test_that("extract matches manual fit", {
  ## first that combine chains works
  mu.0 <- extract(fit, "mu.0")
  expect_equal(mu.0, matrix(aperm(samples.mu.0, c(2L, 1L, 3L)), dim(samples.mu.0)[1L] * dim(samples.mu.0)[2L], dim(samples.mu.0)[3L]))
  
  mu.0   <- extract(fit, "mu.0",   combineChains = FALSE)
  mu.1   <- extract(fit, "mu.1",   combineChains = FALSE)
  cate   <- extract(fit, "cate",   combineChains = FALSE)
  icate  <- extract(fit, "icate",  combineChains = FALSE)
  mu.obs <- extract(fit, "mu.obs", combineChains = FALSE)
  
  expect_equal(mu.0, samples.mu.0) 
  expect_equal(mu.1, samples.mu.1)
  
  expect_equal(icate, samples.mu.1 - samples.mu.0)
  expect_equal(mu.obs, samples.obs)
  
  groups <- levels(as.factor(testData$g))
  expect_equal(length(cate), length(groups))
  for (group in groups)
    expect_equal(cate[[as.character(group)]], apply(icate[,,testData$g == group], c(1L, 2L), mean))
})

test_that("ppd-based estimates match manual", {
  expect_equal(as.numeric((testData$y - fitted(fit, "y.cf")) * (2 * testData$z - 1)),
               fitted(fit, "ite"))
  expect_equal(mean((testData$y - fitted(fit, "y.cf")) * (2 * testData$z - 1)),
               sum(fitted(fit, "sate") * (table(testData$g) / length(testData$y))))
  expect_equal(testData$y[testData$z == 1], fitted(fit, "y.1")[testData$z == 1])
  expect_equal(testData$y[testData$z == 0], fitted(fit, "y.0")[testData$z == 0])
})

test_that("summary object contains correct information", {
  sum <- summary(fit)
  
  testCall <- parse(text = 'bartc(response = y, treatment = z, confounders = x, data = testData, method.trt = "glm", group.by = g, group.effects = TRUE, use.ranef = FALSE, n.samples = 50L,
             n.burn = 25L, n.chains = 4L, n.threads = 1L, verbose = FALSE)')[[1L]]
  
  expect_true(length(testCall) == length(sum$call) && all(names(sum$call) %in% names(testCall)) &&
              all(names(testCall) %in% names(sum$call)) && all(sum$call[order(names(sum$call))] == testCall[order(names(testCall))]))
  expect_equal(sum$method.rsp, "bart")
  expect_equal(sum$method.trt, "glm")
  expect_equal(sum$ci.info$ci.style, eval(formals(bartCause:::summary.bartcFit)$ci.style)[1L])
  expect_equal(sum$ci.info$ci.level, eval(formals(bartCause:::summary.bartcFit)$ci.level))
  expect_equal(tail(sum$n.obs, 1L), length(testData$y))
  expect_equal(sum$n.samples, 50L)
  expect_equal(sum$n.chains, 4L)
  expect_equal(head(sum$estimates$estimate, -1L), unname(fitted(fit, "cate")))
})

test_that("generics work for p.weights", {
  pfit <- bartc(y, z, x, testData, method.trt = "bart", method.rsp = "p.weight", estimand = "att",
                group.by = g, group.effects = TRUE, n.chains = 3L,
                n.samples = 7L, n.burn = 3L, n.threads = 1L, verbose = FALSE)
  pfit.sum <- summary(pfit)
  
  p.weights  <- extract(pfit, "p.weights", sample = "all")
  
  groups <- levels(as.factor(testData$g))
  g.sel <- lapply(groups, function(group) which(testData$g == group))
  
  boundValues <- bartCause:::boundValues
  
  ## match internal bounding
  yBounds <- c(.005, .995)
  p.scoreBounds <- c(0.025, 0.975)
  
  # warnings because "mu.0" isn't meaningful if using p-weights to compute ATT
  mu.0 <- suppressWarnings(extract(pfit, type = "mu.0", sample = "all"))
  mu.1 <- suppressWarnings(extract(pfit, type = "mu.1", sample = "all"))
  p.score <- extract(pfit, sample = "all", type = "p.score")
  
  for (j in seq_along(groups)) {
    m <- min(testData$y[g.sel[[j]]]); M <- max(testData$y[g.sel[[j]]])
    
    mu.hat.0 <- boundValues((boundValues(mu.0[,g.sel[[j]]], c(m, M)) - m) / (M - m), yBounds)
    mu.hat.1 <- boundValues((boundValues(mu.1[,g.sel[[j]]], c(m, M)) - m) / (M - m), yBounds)
    
    icate <- mu.hat.1 - mu.hat.0

    # replicate internal with:
    # temp <- bartCause:::getPWeightEstimates(testData$y[g.sel[[j]]], testData$z[g.sel[[j]]], NULL, "att", mu.hat.0, mu.hat.1,
    #                            extract(pfit, sample = "all", type = "p.score")[g.sel[[j]],], yBounds, p.scoreBounds)
    # f <- bartCause:::getPWeightFunction("att", NULL, icate, boundValues(p.weights[g.sel[[j]],], p.scoreBounds))
    # mean(f(testData$z[g.sel[[j]]], NULL, icate, p.score[g.sel[[j]],])) * (M - m)
    
    est.unscaled <- mean(apply((icate * boundValues(p.score[,g.sel[[j]]], p.scoreBounds)), 2L, mean)) / mean(testData$z[g.sel[[j]]])
    expect_equal(pfit.sum$est$estimate[j], est.unscaled * (M - m))
  }
  
  expect_equal(apply(p.weights, length(dim(p.weights)), mean), fitted(pfit, "p.weights", sample = "all"))
})

test_that("summary works with different styles", {
  expect_is(summary(fit, ci.style = "norm"),  "bartcFit.summary")
  expect_is(summary(fit, ci.style = "quant"), "bartcFit.summary")
  expect_is(summary(fit, ci.style = "hpd"),   "bartcFit.summary")
  expect_is(summary(fit, pate.style = "ppd"), "bartcFit.summary")
  
  oldWarn <- getOption("warn")
  if (!requireNamespace("tmle", quietly = TRUE))
    options(warn = -1)
  
  fit <- bartc(y, z, x, testData, method.trt = "glm", method.rsp = "tmle", verbose = FALSE,
               group.by = g, group.effects = TRUE ,use.ranef = FALSE,
               n.chains = 2L, n.threads = 1L, n.burn = 0L, n.samples = 7L, n.trees = 13L)
  
  options(warn = oldWarn)
  
  expect_is(summary(fit), "bartcFit.summary")
  expect_is(summary(fit, pate.style = "ppd"), "bartcFit.summary")
})

test_that("summary gives consistent answers with grouped data", {
  inGroupFit <- bartc(y, z, x, testData, estimand = "ate",
                      group.by = g, group.effects = TRUE,
                      method.trt = "bart", method.rsp = "bart", verbose = FALSE,
                      n.chains = 2L, n.threads = 1L, n.burn = 0L, n.samples = 7L, n.trees = 13L)
  

  sum.g.cate <- summary(inGroupFit, target = "cate")
  sum.g.sate <- summary(inGroupFit, target = "sate")
  sum.g.pate <- summary(inGroupFit, target = "pate")
  
  expect_equal(nrow(sum.g.cate$estimates), length(unique(testData$g)) + 1L)
  expect_equal(sum.g.cate$estimates$estimate, sum.g.sate$estimates$estimate)
  expect_equal(sum.g.cate$estimates$estimate, sum.g.pate$estimates$estimate)
  expect_true(all(sum.g.pate$estimates$sd > sum.g.sate$estimates$sd &
                  sum.g.sate$estimates$sd > sum.g.cate$estimates$sd))
  expect_equal(unname(sapply(extract(inGroupFit, "cate"), mean)),
               head(sum.g.cate$estimates$estimate, -1L))
})

