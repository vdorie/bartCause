context("generic functions")

source(system.file("common", "groupedData.R", package = "bartCause"))

set.seed(22)
fit <- bartc(y, z, x, data = testData, method.trt = "glm", n.samples = 50L,
             group.by = g, use.ranef = FALSE, group.effects = TRUE,
             n.burn = 25L, n.chains = 4L, n.threads = 1L, verbose = FALSE)

p.score <- fitted(glm(z ~ x + g, family = stats::binomial, data = testData))

set.seed(22)
x.train <- cbind(z = testData$z, testData$x, p.score, testData$g)
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

test_that("combine chains works as expected", {
  combineChains <- bartCause:::combineChains
  
  mu.obs <- extract(fit, "mu.obs")
  expect_equal(as.vector(mu.obs), as.vector(aperm(fit$mu.hat.obs, c(2, 1, 3))))
  sigma <- extract(fit, "sigma")
  expect_equal(sigma, as.vector(t(fit$fit.rsp$sigma)))
})

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
  pfit <- bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "p.weight", estimand = "att",
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
  
  set.seed(22)
  unweighted_fit <- bartc(y, z, x, data = testData, method.trt = "glm", verbose = FALSE,
                          n.chains = 2L, n.threads = 1L, n.burn = 0L, n.samples = 7L, n.trees = 13L)


  expect_is(unweighted_summary <- summary(unweighted_fit, pate.style = "var.exp"), "bartcFit.summary")
  
  set.seed(22)
  weighted_fit <- bartc(y, z, x, data = testData, method.trt = "glm", verbose = FALSE,
                        weights = rep(1, length(testData$y)),
                        n.chains = 2L, n.threads = 1L, n.burn = 0L, n.samples = 7L, n.trees = 13L)
  expect_is(weighted_summary <- summary(weighted_fit, pate.style = "var.exp"), "bartcFit.summary")
  
  expect_equal(unweighted_summary$estimate$est, weighted_summary$estimate$est)
  expect_equal(unweighted_summary$estimate$sd, weighted_summary$estimate$sd)
  
  icates <- extract(unweighted_fit, "icate")
  n.samples <- dim(icates)[1L]
  n.obs <- dim(icates)[2L]
  
  var_tot <- var(as.vector(icates))
  var_w <- var(extract(unweighted_fit, "cate"))
  var_b <- mean(apply((icates - apply(icates, 1, mean))^2, 1, sum) / (n.obs - 1))
  
  
  expect_equal(sqrt(var_w + var_b), unweighted_summary$estimate$sd)
  expect_equal(var_tot, (var_w * (n.samples - 1) / n.samples + var_b * (n.obs - 1) / n.obs) * n.obs * n.samples / (n.obs * n.samples - 1))
})

test_that("summary works with different styles for method tmle", {
  skip_on_cran()

  oldWarn <- getOption("warn")
  if (!requireNamespace("tmle", quietly = TRUE))
    options(warn = -1)
  
  fit <- bartc(y, z, x, data = testData, method.trt = "glm", method.rsp = "tmle", verbose = FALSE,
               group.by = g, group.effects = TRUE ,use.ranef = FALSE,
               n.chains = 2L, n.threads = 1L, n.burn = 0L, n.samples = 7L, n.trees = 13L)
  
  options(warn = oldWarn)
  
  expect_is(summary(fit), "bartcFit.summary")
  expect_is(summary(fit, pate.style = "ppd"), "bartcFit.summary")
})

test_that("summary works with att/atc", {
  fit <- bartc(y, z, x, data = testData, estimand = "att",
               method.trt = "bart", method.rsp = "bart", verbose = FALSE,
               n.chains = 2L, n.threads = 1L, n.burn = 0L, n.samples = 7L, n.trees = 13L)
  expect_is(summary(fit, "pate"), "bartcFit.summary")
  expect_is(summary(fit, "sate"), "bartcFit.summary")
  expect_is(summary(fit, "cate"), "bartcFit.summary")
})

test_that("summary gives consistent answers with grouped data", {
  inGroupFit <- bartc(y, z, x, data = testData, estimand = "ate",
                      group.by = g, group.effects = TRUE,
                      method.trt = "bart", method.rsp = "bart", verbose = FALSE,
                      n.chains = 2L, n.threads = 1L, n.burn = 0L, n.samples = 7L, n.trees = 13L)
  

  sum.g.cate <- summary(inGroupFit, target = "cate")
  sum.g.sate <- summary(inGroupFit, target = "sate")
  sum.g.pate <- summary(inGroupFit, target = "pate")
  
  # test that sub group estimates are actually subgroup estimates
  samples.icate <- extract(inGroupFit, "icate", "all")
  samples.gcate <- lapply(unique(testData$g), function(j) rowMeans(samples.icate[,testData$g == j]))
  names(samples.gcate) <- unique(testData$g)
  expect_equal(sum.g.cate$estimates[names(samples.gcate),]$estimate, unname(sapply(samples.gcate, mean)))
  expect_equal(sum.g.cate$estimates[names(samples.gcate),]$sd, unname(sapply(samples.gcate, sd)))
  expect_equal(sum.g.cate$estimates["total",]$estimate, mean(samples.icate))
  expect_equal(sum.g.cate$estimates["total",]$sd, sd(rowMeans(samples.icate)))

  
  expect_equal(nrow(sum.g.cate$estimates), length(unique(testData$g)) + 1L)
  expect_true(length(unique(sum.g.sate$estimates$estimate)) > 1L)

  expect_equal(sum.g.cate$estimates$estimate, sum.g.pate$estimates$estimate)

  expect_true(all(sum.g.pate$estimates$sd > sum.g.cate$estimates$sd))

  expect_equal(unname(sapply(extract(inGroupFit, "cate"), mean)),
               head(sum.g.cate$estimates$estimate, -1L))
})

test_that("common support cutoffs are being applied consistently", {
  n.chains <- 2L
  n.samples <- 7L
  n.obs <- length(testData$y)
  fit <- bartc(y, z, x, data = testData, estimand = "ate",
               method.trt = "bart", method.rsp = "bart", verbose = FALSE,
               commonSup.rule = "sd", seed = 5,
               n.chains = n.chains, n.threads = 1L, n.burn = 0L, n.samples = n.samples,
               n.trees = 13L)
  
  sum.cate <- summary(fit, target = "cate")
  sum.sate <- summary(fit, target = "sate")
  
  icates <- extract(fit, "icate", "all")[,fit$commonSup.sub]
  y.obs <- as.vector(testData$y)
  
  
  oldSeed <- .GlobalEnv$.Random.seed
  .GlobalEnv$.Random.seed <- fit$seed
  
  mu.cf <- extract(fit, "mu.cf", combineChains = TRUE)
  iscates <- t((y.obs - t(mu.cf)) * (2 * testData$z - 1))
  iscates <- iscates[,fit$commonSup.sub]
  
  mu.cf <- aperm(array(mu.cf, c(n.samples, n.chains, n.obs)), c(2L, 1L, 3L))
  expect_equal(mu.cf, extract(fit, "mu.cf", combineChains = FALSE))
  
  sigma <- rep(extract(fit, "sigma", combineChains = FALSE), times = n.samples)
  
  epsilon <- rnorm(prod(dim(mu.cf)), 0, sigma)
  
  .GlobalEnv$.Random.seed <- oldSeed
  
  y.cf <- mu.cf + epsilon
  expect_equal(y.cf, extract(fit, "y.cf", combineChains = FALSE))
  y.cf <- matrix(aperm(y.cf, c(2L, 1L, 3L)), nrow = n.samples * n.chains)
  expect_equal(y.cf, extract(fit, "y.cf", combineChains = TRUE))
  
  ites <- t((y.obs - t(y.cf)) * ifelse(testData$z == 1, 1, -1))
  ites <- ites[,fit$commonSup.sub]
  
  expect_equal(sd(apply(ites, 1, mean)), sd(extract(fit, "sate")))
  
  expect_true(!is.nan(sum.cate$estimates$estimate) && is.finite(sum.cate$estimates$estimate))
  expect_true(!is.nan(sum.sate$estimates$estimate) && is.finite(sum.sate$estimates$estimate))
  expect_equal(sum.cate$estimates$estimate, mean(icates))
  expect_equal(sum.cate$estimates$sd, sd(apply(icates, 1, mean)))
  expect_equal(sum.sate$estimates$estimate, mean(iscates))
  expect_equal(sum.sate$estimates$sd,
               sqrt(var(rowMeans(iscates)) + mean(extract(fit, "sigma")^2) / sum(fit$commonSup.sub)))
})

source(system.file("common", "friedmanData.R", package = "bartCause"))

test_that("sate summary is correct quantity", {
  data <- generateFriedmanData(n = 100, causal = TRUE)

  fit <- bartc(y, z, x, data = data, verbose = FALSE, samples = 50L,
               n.burn = 25L, n.chains = 4L, n.threads = 1L, seed = 0)

  fit_sum <- summary(fit, "sate")
  samples.sate <- extract(fit, "sate")
  
  expect_true(abs(fit_sum$estimates$estimate - mean(samples.sate)) < 1e-1)
  expect_true(abs(fit_sum$estimates$sd - sd(samples.sate)) < 1e-1)
})

