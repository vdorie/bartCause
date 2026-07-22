context("common support diagnostics")

source(system.file("common", "linearData.R", package = "bartCause"))
testData$g <- sample(3L, nrow(testData$x), replace = TRUE)

test_that("sd common support diagnostic works", {
  expect_is(bartc(y, z, x, data = testData,
                  method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
                  n.burn = 0L, n.samples = 3L, n.trees = 7L, n.chains = 1L, n.threads = 1L,
                  commonSup.rule = "sd", maxIter = 2L), "bartcFit")
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.burn = 0L, n.samples = 3L, n.trees = 7L, n.chains = 1L, n.threads = 1L,
               n.thin = 1L,
               group.by = g,
               commonSup.rule = "sd", maxIter = 2L)
  expect_is(fit, "bartcFit")
  expect_is(summary(fit), "bartcFit.summary")
  
  fit <- bartc(y, z, x, data = testData,
               method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
               n.burn = 0L, n.samples = 3L, n.trees = 7L, n.chains = 1L, n.threads = 1L,
               n.thin = 1L,
               group.by = g, group.effects = TRUE,
               commonSup.rule = "sd", maxIter = 2L)
  expect_is(fit, "bartcFit")
  expect_is(summary(fit), "bartcFit.summary")
})

test_that("chisq common support diagnostic works", {
  expect_is(bartc(y, z, x, data = testData,
                  method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
                  n.burn = 0L, n.samples = 3L, n.trees = 7L, n.chains = 1L, n.threads = 1L,
                  commonSup.rule = "chisq", maxIter = 2L), "bartcFit")
  expect_is(bartc(y, z, x, data = testData,
                  method.rsp = "p.weight", method.trt = "bart", estimand = "att", verbose = FALSE,
                  n.burn = 0L, n.samples = 3L, n.trees = 7L, n.chains = 1L, n.threads = 1L,
                  n.thin = 1L,
                  group.by = g,
                  commonSup.rule = "chisq", maxIter = 2L), "bartcFit")
})

test_that("getCommonSupportSubset validates its arguments", {
  expect_error(bartCause:::getCommonSupportSubset(1:5, 1:5, "bogus", 1, rep(1, 5), rep(FALSE, 5)),
               "commonSup.rule must be one of")
  expect_error(bartCause:::getCommonSupportSubset(1:5, 1:5, "sd", NA_real_, rep(1, 5), rep(FALSE, 5)),
               "commonSup.cut must be a real number")
  expect_error(bartCause:::getCommonSupportSubset(1:5, 1:5, "chisq", 1.5, rep(1, 5), rep(FALSE, 5)),
               "commonSup.cut must be in \\(0, 1\\)")
})

test_that("weighted p.weight estimates work for the atc estimand with grouped effects", {
  testData$w <- runif(length(testData$y), 0.5, 1.5)

  set.seed(22)
  fit <- bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "p.weight", estimand = "atc",
               weights = w, group.by = g, group.effects = TRUE, verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L)
  fit.sum <- summary(fit)

  boundValues <- bartCause:::boundValues
  yBounds <- c(.005, .995)
  p.scoreBounds <- c(0.025, 0.975)

  groups <- levels(as.factor(testData$g))
  mu.0 <- suppressWarnings(extract(fit, type = "mu.0", sample = "all"))
  mu.1 <- suppressWarnings(extract(fit, type = "mu.1", sample = "all"))
  p.score <- extract(fit, sample = "all", type = "p.score")

  for (j in seq_along(groups)) {
    sel <- which(testData$g == groups[j])
    wj <- testData$w[sel]; wj <- wj / sum(wj)
    m <- min(testData$y[sel]); M <- max(testData$y[sel])

    mu.hat.0 <- boundValues((boundValues(mu.0[,sel], c(m, M)) - m) / (M - m), yBounds)
    mu.hat.1 <- boundValues((boundValues(mu.1[,sel], c(m, M)) - m) / (M - m), yBounds)
    icate <- mu.hat.1 - mu.hat.0
    ps.sel <- boundValues(p.score[,sel], p.scoreBounds)

    num <- sapply(seq_len(nrow(icate)), function(i) sum(icate[i,] * (1 - ps.sel[i,]) * wj))
    den <- sapply(seq_len(nrow(ps.sel)), function(i) sum((1 - ps.sel[i,]) * wj))
    est.manual <- mean(num / den) * (M - m)

    expect_equal(fit.sum$estimates$estimate[j], est.manual)
  }
})

