
## tmle::tmle is replaced by a stub that returns the weighted mean of Y as the estimate and
## the weighted mean of A as the standard error, with the obsWeights it was handed. Its answer
## then identifies which rows and which weights reached the tmle call, without the tmle package.
stubTMLE <- function(Y, A, W, Q, g1W, obsWeights = NULL, ...) {
  w <- if (is.null(obsWeights)) rep(1, length(Y)) else obsWeights
  if (length(w) != length(Y)) stop("obsWeights of length ", length(w), " for ", length(Y), " rows")
  est <- list(psi = sum(w * Y) / sum(w), var.psi = (sum(w * A) / sum(w))^2)
  list(estimates = list(ATE = est, ATT = est, ATC = est))
}
environment(stubTMLE) <- baseenv()
## the same with no obsWeights argument, as tmle before 2.0.0
stubTMLE1 <- function(Y, A, W, Q, g1W) stubTMLE(Y, A, W, Q, g1W)
environment(stubTMLE1) <- list2env(list(stubTMLE = stubTMLE), parent = baseenv())

source(system.file("common", "groupedData.R", package = "bartCause"))
wdata <- data.frame(y = as.vector(testData$y), z = testData$z, x1 = testData$x[,1], x2 = testData$x[,2],
                    x3 = testData$x[,3], grp = factor(testData$g))
n <- nrow(wdata)
set.seed(101)
wdata$w <- runif(n, 0.5, 1.5)
keep <- rep(c(TRUE, TRUE, FALSE), length.out = n)

fitStub <- function(stub, data = wdata, estimand = "att", weighted = TRUE, post = TRUE, n.threads = 1L,
                    cs = "none", grouped = FALSE, subset = NULL) {
  args <- list(quote(y), quote(z), quote(x1 + x2 + x3), data = quote(data), method.trt = "bart", method.rsp = "tmle",
               estimand = estimand, commonSup.rule = cs, verbose = FALSE, n.burn = 3L, n.samples = 4L,
               n.trees = 7L, n.chains = 2L, n.threads = n.threads, posteriorOfTMLE = post)
  if (weighted) args$weights <- quote(w)
  if (cs == "sd") args$commonSup.cut <- -0.5
  if (grouped) args <- c(args, list(group.by = quote(grp), group.effects = TRUE, use.ranef = FALSE))
  if (!is.null(subset)) args$subset <- quote(subset)
  local_mocked_bindings(getTMLEFunction = function(...) stub, .package = "bartCause", .env = parent.frame())
  set.seed(8)
  suppressMessages(suppressWarnings(do.call(bartc, args)))
}

## what the stub should have returned for a set of rows of the data, over every draw
expectRows <- function(est, rows, data = wdata, weighted = TRUE, label = NULL) {
  w <- if (weighted) data$w[rows] else rep(1, length(rows))
  want <- c(weighted.mean(data$y[rows], w), weighted.mean(data$z[rows], w))
  got <- if (is.null(dim(est))) rbind(est[c("est", "se")])
         else if (length(dim(est)) == 3L) cbind(as.vector(est[,,"est"]), as.vector(est[,,"se"]))
         else est[, c("est", "se"), drop = FALSE]
  expect_true(all(abs(sweep(got, 2L, want)) < 1e-10), label = label)
}
## the rows a fit's estimate covers: the kept, observed rows that common support leaves in
rowsOf <- function(fit, base, group = NULL, data = wdata) {
  rows <- base[fit$commonSup.sub & !is.na(data$y[base])]
  if (!is.null(group)) rows <- rows[data$grp[rows] == group]
  rows
}

test_that("weights reach tmle for the rows of each response fit", {
  all <- seq_len(n)
  for (estimand in c("att", "ate")) for (post in c(TRUE, FALSE)) {
    fit <- fitStub(stubTMLE, estimand = estimand, post = post)
    expectRows(fit$est, all, label = paste("pooled", estimand, post))
  }

  ## common support drops rows
  for (post in c(TRUE, FALSE)) {
    fit <- fitStub(stubTMLE, post = post, cs = "sd")
    expect_true(any(!fit$commonSup.sub))
    expectRows(fit$est, rowsOf(fit, all), label = paste("support", post))
    ## refit() hands tmle the same weights
    rf <- suppressWarnings(refit(fit, commonSup.rule = "sd", commonSup.cut = -0.5))
    expectRows(rf$est, rowsOf(rf, all), label = paste("refit", post))
    expectRows(suppressWarnings(refit(fit, commonSup.rule = "none"))$est, all, label = paste("refit none", post))
  }

  ## group effects
  for (post in c(TRUE, FALSE)) {
    fit <- fitStub(stubTMLE, post = post, grouped = TRUE)
    for (g in levels(wdata$grp))
      expectRows(fit$est[[g]], rowsOf(fit, all, g), label = paste("group", g, post))
  }
})

test_that("weights reach tmle for the rows of subset", {
  base <- which(keep)
  for (post in c(TRUE, FALSE)) {
    fit <- fitStub(stubTMLE, post = post, subset = keep)
    expectRows(fit$est, base, label = paste("subset", post))
  }
  fit <- fitStub(stubTMLE, grouped = TRUE, subset = keep)
  for (g in levels(wdata$grp)) expectRows(fit$est[[g]], rowsOf(fit, base, g), label = paste("subset group", g))
})

test_that("tmle is handed the response and weights of the rows that are not missing", {
  mdata <- wdata
  mdata$y[c(3L, 17L, 40L, 77L)] <- NA
  all <- seq_len(n)
  for (post in c(TRUE, FALSE)) {
    fit <- fitStub(stubTMLE, data = mdata, post = post)
    expectRows(fit$est, rowsOf(fit, all, data = mdata), data = mdata, label = paste("missing", post))
  }
  fit <- fitStub(stubTMLE, data = mdata, grouped = TRUE)
  for (g in levels(wdata$grp))
    expectRows(fit$est[[g]], rowsOf(fit, all, g, data = mdata), data = mdata, label = paste("missing group", g))
})

test_that("tmle workers are handed the weights", {
  skip_on_cran()
  fit <- fitStub(stubTMLE, n.threads = 2L)
  expectRows(fit$est, seq_len(n), label = "threads")
  fit <- fitStub(stubTMLE, n.threads = 2L, grouped = TRUE)
  for (g in levels(wdata$grp)) expectRows(fit$est[[g]], rowsOf(fit, seq_len(n), g), label = paste("threads group", g))
})

test_that("tmle is called without obsWeights when there are no weights", {
  all <- seq_len(n)
  for (post in c(TRUE, FALSE)) {
    fit <- fitStub(stubTMLE1, weighted = FALSE, post = post)
    expectRows(fit$est, all, weighted = FALSE, label = paste("pooled", post))
    fit <- fitStub(stubTMLE1, weighted = FALSE, post = post, grouped = TRUE)
    for (g in levels(wdata$grp))
      expectRows(fit$est[[g]], rowsOf(fit, all, g), weighted = FALSE, label = paste("group", g, post))
  }
})

test_that("weighted tmle stops before any model is fit when the tmle package is too old", {
  getTMLEFunction <- bartCause:::getTMLEFunction
  local_mocked_bindings(
    getTMLEFunction = function(weighted, ...) getTMLEFunction(weighted, version = package_version("1.5.0")),
    ## errors if the fit gets as far as the treatment model
    getBartTreatmentFit = function(...) stop("the treatment model was fit"),
    .package = "bartCause")
  expect_error(bartc(y, z, x1 + x2 + x3, data = wdata, weights = w, method.trt = "bart", method.rsp = "tmle",
                     n.burn = 1L, n.samples = 1L, n.trees = 1L, n.chains = 1L, n.threads = 1L),
               "version 2.0.0 or later")
})

test_that("weights and responses are paired with their rows in p.weight fits", {
  ## the rows of subset, as if they were the data
  args <- list(quote(y), quote(z), quote(x1 + x2 + x3), method.trt = "bart", method.rsp = "p.weight",
               verbose = FALSE, n.burn = 3L, n.samples = 5L, n.trees = 7L, n.chains = 2L, n.threads = 1L)
  wkeep <- wdata[keep,]
  fit <- function(...) suppressWarnings(suppressMessages(do.call(bartc, c(args, list(...)))))
  for (estimand in c("att", "ate")) {
    set.seed(61)
    fit.subset <- fit(data = quote(wdata), weights = quote(w), subset = quote(keep), estimand = estimand)
    set.seed(61)
    fit.rows <- fit(data = quote(wkeep), weights = quote(w), estimand = estimand)
    expect_equal(fit.subset$est, fit.rows$est, label = estimand)
  }

  ## missing responses leave a defined estimate and standard error
  mdata <- wdata
  mdata$y[c(3L, 17L, 40L, 77L)] <- NA
  for (weighted in c(FALSE, TRUE)) {
    set.seed(62)
    fit.missing <- if (weighted) fit(data = quote(mdata), estimand = "att", weights = quote(w))
                   else fit(data = quote(mdata), estimand = "att")
    expect_true(all(is.finite(fit.missing$est)), label = paste("missing", weighted))
  }

  ## the standard error is the estimator's on the observed responses of the rows that are kept
  for (weighted in c(FALSE, TRUE)) {
    set.seed(63)
    pfit <- if (weighted) fit(data = quote(mdata), estimand = "att", weights = quote(w))
            else fit(data = quote(mdata), estimand = "att")
    mu.hat.0 <- suppressWarnings(aperm(extract(pfit, "mu.0", sample = "all", combineChains = FALSE), c(3L, 1L, 2L)))
    mu.hat.1 <- suppressWarnings(aperm(extract(pfit, "mu.1", sample = "all", combineChains = FALSE), c(3L, 1L, 2L)))
    p.score <- aperm(pfit$samples.p.score, c(3L, 1L, 2L))
    expect_equal(dim(mu.hat.0)[1L], n)
    manual <- bartCause:::getPWeightEstimates(mdata$y, pfit$trt, if (weighted) mdata$w, "att", mu.hat.0, mu.hat.1, p.score,
                                              c(.005, .995), c(0.025, 0.975))
    expect_equal(pfit$est, manual, label = paste("manual", weighted))
  }
})

test_that("the built-in tmle estimator runs on draws and on a single set of means", {
  local_mocked_bindings(getTMLEFunction = function(weighted, ...) NULL, .package = "bartCause")
  set.seed(41)
  mu.hat.0 <- rnorm(n, 0, 0.3)
  mu.hat.1 <- mu.hat.0 + 0.5 + rnorm(n, 0, 0.1)
  p.score <- runif(n, 0.2, 0.8)
  getEst <- function(estimand, mu.hat.0, mu.hat.1, p.score)
    suppressWarnings(bartCause:::getTMLEEstimates(wdata$y, wdata$z, NULL, estimand, mu.hat.0, mu.hat.1, p.score,
                                                  c(.005, .995), c(0.025, 0.975), 0.001, 20L, n.threads = 1L))
  ## recorded from the estimator before weights were routed to the tmle package
  expect_equal(getEst("att", mu.hat.0, mu.hat.1, p.score), c(0.14705762506, 0.07409820776), tolerance = 1e-8, check.attributes = FALSE)
  expect_equal(getEst("atc", mu.hat.0, mu.hat.1, p.score), c(0.31851007328, 0.04977879461), tolerance = 1e-8, check.attributes = FALSE)
  expect_equal(getEst("ate", mu.hat.0, mu.hat.1, p.score), c(0.14439369418, 0.02480104011), tolerance = 1e-8, check.attributes = FALSE)
  ## one set of means is one draw
  est <- getEst("att", cbind(mu.hat.0), cbind(mu.hat.1), cbind(p.score))
  expect_equal(unname(est[1L,]), c(0.14705762506, 0.07409820776), tolerance = 1e-8)
})

test_that("the p.weight standard error divides each draw by its own sum of scores", {
  set.seed(51)
  n.draws <- 6L
  mu.hat.0 <- matrix(rnorm(n * n.draws, 0, 0.3), n, n.draws)
  mu.hat.1 <- mu.hat.0 + 0.5 + matrix(rnorm(n * n.draws, 0, 0.1), n, n.draws)
  p.score <- matrix(runif(n * n.draws, 0.2, 0.8), n, n.draws)
  getEst <- function(weights, estimand, cols = seq_len(n.draws))
    bartCause:::getPWeightEstimates(wdata$y, wdata$z, weights, estimand, mu.hat.0[,cols,drop=FALSE],
                                    mu.hat.1[,cols,drop=FALSE], p.score[,cols,drop=FALSE], c(.005, .995), c(0.025, 0.975))
  for (estimand in c("att", "atc", "ate")) {
    unweighted <- getEst(NULL, estimand)
    ## weights that are all one are no weights
    expect_equal(getEst(rep(1, n), estimand), unweighted, label = estimand)
    for (weights in list(NULL, wdata$w)) {
      ## and a draw's standard error does not depend on the other draws
      for (i in seq_len(n.draws))
        expect_equal(getEst(weights, estimand)[i,], getEst(weights, estimand, i)[1L,], label = estimand)
      expect_equal(unname(getEst(weights, estimand, c(1L, 1L, 1L))[,"se"]),
                   rep(unname(getEst(weights, estimand, 1L)[,"se"]), 3L), label = estimand)
    }
  }
})

test_that("refit reproduces att and atc estimates", {
  for (estimand in c("att", "atc")) {
    set.seed(22)
    pfit <- bartc(y, z, x1 + x2 + x3, data = wdata, method.trt = "bart", method.rsp = "p.weight", estimand = estimand,
                  verbose = FALSE, n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L)
    expect_equal(suppressWarnings(refit(pfit))$est, pfit$est, label = estimand)
  }

  fit <- fitStub(stubTMLE, estimand = "atc")
  expectRows(suppressWarnings(refit(fit))$est, seq_len(n), label = "tmle atc")
})

test_that("the p.weight standard error under missing responses is the influence curve's over the observed rows", {
  set.seed(71)
  n.draws <- 4L
  mu.hat.0 <- matrix(rnorm(n * n.draws, 0, 0.3), n, n.draws)
  mu.hat.1 <- mu.hat.0 + 0.5
  p.score <- matrix(runif(n * n.draws, 0.2, 0.8), n, n.draws)
  y <- wdata$y; y[c(3L, 17L, 40L, 77L)] <- NA
  obs <- !is.na(y)
  got <- bartCause:::getPWeightEstimates(y, wdata$z, NULL, "ate", mu.hat.0, mu.hat.1, p.score, NULL, NULL)
  ## by hand, unweighted ate: bounded and scaled as the estimator does
  m <- min(y, na.rm = TRUE); M <- max(y, na.rm = TRUE)
  r <- c(m - 0.1 * abs(m), M + 0.1 * abs(M))
  y.st <- (pmin(pmax(y, r[1L]), r[2L]) - m) / (M - m)
  for (i in seq_len(n.draws)) {
    q0 <- (pmin(pmax(mu.hat.0[,i], m), M) - m) / (M - m); q1 <- (pmin(pmax(mu.hat.1[,i], m), M) - m) / (M - m)
    e <- p.score[,i]; z <- wdata$z
    psi <- mean(q1 - q0)
    ic <- (z * (1 - (1 - e) / e) + (1 - z) * (1 - e / (1 - e))) * (y.st - ifelse(z == 1, q1, q0)) + (q1 - q0 - psi)
    expect_equal(unname(got[i, "se"]), sd(ic[obs]) / sqrt(sum(obs)))
  }
})

test_that("weights are the rows of subset when the variables are not in a data frame", {
  y <- wdata$y; z <- wdata$z; x <- as.matrix(wdata[, c("x1", "x2", "x3")]); w <- wdata$w
  yk <- y[keep]; zk <- z[keep]; xk <- x[keep,]; wk <- w[keep]
  args <- list(method.trt = "glm", method.rsp = "p.weight", estimand = "att", verbose = FALSE,
               n.burn = 3L, n.samples = 5L, n.trees = 7L, n.chains = 2L, n.threads = 1L)
  set.seed(61)
  fit.subset <- suppressWarnings(suppressMessages(do.call(bartc, c(list(quote(y), quote(z), quote(x), weights = quote(w), subset = quote(keep)), args))))
  set.seed(61)
  fit.rows <- suppressWarnings(suppressMessages(do.call(bartc, c(list(quote(yk), quote(zk), quote(xk), weights = quote(wk)), args))))
  expect_equal(fit.subset$est, fit.rows$est)
})

test_that("a fit with the built-in estimator and posteriorOfTMLE = FALSE carries one named estimate", {
  local_mocked_bindings(getTMLEFunction = function(weighted, ...) NULL, .package = "bartCause")
  set.seed(8)
  fit <- suppressWarnings(suppressMessages(bartc(y, z, x1 + x2 + x3, data = wdata, method.trt = "glm", method.rsp = "tmle", estimand = "att",
    posteriorOfTMLE = FALSE, verbose = FALSE, n.burn = 3L, n.samples = 4L, n.trees = 7L, n.chains = 2L, n.threads = 1L)))
  expect_null(dim(fit$est))
  expect_named(fit$est, c("est", "se"))
  expect_true(all(is.finite(fit$est)))
  expect_true(is.finite(suppressWarnings(summary(fit))$estimates$estimate))
})

test_that("refit of an att fit with method.rsp = bart reproduces its summary", {
  set.seed(23)
  fit <- bartc(y, z, x1 + x2 + x3, data = wdata, method.trt = "glm", method.rsp = "bart", estimand = "att",
               verbose = FALSE, n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L)
  refitted <- refit(fit)
  set.seed(5); s1 <- summary(fit, target = "cate")$estimates
  set.seed(5); s2 <- summary(refitted, target = "cate")$estimates
  expect_equal(s2, s1)
})

test_that("an unweighted tmle fit never asks for the weighted tmle function", {
  local_mocked_bindings(getTMLEFunction = function(weighted, ...) if (weighted) stop("asked for the weighted function") else stubTMLE1,
                        .package = "bartCause")
  expect_error(suppressWarnings(suppressMessages(bartc(y, z, x1 + x2 + x3, data = wdata, method.trt = "glm", method.rsp = "tmle",
    verbose = FALSE, n.burn = 3L, n.samples = 2L, n.trees = 7L, n.chains = 1L, n.threads = 1L))), NA)
})

test_that("grouped fits pair the response with its rows when some are missing", {
  mdata <- wdata
  mdata$y[c(3L, 17L, 40L, 77L)] <- NA
  all <- seq_len(n)
  fit <- fitStub(stubTMLE, data = mdata, grouped = TRUE, post = FALSE)
  for (g in levels(wdata$grp))
    expectRows(fit$est[[g]], rowsOf(fit, all, g, data = mdata), data = mdata, label = paste("missing group nopost", g))
  ## p.weight: each group's estimate has a defined standard error, and that of the group's own rows
  set.seed(64)
  pfit <- suppressWarnings(suppressMessages(bartc(y, z, x1 + x2 + x3, data = mdata, method.trt = "bart", method.rsp = "p.weight", estimand = "att",
    group.by = grp, group.effects = TRUE, use.ranef = FALSE, verbose = FALSE, n.burn = 3L, n.samples = 5L, n.trees = 7L, n.chains = 2L, n.threads = 1L)))
  mu.hat.0 <- suppressWarnings(aperm(extract(pfit, "mu.0", sample = "all", combineChains = FALSE), c(3L, 1L, 2L)))
  mu.hat.1 <- suppressWarnings(aperm(extract(pfit, "mu.1", sample = "all", combineChains = FALSE), c(3L, 1L, 2L)))
  p.score <- aperm(pfit$samples.p.score, c(3L, 1L, 2L))
  for (g in levels(mdata$grp)) {
    rows <- mdata$grp == g
    manual <- bartCause:::getPWeightEstimates(mdata$y[rows], pfit$trt[rows], NULL, "att", mu.hat.0[rows,,,drop = FALSE], mu.hat.1[rows,,,drop = FALSE],
                                              p.score[rows,,,drop = FALSE], c(.005, .995), c(0.025, 0.975))
    expect_equal(pfit$est[[g]], manual, label = paste("manual group", g))
  }
})

test_that("a seeded fit on two tmle workers is reproducible", {
  skip_on_cran()
  drawTMLE <- function(Y, A, W, Q, g1W, obsWeights = NULL, ...) {
    est <- list(psi = stats::runif(1L), var.psi = 1)
    list(estimates = list(ATE = est, ATT = est, ATC = est))
  }
  environment(drawTMLE) <- baseenv()
  expect_equal(fitStub(drawTMLE, n.threads = 2L)$est, fitStub(drawTMLE, n.threads = 2L)$est)
})
