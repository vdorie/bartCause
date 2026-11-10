context("population weights and the response fit's weights")

set.seed(23)
n <- 120L
pwData <- data.frame(x1 = rnorm(n), x2 = runif(n))
pwData$z <- rbinom(n, 1L, plogis(pwData$x1))
pwData$y <- pwData$x1 + pwData$z * (1 + 2 * (pwData$x2 > 0.5)) + rnorm(n)
pwData$w <- ifelse(pwData$x2 > 0.5, 4, 1) * runif(n, 0.8, 1.2)
pwData$cnt <- rep(1:3, length.out = n)

## a bartc call on pwData with the arguments as written, evaluated in the
## caller's frame as a user's call is: the weights expressions are read
## within data first
pwFit <- function(..., data) {
  call <- match.call()
  call <- as.call(c(list(quote(bartCause::bartc), quote(y), quote(z), quote(x1 + x2)), as.list(call)[-1L]))
  if (is.null(call[["data"]])) call[["data"]] <- quote(pwData)
  defaults <- list(verbose = FALSE, seed = 7, n.samples = 20L, n.burn = 10L, n.trees = 10L, n.chains = 2L,
                   n.threads = 1L)
  for (name in names(defaults)) if (is.null(call[[name]])) call[[name]] <- defaults[[name]]
  suppressMessages(eval(call, parent.frame()))
}

## the weighted mean of the individual effects over 'rows', per draw
weightedCATE <- function(fit, w, rows = rep(TRUE, length(fit$trt))) {
  icate <- extract(fit, "icate", sample = "all")
  if (is.null(w)) w <- rep(1, ncol(icate))
  as.vector(icate[, rows, drop = FALSE] %*% (w[rows] / sum(w[rows])))
}

test_that("population weights set the average, and the response model is fit without them", {
  fit <- pwFit(weights = w)
  expect_length(fit$data.rsp@weights, 0L)
  expect_equal(fit$weights, pwData$w)
  expect_null(fit$weights.rsp)
  ## extract() lays its draws out chain by chain within a sample
  expect_equal(sort(extract(fit, "cate")), sort(weightedCATE(fit, pwData$w)))
  expect_equal(fitted(fit, "cate"), mean(weightedCATE(fit, pwData$w)))
  ## weights.rsp = NULL is the same call
  fitNull <- pwFit(weights = w, weights.rsp = NULL)
  expect_equal(fitNull$mu.hat.cf, fit$mu.hat.cf)
  expect_equal(summary(fitNull)$estimates, summary(fit)$estimates)
  ## the weights' variables, or their logarithm, are confounders written by the caller
  expect_is(suppressMessages(bartc(y, z, x1 + x2 + log(w), data = pwData, weights = w, verbose = FALSE,
                                   n.samples = 5L, n.burn = 5L, n.trees = 5L, n.chains = 1L, n.threads = 1L)),
            "bartcFit")
})

test_that("weights.rsp alone weights the response fit and leaves the average equal", {
  fit <- pwFit(weights.rsp = w)
  expect_equal(fit$data.rsp@weights, pwData$w)
  expect_null(fit$weights)
  expect_equal(fit$weights.rsp, pwData$w)
  expect_equal(fitted(fit, "cate"), mean(weightedCATE(fit, NULL)))
  expect_null(summary(fit)$n.eff)
})

test_that("both vectors: the fit takes weights.rsp whole, the average the population weights", {
  p <- rep(c(1, 2), length.out = n)
  fit <- pwFit(weights = w, weights.rsp = w * p)
  expect_equal(fit$data.rsp@weights, pwData$w * p)
  expect_equal(fitted(fit, "cate"), mean(weightedCATE(fit, pwData$w)))
  ## counts of identical rows, read as counts by the fit and as shares by the average
  fit <- pwFit(weights = cnt, weights.rsp = cnt, weights.type = "frequency")
  expect_identical(fit$fit.rsp$weights.type, "frequency")
  expect_equal(fitted(fit, "cate"), mean(weightedCATE(fit, pwData$cnt)))
  ## the reading may ride args.rsp
  fit <- pwFit(weights = cnt, weights.rsp = cnt, args.rsp = list(weights.type = "frequency"))
  expect_identical(fit$fit.rsp$weights.type, "frequency")
})

test_that("weights are read within data, then in the caller's frame, and cut to subset", {
  keep <- pwData$x1 > -1
  fit <- pwFit(weights = w, subset = keep)
  expect_equal(fit$weights, pwData$w[keep])
  expect_equal(fitted(fit, "cate"), mean(weightedCATE(fit, pwData$w[keep])))

  callerWeights <- pwData$w * 2
  fit <- pwFit(weights = callerWeights, weights.rsp = callerWeights, subset = keep)
  expect_equal(fit$weights, callerWeights[keep])
  expect_equal(fit$data.rsp@weights, callerWeights[keep])

  ## no data: vectors
  fit <- suppressMessages(bartc(pwData$y, pwData$z, as.matrix(pwData[c("x1", "x2")]), weights = pwData$w,
                                weights.rsp = pwData$cnt, subset = keep, verbose = FALSE, seed = 7,
                                n.samples = 10L, n.burn = 5L, n.trees = 5L, n.chains = 1L, n.threads = 1L))
  expect_equal(fit$weights, pwData$w[keep])
  expect_equal(fit$weights.rsp, pwData$cnt[keep])
  expect_equal(fit$data.rsp@weights, pwData$cnt[keep])
})

test_that("a weights-only call keeps today's estimand weights", {
  ## before the split the response fit held the weights and the averages read
  ## them back from it; with no missing response those were the same numbers
  fit <- pwFit(weights = w)
  pre <- pwFit(weights = w, weights.rsp = w)
  expect_equal(fit$weights, pre$data.rsp@weights)
  old <- pre
  old[c("weights", "weights.rsp")] <- NULL
  expect_equal(fitted(old, "cate"), fitted(pre, "cate"))
})

test_that("a fit saved before the weights were stored reads them from its response fit", {
  pre <- pwFit(weights = w, weights.rsp = w)
  old <- pre
  old[c("weights", "weights.rsp")] <- NULL
  for (target in c("pate", "sate", "cate"))
    expect_equal(summary(old, target = target)$estimates, summary(pre, target = target)$estimates, label = target)
  expect_equal(extract(old, "y.cf"), extract(pre, "y.cf"))
})

test_that("NULL weights are absent, as written or as a variable's value", {
  unweighted <- pwFit()
  noWeights <- NULL
  for (fit in list(pwFit(weights = NULL), pwFit(weights = noWeights), pwFit(weights = NULL, weights.rsp = NULL),
                   pwFit(weights = NULL, method.trt = "glm"))) {
    expect_null(fit$weights)
    expect_null(fit$weights.rsp)
    expect_equal(length(fit$data.rsp@weights), 0L)
  }
  expect_equal(extract(pwFit(weights = NULL), "pate"), extract(unweighted, "pate"))
})

test_that("population weights are refused by name", {
  ## nothing is fit before the refusal
  local_mocked_bindings(getBartTreatmentFit = function(...) stop("the treatment model was fit"), .package = "bartCause")
  refuse <- function(..., pattern) {
    call <- match.call()
    call[["pattern"]] <- NULL
    call <- as.call(c(list(quote(bartCause::bartc), quote(y), quote(z), quote(x1 + x2), data = quote(pwData)),
                      as.list(call)[-1L], list(verbose = FALSE)))
    expect_error(eval(call, parent.frame()), pattern)
  }
  refuse(weights = w[-1L], pattern = "'weights' must have one value for each row \\(120\\), but its length is 119")
  refuse(weights = as.character(w), pattern = "'weights' must be a numeric vector.*weights\\(design\\)")
  refuse(weights = ifelse(x1 > 1.5, NA, w), pattern = "'weights' is missing at rows? [0-9]")
  refuse(weights = ifelse(x1 > 1.5, -1, w), pattern = "'weights' must be finite and non-negative")
  refuse(weights = ifelse(x1 > 1.5, Inf, w), pattern = "'weights' must be finite and non-negative")
  refuse(weights = 0 * w, pattern = "'weights' must be positive for at least one of the rows modeled")
  refuse(weights.rsp = w[-1L], pattern = "'weights.rsp' must have one value for each row")
  ## a missing weight outside of the subset is not read
  expect_error(bartc(y, z, x1 + x2, data = pwData, weights = ifelse(x1 > 1.5, NA, w), subset = x1 <= 1.5,
                     verbose = FALSE), "the treatment model was fit")

  ## the response fit's reading and weights come by weights.rsp
  refuse(weights = w, weights.type = "frequency", pattern = "weights.type says how the response fit reads weights.rsp")
  refuse(weights = w, args.rsp = list(weights.type = "precision"), pattern = "weights.type says how")
  refuse(weights = w, args.rsp = list(weights = w), pattern = "args.rsp gives 'weights'")
  ## a 0/1 response's probit has no weighted form
  binary <- pwData
  binary$y <- as.numeric(binary$y > 1)
  expect_error(bartc(y, z, x1 + x2, data = binary, weights.rsp = w, verbose = FALSE),
               "a 0/1 response is fit by probit BART, which takes weights.rsp of 0 or 1 only")
})

test_that("a row out of a continuous response fit cannot be in the average", {
  out <- rep(c(0, 1, 1), length.out = n)
  expect_error(pwFit(weights.rsp = out), "have weights.rsp 0 but a positive population weight")
  ## refused before any model is fit, in the singular for one row
  local({
    local_mocked_bindings(getBartTreatmentFit = function(...) stop("the treatment model was fit"),
                          getBartResponseFit  = function(...) stop("the response model was fit"), .package = "bartCause")
    expect_error(pwFit(weights.rsp = ifelse(seq_len(n) == 7L, 0, cnt)),
                 "^row 7 has weights.rsp 0 but a positive population weight")
    expect_error(pwFit(weights = w, weights.rsp = out), "rows 1, 4, 7, 10, 13 and 35 more have weights.rsp 0")
  })
  expect_error(pwFit(weights = w, weights.rsp = out), "rows 1, 4, 7, 10, 13 and 35 more have weights.rsp 0")
  ## with a population weight of 0 the row informs nothing and is drawn as NA
  fit <- pwFit(weights = w * out, weights.rsp = out)
  y.cf <- extract(fit, "y.cf")
  expect_true(all(is.na(y.cf[, out == 0])))
  expect_false(anyNA(y.cf[, out == 1]))
  expect_true(is.finite(fitted(fit, "pate")))
  expect_true(all(is.finite(unlist(summary(fit, target = "sate")$estimates))))
  ## the reverse informs the fit and not the average
  fit <- pwFit(weights = w * out, weights.rsp = w)
  expect_equal(fitted(fit, "cate"), mean(weightedCATE(fit, pwData$w * out)))
  ## a binary response's 0/1 weights are membership, its counterfactual drawn from p
  binary <- pwData
  binary$y <- as.numeric(binary$y > 1)
  fit <- pwFit(weights.rsp = out, data = binary)
  expect_false(anyNA(extract(fit, "y.cf")))
})

test_that("the ATT and ATC average over their own rows' weights", {
  ## the synthesis's probe: weights 5 on x > 0, the weighted fit of bartCause main
  set.seed(1)
  m <- 300
  x <- rnorm(m); z <- rbinom(m, 1, plogis(x)); y <- x + z * (1 + 2 * (x > 0)) + rnorm(m)
  w <- ifelse(x > 0, 5, 1)
  d <- data.frame(y, z, x, w)
  f <- suppressMessages(bartc(y, z, x, data = d, weights = w, weights.rsp = w, estimand = "att", n.samples = 200,
                              n.burn = 100, n.chains = 2, n.threads = 1, verbose = FALSE, seed = 1))
  ic <- extract(f, "icate", sample = "all")
  hand <- mean(ic[, z == 1] %*% (w[z == 1] / sum(w[z == 1])))
  expect_equal(mean(extract(f, "cate")), hand)
  expect_equal(fitted(f, "cate"), hand)
  expect_equal(summary(f, target = "cate")$estimates$estimate, hand)
  expect_equal(summary(f, target = "pate")$estimates$estimate, hand)
  ## 3.324 before, the weights normalized over every row
  expect_equal(hand, 2.734, tolerance = 5e-4)

  fit <- pwFit(weights = w, estimand = "atc")
  expect_equal(fitted(fit, "cate"), mean(weightedCATE(fit, pwData$w, pwData$z == 0)))
})

test_that("summary's conditional average reads the weights", {
  fit <- pwFit(weights = w)
  cate <- weightedCATE(fit, pwData$w)
  est <- summary(fit, target = "cate")$estimates
  expect_equal(est$estimate, mean(cate))
  expect_equal(est$sd, sd(cate))
  expect_false(isTRUE(all.equal(est$estimate, mean(weightedCATE(fit, NULL)))))
})

test_that("common support cuts the weights with the rows", {
  fit <- pwFit(weights = w, estimand = "att", commonSup.rule = "sd", commonSup.cut = -1)
  expect_true(any(!fit$commonSup.sub & pwData$z == 1))
  rows <- fit$commonSup.sub & pwData$z == 1
  cate <- weightedCATE(fit, pwData$w, rows)
  expect_warning(est <- summary(fit, target = "cate")$estimates, NA)
  expect_equal(est$estimate, mean(cate))
  expect_equal(fitted(fit, "cate"), mean(cate))
  ## the summary names the share of the treated rows' weight the cut removes
  treated <- pwData$z == 1
  share <- sum(pwData$w[treated & !fit$commonSup.sub]) / sum(pwData$w[treated])
  expect_equal(summary(fit)$weight.cut, share)
  expect_output(print(summary(fit)), "of the population weight; the estimate describes")
  ## refit reaches the same
  refitted <- refit(pwFit(weights = w, estimand = "att"), commonSup.rule = "sd", commonSup.cut = -1)
  expect_equal(fitted(refitted, "cate"), mean(cate))
})

test_that("missing responses keep their weights and rows", {
  mdata <- pwData
  mdata$y[c(3L, 50L, 77L)] <- NA
  fit <- pwFit(weights = w, data = mdata)
  expect_equal(fit$weights, mdata$w)
  expect_equal(fitted(fit, "cate"), mean(weightedCATE(fit, mdata$w)))
  for (target in c("pate", "sate", "cate"))
    expect_warning(expect_true(all(is.finite(unlist(summary(fit, target = target)$estimates)))), NA)
  ## a missing response is drawn from its posterior predictive
  ite <- extract(fit, "ite")
  expect_equal(dim(ite), c(40L, n))
  expect_false(anyNA(ite))
  expect_false(anyNA(extract(fit, "sate")))
  y.1 <- extract(fit, "y.1", sample = "all")
  observed <- !is.na(mdata$y) & mdata$z == 1
  expect_equal(y.1[1L, observed], mdata$y[observed])
})

test_that("a one-chain weighted fit extracts", {
  fit <- suppressMessages(bartc(y, z, x1 + x2, data = pwData, weights = w, verbose = FALSE, seed = 7,
                                n.samples = 20L, n.burn = 10L, n.trees = 10L, n.chains = 1L, n.threads = 1L))
  expect_equal(dim(extract(fit, "pate", combineChains = FALSE)), c(1L, 20L))
  expect_equal(extract(fit, "cate"), weightedCATE(fit, pwData$w))
})

test_that("group estimates read the weights, and their total each group's share of them", {
  grp <- factor(rep(c("a", "b", "c"), c(20L, 40L, 60L)))
  gdata <- cbind(pwData, grp = grp)
  fit <- pwFit(weights = w, data = gdata, group.by = grp, group.effects = TRUE, use.ranef = FALSE)
  for (target in c("cate", "pate")) {
    est <- summary(fit, target = target)$estimates
    for (level in levels(grp))
      expect_equal(est[level, "estimate"], mean(weightedCATE(fit, gdata$w, grp == level)), label = paste(target, level))
    expect_equal(est["total", "estimate"], mean(weightedCATE(fit, gdata$w)), label = target)
  }
  ## interval ends combine the groups' draws by weight share; extract() lays
  ## its draws out chain by chain within a sample
  est <- summary(fit, target = "cate", ci.style = "quant")$estimates
  expect_equal(unlist(est["total", c("ci.lower", "ci.upper")]),
               setNames(quantile(weightedCATE(fit, gdata$w), c(0.025, 0.975)), c("ci.lower", "ci.upper")))
})

test_that("a group of no population weight is NA, the others and the total estimated", {
  grp <- factor(rep(c("a", "b", "c"), c(20L, 40L, 60L)))
  gdata <- cbind(pwData, grp = grp)
  gdata$w0 <- ifelse(grp == "a", 0, gdata$w)
  fit <- pwFit(weights = w0, data = gdata, group.by = grp, group.effects = TRUE, use.ranef = FALSE)
  expect_true(all(is.na(extract(fit, "cate")[["a"]])))
  expect_true(is.na(fitted(fit, "cate")[["a"]]))
  for (target in c("pate", "sate", "cate")) {
    est <- summary(fit, target = target)$estimates
    expect_true(all(is.na(est["a", c("estimate", "sd", "ci.lower", "ci.upper")])), label = target)
    expect_equal(est["a", "n"], 20)
    expect_true(all(is.finite(unlist(est[c("b", "c", "total"), ]))), label = target)
  }
  est <- summary(fit, target = "cate", ci.style = "quant")$estimates
  expect_equal(est["b", "estimate"], mean(weightedCATE(fit, gdata$w0, grp == "b")))
  expect_equal(est["total", "estimate"], mean(weightedCATE(fit, gdata$w0)))
  expect_equal(unlist(est["total", c("ci.lower", "ci.upper")]),
               setNames(quantile(weightedCATE(fit, gdata$w0), c(0.025, 0.975)), c("ci.lower", "ci.upper")))
  ## the estimand's whole weight 0 stays an error
  gdata$wt <- ifelse(gdata$z == 1, 0, gdata$w)
  fit <- pwFit(weights = wt, data = gdata, estimand = "att", group.by = grp, group.effects = TRUE, use.ranef = FALSE)
  expect_error(summary(fit), "the population weights are 0 on every row the estimate averages over")
})

test_that("a glm treatment model is fit at the weights' mean of 1", {
  fit1 <- pwFit(weights = w, method.trt = "glm")
  expect_warning(fit2 <- pwFit(weights = w * 1e4, method.trt = "glm"), NA)
  expect_equal(coef(fit2$fit.trt), coef(fit1$fit.trt))
  expect_equal(fit1$fit.trt$prior.weights, setNames(pwData$w / mean(pwData$w), rownames(pwData)))
  ## a BART treatment model is fit without them, and says so only where the
  ## estimator reads the scores
  small <- function(...) bartc(y, z, x1 + x2, data = pwData, weights = w, ..., verbose = FALSE,
                               n.samples = 5L, n.burn = 5L, n.trees = 5L, n.chains = 1L, n.threads = 1L)
  expect_message(small(method.rsp = "p.weight"), "propensity score model is fit unweighted")
  expect_message(small(), NA)
})

test_that("no estimate or interval depends on the population weights' scale", {
  fit1 <- pwFit(weights = w)
  fit2 <- pwFit(weights = w * 1000)
  for (target in c("pate", "sate", "cate"))
    expect_equal(summary(fit2, target = target)$estimates, summary(fit1, target = target)$estimates, label = target)
  expect_equal(summary(fit2)$n.eff, summary(fit1)$n.eff)
  fit1 <- pwFit(weights = w, weights.rsp = cnt)
  fit2 <- pwFit(weights = w * 1000, weights.rsp = cnt)
  for (target in c("pate", "sate", "cate"))
    expect_equal(summary(fit2, target = target)$estimates, summary(fit1, target = target)$estimates, label = target)
})

test_that("interval noise terms follow the population weights and the response fit", {
  sigma2 <- function(fit) mean(extract(fit, "sigma")^2)
  ## unweighted fit: Kish's sum of squared shares
  fit <- pwFit(weights = w)
  a <- pwData$w / sum(pwData$w)
  cate <- weightedCATE(fit, pwData$w)
  expect_equal(summary(fit, target = "pate")$estimates$sd, sqrt(var(cate) + 2 * sigma2(fit) * sum(a^2)))
  mu.cf <- extract(fit, "mu.cf", "all")
  scate <- as.vector((rep(pwData$y, each = nrow(mu.cf)) - mu.cf) %*% (a * ifelse(pwData$z == 1, 1, -1)))
  expect_equal(summary(fit, target = "sate")$estimates$sd, sqrt(var(scate) + sigma2(fit) * sum(a^2)))
  ## precision weights: each row's noise at sigma^2 / w
  fit <- pwFit(weights = w, weights.rsp = w)
  cate <- weightedCATE(fit, pwData$w)
  expect_equal(summary(fit, target = "pate")$estimates$sd, sqrt(var(cate) + 2 * sigma2(fit) * sum(a^2 / pwData$w)))
  ## no weights: 1 / n, as before
  fit <- pwFit()
  expect_equal(summary(fit, target = "pate")$estimates$sd,
               sqrt(var(weightedCATE(fit, NULL)) + 2 * sigma2(fit) / n))
  ## a 0/1 response: p (1 - p) in place of sigma^2
  binary <- pwData
  binary$y <- as.numeric(binary$y > 1)
  fit <- pwFit(weights = w, data = binary)
  p.obs <- extract(fit, "mu.obs", "all")
  p.cf  <- extract(fit, "mu.cf", "all")
  noise <- mean((p.obs * (1 - p.obs) + p.cf * (1 - p.cf)) %*% a^2)
  expect_equal(summary(fit, target = "pate")$estimates$sd, sqrt(var(weightedCATE(fit, pwData$w)) + noise))
})

test_that("under missing responses only the SATE adds the noise of the observed draw", {
  sigma2 <- function(fit) mean(extract(fit, "sigma")^2)
  mdata <- pwData
  miss <- c(3L, 18L, 50L, 77L, 101L)
  mdata$y[miss] <- NA
  missing <- is.na(mdata$y)
  sgn <- ifelse(mdata$z == 1, 1, -1)
  sateDraws <- function(fit, a) {
    mu.cf <- extract(fit, "mu.cf", "all")
    iscate <- (rep(mdata$y, each = nrow(mu.cf)) - mu.cf) * rep(sgn, each = nrow(mu.cf))
    iscate[, missing] <- extract(fit, "icate", "all")[, missing]
    as.vector(iscate %*% a)
  }
  ## no weights: the PATE's term is 1 / n, as in 1.0-10, and the SATE's adds
  ## m / n^2 for the m missing rows
  fit <- pwFit(data = mdata)
  a <- rep(1 / n, n)
  expect_equal(summary(fit, target = "pate")$estimates$sd,
               sqrt(var(weightedCATE(fit, NULL)) + 2 * sigma2(fit) / n))
  expect_equal(summary(fit, target = "sate")$estimates$sd,
               sqrt(var(sateDraws(fit, a)) + sigma2(fit) * (1 / n + sum(missing) / n^2)))
  ## population weights and precision weights.rsp
  fit <- pwFit(weights = w, weights.rsp = cnt, data = mdata)
  a <- mdata$w / sum(mdata$w)
  expect_equal(summary(fit, target = "pate")$estimates$sd,
               sqrt(var(weightedCATE(fit, mdata$w)) + 2 * sigma2(fit) * sum(a^2 / mdata$cnt)))
  expect_equal(summary(fit, target = "sate")$estimates$sd,
               sqrt(var(sateDraws(fit, a)) + sigma2(fit) * (sum(a^2 / mdata$cnt) + sum(a[missing]^2 / mdata$cnt[missing]))))
})

test_that("a counterfactual is drawn with its noise as the response fit has it", {
  noiseSD <- function(fit) apply(extract(fit, "y.cf") - extract(fit, "mu.cf"), 2L, sd) / sqrt(mean(extract(fit, "sigma")^2))
  big <- c(1, 16)[1L + (seq_len(n) > n / 2)]
  bigData <- cbind(pwData, big = big)
  fit <- suppressMessages(bartc(y, z, x1 + x2, data = bigData, weights.rsp = big, verbose = FALSE, seed = 7,
                                n.samples = 400L, n.burn = 10L, n.trees = 10L, n.chains = 1L, n.threads = 1L))
  ratio <- mean(noiseSD(fit)[big == 1]) / mean(noiseSD(fit)[big == 16])
  expect_equal(ratio, 4, tolerance = 0.1)
  ## counts: one unit at sigma
  fit <- suppressMessages(bartc(y, z, x1 + x2, data = bigData, weights.rsp = big, weights.type = "frequency",
                                verbose = FALSE, seed = 7, n.samples = 400L, n.burn = 10L, n.trees = 10L,
                                n.chains = 1L, n.threads = 1L))
  ratio <- mean(noiseSD(fit)[big == 1]) / mean(noiseSD(fit)[big == 16])
  expect_equal(ratio, 1, tolerance = 0.1)
  ## new rows are drawn at sigma, as before
  fit <- suppressMessages(bartc(y, z, x1 + x2, data = bigData, weights.rsp = big, method.trt = "none",
                                keepTrees = TRUE, verbose = FALSE, seed = 7, n.samples = 400L, n.burn = 10L,
                                n.trees = 10L, n.chains = 1L, n.threads = 1L))
  newData <- bigData[big == 16, ]
  noise <- predict(fit, newData, type = "y") - predict(fit, newData, type = "mu")
  expect_equal(mean(apply(noise, 2L, sd)) / sqrt(mean(extract(fit, "sigma")^2)), 1, tolerance = 0.1)
})

test_that("summary prints the weights' effective sample size", {
  fit <- pwFit(weights = w)
  expect_equal(summary(fit)$n.eff, sum(pwData$w)^2 / sum(pwData$w^2))
  ## over the rows the estimand averages, those a common support rule cut out
  fit <- pwFit(weights = w, estimand = "att", commonSup.rule = "sd", commonSup.cut = -1)
  rows <- pwData$z == 1 & fit$commonSup.sub
  expect_true(any(pwData$z == 1 & !fit$commonSup.sub))
  expect_equal(summary(fit)$n.eff, sum(pwData$w[rows])^2 / sum(pwData$w[rows]^2))
  expect_output(print(summary(fit)), "Effective sample size of the population weights over the rows averaged")
  expect_failure(expect_output(print(summary(pwFit())), "Effective sample size"))
})

test_that("p.weight and tmle average with the population weights and fit with weights.rsp", {
  ## the propensity-weighted estimate is the scores-and-weights weighted mean
  check <- function(fit, w) {
    mu.0 <- suppressWarnings(aperm(extract(fit, "mu.0", sample = "all", combineChains = FALSE), c(3L, 1L, 2L)))
    mu.1 <- suppressWarnings(aperm(extract(fit, "mu.1", sample = "all", combineChains = FALSE), c(3L, 1L, 2L)))
    p.score <- aperm(fit$samples.p.score, c(3L, 1L, 2L))
    manual <- bartCause:::getPWeightEstimates(pwData$y, fit$trt, w, "att", mu.0, mu.1, p.score, c(.005, .995), c(0.025, 0.975))
    expect_equal(fit$est, manual)
  }
  fit <- pwFit(weights = w, method.rsp = "p.weight", estimand = "att")
  expect_length(fit$data.rsp@weights, 0L)
  check(fit, pwData$w)
  fit <- pwFit(weights.rsp = w, method.rsp = "p.weight", estimand = "att")
  expect_equal(fit$data.rsp@weights, pwData$w)
  check(fit, NULL)
  fit <- pwFit(weights = w, weights.rsp = cnt, method.rsp = "p.weight", estimand = "att")
  expect_equal(fit$data.rsp@weights, pwData$cnt)
  check(fit, pwData$w)

  ## tmle is handed the population weights as obsWeights
  handed <- NULL
  stub <- function(Y, A, W, Q, g1W, obsWeights = NULL, ...) {
    handed <<- obsWeights
    est <- list(psi = 0, var.psi = 1)
    list(estimates = list(ATE = est, ATT = est, ATC = est))
  }
  local_mocked_bindings(getTMLEFunction = function(...) stub, .package = "bartCause")
  tfit <- pwFit(weights = w, weights.rsp = cnt, method.rsp = "tmle", estimand = "att")
  expect_equal(handed, pwData$w)
  expect_equal(tfit$data.rsp@weights, pwData$cnt)
  handed <- NULL
  pwFit(weights.rsp = cnt, method.rsp = "tmle", estimand = "att")
  expect_null(handed)
})

test_that("bcf takes weights.rsp and averages with the population weights", {
  fit <- suppressMessages(bartc(y, z, x1 + x2, data = pwData, weights = w, weights.rsp = cnt, method.rsp = "bcf",
                                verbose = FALSE, seed = 7, n.samples = 10L, n.burn = 5L, n.trees = 5L,
                                n.chains = 1L, n.threads = 1L))
  expect_equal(fit$data.rsp@weights, pwData$cnt)
  expect_equal(fitted(fit, "cate"), mean(weightedCATE(fit, pwData$w)))
  fit <- suppressMessages(bartc(y, z, x1 + x2, data = pwData, weights.rsp = cnt, weights.type = "frequency",
                                method.rsp = "bcf", verbose = FALSE, seed = 7, n.samples = 10L, n.burn = 5L,
                                n.trees = 5L, n.chains = 1L, n.threads = 1L))
  expect_identical(fit$fit.rsp$weights.type, "frequency")
})
