context("bcf and bartBCF")

source(system.file("common", "linearData.R", package = "bartCause"))

linearFrame <- with(testData, data.frame(y = y, z = z, x1 = x[,1], x2 = x[,2], x3 = x[,3]))
n.obs <- nrow(linearFrame)

## The sequence bcf() drives dbarts through, written out by hand. Every
## structural test compares against this rather than against itself.
handBCFSampler <- function(frame, n.trees = 20L, n.trees.treatment = 10L,
                           n.chains = 2L, rngSeed = 11L, extraColumns = character(),
                           subset = seq_len(nrow(frame)), muVars = NULL, tauVars = NULL)
{
  rhs <- paste0(c("x1", "x2", "x3", extraColumns, "z"), collapse = " + ")
  formula <- eval(str2lang(paste0("y ~ ", rhs)))
  ## the basis covers the frame before 'subset', which dbarts restricts it with
  data <- dbarts::dbartsData(formula, data = frame, subset = subset,
                             bases = list(NULL, cbind(1 - frame$z, frame$z)))
  if (is.null(muVars))  muVars  <- setdiff(colnames(data@x), "z")
  if (is.null(tauVars)) tauVars <- setdiff(colnames(data@x), "z")
  control <- dbarts::dbartsControl(n.chains = n.chains, n.threads = 1L, n.trees = n.trees,
                                   n.burn = 5L, n.samples = 7L, verbose = FALSE,
                                   updateState = FALSE, seed = rngSeed)
  dbarts::dbarts(data, control = control, tree.prior = dbarts::dbartsPriors$cgm(2.0, 0.95),
                 forests = list(dbarts::forest(vars = muVars),
                                dbarts::forest(vars = tauVars, n.trees = n.trees.treatment,
                                               base = 0.25, power = 3, sd = 1,
                                               amplitude.prior.variance = 0.5,
                                               update.amplitude = TRUE)))
}

toBartCause <- function(x) aperm(x, c(3L, 2L, 1L))

test_that("bcf reproduces a hand-written dbarts driver bitwise (FB0)", {
  set.seed(101)
  fit <- bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z,
             n.trees = 20L, n.trees.treatment = 10L,
             n.samples = 7L, n.burn = 5L, n.chains = 2L, n.threads = 1L,
             verbose = FALSE, seed = 11L)

  set.seed(101)
  sampler <- handBCFSampler(linearFrame)
  sampler$sampleTreesFromPrior(updateState = FALSE)
  burn    <- sampler$run(0L, 5L, updateState = FALSE)
  samples <- sampler$run(0L, 7L, updateState = FALSE)

  expect_identical(fit$forests$mu,  toBartCause(samples$forestFits[,1L,,]))
  expect_identical(fit$forests$tau, toBartCause(samples$forestFits[,2L,,]))
  expect_identical(unname(fit$glue), unname(toBartCause(samples$glue)))
  expect_identical(fit$mu.hat.obs, toBartCause(samples$train))
  expect_identical(fit$sigma, t(samples$sigma))
  expect_identical(fit$first.sigma, t(burn$sigma))
  expect_identical(unname(fit$varcount$tau), unname(toBartCause(samples$varcount[,2L,,])))

  ## negative half: perturb the sequence and the draws must differ
  set.seed(101)
  perturbed <- handBCFSampler(linearFrame)
  perturbed$sampleTreesFromPrior(updateState = FALSE)
  invisible(perturbed$run(0L, 6L, updateState = FALSE))
  expect_false(identical(fit$forests$mu,
                         toBartCause(perturbed$run(0L, 7L, updateState = FALSE)$forestFits[,1L,,])))
})

test_that("the fitted surface reconstructs from mu, tau, glue and the offset (FB1)", {
  offsetFrame <- linearFrame
  set.seed(3)
  offsetFrame$off <- rnorm(n.obs)

  fit <- bcf(y ~ x1 + x2 + x3, data = offsetFrame, treatment = z, offset = off,
             n.trees = 20L, n.samples = 5L, n.burn = 3L, n.chains = 2L,
             n.threads = 1L, verbose = FALSE, seed = 4L)

  a  <- fit$glue[,,"a"]
  b.0 <- fit$glue[,,"b.0"]
  b.1 <- fit$glue[,,"b.1"]
  b.z <- b.cf <- array(0, dim(fit$forests$mu))
  for (i in seq_along(fit$trt)) {
    b.z[,,i]  <- if (fit$trt[i] == 1) b.1 else b.0
    b.cf[,,i] <- if (fit$trt[i] == 1) b.0 else b.1
  }
  combined <- fit$response.scale * (as.vector(a) * fit$forests$mu + b.z * fit$forests$tau) +
    fit$response.shift
  offsetArray <- array(rep(fit$data@offset, each = prod(dim(combined)[1:2])), dim(combined))

  expect_lt(max(abs(combined + offsetArray - fit$mu.hat.obs)), 1e-14)
  expect_gt(max(abs(combined - fit$mu.hat.obs)), 1e-2)

  ## flipping the b index moves the surface by exactly response.scale * (b_1-z - b_z) * tau
  expect_equal(fit$mu.hat.cf,
               fit$mu.hat.obs + fit$response.scale * (b.cf - b.z) * fit$forests$tau)
  treated <- fit$trt == 1
  delta <- fit$mu.hat.cf - fit$mu.hat.obs
  expect_equal(delta[,,treated],
               (fit$response.scale * as.vector(b.0 - b.1) * fit$forests$tau)[,,treated])
  expect_equal(delta[,,!treated],
               (fit$response.scale * as.vector(b.1 - b.0) * fit$forests$tau)[,,!treated])
})

test_that("bartBCF carries the documented structure", {
  set.seed(11)
  fit <- bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z,
             n.trees = 20L, n.trees.treatment = 10L,
             n.samples = 7L, n.burn = 5L, n.chains = 2L, n.threads = 1L,
             verbose = FALSE, seed = 11L)

  expect_is(fit, "bartBCF")
  expect_equal(names(fit$forests), c("mu", "tau"))
  expect_equal(names(fit$varcount), names(fit$forests))
  expect_equal(dim(fit$forests$mu), c(2L, 7L, n.obs))
  expect_equal(dim(fit$forests$tau), c(2L, 7L, n.obs))
  expect_equal(dim(fit$mu.hat.obs), c(2L, 7L, n.obs))
  expect_equal(dim(fit$mu.hat.cf), c(2L, 7L, n.obs))
  expect_equal(dim(fit$glue), c(2L, 7L, 3L))
  expect_equal(dimnames(fit$glue)[[3L]], c("a", "b.0", "b.1"))
  expect_equal(dim(fit$varcount$tau), c(2L, 7L, 4L))
  expect_equal(dimnames(fit$varcount$tau)[[3L]], c("x1", "x2", "x3", "z"))
  expect_equal(dim(fit$sigma), c(2L, 7L))
  expect_equal(dim(fit$first.sigma), c(2L, 5L))
  expect_equal(as.vector(fit$y), linearFrame$y)
  expect_equal(as.vector(fit$trt), linearFrame$z)
  expect_equal(fit$name.trt, "z")
  expect_null(fit$name.p.score)
  expect_equal(fit$family, "gaussian")
  expect_null(fit$fit)
  expect_equal(fit$n.trees, c(mu = 20L, tau = 10L))

  ## one chain drops the leading margin, but sigma stays a matrix
  set.seed(11)
  fit.1 <- bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z, n.trees = 20L,
               n.samples = 7L, n.burn = 5L, n.chains = 1L, n.threads = 1L,
               verbose = FALSE, seed = 11L)
  expect_equal(dim(fit.1$mu.hat.obs), c(7L, n.obs))
  expect_equal(dim(fit.1$glue), c(7L, 3L))
  expect_equal(dim(fit.1$sigma), c(1L, 7L))
  expect_equal(dim(fit.1$first.sigma), c(1L, 5L))

  ## a zero-burn fit has a zero-column first.sigma rather than none
  fit.0 <- bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z, n.trees = 10L,
               n.samples = 4L, n.burn = 0L, n.chains = 2L, n.threads = 1L, verbose = FALSE)
  expect_equal(dim(fit.0$first.sigma), c(2L, 0L))

  expect_output(print(fit), "Bayesian causal forest")
  expect_output(print(fit), "mu: 20 trees, tau: 10 trees")
})

test_that("the treatment and the propensity score are masked out by forest (FB2)", {
  set.seed(5)
  fit <- bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z,
             p.score = as.vector(testData$p.score),
             n.trees = 20L, n.samples = 5L, n.burn = 3L, n.chains = 2L,
             n.threads = 1L, verbose = FALSE, keepSampler = TRUE, seed = 6L)

  expect_equal(fit$name.p.score, "ps")
  expect_equal(dimnames(fit$varcount$mu)[[3L]], c("x1", "x2", "x3", "z", "ps"))

  ## the per-draw container
  expect_equal(sum(fit$varcount$tau[,,fit$name.p.score]), 0)
  expect_equal(sum(fit$varcount$tau[,,fit$name.trt]), 0)
  expect_equal(sum(fit$varcount$mu[,,fit$name.trt]), 0)
  expect_gt(sum(fit$varcount$mu[,,fit$name.p.score]), 0)
  ## and the live read, on the kept sampler; the design order is
  ## (x1, x2, x3, z, ps), so the treatment is column 4 and the score column 5
  expect_true(all(fit$fit$getForestVariableCounts(2L)[c(4L, 5L),] == 0))
  expect_true(all(fit$fit$getForestVariableCounts(1L)[4L,] == 0))
  expect_gt(sum(fit$fit$getForestVariableCounts(2L)[1:3,]), 0)

  ## negative half: the same seed with no masks puts counts on both columns
  frame <- linearFrame
  frame$ps <- as.vector(testData$p.score)
  sampler <- handBCFSampler(frame, n.trees = 20L, n.trees.treatment = 50L, rngSeed = 6L,
                            extraColumns = "ps",
                            muVars = c("x1", "x2", "x3", "ps", "z"),
                            tauVars = c("x1", "x2", "x3", "ps", "z"))
  sampler$sampleTreesFromPrior(updateState = FALSE)
  invisible(sampler$run(0L, 3L, updateState = FALSE))
  unmasked <- sampler$run(0L, 5L, updateState = FALSE)
  expect_gt(sum(unmasked$varcount[4L,2L,,]), 0)
  expect_gt(sum(unmasked$varcount[5L,2L,,]), 0)
  expect_gt(sum(unmasked$varcount[5L,1L,,]), 0)
})

test_that("every subset kind reaches the same fit, and the basis covers the full data (FB11)", {
  namedFrame <- linearFrame
  rownames(namedFrame) <- paste0("r", seq_len(n.obs))
  keep <- seq_len(60L)

  fitBy <- function(subset)
    bcf(y ~ x1 + x2 + x3, data = namedFrame, treatment = z, subset = subset,
        n.trees = 10L, n.samples = 3L, n.burn = 2L, n.chains = 1L, n.threads = 1L,
        verbose = FALSE, seed = 12L)

  fit.int <- fitBy(keep)
  expect_equal(fit.int$n.obs, 60L)
  expect_equal(as.vector(fit.int$trt), namedFrame$z[keep])

  fit.neg <- fitBy(-seq.int(61L, n.obs))
  expect_equal(fit.neg$n.obs, 60L)
  expect_identical(fit.neg$mu.hat.obs, fit.int$mu.hat.obs)

  fit.chr <- fitBy(paste0("r", keep))
  expect_equal(fit.chr$n.obs, 60L)
  expect_identical(fit.chr$mu.hat.obs, fit.int$mu.hat.obs)

  mask <- rep(FALSE, n.obs); mask[keep] <- TRUE
  fit.lgl <- fitBy(mask)
  expect_equal(fit.lgl$n.obs, 60L)
  expect_identical(fit.lgl$mu.hat.obs, fit.int$mu.hat.obs)

  ## man/bartc.Rd's canonical workflow: refit under a common support rule, then
  ## refit on the logical subset it produces
  base <- bartc(y, z, x, data = testData, method.rsp = "bcf", method.trt = "glm",
                n.trees = 10L, n.samples = 5L, n.burn = 3L, n.chains = 2L,
                n.threads = 1L, verbose = FALSE)
  cut <- refit(base, commonSup.rule = "sd")
  onSupport <- bartc(y, z, x, data = testData, subset = cut$commonSup.sub,
                     method.rsp = "bcf", method.trt = "glm",
                     n.trees = 10L, n.samples = 5L, n.burn = 3L, n.chains = 2L,
                     n.threads = 1L, verbose = FALSE)
  expect_equal(length(onSupport$trt), sum(cut$commonSup.sub))

  ## the basis contract bcf() hands its treatment basis over under: one covering
  ## the data before 'subset' is restricted to the kept rows, and one already at
  ## the kept-row count is refused loudly rather than aligned by position
  aligned <- dbarts::dbartsData(y ~ x1 + x2 + x3 + z, data = namedFrame, subset = keep,
                                bases = list(NULL, cbind(1 - namedFrame$z, namedFrame$z)))
  expect_equal(as.vector(aligned@bases[[2L]][, 2L]), as.numeric(namedFrame$z[keep]))
  expect_error(
    dbarts::dbartsData(y ~ x1 + x2 + x3 + z, data = namedFrame, subset = keep,
                       bases = list(NULL, cbind(1 - namedFrame$z[keep], namedFrame$z[keep]))),
    paste0("matching 'subset' (", length(keep), ") but not the full data (", n.obs, " rows)"),
    fixed = TRUE)
  ## a subset that loses an arm is refused by name
  expect_error(fitBy(which(namedFrame$z == 1)), "treatment arm with no observations")
  expect_error(fitBy(c("r1", "nope")), "names rows not present in the data")
})

test_that("a missing response is refused rather than dropped (FB4)", {
  missFrame <- linearFrame
  missFrame$y[seq_len(5L)] <- NA

  expect_error(
    bcf(y ~ x1 + x2 + x3, data = missFrame, treatment = z, n.trees = 5L,
        n.samples = 3L, n.burn = 2L, n.chains = 1L, verbose = FALSE),
    "cannot be fit with missing response values")

  missData <- testData
  missData$y[seq_len(5L)] <- NA
  expect_error(
    bartc(y, z, x, data = missData, method.rsp = "bcf", method.trt = "glm",
          n.samples = 3L, n.burn = 2L, n.chains = 1L, verbose = FALSE),
    "cannot fit with missing response values")

  ## the same data under 'bart' fits
  expect_is(bartc(y, z, x, data = missData, method.rsp = "bart", method.trt = "glm",
                  n.trees = 7L, n.samples = 3L, n.burn = 2L, n.chains = 1L,
                  n.threads = 1L, verbose = FALSE),
            "bartcFit")
})

test_that("predict refuses on both surfaces, naming the replay door (FB5)", {
  fit <- bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z, n.trees = 10L,
             n.samples = 3L, n.burn = 2L, n.chains = 2L, n.threads = 1L, verbose = FALSE)
  expect_error(predict(fit, linearFrame), "per-forest saved-tree replay")

  bartcBCF <- bartc(y, z, x, data = testData, method.rsp = "bcf", method.trt = "glm",
                    n.trees = 10L, n.samples = 3L, n.burn = 2L, n.chains = 2L,
                    n.threads = 1L, verbose = FALSE)
  expect_error(predict(bartcBCF, testData$x), "per-forest saved-tree replay")
})

test_that("a binary response is linked last and reports no sigma (FB7)", {
  source(system.file("common", "binaryData.R", package = "bartCause"), local = TRUE)
  binaryFrame <- with(testData, data.frame(y = y, z = z, x1 = x[,1], x2 = x[,2], x3 = x[,3]))

  fit <- bcf(y ~ x1 + x2 + x3, data = binaryFrame, treatment = z, n.trees = 20L,
             n.samples = 5L, n.burn = 3L, n.chains = 2L, n.threads = 1L, verbose = FALSE)

  expect_equal(fit$family, "probit")
  expect_null(fit$sigma)
  expect_null(fit$first.sigma)
  expect_true(all(fit$mu.hat.obs > 0 & fit$mu.hat.obs < 1))
  expect_true(all(fit$mu.hat.cf > 0 & fit$mu.hat.cf < 1))
  expect_error(extract(fit, "sigma"), "does not have a residual standard deviation")

  bartcBCF <- bartc(y, z, x, data = testData, method.rsp = "bcf", method.trt = "glm",
                    n.trees = 20L, n.samples = 5L, n.burn = 3L, n.chains = 2L,
                    n.threads = 1L, verbose = FALSE)
  expect_true(bartCause:::responseIsBinary(bartcBCF))
  expect_error(extract(bartcBCF, "sigma"), "binary response model does not have")
})

test_that("plot_sigma reads a bcf fit's per-chain sigma (FB10)", {
  fit <- bartc(y, z, x, data = testData, method.rsp = "bcf", method.trt = "glm",
               n.trees = 15L, n.samples = 5L, n.burn = 3L, n.chains = 2L,
               n.threads = 1L, verbose = FALSE)
  expect_equal(dim(fit$fit.rsp$first.sigma), c(2L, 3L))
  expect_equal(dim(fit$fit.rsp$sigma), c(2L, 5L))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot_sigma(fit))

  ## negative half: without first.sigma the call must fail
  broken <- fit
  broken$fit.rsp$first.sigma <- NULL
  expect_error(suppressWarnings(plot_sigma(broken)))

  ## third leg: a zero-burn fit still plots
  zeroBurn <- bartc(y, z, x, data = testData, method.rsp = "bcf", method.trt = "glm",
                    n.trees = 15L, n.samples = 5L, n.burn = 0L, n.chains = 2L,
                    n.threads = 1L, verbose = FALSE)
  expect_equal(dim(zeroBurn$fit.rsp$first.sigma), c(2L, 0L))
  expect_silent(plot_sigma(zeroBurn))
})

test_that("the bartBCF accessors return the documented quantities", {
  set.seed(7)
  fit <- bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z, n.trees = 15L,
             n.samples = 5L, n.burn = 3L, n.chains = 2L, n.threads = 1L,
             verbose = FALSE, seed = 8L)

  expect_equal(dim(extract(fit, "mu.obs")), c(10L, n.obs))
  expect_identical(extract(fit, "mu.obs", combineChains = FALSE), fit$mu.hat.obs)
  expect_equal(dim(extract(fit, "glue")), c(10L, 3L))
  expect_equal(length(extract(fit, "sigma")), 10L)
  expect_identical(extract(fit, "mu", combineChains = FALSE), fit$forests$mu)
  expect_identical(extract(fit, "tau", combineChains = FALSE), fit$forests$tau)

  ## icate is response.scale * (b.1 - b.0) * tau, which on a gaussian fit is
  ## exactly mu.1 - mu.0
  expect_equal(extract(fit, "icate"), extract(fit, "mu.1") - extract(fit, "mu.0"))
  ## and mu.1/mu.0 recover the observed surface where the condition was observed
  mu.1 <- extract(fit, "mu.1", combineChains = FALSE)
  expect_equal(mu.1[,,fit$trt == 1], fit$mu.hat.obs[,,fit$trt == 1])
  expect_equal(mu.1[,,fit$trt == 0], fit$mu.hat.cf[,,fit$trt == 0])

  expect_true(is.list(extract(fit, "varcount")))
  expect_equal(names(extract(fit, "varcount")), c("mu", "tau"))
  expect_identical(extract(fit, "varcount", forest = "tau", combineChains = FALSE),
                   fit$varcount$tau)
  expect_error(extract(fit, "mu", forest = "tau"), "applies only to type = 'varcount'")
  expect_error(extract(fit, "not-a-type"), "type must be in")

  expect_equal(fitted(fit, "mu.obs"), apply(fit$mu.hat.obs, 3L, mean))
  expect_equal(fitted(fit, "sigma"), mean(fit$sigma))
  expect_equal(fitted(fit, "glue"), apply(fit$glue, 3L, mean))
  expect_equal(residuals(fit), as.vector(fit$y) - fitted(fit, "mu.obs"))
  expect_warning(extract(fit, value = "mu.obs"), "called with unknown argument")
})

test_that("bcf refuses what it cannot express", {
  expect_error(bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z,
                   mu.blocks = dbarts::blocks(c("x1", "x2"))),
               "prognostic-forest block partition is not supported")
  expect_error(bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z,
                   blocks = dbarts::blocks(c("x1", "x2"))),
               "prognostic-forest block partition is not supported")
  expect_error(bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z,
                   moderators = c("x1", "z"), n.trees = 5L, n.samples = 2L,
                   n.burn = 1L, n.chains = 1L, verbose = FALSE),
               "'moderators' cannot include 'z'")
  expect_error(bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z,
                   moderators = "not-a-column", n.trees = 5L, n.samples = 2L,
                   n.burn = 1L, n.chains = 1L, verbose = FALSE),
               "names columns not in the design")
  expect_error(bcf(y ~ x1 + x2 + x3, data = linearFrame),
               "'treatment' variable must be specified")
  expect_error(bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = x1,
                   n.trees = 5L, n.samples = 2L, n.burn = 1L, n.chains = 1L,
                   verbose = FALSE),
               "must be a binary treatment coded 0 and 1")
  expect_error(bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z,
                   forests = list()),
               "is set by bcf\\(\\) itself")
})

test_that("moderators restrict the treatment forest and nothing else", {
  set.seed(9)
  fit <- bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z, moderators = "x1",
             n.trees = 15L, n.samples = 5L, n.burn = 3L, n.chains = 2L,
             n.threads = 1L, verbose = FALSE, seed = 10L)
  tauCounts <- apply(fit$varcount$tau, 3L, sum)
  expect_gt(tauCounts[["x1"]], 0)
  expect_equal(unname(tauCounts[c("x2", "x3", "z")]), c(0, 0, 0))
  muCounts <- apply(fit$varcount$mu, 3L, sum)
  expect_gt(sum(muCounts[c("x1", "x2", "x3")]), 0)
  expect_equal(muCounts[["z"]], 0)
})

test_that("bcf takes an x/y interface as bart2 does", {
  fit <- bcf(testData$x, testData$y, treatment = testData$z,
             p.score = as.vector(testData$p.score),
             n.trees = 15L, n.samples = 4L, n.burn = 3L, n.chains = 2L,
             n.threads = 1L, verbose = FALSE)
  expect_equal(dimnames(fit$varcount$mu)[[3L]], c("V1", "V2", "V3", "z", "ps"))
  expect_equal(fit$name.trt, "z")
  expect_equal(fit$name.p.score, "ps")
  expect_equal(sum(fit$varcount$tau[,,"ps"]), 0)
  expect_equal(dim(fit$mu.hat.obs), c(2L, 4L, n.obs))
})

test_that("bcf's default n.threads is capped at n.chains (dbarts dec-B115)", {
  # both defaults resolve dbarts::guessNumCores(); on a machine with more
  # cores than chains, an uncapped default warns "n.threads (N) exceeds
  # n.chains (M)" out of dbartsControl
  expect_no_warning(
    bcf(y ~ x1 + x2 + x3, data = linearFrame, treatment = z,
        n.trees = 5L, n.samples = 3L, n.burn = 2L, verbose = FALSE),
    message = "n.threads.*exceeds n.chains"
  )
})
