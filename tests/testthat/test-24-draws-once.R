context("a fit keeps its posterior draws once")

## The draws as the response fitters read them out of a dbarts fit before they
## were read out on demand: one extract call each for the training and the test
## rows, the one-chain margin dropped, and the rows with a missing response
## taken from the test rows after the counterfactual ones.
referenceDraws <- function(bartFit, combineChains, missingRows = FALSE) {
  train <- dbarts::extract(bartFit, sample = "train", combineChains = combineChains)
  test  <- dbarts::extract(bartFit, sample = "test",  combineChains = combineChains)
  if (!combineChains && bartFit$n.chains == 1L && length(dim(train)) == 3L) {
    train <- array(train, dim(train)[-1L], dimnames(train)[-1L])
    test  <- array(test,  dim(test)[-1L],  dimnames(test)[-1L])
  }
  if (!any(missingRows)) {
    dimnames(test) <- dimnames(train)
    return(list(obs = train, cf = test))
  }
  n <- length(missingRows)
  n.mis <- sum(missingRows)
  three <- length(dim(train)) > 2L
  obs <- cf <- if (three) array(0, c(dim(train)[1L:2L], n)) else matrix(0, nrow(train), n)
  pick <- function(x, rows) if (three) x[,,rows, drop = FALSE] else x[,rows, drop = FALSE]
  put <- function(x, rows, value) { if (three) x[,,rows] <- value else x[,rows] <- value; x }
  obs <- put(obs, !missingRows, train)
  obs <- put(obs,  missingRows, pick(test, n + seq_len(n.mis)))
  cf  <- pick(test, seq_len(n))
  list(obs = obs, cf = cf)
}

n <- 90L
set.seed(24)
df <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), g = factor(rep(1:3, length.out = n)))
df$z <- rbinom(n, 1, plogis(0.5 * df$x1))
df$y <- df$x1 + df$z * df$x2 + rnorm(n)
df$yb <- rbinom(n, 1, pnorm(0.3 * df$x1 + 0.5 * df$z))
df$w <- runif(n, 0.5, 1.5)
dfm <- df
dfm$y[c(4L, 19L, 41L)] <- NA
dfm$yb[c(6L, 22L)] <- NA

## the arguments after the first few are evaluated where the fit is made, as a
## call to bartc() has them, so that 'weights = w' names a column of the data
fitDraws <- function(..., data = df, response = y, method.rsp = "bart", method.trt = "bart", seed = 5L, envir = parent.frame()) {
  dots <- as.list(substitute(list(...)))[-1L]
  env <- new.env(parent = envir)
  env$dataArg <- data
  env$method.trtArg <- method.trt
  args <- list(response = substitute(response), treatment = quote(z), confounders = quote(x1 + x2 + x3),
               data = quote(dataArg), method.rsp = method.rsp, method.trt = quote(method.trtArg), verbose = FALSE,
               n.burn = 4L, n.samples = 8L, n.trees = 6L, n.chains = 2L, n.threads = 1L, seed = seed)
  args[names(dots)] <- dots
  suppressWarnings(suppressMessages(do.call(bartc, args, envir = env)))
}

drawNames <- c("mu.hat.obs", "mu.hat.cf", "samples.p.score")

## what the propensity fit's draws are, read the way bartc read them
referencePScore <- function(fit, combineChains) {
  samples <- dbarts::extract(fit$fit.trt, combineChains = combineChains)
  if (!combineChains && fit$fit.trt$n.chains == 1L && length(dim(samples)) == 3L)
    samples <- array(samples, dim(samples)[-1L], dimnames(samples)[-1L])
  samples
}

test_that("the three arrays are not held, and read as the response fit's own draws", {
  fit <- fitDraws()
  expect_false(any(drawNames %in% names(fit)))
  expect_false(any(drawNames %in% names(unclass(fit))))
  ref <- referenceDraws(fit$fit.rsp, FALSE)
  expect_identical(fit$mu.hat.obs, ref$obs)
  expect_identical(fit$mu.hat.cf, ref$cf)
  expect_identical(fit$samples.p.score, referencePScore(fit, FALSE))
  expect_identical(fit[["mu.hat.obs"]], fit$mu.hat.obs)
  expect_identical(fit[["mu.hat.cf"]], fit$mu.hat.cf)
  expect_identical(fit[["samples.p.score"]], fit$samples.p.score)
  ## the response fit's own draws are those same numbers
  expect_identical(unname(dbarts::extract(fit$fit.rsp, sample = "train", combineChains = FALSE)), unname(fit$mu.hat.obs))
})

test_that("every other name reads as a plain list reads it", {
  fit <- fitDraws()
  plain <- unclass(fit)
  for (name in names(plain)) {
    expect_identical(fit[[name]], plain[[name]])
    expect_identical(fit[[name, exact = TRUE]], plain[[name]])
  }
  expect_identical(fit$n.chain, plain$n.chain)
  expect_identical(fit$n.chain, fit$n.chains)
  expect_null(fit$fit)
  expect_null(fit[["mu"]])
  expect_null(fit$notAnElement)
  expect_identical(fit[c("trt", "p.score")], plain[c("trt", "p.score")])
  expect_identical(fit[[c("fit.rsp", "n.chains")]], plain[[c("fit.rsp", "n.chains")]])
  expect_identical(fit[["trt"]], fit$trt)
  expect_false(any(drawNames %in% names(with(fit, list(trt = trt)))))
  expect_false(any(grepl("mu.hat", capture.output(str(fit, max.level = 1L)))))
  ## an assignment still works as on a list
  fit$extra <- 1
  expect_identical(fit$extra, 1)
})

test_that("the arrays are those the response fitters made, for every case", {
  checkCase <- function(fit, combineChains, missingRows = FALSE, ...) {
    ref <- referenceDraws(fit$fit.rsp, combineChains, missingRows)
    expect_identical(fit$mu.hat.obs, ref$obs, ...)
    expect_identical(fit$mu.hat.cf, ref$cf, ...)
  }
  ## continuous and binary, with and without a propensity model
  checkCase(fitDraws(), FALSE, info = "continuous")
  checkCase(fitDraws(response = yb), FALSE, info = "binary")
  checkCase(fitDraws(n.chains = 1L), FALSE, info = "one chain")
  checkCase(fitDraws(response = yb, n.chains = 1L), FALSE, info = "one chain, binary")
  fit <- fitDraws(combineChains = TRUE)
  checkCase(fit, TRUE, info = "combined")
  expect_identical(length(dim(fit$mu.hat.obs)), 2L)
  expect_identical(fit$samples.p.score, referencePScore(fit, TRUE))
  checkCase(fitDraws(combineChains = TRUE, n.chains = 1L), TRUE, info = "combined, one chain")
  checkCase(fitDraws(method.trt = "glm"), FALSE, info = "glm propensity")
  checkCase(fitDraws(method.trt = "none"), FALSE, info = "no propensity")
  ## common support, weights and a subset change which rows are averaged, not the arrays
  checkCase(fitDraws(commonSup.rule = "sd", commonSup.cut = 0.5), FALSE, info = "common support")
  checkCase(fitDraws(weights = w), FALSE, info = "weights")
  checkCase(fitDraws(weights.rsp = w), FALSE, info = "response weights")
  checkCase(fitDraws(subset = x1 > -1), FALSE, info = "subset")
  checkCase(fitDraws(estimand = "att"), FALSE, info = "att")
  ## the propensity fit's own arguments
  fit <- fitDraws(args.trt = list(combineChains = TRUE))
  expect_identical(fit$samples.p.score, referencePScore(fit, TRUE))
  expect_identical(length(dim(fit$mu.hat.obs)), 3L)
  ## grouped fits
  fit <- fitDraws(group.by = g, use.ranef = FALSE, group.effects = TRUE)
  checkCase(fit, FALSE, info = "group.by")
})

test_that("a missing response has its row's expected outcome from the test rows", {
  for (case in list(list(), list(n.chains = 1L), list(combineChains = TRUE))) {
    combine <- isTRUE(case$combineChains)
    fit <- do.call(fitDraws, c(case, list(data = dfm)))
    ref <- referenceDraws(fit$fit.rsp, combine, fit$missingRows)
    expect_identical(fit$mu.hat.obs, ref$obs)
    expect_identical(fit$mu.hat.cf, ref$cf)
    expect_identical(dim(fit$mu.hat.obs)[length(dim(fit$mu.hat.obs))], n)
  }
  fit <- fitDraws(data = dfm, response = yb)
  ref <- referenceDraws(fit$fit.rsp, FALSE, fit$missingRows)
  expect_identical(fit$mu.hat.obs, ref$obs)
  expect_identical(fit$mu.hat.cf, ref$cf)
})

test_that("the propensity-weighted routes and the bcf route read the same way", {
  fit <- fitDraws(method.rsp = "p.weight")
  ref <- referenceDraws(fit$fit.rsp, FALSE)
  expect_identical(fit$mu.hat.obs, ref$obs)
  expect_identical(fit$mu.hat.cf, ref$cf)
  expect_identical(fit$samples.p.score, referencePScore(fit, FALSE))
  fit <- fitDraws(method.rsp = "p.weight", data = dfm)
  ref <- referenceDraws(fit$fit.rsp, FALSE, fit$missingRows)
  expect_identical(fit$mu.hat.obs, ref$obs)
  
  fit <- fitDraws(method.rsp = "bcf")
  expect_false(any(drawNames %in% names(fit)))
  expect_identical(fit$mu.hat.obs, fit$fit.rsp$mu.hat.obs)
  expect_identical(fit$mu.hat.cf, fit$fit.rsp$mu.hat.cf)
  fit <- fitDraws(method.rsp = "bcf", response = yb, n.chains = 1L)
  expect_identical(fit$mu.hat.obs, fit$fit.rsp$mu.hat.obs)
  
  skip_if_not_installed("tmle")
  fit <- fitDraws(method.rsp = "tmle", posteriorOfTMLE = FALSE)
  ref <- referenceDraws(fit$fit.rsp, FALSE)
  expect_identical(fit$mu.hat.obs, ref$obs)
  expect_identical(fit$samples.p.score, referencePScore(fit, FALSE))
})

test_that("a propensity matrix that was given is held, as no fit has it", {
  scores <- matrix(rep(seq(0.2, 0.8, length.out = n), 16L), n, 16L)
  fit <- fitDraws(method.trt = scores)
  expect_identical(names(fit)[names(fit) %in% drawNames], "samples.p.score")
  expect_identical(fit$samples.p.score, scores)
  expect_identical(fit[["samples.p.score"]], scores)
  expect_false(any(c("mu.hat.obs", "mu.hat.cf") %in% names(fit)))
  ref <- referenceDraws(fit$fit.rsp, FALSE)
  expect_identical(fit$mu.hat.obs, ref$obs)
  ## a propensity vector has no draws
  fit <- fitDraws(method.trt = seq(0.2, 0.8, length.out = n))
  expect_null(fit$samples.p.score)
  expect_null(fit[["samples.p.score"]])
})

test_that("a stan4bart fit's draws are laid out as bartc laid them out", {
  ## the sampler's draws are observations by samples by chains, and combined
  ## observations by (samples * chains)
  fake <- structure(list(), class = "stan4bartFit")
  train <- array(seq_len(5 * 6 * 2), c(5L, 6L, 2L))
  test  <- train + 1000
  mockExtract <- function(object, type, sample = "train", combine_chains = TRUE, ...) {
    x <- if (sample == "test") test else train
    if (combine_chains) matrix(x, dim(x)[1L], prod(dim(x)[2L:3L])) else x
  }
  testthat::local_mocked_bindings(extract = mockExtract, .package = "bartCause")
  object <- structure(list(fit.rsp = fake, fit.trt = NULL, combineChains = c(rsp = FALSE, trt = NA)), class = "bartcFit")
  draws <- bartCause:::bartcDraws(object, "mu")
  expect_identical(draws$mu.hat.obs, aperm(train, c(3L, 2L, 1L)))
  expect_identical(draws$mu.hat.cf, aperm(test, c(3L, 2L, 1L)))
  object$combineChains <- c(rsp = TRUE, trt = NA)
  draws <- bartCause:::bartcDraws(object, "mu")
  expect_identical(draws$mu.hat.obs, t(matrix(train, 5L, 12L)))
  expect_identical(draws$mu.hat.cf, t(matrix(test, 5L, 12L)))
  expect_identical(object$mu.hat.cf, draws$mu.hat.cf)
})

test_that("a fit saved by an earlier version, holding the arrays, reads them as stored", {
  fit <- fitDraws()
  held <- unclass(fit)
  held$mu.hat.obs <- array(1, c(2L, 3L, 4L))
  held$mu.hat.cf  <- array(2, c(2L, 3L, 4L))
  held$samples.p.score <- matrix(3, 2L, 2L)
  held <- structure(held, class = "bartcFit")
  expect_true(all(drawNames %in% names(held)))
  expect_identical(held$mu.hat.obs, array(1, c(2L, 3L, 4L)))
  expect_identical(held[["mu.hat.cf"]], array(2, c(2L, 3L, 4L)))
  expect_identical(held$samples.p.score, matrix(3, 2L, 2L))
  ## and one that holds them but no record of the layout
  held$combineChains <- NULL
  expect_identical(held$mu.hat.obs, array(1, c(2L, 3L, 4L)))
  ## a held array is what the readers use
  expect_identical(extract(held, "mu.obs", combineChains = FALSE, sample = "all"), array(1, c(2L, 3L, 4L)))
  
  ## a fit that holds the arrays the fit would compute answers the readers as the fit that does not
  stored <- unclass(fit)
  for (name in drawNames) stored[[name]] <- fit[[name]]
  stored <- structure(stored, class = "bartcFit")
  for (type in c("mu.obs", "mu.cf", "mu.0", "mu.1", "icate", "p.score", "p.weights", "sigma", "cate"))
    expect_identical(extract(stored, type, combineChains = FALSE), extract(fit, type, combineChains = FALSE), info = type)
})

test_that("fitted, extract, refit, summary and the plots read the same from a fit that holds the arrays", {
  for (arguments in list(list(), list(response = quote(yb), n.chains = 1L), list(data = dfm), list(method.rsp = "p.weight"),
                         list(method.rsp = "bcf"), list(method.trt = "glm", commonSup.rule = "sd", commonSup.cut = 0.5),
                         list(combineChains = TRUE))) {
    fit <- suppressWarnings(do.call(fitDraws, arguments))
    stored <- unclass(fit)
    for (name in drawNames) stored[name] <- list(fit[[name]])
    stored <- structure(stored, class = "bartcFit")
    label <- paste(deparse(arguments), collapse = "")
    
    expect_identical(fitted(fit), fitted(stored), info = label)
    types <- c("mu.obs", "mu.cf", "icate", "p.score", "ev", "y.cf")
    types <- if (fit$method.rsp %in% c("p.weight", "tmle")) c("mu.obs", "mu.cf", "p.score", "pate") else
      c("mu.obs", "mu.cf", "mu.0", "mu.1", "icate", "p.score", "p.weights", "cate", "pate", "sate", "y.0", "y.1", "ite", "y.cf")
    for (type in types) {
      set.seed(11); a <- extract(fit, type, combineChains = FALSE)
      set.seed(11); b <- extract(stored, type, combineChains = FALSE)
      expect_identical(a, b, info = paste(label, type))
    }
    if (!isTRUE(arguments$combineChains)) {
      expect_identical(summary(fit)$estimates, summary(stored)$estimates, info = label)
      expect_identical(unclass(summary(fit))[c("n.obs", "n.samples")], unclass(summary(stored))[c("n.obs", "n.samples")], info = label)
    }
    
    refitted <- suppressWarnings(refit(fit, commonSup.rule = "sd", commonSup.cut = 0.5))
    refittedStored <- suppressWarnings(refit(stored, commonSup.rule = "sd", commonSup.cut = 0.5))
    expect_identical(refitted$est, refittedStored$est, info = label)
    expect_identical(refitted$commonSup.sub, refittedStored$commonSup.sub, info = label)
    expect_false(any(drawNames %in% names(refitted)))
    
    if (!isTRUE(arguments$combineChains)) {
      pdf(NULL)
      on.exit(dev.off(), add = TRUE)
      if (!is.null(fit$fit.rsp$sigma) || fit$method.rsp == "bcf") {
        if (is.null(fit$fit.rsp[["sigma"]])) NULL else expect_silent(plot_sigma(fit))
      }
      expect_silent(plot_indiv(fit))
    }
  }
})

test_that("each read-out of the draws is made once for each call", {
  counter <- new.env()
  counter$n <- 0L
  readOut <- bartCause:::bartResponseDraws
  testthat::local_mocked_bindings(
    bartResponseDraws = function(...) { counter$n <- counter$n + 1L; readOut(...) },
    .package = "bartCause")
  
  fit <- fitDraws()
  counter$n <- 0L
  for (type in c("icate", "mu.0", "mu.1", "cate", "pate", "sate", "y.1", "ite", "y.cf"))
    { counter$n <- 0L; extract(fit, type); expect_identical(counter$n, 1L, info = type) }
  counter$n <- 0L; extract(fit, "p.score"); expect_identical(counter$n, 0L)
  counter$n <- 0L; extract(fit, "sigma"); expect_identical(counter$n, 0L)
  counter$n <- 0L; fitted(fit); expect_identical(counter$n, 1L)
  counter$n <- 0L; summary(fit); expect_identical(counter$n, 1L)
  counter$n <- 0L; refit(fit, commonSup.rule = "sd", commonSup.cut = 0.5); expect_identical(counter$n, 1L)
  
  fit <- fitDraws(method.rsp = "p.weight")
  counter$n <- 0L; suppressWarnings(refit(fit, commonSup.rule = "sd", commonSup.cut = 0.5)); expect_identical(counter$n, 1L)
  
  ## a held sigma reads the draws' shape from the draws at hand
  fit <- fitDraws(args.rsp = list(family = dbarts::dbartsFamilies$gaussian(sigma = dbarts::dbartsPriors$fixed(0.5))))
  counter$n <- 0L; extract(fit, "y.1"); expect_identical(counter$n, 1L)
})

test_that("a fit saves smaller than one that holds the arrays too", {
  ## made where the test's own objects are out of reach of the fit's environments
  fit <- fitDraws(n.samples = 100L, n.chains = 4L, envir = globalenv())
  held <- unclass(fit)
  for (name in drawNames) held[name] <- list(fit[[name]])
  held <- structure(held, class = "bartcFit")
  
  pathNew <- tempfile(fileext = ".rds")
  pathHeld <- tempfile(fileext = ".rds")
  on.exit(unlink(c(pathNew, pathHeld)), add = TRUE)
  saveRDS(fit, pathNew, compress = FALSE)
  saveRDS(held, pathHeld, compress = FALSE)
  ## the six arrays are most of the file
  expect_lt(file.size(pathNew), 0.7 * file.size(pathHeld))
  
  read <- readRDS(pathNew)
  expect_identical(read$mu.hat.obs, fit$mu.hat.obs)
  expect_identical(read$samples.p.score, fit$samples.p.score)
})
