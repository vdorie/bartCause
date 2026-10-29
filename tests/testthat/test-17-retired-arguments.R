context("arguments of bartCause 1.0-10 that dbarts 1.0-0 spells another way")

source(system.file("common", "linearData.R", package = "bartCause"), local = TRUE)

x <- testData$x
y <- testData$y
z <- testData$z
y.binary <- as.integer(y > median(y))

## the call is built with the arguments' values in place: an argument given
## through a wrapper's '...' would reach bartc as '..1', which bartc's own
## argument handling does not resolve; a list for args.rsp or args.trt, and an
## expression for tree.prior or family, are passed quoted
smallFit <- function(response = y, ...)
  eval(as.call(c(list(quote(bartCause::bartc), response, quote(z), quote(x),
                      n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 6L,
                      n.trees = 10L, keepTrees = TRUE, verbose = FALSE),
                 list(...))))

treePrior <- function(fit) {
  model <- fit$fit$model
  c(power = model@tree.prior@power, base = model@tree.prior@base)
}
residPrior <- function(fit) {
  model <- fit$fit$model
  c(df = model@resid.prior@df, quantile = model@resid.prior@quantile)
}

test_that("power, base, sigdf and sigquant reach the response fit as 1.0-10 took them", {
  default <- smallFit(method.trt = "none")
  expect_equal(treePrior(default$fit.rsp), c(power = 2, base = 0.95))
  expect_equal(residPrior(default$fit.rsp), c(df = 3, quantile = 0.9))
  
  fit <- smallFit(method.trt = "none", power = 1.5, base = 0.7, sigdf = 7, sigquant = 0.5)
  expect_equal(treePrior(fit$fit.rsp), c(power = 1.5, base = 0.7))
  expect_equal(residPrior(fit$fit.rsp), c(df = 7, quantile = 0.5))
  
  fit <- smallFit(method.trt = "none",
                  args.rsp = quote(list(power = 1.5, base = 0.7, sigdf = 7, sigquant = 0.5)))
  expect_equal(treePrior(fit$fit.rsp), c(power = 1.5, base = 0.7))
  expect_equal(residPrior(fit$fit.rsp), c(df = 7, quantile = 0.5))
  
  ## each is taken alone, the other at 0.9-x's default
  expect_equal(treePrior(smallFit(method.trt = "none", base = 0.7)$fit.rsp), c(power = 2, base = 0.7))
  expect_equal(residPrior(smallFit(method.trt = "none", sigdf = 7)$fit.rsp), c(df = 7, quantile = 0.9))
  expect_equal(residPrior(smallFit(method.trt = "none", sigquant = 0.5)$fit.rsp), c(df = 3, quantile = 0.5))
  
  ## and through crossvalidation, which refits the same call with xbart
  fit <- smallFit(method.trt = "none", crossvalidate = TRUE, power = 1.5, sigdf = 7)
  expect_equal(treePrior(fit$fit.rsp)[["power"]], 1.5)
  expect_equal(residPrior(fit$fit.rsp)[["df"]], 7)
})

test_that("power and base reach the propensity fit through args.trt", {
  fit <- smallFit(args.trt = quote(list(power = 1.2, base = 0.5)), power = 1.5, base = 0.7)
  expect_equal(treePrior(fit$fit.trt), c(power = 1.2, base = 0.5))
  expect_equal(treePrior(fit$fit.rsp), c(power = 1.5, base = 0.7))
  
  ## from the dots, both fits take them
  fit <- smallFit(power = 1.5, base = 0.7)
  expect_equal(treePrior(fit$fit.trt), c(power = 1.5, base = 0.7))
})

test_that("a retired name and its successor together are refused, naming both", {
  expect_error(smallFit(method.trt = "none", tree.prior = quote(cgm(2, 0.9)), power = 1.5),
               "'power' and 'tree.prior' cannot both be given")
  expect_error(smallFit(method.trt = "none", tree.prior = quote(cgm(2, 0.9)), power = 1.5, base = 0.7),
               "'power', 'base' and 'tree.prior' cannot both be given")
  expect_error(smallFit(method.trt = "none", family = quote(gaussian(sigma = fixed(1))), sigdf = 4),
               "'sigdf' and 'family' cannot both be given")
  expect_error(smallFit(method.trt = "none", seed = 3L, rngSeed = 4L),
               "'rngSeed' and 'seed' cannot both be given")
})

test_that("rngSeed is the seed of the fit", {
  draws <- function(...) smallFit(method.trt = "none", ...)$mu.hat.obs
  expect_equal(draws(rngSeed = 3L), draws(rngSeed = 3L))
  expect_false(isTRUE(all.equal(draws(rngSeed = 3L), draws(rngSeed = 4L))))
})

test_that("an argument that does not reach the sampler is named in one warning", {
  expect_warning(smallFit(method.trt = "none", foo = 1, bar = 2),
                 "arguments not used: 'foo', 'bar' \\(not an argument of dbarts::bart or dbartsControl\\)")
  expect_warning(smallFit(method.trt = "none", rngKind = "default"),
                 "'rngKind' \\(dbarts 1.0-0 has no choice of random number generator\\)")
  expect_warning(smallFit(method.trt = "none", foo = 1, rngKind = "default"),
                 "'rngKind' \\(dbarts 1.0-0.*; 'foo' \\(not an argument")
  
  ## the propensity model says the same of the arguments it is given
  expect_warning(smallFit(args.trt = quote(list(foo = 1))),
                 "'foo' \\(not an argument of dbarts::bart or dbartsControl\\)")
  
  ## a binary response has no residual scale to give a prior
  expect_warning(smallFit(y.binary, method.trt = "none", sigdf = 4),
                 "'sigdf' \\(a binary response has no residual scale")
  
  ## nothing the caller is entitled to give is reported
  expect_warning(smallFit(power = 1.5, base = 0.7, sigdf = 7, sigquant = 0.5, rngSeed = 3L,
                          args.trt = quote(list(power = 1.2))), NA)
  expect_warning(smallFit(), NA)
})

test_that("sigdf and sigquant are named as unused where stan4bart fits the response", {
  skip_if_not_installed("stan4bart")
  
  test.df <- data.frame(y, z, x, g = factor(rep_len(letters[1:4], length(y))))
  fit <- function(...)
    eval(as.call(c(list(quote(bartCause::bartc), quote(y), quote(z), quote(X1 + X2 + X3),
                        data = quote(test.df), group.by = quote(g), method.trt = "none",
                        n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 6L,
                        n.trees = 10L, verbose = FALSE), list(...))))
  
  expect_warning(fit(sigdf = 7, sigquant = 0.5),
                 "arguments not used: 'sigdf', 'sigquant' \\(the parametric component draws the residual standard deviation")
  expect_warning(fit(foo = 1), "'foo' \\(not an argument of stan4bart\\)")
  expect_warning(fit(power = 1.5, base = 0.7, n.trees = 10L, k = 2), NA)
  expect_warning(fit(), NA)
})

test_that("arguments for the response method are not reported by the propensity fit", {
  expect_warning(bartc(y, z, x, method.rsp = "tmle", method.trt = "bart", maxIter = 2L,
                       n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 6L, n.trees = 10L,
                       verbose = FALSE),
                 NA)
  expect_warning(bartc(y, z, x, method.rsp = "p.weight", method.trt = "bart", maxIter = 2L,
                       n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 6L, n.trees = 10L,
                       verbose = FALSE),
                 NA)
  
  skip_if_not_installed("stan4bart")
  test.df <- data.frame(y, z, x, g = factor(rep_len(letters[1:4], length(y))))
  expect_warning(bartc(y, z, X1 + X2 + X3, data = test.df, group.by = g, method.rsp = "p.weight",
                       method.trt = "bart", maxIter = 2L, chains = 1L, iter = 8L, warmup = 4L,
                       bart_args = list(n.trees = 7L), verbose = FALSE),
                 NA)
})
test_that("an argument of dbartsControl is not reported as unused by the propensity fit", {
  ## proposal.probs, which dbarts takes through bart's '...', was named as unused
  ## for the propensity model alone
  warned <- character()
  withCallingHandlers(
    smallFit(proposal.probs = c(birth_death = 0.6, swap = 0, change = 0.4, birth = 0.5)),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("arguments not used", warned, fixed = TRUE)))
  
  ## one that neither takes is still named
  expect_warning(smallFit(method.rsp = "p.weight", foo = 1), "arguments not used: 'foo'")
})

