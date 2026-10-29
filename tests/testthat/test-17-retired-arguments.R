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
  ## both in one list reach the fit together
  expect_error(smallFit(method.trt = "none", args.rsp = quote(list(seed = 3L, rngSeed = 4L))),
               "'rngSeed' and 'seed' cannot both be given")
})

test_that("a retired name beside its successor in args.rsp or args.trt is overridden by the list entry, and named", {
  ## the fit and the warnings of a call, of which there is one for the call
  fitAndWarn <- function(...) {
    warned <- character()
    fit <- withCallingHandlers(smallFit(...), warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(fit = fit, warned = warned)
  }
  listOf <- function(...) str2lang(sprintf("list(%s)", paste0(c(...), collapse = ", ")))
  prior <- "tree.prior = dbarts::dbartsPriors$cgm(2, 0.9)"
  
  ## the list governs, as it does for any other argument, and the argument it overrides is named
  listOnly <- smallFit(method.trt = "none", seed = 7L, args.rsp = listOf(prior))
  both <- fitAndWarn(method.trt = "none", seed = 7L, power = 1.5, base = 0.7, args.rsp = listOf(prior))
  expect_equal(treePrior(both$fit$fit.rsp), c(power = 2, base = 0.9))
  expect_equal(both$fit$mu.hat.obs, listOnly$mu.hat.obs)
  expect_equal(both$warned, "arguments not used: 'power', 'base' (args.rsp gives 'tree.prior', which governs)")
  
  ## and the other way about: the retired name in the list, the successor beside it
  both <- fitAndWarn(method.trt = "none", tree.prior = quote(cgm(2, 0.9)), args.rsp = quote(list(power = 1.5)))
  expect_equal(treePrior(both$fit$fit.rsp)[["power"]], 1.5)
  expect_equal(both$warned, "arguments not used: 'tree.prior' (args.rsp gives 'power', which governs)")
  
  ## a tree prior beside one of its parts in the list is not kept in part
  both <- fitAndWarn(method.trt = "none", tree.prior = quote(cgm(3, 0.8)), args.rsp = quote(list(base = 0.7)))
  expect_equal(treePrior(both$fit$fit.rsp), c(power = 2, base = 0.7))
  expect_equal(both$warned, "arguments not used: 'tree.prior' (args.rsp gives 'base', which governs)")
  
  both <- fitAndWarn(method.trt = "none", sigdf = 4,
                     args.rsp = listOf("family = dbarts::dbartsFamilies$gaussian(sigma = dbarts::dbartsPriors$fixed(1))"))
  expect_equal(both$warned, "arguments not used: 'sigdf' (args.rsp gives 'family', which governs)")
  
  ## the propensity fit, and both lists in one call, are named in the one warning
  both <- fitAndWarn(args.trt = listOf(prior), power = 1.5)
  expect_equal(treePrior(both$fit$fit.trt), c(power = 2, base = 0.9))
  expect_equal(treePrior(both$fit$fit.rsp)[["power"]], 1.5)
  expect_equal(both$warned, "arguments not used: 'power' (args.trt gives 'tree.prior', which governs)")
  both <- fitAndWarn(args.trt = listOf(prior), args.rsp = listOf(prior), power = 1.5)
  expect_equal(both$warned, paste0("arguments not used: 'power' (args.trt gives 'tree.prior', which governs); ",
                                   "'power' (args.rsp gives 'tree.prior', which governs)"))
  
  ## nothing is said where nothing is overridden
  expect_equal(fitAndWarn(method.trt = "none", power = 1.5, args.rsp = quote(list(n.trees = 11L)))$warned, character())
  
  ## each set of retired names has its successor, in either direction
  drop <- bartCause:::dropOverriddenSpellings
  call <- quote(f(a = 1, power = 1, base = 2, split.probs = 3, sigdf = 4, sigquant = 5, rngSeed = 6, tree.prior = 7))
  expect_equal(names(drop(call, "tree.prior")), c("", "a", "sigdf", "sigquant", "rngSeed", "tree.prior"))
  expect_equal(names(drop(call, "family")), c("", "a", "power", "base", "split.probs", "rngSeed", "tree.prior"))
  expect_equal(names(drop(call, "seed")), c("", "a", "power", "base", "split.probs", "sigdf", "sigquant", "tree.prior"))
  expect_equal(names(drop(call, "base")), c("", "a", "power", "base", "split.probs", "sigdf", "sigquant", "rngSeed"))
  expect_equal(names(drop(quote(f(sigdf = 1, family = 2)), "sigquant")), c("", "sigdf"))
  expect_equal(names(drop(quote(f(sigdf = 1, family = 2)), "n.trees")), c("", "sigdf", "family"))
  expect_equal(attr(drop(call, c("power", "base")), "overridden"),
               c(tree.prior = "the list gives 'power', 'base', which govern"))
  
  ## each of the pair given beside the other in the call is still refused
  expect_error(smallFit(method.trt = "none", power = 1.5, tree.prior = quote(cgm(2, 0.9))),
               "'power' and 'tree.prior' cannot both be given")
})

test_that("bartc's own seed governs where rngSeed is also given, as in 1.0-10", {
  draws <- function(...) smallFit(method.trt = "none", ...)$mu.hat.obs
  expect_warning(both <- draws(seed = 3L, rngSeed = 4L),
                 "arguments not used: 'rngSeed' \\(seed gives each fit its own seed, which governs\\)")
  expect_equal(both, draws(seed = 3L))
  expect_false(isTRUE(all.equal(both, suppressWarnings(draws(seed = 5L, rngSeed = 4L)))))
  
  ## also when the rngSeed is in args.rsp or args.trt
  expect_warning(inList <- draws(seed = 3L, args.rsp = quote(list(rngSeed = 4L))), "'rngSeed'")
  expect_equal(inList, draws(seed = 3L))
  expect_warning(smallFit(seed = 3L, args.trt = quote(list(rngSeed = 4L))), "'rngSeed'")
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
