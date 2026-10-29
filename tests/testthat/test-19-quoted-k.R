context("a quoted prior for k with group.by")

test_that("a bare name for k is refused with the forms that work, where stan4bart fits the response", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")
  
  pattern <- paste0("'k' is a bare name, which is not taken as a prior; ",
                    "write k = chi(1.25, Inf) directly in the call, ",
                    "or give dbartsPriors$chi(1.25, Inf) or the string \"chi(1.25, Inf)\" instead")
  
  ## a name, which stan4bart would look up
  name <- quote(quote(chi))
  expect_error(fit(group.by = quote(country), k = name), pattern, fixed = TRUE)
  expect_error(fit(group.by = quote(country), args.rsp = call("list", k = name)), pattern, fixed = TRUE)
  expect_error(fit(group.by = quote(country), args.trt = call("list", k = name)), pattern, fixed = TRUE)
  ## and held in a variable
  held <- list(k = quote(chi))
  expect_error(fit(group.by = quote(country), args.rsp = quote(held)), pattern, fixed = TRUE)
  
  ## the forms the message names are the ones that fit
  for (k in list(dbarts::dbartsPriors$chi(1.25, Inf), "chi(1.25, Inf)")) {
    held <- list(k = k)
    fitted <- fit(group.by = quote(country), args.rsp = quote(held), args.trt = quote(held))
    expect_true(inherits(fitted$fit.rsp, "stan4bartFit"))
    expect_true(all(is.finite(fitted(fitted))))
  }
  fitted <- fit(group.by = quote(country), args.rsp = quote(list(k = "chi(1.25, Inf)")))
  expect_true(inherits(fitted$fit.rsp, "stan4bartFit"))
  expect_true(inherits(fit(group.by = quote(country), k = 2)$fit.rsp, "stan4bartFit"))
})

test_that("a quoted call for k is the prior it holds, where stan4bart fits the response, as it was in 1.0-10", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")
  
  ## k = quote(chi(1.25, Inf)), as an argument, in args.rsp or args.trt, and held in a list
  quoted <- quote(quote(chi(1.25, Inf)))
  held <- list(k = quote(chi(1.25, Inf)))
  fits <- list(fit(group.by = quote(country), k = quoted),
               fit(group.by = quote(country), args.rsp = call("list", k = quoted)),
               fit(group.by = quote(country), args.rsp = quote(held)))
  for (fitted in fits) {
    expect_true(inherits(fitted$fit.rsp, "stan4bartFit"))
    ## a prior on k has draws of k; a fixed k has none
    expect_false(is.null(fitted$fit.rsp$k))
    expect_true(all(is.finite(fitted(fitted))))
  }
  expect_true(inherits(fit(group.by = quote(country), args.trt = call("list", k = quoted))$fit.trt, "stan4bartFit"))
  expect_true(inherits(fit(group.by = quote(country), args.trt = quote(held))$fit.trt, "stan4bartFit"))
})

test_that("a quoted prior for k is still refused by dbarts where it fits the response", {
  expect_error(fit(args.rsp = quote(list(k = quote(chi)))), "'k' is a quoted expression")
  expect_error(fit(args.rsp = quote(list(k = quote(chi(1.25, Inf))))), "'k' is a quoted expression")
})

## 'k' is a partial name for bartc's keepCall, which R matches it to; the
## value is a prior written in the call, chi(1.25, Inf), that only the fits can
## evaluate, so it must reach them unevaluated

test_that("a prior written as k reaches the fits of both paths, where R matched it to keepCall", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")
  
  for (group.by in list(NULL, quote(country))) {
    label <- if (is.null(group.by)) "dbarts" else "stan4bart"
    args <- if (is.null(group.by)) list() else list(group.by = group.by)
    
    hyper <- do.call(fit, c(args, list(k = quote(chi(1.25, Inf)))), quote = TRUE)
    expect_true(inherits(hyper$fit.rsp, if (is.null(group.by)) "bart" else "stan4bartFit"), label = label)
    expect_true(all(is.finite(fitted(hyper))), label = label)
    ## a prior on k has draws of k; a fixed k has none
    expect_false(is.null(hyper$fit.rsp$k), label = label)
    expect_identical(hyper$call$k, quote(chi(1.25, Inf)), label = label)
    expect_null(hyper$call[["keepCall"]], label = label)
    
    fixed <- do.call(fit, c(args, list(k = 2)), quote = TRUE)
    expect_null(fixed$fit.rsp$k, label = label)
    expect_equal(fixed$call$k, 2, label = label)
    expect_null(fixed$call[["keepCall"]], label = label)
  }
})

test_that("keepCall, by its name or a partial name, still decides whether the call is kept", {
  expect_identical(fit(keepCall = FALSE)$call, call("NULL"))
  expect_identical(fit(keepC = FALSE)$call, call("NULL"))
  expect_identical(fit(keepCall = TRUE)$call[[1L]], quote(bartCause::bartc))
  
  ## with k beside it, k is an argument for the fits and keepCall is the formal
  expect_identical(fit(keepCall = FALSE, k = 2)$call, call("NULL"))
  expect_equal(fit(k = 2)$fit.rsp$leaf.prior$leaf.prior@k, 2)
})

test_that("k is told from other arguments that happen to hold the same value", {
  ## the partial match is found by name; two arguments with equal values were
  ## both renamed to k
  args <- c(list(quote(bartCause::bartc), quote(y), quote(z), quote(x1 + x2), data = quote(formData)),
            modifyList(small, list(n.trees = 5L)), list(k = 5))
  fitted <- eval(as.call(args))
  expect_equal(sum(names(fitted$call) == "k"), 1L)
  expect_equal(sum(names(fitted$call) == "n.trees"), 1L)
  expect_equal(fitted$fit.rsp$leaf.prior$leaf.prior@k, 5)
  expect_equal(fitted$fit.rsp$call$n.trees, 5)
})

## the fitters read some arguments of their '...' themselves, which must not
## evaluate the rest of them

test_that("one argument of the dots is read without evaluating the others", {
  call <- quote(f(a = 1, k = chi(1.25, Inf), b = stop("not to be evaluated")))
  expect_equal(bartCause:::getDotsArgument(call, "a", environment(), 0), 1)
  expect_equal(bartCause:::getDotsArgument(call, "missing", environment(), 0), 0)
  held <- 3
  expect_equal(bartCause:::getDotsArgument(quote(f(a = held)), "a", environment(), 0), 3)
})

test_that("a prior written as k reaches the fits beside combineChains and with the tmle response method", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")
  
  prior <- quote(chi(1.25, Inf))
  for (group.by in list(NULL, quote(country))) {
    args <- if (is.null(group.by)) list() else list(group.by = group.by)
    label <- if (is.null(group.by)) "dbarts" else "stan4bart"
    
    combined <- do.call(fit, c(args, list(k = prior, combineChains = TRUE)), quote = TRUE)
    expect_false(is.null(combined$fit.rsp$k), label = label)
    expect_true(all(is.finite(fitted(combined))), label = label)
    
    tmle <- suppressWarnings(do.call(fit, c(args, list(k = prior, method.rsp = "tmle", maxIter = 2L)), quote = TRUE))
    expect_false(is.null(tmle$fit.rsp$k), label = label)
    expect_true(all(is.finite(tmle$est)), label = label)
  }
  
  ## the tmle estimates of each group read n.threads from the call as well
  for (use.ranef in c(TRUE, FALSE)) {
    label <- paste("group.effects, use.ranef =", use.ranef)
    grouped <- suppressWarnings(fit(k = prior, method.rsp = "tmle", maxIter = 2L, group.by = quote(country),
                                    group.effects = TRUE, use.ranef = use.ranef))
    expect_false(is.null(grouped$fit.rsp$k), label = label)
    expect_equal(names(grouped$est), levels(formData$country), label = label)
    expect_true(all(vapply(grouped$est, function(e) all(is.finite(e)), NA)), label = label)
  }
})
