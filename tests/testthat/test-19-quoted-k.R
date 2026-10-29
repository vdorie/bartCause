context("a quoted prior for k with group.by")

test_that("a quoted prior for k is refused with the forms that work, where stan4bart fits the response", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")
  
  pattern <- paste0("'k' is a quoted expression, which is not taken as a prior; ",
                    "give dbartsPriors$chi(1.25, Inf) or the string \"chi(1.25, Inf)\" instead")
  
  ## a symbol, which stan4bart would look up, and a call, which it would evaluate
  ## to a prior that dbarts refuses to take quoted
  for (k in list(quote(quote(chi)), quote(quote(chi(1.25, Inf))))) {
    expect_error(fit(group.by = quote(country), args.rsp = call("list", k = k)), pattern, fixed = TRUE)
    expect_error(fit(group.by = quote(country), args.trt = call("list", k = k)), pattern, fixed = TRUE)
  }
  ## and held in a variable
  held <- list(k = quote(chi(1.25, Inf)))
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
