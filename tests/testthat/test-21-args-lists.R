context("args.rsp and args.trt as variables")

nTrees <- function(sampler) sampler$fit$control@n.trees

test_that("args.rsp and args.trt take a variable or any expression that evaluates to a list", {
  ar <- list(n.trees = 7L)
  expect_equal(nTrees(fit(args.rsp = quote(ar))$fit.rsp), 7L)
  expect_equal(nTrees(fit(args.rsp = quote(list(n.trees = 7L)))$fit.rsp), 7L)
  expect_equal(nTrees(fit(args.rsp = quote(c(ar, list(k = 2))))$fit.rsp), 7L)
  expect_equal(nTrees(fit(args.rsp = quote(base::list(n.trees = 7L)))$fit.rsp), 7L)
  expect_equal(nTrees(fit(args.rsp = ar)$fit.rsp), 7L)
  expect_equal(nTrees(fit()$fit.rsp), 10L)
  
  at <- list(n.trees = 5L)
  expect_equal(nTrees(fit(args.trt = quote(at))$fit.trt), 5L)
  expect_equal(nTrees(fit(args.trt = quote(at))$fit.rsp), 10L)
  expect_equal(nTrees(fit(args.rsp = quote(ar), args.trt = quote(at))$fit.trt), 5L)
  expect_equal(nTrees(fit(args.rsp = quote(ar), args.trt = quote(at))$fit.rsp), 7L)
  
  ## the same list read from a variable and written out is the same fit
  prior <- list(power = 1.5, base = 0.7)
  set.seed(1)
  written <- fit(args.rsp = quote(list(power = 1.5, base = 0.7)))
  set.seed(1)
  variable <- fit(args.rsp = quote(prior))
  expect_equal(variable$mu.hat.obs, written$mu.hat.obs)
  expect_equal(variable$fit.rsp$fit$model@tree.prior@power, 1.5)
  expect_equal(variable$fit.rsp$fit$model@tree.prior@base, 0.7)
  
  ## a value that is language reaches the fit as a value rather than being evaluated
  quoted <- list(k = quote(chi))
  expect_error(fit(args.rsp = quote(quoted)), "'k' is a quoted expression")
  ## and a NULL value is a value
  nothing <- list(n.trees = 7L, k = NULL)
  expect_equal(nTrees(fit(args.rsp = quote(nothing))$fit.rsp), 7L)
  
  ## an empty list is no arguments
  none <- list()
  expect_equal(nTrees(fit(args.rsp = quote(none))$fit.rsp), 10L)
  
  expect_error(fit(args.rsp = 3), "'args.rsp' must be a list of arguments")
  unnamed <- list(7L)
  expect_error(fit(args.rsp = quote(unnamed)), "every element of 'args.rsp' must be named")
  expect_error(fit(args.trt = quote(unnamed)), "every element of 'args.trt' must be named")
})

test_that("a prior written in a list for args.rsp or args.trt reaches the fit as written", {
  ## chi is not exported by dbarts, so a list that was evaluated first would fail on it
  set.seed(3)
  rsp <- fit(method.trt = "none", k = quote(chi(1.25, Inf)))
  set.seed(3)
  listed <- fit(method.trt = "none", args.rsp = quote(list(k = chi(1.25, Inf))))
  expect_identical(listed$mu.hat.obs, rsp$mu.hat.obs)
  expect_identical(listed$mu.hat.cf, rsp$mu.hat.cf)
  
  ## the propensity fit is the first drawn, so it matches the top-level k as well
  set.seed(3)
  top <- fit(k = quote(chi(1.25, Inf)))
  set.seed(3)
  listed <- fit(args.trt = quote(list(k = chi(1.25, Inf))))
  expect_identical(listed$p.score, top$p.score)
  
  ## the forms that were already read: held in a variable, and the string
  held <- list(k = dbarts::dbartsPriors$chi(1.25, Inf))
  set.seed(3)
  variable <- fit(method.trt = "none", args.rsp = quote(held))
  expect_identical(variable$mu.hat.obs, rsp$mu.hat.obs)
  set.seed(3)
  string <- fit(method.trt = "none", args.rsp = quote(list(k = "chi(1.25, Inf)")))
  expect_identical(string$mu.hat.obs, variable$mu.hat.obs)
})

test_that("each entry of a list written out is evaluated once", {
  counter <- new.env()
  counter$rsp <- 0L
  counter$trt <- 0L
  fitted <- fit(args.rsp = quote(list(n.trees = {counter$rsp <- counter$rsp + 1L; 7L})),
                args.trt = quote(list(n.trees = {counter$trt <- counter$trt + 1L; 5L})))
  expect_equal(counter$rsp, 1L)
  expect_equal(counter$trt, 1L)
  expect_equal(nTrees(fitted$fit.rsp), 7L)
  expect_equal(nTrees(fitted$fit.trt), 5L)
})

test_that("an empty list, NULL and an unnamed entry are as before", {
  expect_equal(nTrees(fit(args.rsp = quote(list()))$fit.rsp), 10L)
  expect_equal(nTrees(fit(args.trt = quote(list()))$fit.trt), 10L)
  expect_equal(nTrees(fit(args.rsp = quote(NULL), args.trt = quote(NULL))$fit.rsp), 10L)
  expect_error(fit(args.rsp = quote(list(7L))), "every element of 'args.rsp' must be named")
  expect_error(fit(args.trt = quote(list(7L))), "every element of 'args.trt' must be named")
})

test_that("a forest constraint in a list for args.rsp with method.rsp = 'bcf' is read as before", {
  ## an entry is evaluated by the response fitter, so the constraint is built beforehand
  fitted <- fit(method.rsp = "bcf",
                args.rsp = quote(list(mu.interactions = dbarts::dbartsForests$interactions(max.order = 1))))
  expect_s3_class(fitted, "bartcFit")
  expect_error(fit(method.rsp = "bcf", args.rsp = quote(list(mu.interactions = interactions(max.order = 1)))),
               "interactions")
})
