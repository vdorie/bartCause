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
