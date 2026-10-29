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
