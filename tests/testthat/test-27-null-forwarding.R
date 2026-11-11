context("an argument written NULL reaches the fitter")

source(system.file("common", "groupedData.R", package = "bartCause"))

df <- with(testData, data.frame(y = y, z = z, V1 = x[,1L], V2 = x[,2L], V3 = x[,3L], g = g))

s4bResponse <- function(...) {
  set.seed(22)
  fit <- stan4bart::stan4bart(y ~ bart(z + V1 + V2 + V3) + (1 | g), df, treatment = z,
                              chains = 1L, iter = 16L, warmup = 8L, verbose = -1L,
                              bart_args = list(n.trees = 7L), ...)
  aperm(dbarts::extract(fit, sample = "train", combine_chains = FALSE), c(3L, 2L, 1L))
}

testQuietMixing("a NULL prior reaches stan4bart's response fit as a flat prior", {
  skip_if_not_installed("stan4bart")
  ## stan4bart takes prior as an argument from 0.0-14 (stan_args before)
  skip_if_not("prior" %in% names(formals(stan4bart::stan4bart)))
  
  set.seed(22)
  res <- bartCause:::getBartResponseFit(y, z, x, data = testData, estimand = "ate", group.by = g,
                                        commonSup.rule = "none", commonSup.cut = NA,
                                        chains = 1L, iter = 16L, warmup = 8L,
                                        bart_args = list(n.trees = 7L), prior = NULL)
  expect_true("prior" %in% names(res$fit$call))
  expect_null(res$fit$call$prior)
  expect_equal(res$mu.hat.obs, s4bResponse(prior = NULL))
  
  set.seed(22)
  res <- bartCause:::getBartResponseFit(y, z, x, data = testData, estimand = "ate", group.by = g,
                                        commonSup.rule = "none", commonSup.cut = NA,
                                        chains = 1L, iter = 16L, warmup = 8L,
                                        bart_args = list(n.trees = 7L), prior = normal(0, 1))
  expect_equal(res$fit$call$prior, quote(normal(0, 1)))
  expect_equal(res$mu.hat.obs, s4bResponse(prior = normal(0, 1)))
  
  # the other prior argument, and an argument that follows a NULL one
  set.seed(22)
  res <- bartCause:::getBartResponseFit(y, z, x, data = testData, estimand = "ate", group.by = g,
                                        commonSup.rule = "none", commonSup.cut = NA,
                                        chains = 1L, iter = 16L, warmup = 8L,
                                        bart_args = list(n.trees = 7L), prior_aux = NULL, QR = FALSE)
  expect_true("prior_aux" %in% names(res$fit$call))
  expect_false(isTRUE(res$fit$call$QR))
  expect_equal(res$mu.hat.obs, s4bResponse(prior_aux = NULL, QR = FALSE))
})

testQuietMixing("a NULL prior reaches stan4bart's treatment fit", {
  skip_if_not_installed("stan4bart")
  ## stan4bart takes prior as an argument from 0.0-14 (stan_args before)
  skip_if_not("prior" %in% names(formals(stan4bart::stan4bart)))
  
  s4bTreatment <- function(...) {
    set.seed(22)
    fit <- stan4bart::stan4bart(z ~ bart(V1 + V2 + V3) + (1 | g), df, chains = 1L, iter = 16L,
                                warmup = 8L, verbose = -1L, bart_args = list(n.trees = 7L), ...)
    apply(dbarts::extract(fit, combine_chains = FALSE), 1L, mean)
  }
  set.seed(22)
  res <- bartCause:::getBartTreatmentFit(y, z, x, data = testData, chains = 1L, iter = 16L, warmup = 8L,
                                         bart_args = list(n.trees = 7L), group.by = g, use.ranef = TRUE,
                                         prior = NULL)
  expect_true("prior" %in% names(res$fit$call))
  expect_equal(res$p.score, s4bTreatment(prior = NULL))
})

testQuietMixing("bartc passes a NULL prior from its dots and from args.rsp and args.trt", {
  skip_if_not_installed("stan4bart")
  ## stan4bart takes prior as an argument from 0.0-14 (stan_args before)
  skip_if_not("prior" %in% names(formals(stan4bart::stan4bart)))
  
  ## do.call, so that a written prior is the expression the fit resolves
  fitBartc <- function(...)
    suppressMessages(do.call(bartc, c(list(quote(y), quote(z), quote(V1 + V2 + V3), data = quote(df), group.by = quote(g),
                                           verbose = FALSE, seed = 1L, chains = 1L, iter = 16L, warmup = 8L,
                                           bart_args = list(n.trees = 7L)), list(...))))
  
  fit <- fitBartc(prior = NULL)
  expect_true("prior" %in% names(fit$fit.rsp$call))
  expect_true("prior" %in% names(fit$fit.trt$call))
  
  fit <- fitBartc(args.rsp = list(prior = NULL))
  expect_true("prior" %in% names(fit$fit.rsp$call))
  expect_false("prior" %in% names(fit$fit.trt$call))
  
  fit <- fitBartc(args.trt = list(prior = NULL))
  expect_false("prior" %in% names(fit$fit.rsp$call))
  expect_true("prior" %in% names(fit$fit.trt$call))
  
  fit <- fitBartc(prior = quote(normal(0, 1)))
  expect_equal(fit$fit.rsp$call$prior, quote(normal(0, 1)))
  expect_equal(fit$fit.trt$call$prior, quote(normal(0, 1)))
  
})

test_that("an argument written NULL that a dbarts or glm fit leaves out is left out", {
  args <- list(data = testData, n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L)
  
  set.seed(22)
  base <- do.call(bartCause:::getBartTreatmentFit, c(list(quote(y), quote(z), quote(x)), args))
  set.seed(22)
  res <- do.call(bartCause:::getBartTreatmentFit, c(list(quote(y), quote(z), quote(x)), args, list(monotone = NULL)))
  expect_equal(res$p.score, base$p.score)
  
  set.seed(22)
  base <- do.call(bartCause:::getBartResponseFit, c(list(quote(y), quote(z), quote(x)), args,
                                                   list(estimand = "ate", group.by = NULL, commonSup.rule = "none", commonSup.cut = NA)))
  set.seed(22)
  res <- do.call(bartCause:::getBartResponseFit, c(list(quote(y), quote(z), quote(x)), args,
                                                  list(estimand = "ate", group.by = NULL, commonSup.rule = "none", commonSup.cut = NA,
                                                       monotone = NULL)))
  expect_equal(res$mu.hat.obs, base$mu.hat.obs)
  
  res <- bartCause:::getGLMTreatmentFit(y, z, x, data = testData, start = NULL)
  expect_equal(res$p.score, bartCause:::getGLMTreatmentFit(y, z, x, data = testData)$p.score)
  
  set.seed(22)
  fit <- suppressMessages(bartc(y, z, x, data = testData, verbose = FALSE, seed = 1L, n.chains = 1L, n.threads = 1L,
                                n.burn = 3L, n.samples = 13L, n.trees = 7L, monotone = NULL))
  set.seed(22)
  base <- suppressMessages(bartc(y, z, x, data = testData, verbose = FALSE, seed = 1L, n.chains = 1L, n.threads = 1L,
                                 n.burn = 3L, n.samples = 13L, n.trees = 7L))
  expect_equal(fitted(fit), fitted(base))
})
