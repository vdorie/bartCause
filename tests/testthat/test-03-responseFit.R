context("bartc response fits")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("bart fit matches manual call", {
  set.seed(22)
  bartcFit <- bartCause:::getBartResponseFit(y, z, x, data = testData, estimand = "ate",
                                             group.by = NULL, commonSup.rule = "none", commonSup.cut = NA,
                                             n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L)
  x.train <- with(testData, cbind(z, x))
  # colnames(x.train) <- c("x1", "x2", "x3", "z")
  x.test <- x.train
  x.test[,"z"] <- 1 - x.test[,"z"]
  y <- testData$y
  set.seed(22)
  bartFit <- dbarts::bart2(x.train, y, x.test, n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L, verbose = FALSE)
      
  expect_equal(bartFit$yhat.train, bartcFit$fit$yhat.train)
  expect_equal(bartFit$yhat.test,  bartcFit$fit$yhat.test)
})

test_that("bcf fit matches manual call", {
  n <- length(testData$y)
  set.seed(22)
  bcfFit <- bartCause:::getBCFResponseFit(y, z, x, data = testData, estimand = "ate",
                                          group.by = NULL, commonSup.rule = "none", commonSup.cut = NA,
                                          n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L,
                                          n.trees = 7L, seed = 5L)
  # the literal builder assembles y ~ z + V1 + V2 + V3, so the design is
  # (z, V1, V2, V3) and both forests are masked out of the treatment column
  df <- with(testData, data.frame(y = y, z = z, V1 = x[,1], V2 = x[,2], V3 = x[,3]))
  data <- dbarts::dbartsData(y ~ z + V1 + V2 + V3, data = df, subset = seq_len(n),
                             bases = list(NULL, cbind(1 - df$z, df$z)))
  control <- dbarts::dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 7L, n.burn = 3L,
                                   n.samples = 13L, verbose = FALSE, updateState = FALSE,
                                   seed = 5L)
  set.seed(22)
  sampler <- dbarts::dbarts(data, control = control,
                            tree.prior = dbarts::dbartsPriors$cgm(2.0, 0.95),
                            forests = list(dbarts::forest(vars = c("V1", "V2", "V3")),
                                           dbarts::forest(vars = c("V1", "V2", "V3"), n.trees = 50L,
                                                          base = 0.25, power = 3, sd = 1,
                                                          amplitude.prior.variance = 0.5,
                                                          update.amplitude = TRUE)))
  sampler$sampleTreesFromPrior(updateState = FALSE)
  burn    <- sampler$run(0L, 3L, updateState = FALSE)
  samples <- sampler$run(0L, 13L, updateState = FALSE)

  expect_identical(bcfFit$fit$forests$mu, t(samples$forestFits[,1L,]))
  expect_identical(bcfFit$fit$forests$tau, t(samples$forestFits[,2L,]))
  expect_identical(bcfFit$mu.hat.obs, t(samples$train))
  expect_identical(bcfFit$fit$sigma, matrix(samples$sigma, nrow = 1L))
  expect_identical(bcfFit$fit$first.sigma, matrix(burn$sigma, nrow = 1L))
  expect_equal(bcfFit$name.trt, "z")
  expect_equal(as.vector(bcfFit$trt), testData$z)
  expect_equal(bcfFit$missingRows, rep_len(FALSE, n))
})

test_that("getBCFResponseFit defaults n.chains to 10 and refuses what bcf cannot express", {
  set.seed(22)
  res <- bartCause:::getBCFResponseFit(y, z, x, data = testData, estimand = "ate", group.by = NULL,
                                       commonSup.rule = "none", commonSup.cut = NA,
                                       n.threads = 1L, n.burn = 2L, n.samples = 3L, n.trees = 3L)
  expect_equal(dim(res$mu.hat.obs), c(10L, 3L, length(testData$y)))

  expect_error(
    bartCause:::getBCFResponseFit(y, z, x, data = testData, estimand = "ate", group.by = NULL,
                                  commonSup.rule = "none", commonSup.cut = NA, crossvalidate = TRUE,
                                  n.chains = 1L, n.threads = 1L, n.burn = 2L, n.samples = 3L),
    "crossvalidation is not supported for response method 'bcf'")
  expect_error(
    bartCause:::getBCFResponseFit(y, z, x, parametric = x[,1L], data = testData, estimand = "ate",
                                  commonSup.rule = "none", commonSup.cut = NA,
                                  n.chains = 1L, n.threads = 1L, n.burn = 2L, n.samples = 3L),
    "does not support 'parametric'")
  expect_error(bartCause:::getBCFResponseFit(treatment = z, confounders = x, data = testData),
               "'response' variable must be specified")
  expect_error(bartCause:::getBCFResponseFit(y, z, x, data = testData, estimand = "bogus",
                                             commonSup.rule = "none", commonSup.cut = NA),
               "estimand must be one of")
})

test_that("the response builders return the resolved propensity score name", {
  built <- bartCause:::getResponseLiteralCall(dbarts::dbartsData, testData$y, testData$z, testData$x,
                                              p.score = testData$p.score)
  expect_equal(built$p.score, "ps")
  expect_equal(length(built), 5L)

  withoutScore <- bartCause:::getResponseLiteralCall(dbarts::dbartsData, testData$y, testData$z,
                                                     testData$x)
  expect_null(withoutScore$p.score)
})

test_that("p.weight fits", {
  set.seed(22)
  testData$w <- 1 + rpois(length(testData$y), 0.5)
  
  testCall <- quote(bartc(y, z, x, data = testData, method.trt = "glm", method.rsp = "p.weight",
                          n.chains = 1L, n.threads = 1L, n.samples = 13L, n.burn = 3L, n.trees = 7L,
                          verbose = FALSE))
  
  expect_is(eval(testCall), "bartcFit")
  
  testCall$method.trt <- "bart"
  expect_is(eval(testCall), "bartcFit")
  
  testCall$method.trt <- "glm"
  testCall$weights <- quote(w)
  expect_is(eval(testCall), "bartcFit")

  ## As of dbarts 1.0-0 a weighted probit BART has no tractable latent form and
  ## is refused, so a BART propensity model is fit unweighted and says so once;
  ## the weights still enter the treatment-effect estimators.
  testCall$method.trt <- "bart"
  expect_message(expect_is(eval(testCall), "bartcFit"),
                 "propensity score model is fit unweighted")

  ## multiple chains
  testCall$n.chains  <- 4L
  testCall$method.trt <- "glm"
  testCall$weights <- NULL

  expect_is(eval(testCall), "bartcFit")

  testCall$method.trt <- "bart"
  expect_is(eval(testCall), "bartcFit")

  testCall$method.trt <- "glm"
  testCall$weights <- quote(w)
  expect_is(eval(testCall), "bartcFit")

  testCall$method.trt <- "bart"
  expect_message(expect_is(eval(testCall), "bartcFit"),
                 "propensity score model is fit unweighted")
})

source(system.file("common", "groupedData.R", package = "bartCause"))

test_that("varying intercept fit matches manual call", {
  skip_if_not_installed("stan4bart")
  df <- with(testData, data.frame(y = y, z = z, V1 = x[,1L], V2 = x[,2L], V3 = x[,3L], g = g))
  
  set.seed(22)
  bartcFit <- bartCause:::getBartResponseFit(y, z, x, data = testData, estimand = "ate", group.by = g,
                                             commonSup.rule = "none", commonSup.cut = NA,
                                             chains = 1L, iter = 16L, warmup = 8L,
                                             bart_args = list(n.trees = 7L))
  expect_true(inherits(bartcFit$fit, "stan4bartFit"))
  expect_equal(bartcFit$fit$call$formula, str2lang("y ~ bart(z + V1 + V2 + V3) + (1 | g)"))
  
  set.seed(22)
  s4bFit <- stan4bart::stan4bart(y ~ bart(z + V1 + V2 + V3) + (1 | g), df, treatment = z,
                                 chains = 1L, iter = 16L, warmup = 8L, verbose = -1L,
                                 bart_args = list(n.trees = 7L))
  
  expect_equal(bartcFit$mu.hat.obs,
               aperm(dbarts::extract(s4bFit, sample = "train", combine_chains = FALSE), c(3L, 2L, 1L)))
  expect_equal(bartcFit$mu.hat.cf,
               aperm(dbarts::extract(s4bFit, sample = "test", combine_chains = FALSE), c(3L, 2L, 1L)))
})

# commenting this out until crossvalidation calls have more control over run time
if (FALSE) test_that("xbart fit matches manual call", {
  set.seed(22)
  res <- bartCause:::getBartResponseFit(y, z, x, data = testData,
                                        estimand = "ate", group.by = NULL, commonSup.rule = "none", commonSup.cut = NA,
                                        n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L,
                                        crossvalidate = TRUE)
})

test_that("getBartResponseFit requires response, treatment, confounders, and a valid estimand", {
  expect_error(bartCause:::getBartResponseFit(treatment = z, confounders = x, data = testData),
               "'response' variable must be specified")
  expect_error(bartCause:::getBartResponseFit(y, confounders = x, data = testData),
               "'treatment' variable must be specified")
  expect_error(bartCause:::getBartResponseFit(y, z, data = testData),
               "'confounders' variable must be specified")
  expect_error(bartCause:::getBartResponseFit(y, z, x, data = testData, estimand = "bogus"),
               "estimand must be one of")
})

test_that("getBartResponseFit rejects crossvalidate with a varying-intercept model", {
  expect_error(
    bartCause:::getBartResponseFit(y, z, x, data = testData, estimand = "ate", group.by = g,
                                   commonSup.rule = "none", commonSup.cut = NA, crossvalidate = TRUE,
                                   n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L),
    "crossvalidation not yet supported"
  )
})

test_that("getBartResponseFit rejects a missing response with a varying-intercept model", {
  skip_if_not_installed("stan4bart")
  missData <- testData
  missData$y[seq_len(5L)] <- NA
  expect_error(
    bartCause:::getBartResponseFit(y, z, x, data = missData, estimand = "ate", group.by = g,
                                   commonSup.rule = "none", commonSup.cut = NA,
                                   chains = 1L, iter = 16L, warmup = 8L, bart_args = list(n.trees = 7L)),
    "cannot fit with missing response values"
  )
})

test_that("getBartResponseFit defaults n.chains to 10 when unspecified", {
  set.seed(22)
  res <- bartCause:::getBartResponseFit(y, z, x, data = testData, estimand = "ate", group.by = NULL,
                                        commonSup.rule = "none", commonSup.cut = NA,
                                        n.threads = 1L, n.burn = 2L, n.samples = 3L, n.trees = 3L)
  expect_equal(dim(res$mu.hat.obs), c(10L, 3L, length(testData$y)))
})

test_that("bartc handles missing response data with a single chain (no combined-chains reshape)", {
  # regression coverage: the missing-data reshape path in getBartResponseFit has
  # a 2-D (single chain) branch that combined-chains-only tests never exercised
  missData <- testData
  missData$y[seq_len(10L)] <- NA

  fit <- bartc(y, z, x, data = missData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L)
  expect_is(fit, "bartcFit")
  expect_equal(fit$n.chains, 1L)
  expect_equal(dim(fit$mu.hat.obs), c(13L, length(missData$y)))
  expect_equal(sum(fit$missingRows), 10L)

  icate <- extract(fit, "icate", combineChains = FALSE)
  expect_false(anyNA(icate))
  expect_equal(dim(icate), c(13L, length(missData$y)))
})

