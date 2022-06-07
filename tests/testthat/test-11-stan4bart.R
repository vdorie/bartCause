context("semiparametric models fit with stan4bart")

source(system.file("common", "friedmanData.R", package = "bartCause"))

testData <- generateFriedmanData(100, ranef = TRUE, causal = TRUE)

test.df <- with(testData, data.frame(x, y, z, g.1))

test_that("semiparametric models are consistent with each other", {
  # Because this is not documented, to enable this test execute from R
  #   Sys.setenv(NOT_CRAN = "true")
  # or from shell
  #   export NOT_CRAN=true
  skip_on_cran()
  skip_if_not_installed("stan4bart")
  
  seed <- 0
  fit1 <- bartc(y, z, X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                group.by = g.1,
                seed = seed,
                data = test.df,
                verbose = FALSE)

  summary1 <- summary(fit1)
  
  fit2 <- bartc(y, z, X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                parametric = (1 | g.1),
                seed = seed,
                data = test.df,
                verbose = FALSE)
  
  summary2 <- summary(fit2)
  expect_true(inherits(fit2$fit.trt, "stan4bartFit"))
  expect_equal(fit2$fit.trt$call$formula, str2lang("z ~ (1 | g.1) + bart(X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + 
    X9 + X10)"))
  data_frame <- test.df[setdiff(names(test.df), "y")]
  fit_frame <- fit2$fit.trt$frame[names(data_frame)]
  expect_equal(fit_frame, data_frame)
  
  expect_true(inherits(fit2$fit.rsp, "stan4bartFit"))
  expect_equal(fit2$fit.rsp$call$formula, str2lang("y ~ z + ps + bart(X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + 
    X9 + X10) + (1 | g.1)"))
  expect_equal(fit2$fit.rsp$frame[names(test.df)], test.df)
  
  fit3 <- bartc(y, z, X1 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10,
                parametric = X4 + (1 | g.1),
                seed = seed,
                data = test.df,
                verbose = FALSE)
  
  summary3 <- summary(fit3)
  expect_true(inherits(fit3$fit.trt, "stan4bartFit"))
  expect_equal(fit3$fit.trt$call$formula, str2lang("z ~ X4 + (1 | g.1) + bart(X1 + X2 + X3 + X5 + X6 + X7 + X8 + 
    X9 + X10)"))
  data_frame <- test.df[setdiff(names(test.df), "y")]
  fit_frame <- fit3$fit.trt$frame[names(data_frame)]
  expect_equal(fit_frame, data_frame)
  
  expect_true(inherits(fit3$fit.rsp, "stan4bartFit"))
  expect_equal(fit3$fit.rsp$call$formula, str2lang("y ~ z + ps + bart(X1 + X2 + X3 + X5 + X6 + X7 + X8 + 
    X9 + X10) + (X4 + (1 | g.1))"))
  expect_equal(fit3$fit.rsp$frame[names(test.df)], test.df)

  
  expect_in_range <- function(x, r) expect_true(all(x >= r[1L] & x <= r[2L]))
  
  expect_in_range(summary1$estimates$estimate, c(5.0, 5.6))
  expect_in_range(summary2$estimates$estimate, c(5.0, 5.6))
  expect_in_range(summary3$estimates$estimate, c(5.0, 5.6))
  
  expect_in_range(summary1$estimates$sd, c(0.55, 0.72))
  expect_in_range(summary2$estimates$sd, c(0.55, 0.72))
  expect_in_range(summary3$estimates$sd, c(0.55, 0.72))
})


test_that("semiparametric predict works", {
  skip_if_not_installed("stan4bart")
  
  seed <- 0
  fit <- bartc(y, z, X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
               parametric = (1 | g.1),
               seed = seed,
               data = test.df,
               verbose = FALSE,
               iter = 10, warmup = 3, chains = 2,
               bart_args = list(keepTrees = TRUE, n.trees = 7))
  
  expect_true(is.finite(fitted(fit)))
  
  samples.train <- extract(fit, "icate")
  samples.test  <- predict(fit, test.df, type = "icate")
  
  expect_true(sqrt(mean((samples.train - samples.test)^2)) <= 1e-10)
  
  expect_is(sum <- summary(fit), "bartcFit.summary")
  expect_true(is.finite(sum$estimates$estimate))
})


