context("bartc response fits")

source(system.file("common", "linearData.R", package = "BartCause"))

test_that("bart fit matches manual call", {
  set.seed(22)
  bartcFit <- BartCause:::getBartResponseFit(y, z, x, testData, estimand = "ate", group.by = NULL, commonSup.rule = "none", commonSup.cut = NA,
                                         n.chains = 1L, n.threads = 1L, n.burn = 50L, n.samples = 150L, n.trees = 75L)
  x.train <- with(testData, cbind(x, z))
  # colnames(x.train) <- c("x1", "x2", "x3", "z")
  x.test <- x.train
  x.test[,"z"] <- 1 - x.test[,"z"]
  y <- testData$y
  set.seed(22)
  bartFit <- dbarts::bart2(x.train, y, x.test, n.chains = 1L, n.threads = 1L, n.burn = 50L, n.samples = 150L, n.trees = 75L, verbose = FALSE)
      
  expect_equal(bartFit$yhat.train, bartcFit$fit$yhat.train)
  expect_equal(bartFit$yhat.test,  bartcFit$fit$yhat.test)
})

test_that("pweight fits", {
  set.seed(22)
  testData$w <- 1 + rpois(length(testData$y), 0.5)
  
  testCall <- quote(bartc(y, z, x, testData, method.trt = "glm", method.rsp = "pweight",
                          n.chains = 1L, n.threads = 1L, n.samples = 200L, n.burn = 40L,
                          verbose = FALSE))
  
  expect_is(eval(testCall), "bartcFit")
  
  testCall$method.trt <- "bart"
  expect_is(eval(testCall), "bartcFit")
  
  testCall$method.trt <- "glm"
  testCall$weights <- quote(w)
  expect_is(eval(testCall), "bartcFit")
  
  testCall$method.trt <- "bart"
  expect_is(eval(testCall), "bartcFit")
  
  ## multiple chains
  testCall$n.chains  <- 4L
  testCall$n.samples <- 50L
  testCall$n.burn    <- 10L
  testCall$method.trt <- "glm"
  testCall$weights <- NULL
  
  expect_is(eval(testCall), "bartcFit")
  
  testCall$method.trt <- "bart"
  expect_is(eval(testCall), "bartcFit")
  
  testCall$method.trt <- "glm"
  testCall$weights <- quote(w)
  expect_is(eval(testCall), "bartcFit")
  
  testCall$method.trt <- "bart"
  expect_is(eval(testCall), "bartcFit")
})

