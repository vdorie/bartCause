context("bartc response fits")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("bart fit matches manual call", {
  set.seed(22)
  bartcFit <- bartCause:::getBartResponseFit(y, z, x, testData, estimand = "ate", group.by = NULL, commonSup.rule = "none", commonSup.cut = NA,
                                         n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L)
  x.train <- with(testData, cbind(x, z))
  # colnames(x.train) <- c("x1", "x2", "x3", "z")
  x.test <- x.train
  x.test[,"z"] <- 1 - x.test[,"z"]
  y <- testData$y
  set.seed(22)
  bartFit <- dbarts::bart2(x.train, y, x.test, n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L, verbose = FALSE)
      
  expect_equal(bartFit$yhat.train, bartcFit$fit$yhat.train)
  expect_equal(bartFit$yhat.test,  bartcFit$fit$yhat.test)
})

test_that("p.weight fits", {
  set.seed(22)
  testData$w <- 1 + rpois(length(testData$y), 0.5)
  
  testCall <- quote(bartc(y, z, x, testData, method.trt = "glm", method.rsp = "p.weight",
                          n.chains = 1L, n.threads = 1L, n.samples = 13L, n.burn = 3L, n.trees = 7L,
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

source(system.file("common", "groupedData.R", package = "bartCause"))

test_that("rbart_vi fit matches manual call", {
  set.seed(22)
  bartcFit <- bartCause:::getBartResponseFit(y, z, x, testData, estimand = "ate", group.by = g, commonSup.rule = "none", commonSup.cut = NA,
                                             n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L)
  x.train <- with(testData, cbind(x, z))
  # colnames(x.train) <- c("x1", "x2", "x3", "z")
  x.test <- x.train
  x.test[,"z"] <- 1 - x.test[,"z"]
  y <- testData$y
  set.seed(22)
  bartFit <- dbarts::rbart_vi(x.train, y, x.test, group.by = testData$g, group.by.test = testData$g,
                              n.chains = 1L, n.threads = 1L, n.burn = 3L, n.samples = 13L, n.trees = 7L, verbose = FALSE)
      
  expect_equal(bartFit$yhat.train, bartcFit$fit$yhat.train)
  expect_equal(bartFit$yhat.test,  bartcFit$fit$yhat.test)
})

