context("binary outcomes")

source(system.file("common", "binaryData.R", package = "bartCause"))

set.seed(22)
bartcFit <- bartc(y, z, x, data = testData,
                  method.rsp = "bart", method.trt = "bart", verbose = FALSE,
                  n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L)

test_that("binary outcome model matches manual", {
  set.seed(22)
  fit.trt <- dbarts::bart2(z ~ x, testData, verbose = FALSE,
                           n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L)
  p.score <- apply(pnorm(fit.trt$yhat.train), 2L, mean)
  expect_equal(p.score, fitted(bartcFit, type = "p.score"))
  
  x.train <- cbind(z = testData$z, testData$x, ps = p.score)
  x.test  <- cbind(z = 1, testData$x, ps = p.score)
  x.test <- rbind(x.test, x.test)
  x.test[seq.int(nrow(testData$x) + 1L, nrow(x.test)),"z"] <- 0
  
  fit.rsp <- dbarts::bart2(x.train, testData$y, x.test, verbose = FALSE,
                           n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L)
  expect_equal(extract(bartcFit, type = "mu.0"),
               pnorm(fit.rsp$yhat.test)[,seq.int(nrow(testData$x) + 1L, nrow(x.test))])
})

test_that("summary works for binary outcomes", {
  expect_is(summary(bartcFit), "bartcFit.summary")
  expect_is(summary(bartcFit, target = "sate"), "bartcFit.summary")
  expect_is(summary(bartcFit, target = "cate"), "bartcFit.summary")
})


test_that("binary outcome runs with tmle", {
  oldWarn <- getOption("warn")
  if (!requireNamespace("tmle", quietly = TRUE))
    options(warn = -1)
  
  expect_is(bartc(y, z, x, data = testData,
                  method.rsp = "tmle", method.trt = "bart", verbose = FALSE,
                  n.samples = 5L, n.burn = 5L, n.chains = 1L, n.threads = 1L), "bartcFit")
  
  options(warn = oldWarn)
})

