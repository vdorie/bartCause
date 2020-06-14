context("bartc input")

source(system.file("common", "groupedData.R", package = "bartCause"))

test_that("treatment call with data as list argument returns valid output", {
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, testData)
  ## check that it generates the correct call and that the environments point here
  expect_equal(res$call, parse(text = "stats::glm(z ~ x, data = testData)")[[1L]])
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, testData, group.by = g, use.ranef = FALSE)
  expect_equal(res$call, parse(text = "stats::glm(z ~ x + g, data = testData)")[[1L]])
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, testData, group.by = g, use.ranef = TRUE, use.lmer = FALSE)
  expect_equal(res$call, parse(text = "stats::glm(z ~ x, data = testData)")[[1L]])
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, testData, group.by = g, use.ranef = TRUE, use.lmer = TRUE)
  expect_equal(res$call, parse(text = "stats::glm(z ~ x + (1 | g), data = testData)")[[1L]])
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
})

test_that("treatment call with data as data.frame argument returns valid output", {
  df <- with(testData, data.frame(x, z, y, p.score, g))
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, X1 + X2 + X3, df)
  expect_true(res$call == parse(text = "stats::glm(z ~ X1 + X2 + X3, data = df)")[[1L]])
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, X1 + X2 + X3, df, group.by = g, use.ranef = FALSE)
  expect_equal(res$call, parse(text = "stats::glm(z ~ X1 + X2 + X3 + g, data = df)")[[1L]])
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, X1 + X2 + X3, df, group.by = g, use.ranef = TRUE, use.lmer = FALSE)
  expect_equal(res$call, parse(text = "stats::glm(z ~ X1 + X2 + X3, data = df)")[[1L]])
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, X1 + X2 + X3, df, group.by = g, use.ranef = TRUE, use.lmer = TRUE)
  expect_equal(res$call, parse(text = "stats::glm(z ~ X1 + X2 + X3 + (1 | g), data = df)")[[1L]])
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  
  confounders <- "X1 + X2 + X3"
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, confounders, df, group.by = g, use.ranef = FALSE)
  expect_equal(res$call, parse(text = "stats::glm(z ~ X1 + X2 + X3 + g, data = df)")[[1L]])
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, confounders, df, group.by = g, use.ranef = TRUE, use.lmer = FALSE)
  expect_equal(res$call, parse(text = "stats::glm(z ~ X1 + X2 + X3, data = df)")[[1L]])
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, confounders, df, group.by = g, use.ranef = TRUE, use.lmer = TRUE)
  expect_equal(res$call, parse(text = "stats::glm(z ~ X1 + X2 + X3 + (1 | g), data = df)")[[1L]])
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
})

test_that("treatment call with literal arguments retuns valid output", {
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x)
  expect_equal(res$call, parse(text = "stats::glm(z ~ V1 + V2 + V3, data = df)")[[1L]])
  expect_true(all(colnames(res$df) %in% c("z", "V1", "V2", "V3")))
  
  colnames(testData$x) <- c("x1", "x2", "z")
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x)
  expect_equal(res$call, parse(text = "stats::glm(zz ~ x1 + x2 + z, data = df)")[[1L]])
  expect_true(all(colnames(res$df) %in% c("zz", "x1", "x2", "z")))
  
  x <- as.data.frame(testData$x)
  z <- testData$z
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, z, x)
  expect_true(res$call == parse(text = "stats::glm(zz ~ x1 + x2 + z, data = df)")[[1L]])
  
  colnames(testData$x) <- NULL
  
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x, group.by = testData$g, use.ranef = FALSE)
  expect_equal(res$call, parse(text = "stats::glm(z ~ V1 + V2 + V3 + g, data = df)")[[1L]])
  expect_true(all(colnames(res$df) %in% c("z", "V1", "V2", "V3", "g")))
  
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x, group.by = testData$g, use.ranef = TRUE, use.lmer = FALSE)
  expect_equal(res$call, parse(text = "stats::glm(z ~ V1 + V2 + V3, data = df)")[[1L]])
  expect_true(all(colnames(res$df) %in% c("z", "V1", "V2", "V3", "g")))
  
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x, group.by = testData$g, use.ranef = TRUE, use.lmer = TRUE)
  expect_equal(res$call, parse(text = "stats::glm(z ~ V1 + V2 + V3 + (1 | g), data = df)")[[1L]])
  expect_true(all(colnames(res$df) %in% c("z", "V1", "V2", "V3", "g")))
})

test_that("response call with data as list argument returns valid output", {
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, x, testData)
  expect_equal(res$call, parse(text = "stats::lm(y ~ x + z, data = testData)")[[1L]])
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  
  ## p.score in testData
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, x, testData, p.score = p.score)
  expect_equal(res$call, parse(text = "stats::lm(y ~ x + p.score + z, data = testData)")[[1L]])
  
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  
  ## p.score supplied as literal
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, x, testData, p.score = testData$p.score)
  expect_equal(res$call, parse(text = "stats::lm(y ~ x + ps + z, data = data)")[[1L]])
  
  expect_false(identical(res$env, currentEnv))
  expect_false(identical(environment(res$call[[2L]]),currentEnv))
  expect_equal(res$env$data$ps, testData$p.score)
  
  expect_equal(res$trt, "z")
})

test_that("response call with data as data.frame argument returns valid output", {
  df <- with(testData, data.frame(x, y, z, p.score))
  
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, X1 + X2 + X3, df)
  expect_equal(res$call, parse(text = "stats::lm(y ~ X1 + X2 + X3 + z, data = df)")[[1L]])
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  
  ## p.score in testData
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, X1 + X2 + X3, df, p.score = p.score)
  expect_equal(res$call, parse(text = "stats::lm(y ~ X1 + X2 + X3 + p.score + z, data = df)")[[1L]])
  
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  
  ## p.score supplied as literal
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, X1 + X2 + X3, df, p.score = df$p.score)
  expect_equal(res$call, parse(text = "stats::lm(y ~ X1 + X2 + X3 + ps + z, data = data)")[[1L]])
  
  expect_false(identical(res$env, currentEnv))
  expect_false(identical(environment(res$call[[2L]]), currentEnv))
  expect_equal(res$env$data$ps, df$p.score)
  
  expect_equal(res$trt, "z")
  
  
  confounders <- "X1 + X2 + X3"
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, confounders, df)
  expect_equal(res$call, parse(text = "stats::lm(y ~ X1 + X2 + X3 + z, data = df)")[[1L]])
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  
  ## p.score in testData
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, confounders, df, p.score = p.score)
  expect_equal(res$call, parse(text = "stats::lm(y ~ X1 + X2 + X3 + p.score + z, data = df)")[[1L]])
  
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  
  ## p.score supplied as literal
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, confounders, df, p.score = df$p.score)
  expect_equal(res$call, parse(text = "stats::lm(y ~ X1 + X2 + X3 + ps + z, data = data)")[[1L]])
  
  expect_false(identical(res$env, currentEnv))
  expect_false(identical(environment(res$call[[2L]]), currentEnv))
  expect_equal(res$env$data$ps, df$p.score)
  
  expect_equal(res$trt, "z")
})

test_that("response call with data as data.frame argument returns valid output", {
  res <- bartCause:::getResponseLiteralCall(stats::lm, testData$y, testData$z, testData$x)
  expect_equal(res$call, parse(text = "stats::lm(y ~ V1 + V2 + V3 + z, data = df)")[[1L]])
  expect_true(all(colnames(res$df) %in% c("z", "y", "V1", "V2", "V3")))
  
  expect_equal(res$trt, "z")
  
  colnames(testData$x) <- c("x1", "x2", "z")
  res <- bartCause:::getResponseLiteralCall(stats::lm, testData$y, testData$z, testData$x)
  expect_equal(res$call, parse(text = "stats::lm(y ~ x1 + x2 + z + zz, data = df)")[[1L]])
  expect_true(all(colnames(res$df) %in% c("zz", "y", "x1", "x2", "z")))
  
  colnames(testData$x) <- NULL
    
  res <- bartCause:::getResponseLiteralCall(stats::lm, testData$y, testData$z, testData$x, p.score = testData$p.score)
  expect_equal(res$call, parse(text = "stats::lm(y ~ V1 + V2 + V3 + z + ps, data = df)")[[1L]])
  expect_true(all(colnames(res$df) %in% c("z", "y", "V1", "V2", "V3", "ps")))
  
  expect_equal(res$df$ps, testData$p.score)
   
  expect_equal(res$trt, "z")
})

