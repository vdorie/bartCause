context("bartc input")

source(system.file("common", "groupedData.R", package = "bartCause"))

test_that("treatment call with data as list argument returns valid output", {
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, data = testData)
  ## check that it generates the correct call and that the environments point here
  expect_equal(res$call, str2lang("stats::glm(z ~ x, data = testData)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, data = testData, group.by = g, use.ranef = FALSE)
  expect_equal(res$call, str2lang("stats::glm(z ~ x + g, data = testData)"))
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, data = testData, group.by = g, use.ranef = TRUE, use.lmer = FALSE)
  expect_equal(res$call, str2lang("stats::glm(z ~ x, data = testData)"))
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, data = testData, group.by = g, use.ranef = TRUE, use.lmer = TRUE)
  expect_equal(res$call, str2lang("stats::glm(z ~ x + (1 | g), data = testData)"))
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
})

test_that("treatment call with data as data.frame argument returns valid output", {
  df <- with(testData, data.frame(x, z, y, p.score, g))
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, X1 + X2 + X3, data = df)
  expect_true(res$call == str2lang("stats::glm(z ~ X1 + X2 + X3, data = df)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, X1 + X2 + X3, data = df, group.by = g, use.ranef = FALSE)
  expect_equal(res$call, str2lang("stats::glm(z ~ X1 + X2 + X3 + g, data = df)"))
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, X1 + X2 + X3, data = df, group.by = g, use.ranef = TRUE, use.lmer = FALSE)
  expect_equal(res$call, str2lang("stats::glm(z ~ X1 + X2 + X3, data = df)"))
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, X1 + X2 + X3, data = df, group.by = g, use.ranef = TRUE, use.lmer = TRUE)
  expect_equal(res$call, str2lang("stats::glm(z ~ X1 + X2 + X3 + (1 | g), data = df)"))
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  confounders <- "X1 + X2 + X3"
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, confounders, data = df, group.by = g, use.ranef = FALSE)
  expect_equal(res$call, str2lang("stats::glm(z ~ X1 + X2 + X3 + g, data = df)"))
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, confounders, data = df, group.by = g, use.ranef = TRUE, use.lmer = FALSE)
  expect_equal(res$call, str2lang("stats::glm(z ~ X1 + X2 + X3, data = df)"))
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, confounders, data = df, group.by = g, use.ranef = TRUE, use.lmer = TRUE)
  expect_equal(res$call, str2lang("stats::glm(z ~ X1 + X2 + X3 + (1 | g), data = df)"))
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
})

test_that("treatment call with literal arguments retuns valid output", {
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x)
  expect_equal(res$call, str2lang("stats::glm(z ~ V1 + V2 + V3, data = df)"))
  expect_true(all(colnames(res$df) %in% c("z", "V1", "V2", "V3")))
  
  colnames(testData$x) <- c("x1", "x2", "z")
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x)
  expect_equal(res$call, str2lang("stats::glm(zz ~ x1 + x2 + z, data = df)"))
  expect_true(all(colnames(res$df) %in% c("zz", "x1", "x2", "z")))
  
  x <- as.data.frame(testData$x)
  z <- testData$z
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, z, x)
  expect_true(res$call == str2lang("stats::glm(zz ~ x1 + x2 + z, data = df)"))
  
  colnames(testData$x) <- NULL
  
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x, group.by = testData$g, use.ranef = FALSE)
  expect_equal(res$call, str2lang("stats::glm(z ~ V1 + V2 + V3 + g, data = df)"))
  expect_true(all(colnames(res$df) %in% c("z", "V1", "V2", "V3", "g")))
  
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x, group.by = testData$g, use.ranef = TRUE, use.lmer = FALSE)
  expect_equal(res$call, str2lang("stats::glm(z ~ V1 + V2 + V3, data = df)"))
  expect_true(all(colnames(res$df) %in% c("z", "V1", "V2", "V3", "g")))
  
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x, group.by = testData$g, use.ranef = TRUE, use.lmer = TRUE)
  expect_equal(res$call, str2lang("stats::glm(z ~ V1 + V2 + V3 + (1 | g), data = df)"))
  expect_true(all(colnames(res$df) %in% c("z", "V1", "V2", "V3", "g")))
})

test_that("response call with data as list argument returns valid output", {
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, x, data = testData)
  expect_equal(res$call, str2lang("stats::lm(y ~ x + z, data = testData)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  
  ## p.score in testData
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, x, data = testData, p.score = p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ x + p.score + z, data = testData)"))
  
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  
  ## p.score supplied as literal
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, x, data = testData, p.score = testData$p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ x + ps + z, data = data)"))
  
  expect_false(identical(res$env, currentEnv))
  expect_false(identical(environment(res$call[[2L]]),currentEnv))
  expect_equal(res$env$data$ps, testData$p.score)
  
  expect_equal(res$trt, "z")
})

test_that("response call with data as data.frame argument returns valid output", {
  df <- with(testData, data.frame(x, y, z, p.score))
  
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, X1 + X2 + X3, data = df)
  expect_equal(res$call, str2lang("stats::lm(y ~ X1 + X2 + X3 + z, data = df)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  
  ## p.score in testData
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, X1 + X2 + X3, data = df, p.score = p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ X1 + X2 + X3 + p.score + z, data = df)"))
  
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  
  ## p.score supplied as literal
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, X1 + X2 + X3, data = df, p.score = df$p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ X1 + X2 + X3 + ps + z, data = data)"))
  
  expect_false(identical(res$env, currentEnv))
  expect_false(identical(environment(res$call[[2L]]), currentEnv))
  expect_equal(res$env$data$ps, df$p.score)
  
  expect_equal(res$trt, "z")
  
  
  confounders <- "X1 + X2 + X3"
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, confounders, data = df)
  expect_equal(res$call, str2lang("stats::lm(y ~ X1 + X2 + X3 + z, data = df)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  
  ## p.score in testData
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, confounders, data = df, p.score = p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ X1 + X2 + X3 + p.score + z, data = df)"))
  
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  
  ## p.score supplied as literal
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, confounders, data = df, p.score = df$p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ X1 + X2 + X3 + ps + z, data = data)"))
  
  expect_false(identical(res$env, currentEnv))
  expect_false(identical(environment(res$call[[2L]]), currentEnv))
  expect_equal(res$env$data$ps, df$p.score)
  
  expect_equal(res$trt, "z")
})

test_that("response call with data as data.frame argument returns valid output", {
  res <- bartCause:::getResponseLiteralCall(stats::lm, testData$y, testData$z, testData$x)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + V1 + V2 + V3, data = df)"))
  expect_true(all(colnames(res$df) %in% c("z", "y", "V1", "V2", "V3")))
  
  expect_equal(res$trt, "z")
  
  colnames(testData$x) <- c("x1", "x2", "z")
  res <- bartCause:::getResponseLiteralCall(stats::lm, testData$y, testData$z, testData$x)
  expect_equal(res$call, str2lang("stats::lm(y ~ zz + x1 + x2 + z, data = df)"))
  expect_true(all(colnames(res$df) %in% c("zz", "y", "x1", "x2", "z")))
  
  colnames(testData$x) <- NULL
    
  res <- bartCause:::getResponseLiteralCall(stats::lm, testData$y, testData$z, testData$x, p.score = testData$p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + V1 + V2 + V3 + ps, data = df)"))
  expect_true(all(colnames(res$df) %in% c("z", "y", "V1", "V2", "V3", "ps")))
  
  expect_equal(res$df$ps, testData$p.score)
   
  expect_equal(res$trt, "z")
})

