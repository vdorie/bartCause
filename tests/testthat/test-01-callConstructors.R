context("bartc input")

source(system.file("common", "groupedData.R", package = "bartCause"))

test_that("treatment call with data as list argument returns valid output", {
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, data = testData)
  # check that it generates the correct call and that the environments point here
  expect_equal(res$call, str2lang("stats::glm(z ~ x, data = testData)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, data = testData, group.by = g, use.ranef = FALSE)
  expect_equal(res$call, str2lang("stats::glm(z ~ x + g, data = testData)"))
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  # function doesn't make sense, but allows us to run test without relevant package
  # installed 
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, data = testData, parametric = (1 | g), use.lmer = FALSE)
  expect_equal(res$call, str2lang("stats::glm(z ~ (1 | g) + bart(x), data = testData)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, data = testData, parametric = (1 | g), use.lmer = TRUE)
  expect_equal(res$call, str2lang("stats::glm(z ~ (1 | g) + x, data = testData)"))
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, data = testData, group.by = g, use.ranef = TRUE, use.lmer = FALSE)
  expect_equal(res$call, str2lang("stats::glm(z ~ x, data = testData)"))
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, x, data = testData, group.by = g, use.ranef = TRUE, use.lmer = TRUE)
  expect_equal(res$call, str2lang("stats::glm(z ~ x + (1 | g), data = testData)"))
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_error(bartCause:::getTreatmentDataCall(stats::glm, z, x, data = testData, group.by = g, parametric = (1 | g)))
})

test_that("treatment call with data as data.frame argument returns valid output", {
  df <- with(testData, data.frame(x, z, y, p.score, g))
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, X1 + X2 + X3, data = df)
  expect_true(res$call == str2lang("stats::glm(z ~ X1 + X2 + X3, data = df)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, X1, data = df, parametric = X2 + X3, use.lmer = FALSE)
  expect_true(res$call == str2lang("stats::glm(z ~ X2 + X3 + bart(X1), data = df)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, X1, data = df, parametric = X2 + X3, use.lmer = TRUE)
  expect_true(res$call == str2lang("stats::glm(z ~ X2 + X3 + X1, data = df)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, X1 + X2 + X3, data = df, parametric = (1 | g), use.lmer = FALSE)
  expect_true(res$call == str2lang("stats::glm(z ~ (1 | g) + bart(X1 + X2 + X3), data = df)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, X1 + X2 + X3, data = df, parametric = (1 | g), use.lmer = TRUE)
  expect_true(res$call == str2lang("stats::glm(z ~ (1 | g) + (X1 + X2 + X3), data = df)"))
  
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
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, confounders, data = df, parametric = (1 | g), use.lmer = FALSE)
  expect_true(res$call == str2lang("stats::glm(z ~ (1 | g) + bart(X1 + X2 + X3), data = df)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  res <- bartCause:::getTreatmentDataCall(stats::glm, z, confounders, data = df, parametric = (1 | g), use.lmer = TRUE)
  expect_true(res$call == str2lang("stats::glm(z ~ (1 | g) + (X1 + X2 + X3), data = df)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
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
  testData$Z <- model.matrix(~ -1 + as.factor(g), testData)
  colnames(testData$Z) <- NULL
  attr(testData$Z, "assign") <- NULL
  attr(testData$Z, "contrasts") <- NULL
  
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x)
  expect_equal(res$call, str2lang("stats::glm(z ~ V1 + V2 + V3, data = df)"))
  expect_true(all(colnames(res$df) %in% c("z", "V1", "V2", "V3")))
  
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x, parametric = testData$Z, use.lmer = FALSE)
  expect_equal(res$call, str2lang("stats::glm(z ~ V1 + V2 + V3 + bart(V1_bart + V2_bart + V3_bart), data = df)"))
  expect_true(all(colnames(res$df) %in% c("z", "V1", "V2", "V3", "V1_bart", "V2_bart", "V3_bart")))
  
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x, parametric = testData$Z, use.lmer = TRUE)
  expect_equal(res$call, str2lang("stats::glm(z ~ V1 + V2 + V3 + V1_bart + V2_bart + V3_bart, data = df)"))
  expect_true(all(colnames(res$df) %in% c("z", "V1", "V2", "V3", "V1_bart", "V2_bart", "V3_bart")))
  
  colnames(testData$x) <- c("x1", "x2", "z")
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x)
  expect_equal(res$call, str2lang("stats::glm(zz ~ x1 + x2 + z, data = df)"))
  expect_true(all(colnames(res$df) %in% c("zz", "x1", "x2", "z")))
  
  colnames(testData$Z) <- c(paste0("g_", seq_len(ncol(testData$Z) - 1L)), "x2")
  res <- bartCause:::getTreatmentLiteralCall(stats::glm, testData$z, testData$x, parametric = testData$Z, use.lmer = FALSE)
  expect_equal(res$call, str2lang("stats::glm(zz ~ g_1 + g_2 + x2 + bart(x1 + x2_bart + z), data = df)"))
  expect_true(all(colnames(res$df) %in% c("zz", "x1", "x2_bart", "z", "g_1", "g_2", "x2")))
  
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
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, x, parametric = (1 | g), data = testData)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + bart(x + z) + (1 | g), data = testData)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  ## p.score in testData
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, x, data = testData, p.score = p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ x + p.score + z, data = testData)"))
  
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, x, parametric = (1 | g),
                                         data = testData, p.score = p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + p.score + bart(x + z + p.score) + (1 | g), data = testData)"))
  
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  ## p.score supplied as literal
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, x, data = testData, p.score = testData$p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ x + ps + z, data = data)"))
  
  expect_false(identical(res$env, currentEnv))
  expect_false(identical(environment(res$call[[2L]]),currentEnv))
  expect_equal(res$env$data$ps, testData$p.score)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, x, parametric = (1 | g),
                                         data = testData, p.score = testData$p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + ps + bart(x + z + ps) + (1 | g), data = data)"))
  
  expect_false(identical(res$env, currentEnv))
  expect_false(identical(environment(res$call[[2L]]),currentEnv))
  expect_equal(res$env$data$ps, testData$p.score)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
})

test_that("response call with data as data.frame argument returns valid output", {
  df <- with(testData, data.frame(x, y, z, p.score, g))
  
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, X1 + X2 + X3, data = df)
  expect_equal(res$call, str2lang("stats::lm(y ~ X1 + X2 + X3 + z, data = df)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, x, parametric = (1 | g), data = df)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + bart(x + z) + (1 | g), data = df)"))
  
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  ## p.score in testData
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, X1 + X2 + X3, data = df, p.score = p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ X1 + X2 + X3 + p.score + z, data = df)"))
  
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, X1 + X2 + X3, parametric = (1 | g),
                                         data = df, p.score = p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + p.score + bart(X1 + X2 + X3 + z + p.score) + (1 | g), data = df)"))
  
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  ## p.score supplied as literal
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, X1 + X2 + X3, data = df, p.score = df$p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ X1 + X2 + X3 + ps + z, data = data)"))
  
  expect_false(identical(res$env, currentEnv))
  expect_false(identical(environment(res$call[[2L]]), currentEnv))
  expect_equal(res$env$data$ps, df$p.score)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, X1 + X2 + X3, parametric = (1 | g),
                                         data = df, p.score = df$p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + ps + bart(X1 + X2 + X3 + z + ps) + (1 | g), data = data)"))
  
  expect_false(identical(res$env, currentEnv))
  expect_false(identical(environment(res$call[[2L]]), currentEnv))
  expect_equal(res$env$data$ps, df$p.score)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  
  confounders <- "X1 + X2 + X3"
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, confounders, data = df)
  expect_equal(res$call, str2lang("stats::lm(y ~ X1 + X2 + X3 + z, data = df)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, confounders, parametric = (1 | g), data = df)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + bart(X1 + X2 + X3 + z) + (1 | g), data = df)"))
  
  currentEnv <- sys.frame(sys.nframe())
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  
  ## p.score in testData
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, confounders, data = df, p.score = p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ X1 + X2 + X3 + p.score + z, data = df)"))
  
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, confounders, parametric = (1 | g),
                                         data = df, p.score = p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + p.score + bart(X1 + X2 + X3 + z + p.score) + (1 | g), data = df)"))
  
  expect_identical(res$env, currentEnv)
  expect_identical(environment(res$call[[2L]]), currentEnv)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  
  ## p.score supplied as literal
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, confounders, data = df, p.score = df$p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ X1 + X2 + X3 + ps + z, data = data)"))
  
  expect_false(identical(res$env, currentEnv))
  expect_false(identical(environment(res$call[[2L]]), currentEnv))
  expect_equal(res$env$data$ps, df$p.score)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  res <- bartCause:::getResponseDataCall(stats::lm, y, z, confounders, parametric = (1 | g),
                                         data = df, p.score = df$p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + ps + bart(X1 + X2 + X3 + z + ps) + (1 | g), data = data)"))
  
  expect_false(identical(res$env, currentEnv))
  expect_false(identical(environment(res$call[[2L]]), currentEnv))
  expect_equal(res$env$data$ps, df$p.score)
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
})

test_that("response call with literal arguments retuns valid output", {
  res <- bartCause:::getResponseLiteralCall(stats::lm, testData$y, testData$z, testData$x)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + V1 + V2 + V3, data = df)"))
  expect_true(all(colnames(res$df) %in% c("z", "y", "V1", "V2", "V3")))
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  
  testData$Z <- model.matrix(~ -1 + as.factor(g), testData)
  colnames(testData$Z) <- NULL
  attr(testData$Z, "assign") <- NULL
  attr(testData$Z, "contrasts") <- NULL

  res <- bartCause:::getResponseLiteralCall(stats::lm, testData$y, testData$z, testData$x, parametric = testData$Z)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + V1 + V2 + V3 + bart(V1_bart + V2_bart + V3_bart + z), data = df)"))
  expect_true(all(colnames(res$df) %in% c("z", "y", "V1", "V2", "V3", "V1_bart", "V2_bart", "V3_bart")))
  
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  
  colnames(testData$x) <- c("x1", "x2", "z")
  res <- bartCause:::getResponseLiteralCall(stats::lm, testData$y, testData$z, testData$x)
  expect_equal(res$call, str2lang("stats::lm(y ~ zz + x1 + x2 + z, data = df)"))
  expect_true(all(colnames(res$df) %in% c("zz", "y", "x1", "x2", "z")))
  
  expect_equal(res$trt, "zz")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  res <- bartCause:::getResponseLiteralCall(stats::lm, testData$y, testData$z, testData$x, parametric = testData$Z)
  expect_equal(res$call, str2lang("stats::lm(y ~ zz + V1 + V2 + V3 + bart(x1 + x2 + z + zz), data = df)"))
  expect_true(all(colnames(res$df) %in% c("zz", "y", "V1", "V2", "V3", "x1", "x2", "z")))
  
  expect_equal(res$trt, "zz")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  
  colnames(testData$x) <- NULL
  
  res <- bartCause:::getResponseLiteralCall(stats::lm, testData$y, testData$z, testData$x, p.score = testData$p.score)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + V1 + V2 + V3 + ps, data = df)"))
  expect_true(all(colnames(res$df) %in% c("z", "y", "V1", "V2", "V3", "ps")))
  
  expect_equal(res$df$ps, testData$p.score)
   
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
  
  res <- bartCause:::getResponseLiteralCall(stats::lm, testData$y, testData$z, testData$x, p.score = testData$p.score, parametric = testData$Z)
  expect_equal(res$call, str2lang("stats::lm(y ~ z + ps + V1 + V2 + V3 + bart(V1_bart + V2_bart + 
    V3_bart + z + ps), data = df)"))
  expect_true(all(colnames(res$df) %in% c("z", "y", "V1", "V2", "V3", "ps", "V1_bart", "V2_bart", "V3_bart")))
  
  expect_equal(as.vector(res$df[,"ps"]), as.vector(testData$p.score))
   
  expect_equal(res$trt, "z")
  expect_true(length(res$missingRows) > 0L && !any(res$missingRows))
})

