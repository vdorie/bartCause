context("print methods")

source(system.file("common", "linearData.R", package = "bartCause"))

test_that("print.bartcFit shows the call and matches the conditional average estimate", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData, method.trt = "glm", method.rsp = "bart", verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L)

  out <- capture.output(print(fit))
  expect_true(any(grepl("^Call:", out)))
  expect_true(any(grepl("bartc\\(", paste(out, collapse = " "))))

  digits <- max(3L, getOption("digits") - 3L)
  expected <- format(fitted(fit, "cate"), digits = digits)
  expect_true(any(grepl(paste0("Treatment effect \\(ate, conditional average\\): ", expected), out, fixed = FALSE)))
})

test_that("print.bartcFit uses population average for p.weight/tmle methods", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "p.weight", verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L)

  out <- capture.output(print(fit))
  digits <- max(3L, getOption("digits") - 3L)
  expected <- format(fitted(fit, "pate"), digits = digits)
  expect_true(any(grepl("population average", out)))
  expect_true(any(grepl(expected, out, fixed = TRUE)))
})

test_that("print.bartcFit.summary reports call, method, estimates table, and sample counts", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", estimand = "att", verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 2L, n.threads = 1L)
  sum <- summary(fit, target = "cate")
  out <- capture.output(print(sum))

  expect_true(any(grepl("^Call:", out)))
  expect_true(any(grepl("model.rsp:\\s*bart", out)))
  expect_true(any(grepl("model.trt:\\s*bart", out)))
  expect_true(any(grepl("conditional average", out)))

  digits <- max(3L, getOption("digits") - 3L)
  expect_true(any(grepl(format(round(sum$estimates$estimate, digits), nsmall = 0), out)))

  expect_true(any(grepl(paste0("Estimates fit from ", sum$n.obs, " total observations"), out)))
  expect_true(any(grepl("normal approximation", out)))
  expect_true(any(grepl(paste0(sum$n.samples, " posterior samples"), out)))
  expect_true(any(grepl(paste0("times ", sum$n.chains, " chains"), out)))
})

test_that("print.bartcFit.summary omits the population-TE line for non-norm ci styles", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L)

  out.norm  <- capture.output(print(summary(fit, ci.style = "norm")))
  out.quant <- capture.output(print(summary(fit, ci.style = "quant")))

  expect_true(any(grepl("population TE approximated by", out.norm)))
  expect_false(any(grepl("population TE approximated by", out.quant)))
  expect_true(any(grepl("empirical quantiles", out.quant)))
})

test_that("print.bartcFit.summary reports common support cutoffs and suppressed counts", {
  set.seed(22)
  fit <- bartc(y, z, x, data = testData, method.trt = "bart", method.rsp = "bart", verbose = FALSE,
               commonSup.rule = "sd", n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L)
  sum <- summary(fit)
  out <- capture.output(print(sum))

  expect_true(any(grepl("Common support enforced by cutting using 'sd' rule", out)))
  expect_true(any(grepl("Suppressed observations:", out)))
  expect_equal(sum(sum$n.cut), sum(fit$trt & !fit$commonSup.sub) + sum(!fit$trt & !fit$commonSup.sub))
})

test_that("print.bartcFit.summary warns about unstable small-group estimates", {
  set.seed(1)
  n <- 60L
  xx <- matrix(rnorm(3 * n), n, 3)
  zz <- rbinom(n, 1, 0.5)
  yy <- 2 * zz + xx[,1] + rnorm(n)
  gg <- factor(c(rep("small", 5L), rep("big", n - 5L)))

  fit <- bartc(yy, zz, xx, group.by = gg, group.effects = TRUE, use.ranef = FALSE,
               method.trt = "bart", method.rsp = "bart", verbose = FALSE,
               n.burn = 3L, n.samples = 13L, n.trees = 7L, n.chains = 1L, n.threads = 1L)
  sum <- summary(fit)
  expect_true(any(sum$estimates$n <= 10L))

  out <- capture.output(print(sum))
  expect_true(any(grepl("group-size estimates may be unstable", out)))
})
