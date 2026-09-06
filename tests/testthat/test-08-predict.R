context("predict")

source(system.file("common", "linearData.R", package = "bartCause"), local = TRUE)

n.train <- 80L
x <- testData$x[seq_len(n.train),]
y <- testData$y[seq_len(n.train)]
z <- testData$z[seq_len(n.train)]

x.new <- testData$x[seq.int(n.train + 1L, nrow(testData$x)),]
n.test <- nrow(x.new)

test_that("predict gives sane results", {
  n.samples <- 7L
  n.chains  <- 2L
  fit <- bartc(y, z, x, method.trt = "glm", method.rsp = "bart",
               n.chains = n.chains, n.threads = 1L, n.burn = 0L, n.samples = n.samples, n.trees = 13L,
               keepTrees = TRUE,
               verbose = FALSE)
  
  # check predict for single row
  expect_equal(length(predict(fit, x.new[1,], type = "mu.0")), n.samples * n.chains)
  
  p.score <- predict(fit, x.new, type = "p.score")
  mu.1    <- predict(fit, x.new, type = "mu.1", combineChains = FALSE)
  mu.0    <- predict(fit, x.new, type = "mu.0", combineChains = TRUE)
  icate   <- predict(fit, x.new, type = "icate", combineChains = TRUE)
  
  expect_true(is.null(dim(p.score)))
  expect_equal(dim(mu.1), c(n.chains, n.samples, n.test))
  expect_equal(dim(mu.0), c(n.chains * n.samples, n.test))
  expect_equal(as.vector(icate), as.vector(matrix(aperm(mu.1, c(2L, 1L, 3L)), n.samples * n.chains)) - as.vector(mu.0))
})

test_that("predict results matches training data", {
  n.samples <- 7L
  n.chains  <- 2L
  fit <- bartc(y, z, x, method.trt = "bart", method.rsp = "bart",
               n.chains = n.chains, n.threads = 1L, n.burn = 0L, n.samples = n.samples, n.trees = 13L,
               keepTrees = TRUE,
               args.trt = list(k = 1.5), verbose = FALSE)
  
  p.score <- extract(fit, type = "p.score")
  mu.1    <- extract(fit, type = "mu.1")
  mu.0    <- extract(fit, type = "mu.0")
  icate   <- extract(fit, type = "icate")
  mu      <- extract(fit, type = "mu.obs")
  
  p.score.new <- predict(fit, x, type = "p.score")
  mu.1.new    <- predict(fit, x, type = "mu.1")
  mu.0.new    <- predict(fit, x, type = "mu.0")
  icate.new   <- predict(fit, x, type = "icate")
  mu.new      <- predict(fit, cbind(x, z), type = "mu")
  
  expect_equal(p.score, p.score.new)
  expect_equal(mu.0, mu.0.new)
  expect_equal(mu.1, mu.1.new)
  expect_equal(icate, icate.new)
  expect_equal(mu, mu.new)
})

set.seed(22)
g <- sample(3L, nrow(x), replace = TRUE)

n.samples <- 7L
n.chains  <- 2L

test_that("predict works with grouped data, glm trt model", {
 
  fit <- bartc(y, z, x, method.trt = "glm", method.rsp = "bart", group.by = g,
               n.chains = n.chains, n.threads = 1L, n.burn = 0L, n.samples = n.samples, n.trees = 13L,
               keepTrees = TRUE, use.ranef = FALSE,
               args.trt = list(k = 1.5), verbose = FALSE)
  
  p.score <- fitted(fit, type = "p.score")
  mu.1    <- extract(fit, type = "mu.1")
  mu.0    <- extract(fit, type = "mu.0")
  icate   <- extract(fit, type = "icate")
  
  p.score.new <- predict(fit, x, group.by = g, type = "p.score")
  mu.1.new    <- predict(fit, x, group.by = g, type = "mu.1")
  mu.0.new    <- predict(fit, x, group.by = g, type = "mu.0")
  icate.new   <- predict(fit, x, group.by = g, type = "icate")
  
  expect_equal(p.score, p.score.new)
  expect_equal(mu.0, mu.0.new)
  expect_equal(mu.1, mu.1.new)
  expect_equal(icate, icate.new)
})

test_that("predict works with grouped data, glmer trt model", {
  skip_if_not_installed("lme4")

  suppressWarnings(
    fit <- bartc(y, z, x, method.trt = "glm", method.rsp = "bart", group.by = g,
                 n.chains = n.chains, n.threads = 1L, n.burn = 0L, n.samples = n.samples, n.trees = 13L,
                 keepTrees = TRUE, use.ranef = FALSE,
                 args.trt = list(k = 1.5), verbose = FALSE)
  )
  
  p.score <- fitted(fit, type = "p.score")
  mu.1    <- extract(fit, type = "mu.1")
  mu.0    <- extract(fit, type = "mu.0")
  icate   <- extract(fit, type = "icate")
  
  p.score.new <- predict(fit, x, group.by = g, type = "p.score")
  mu.1.new    <- predict(fit, x, group.by = g, type = "mu.1")
  mu.0.new    <- predict(fit, x, group.by = g, type = "mu.0")
  icate.new   <- predict(fit, x, group.by = g, type = "icate")
  
  expect_equal(p.score, p.score.new)
  expect_equal(mu.0, mu.0.new)
  expect_equal(mu.1, mu.1.new)
  expect_equal(icate, icate.new)
})

test_that("predict works with grouped data, bart trt model", {
  skip_if_not_installed("stan4bart")
  fit <- bartc(y, z, x, method.trt = "bart", method.rsp = "bart", group.by = g,
               chains = 2L, iter = 14L, warmup = 7L,
               bart_args = list(keepTrees = TRUE, n.trees = 13L),
               verbose = FALSE)
  
  p.score <- extract(fit, type = "p.score")
  mu.1    <- extract(fit, type = "mu.1")
  mu.0    <- extract(fit, type = "mu.0")
  icate   <- extract(fit, type = "icate")
  
  p.score.new <- predict(fit, x, group.by = g, type = "p.score")
  mu.1.new    <- predict(fit, x, group.by = g, type = "mu.1")
  mu.0.new    <- predict(fit, x, group.by = g, type = "mu.0")
  icate.new   <- predict(fit, x, group.by = g, type = "icate")
  
  expect_equal(p.score, p.score.new)
  expect_equal(mu.0, mu.0.new)
  expect_equal(mu.1, mu.1.new)
  expect_equal(icate, icate.new)
  
  fit <- bartc(y, z, x, method.trt = "bart", method.rsp = "bart", group.by = g,
               n.chains = n.chains, n.threads = 1L, n.burn = 0L, n.samples = n.samples, n.trees = 13L,
               keepTrees = TRUE, use.ranef = FALSE,
               args.trt = list(k = 1.5), verbose = FALSE)
  
  p.score <- extract(fit, type = "p.score")
  mu.1    <- extract(fit, type = "mu.1")
  mu.0    <- extract(fit, type = "mu.0")
  icate   <- extract(fit, type = "icate")
  
  p.score.new <- predict(fit, x, group.by = g, type = "p.score")
  mu.1.new    <- predict(fit, x, group.by = g, type = "mu.1")
  mu.0.new    <- predict(fit, x, group.by = g, type = "mu.0")
  icate.new   <- predict(fit, x, group.by = g, type = "icate")
  
  expect_equal(p.score, p.score.new)
  expect_equal(mu.0, mu.0.new)
  expect_equal(mu.1, mu.1.new)
  expect_equal(icate, icate.new)
})

test_that("predict handles response types y, y.0, y.1, and ite", {
  n.samples <- 7L
  n.chains  <- 2L
  fit <- bartc(y, z, x, method.trt = "glm", method.rsp = "bart",
               n.chains = n.chains, n.threads = 1L, n.burn = 0L, n.samples = n.samples, n.trees = 13L,
               keepTrees = TRUE, verbose = FALSE)

  mu   <- predict(fit, cbind(x.new, z = 1), type = "mu",   combineChains = FALSE)
  mu.0 <- predict(fit, x.new, type = "mu.0", combineChains = FALSE)
  mu.1 <- predict(fit, x.new, type = "mu.1", combineChains = FALSE)
  sigma <- extract(fit, "sigma", combineChains = FALSE)

  set.seed(101)
  y.pred <- predict(fit, cbind(x.new, z = 1), type = "y", combineChains = FALSE)
  expect_equal(dim(y.pred), c(n.chains, n.samples, n.test))

  # regression check: predict(type = "y") used to reference an undefined/wrong
  # 'y' instead of the just-computed 'mu' and would error or silently reuse
  # whatever 'y' happened to be visible in the calling frame
  set.seed(101)
  sigma.rep <- rep_len(sigma, length(sigma) * n.test)
  epsilon <- rnorm(length(sigma.rep), 0, sigma.rep)
  dim(epsilon) <- dim(mu)
  expect_equal(y.pred, mu + epsilon)

  y.0 <- predict(fit, x.new, type = "y.0", combineChains = FALSE)
  y.1 <- predict(fit, x.new, type = "y.1", combineChains = FALSE)
  expect_equal(dim(y.0), dim(mu.0))
  expect_equal(dim(y.1), dim(mu.1))

  # a single "ite" call draws y.0 then y.1 from the ppd in that order; splitting
  # the same draws across two calls with a continued RNG stream must reproduce it
  set.seed(303)
  ite <- predict(fit, x.new, type = "ite", combineChains = FALSE)
  set.seed(303)
  y.0.split <- predict(fit, x.new, type = "y.0", combineChains = FALSE)
  y.1.split <- predict(fit, x.new, type = "y.1", combineChains = FALSE)
  expect_equal(ite, y.1.split - y.0.split)
})

test_that("predict enforces method/keepTrees preconditions", {
  fit.pweight <- bartc(y, z, x, method.trt = "bart", method.rsp = "p.weight",
                       n.chains = 1L, n.threads = 1L, n.burn = 0L, n.samples = 7L, n.trees = 13L,
                       verbose = FALSE)
  expect_error(predict(fit.pweight, x.new, type = "mu.0"), "requires method.rsp == 'bart'")

  fit.nokeep <- bartc(y, z, x, method.trt = "glm", method.rsp = "bart",
                      n.chains = 1L, n.threads = 1L, n.burn = 0L, n.samples = 7L, n.trees = 13L,
                      verbose = FALSE)
  expect_error(predict(fit.nokeep, x.new, type = "mu.0"), "keepTrees == TRUE")

  fit.none <- bartc(y, z, x, method.trt = "none", method.rsp = "bart", estimand = "att",
                    n.chains = 1L, n.threads = 1L, n.burn = 0L, n.samples = 7L, n.trees = 13L,
                    keepTrees = TRUE, verbose = FALSE)
  expect_error(predict(fit.none, x.new, type = "p.score"), "requires method.trt to specify a model")

  expect_error(predict(fit.none, x.new, type = "not-a-type"), "type must be in")
})

test_that("predict resolves the propensity score column from the fit, not by name", {
  # A confounder named to defeat the old name ladder: the response builders take
  # "ps" for the score, since "ps" is free, while a ladder that walks the design
  # names by stem lands on the confounder "psps" instead. Without the recorded
  # name, predict wrote the predicted scores over the confounder's column, and
  # the real score column never received them.
  x.collide <- x
  colnames(x.collide) <- c("psps", "x2", "x3")

  set.seed(22)
  fit <- bartc(y, z, x.collide, method.trt = "bart", method.rsp = "bart",
               n.chains = 2L, n.threads = 1L, n.burn = 3L, n.samples = 7L, n.trees = 13L,
               keepTrees = TRUE, verbose = FALSE)

  # the design carries both columns, and only one of them is the score
  expect_true(all(c("psps", "ps") %in% colnames(fit$data.rsp@x)))
  expect_equal(fit$name.p.score, "ps")

  # the ladder predict used to run, verbatim, resolves the confounder instead
  ladderName <- "ps"
  predictors.rsp <- colnames(fit$data.rsp@x)
  while (any(startsWith(predictors.rsp, ladderName)) &&
         paste0(ladderName, "ps") %in% predictors.rsp) ladderName <- paste0(ladderName, "ps")
  expect_equal(ladderName, "psps")
  expect_false(identical(ladderName, fit$name.p.score))

  newdata <- as.data.frame(x.collide[seq_len(5L),])
  scores  <- apply(predict(fit, newdata, type = "p.score", combineChains = FALSE), 3L, mean)

  # what predict must be doing: install the scores in the score's own column and
  # leave the confounder alone
  correct <- newdata
  correct[[fit$name.p.score]] <- scores
  correct[[fit$name.trt]] <- 1
  expect_equal(predict(fit, newdata, type = "mu.1", combineChains = FALSE),
               predict(fit$fit.rsp, correct, combineChains = FALSE))

  # and the placement is observable: writing the scores over the confounder as
  # well moves the prediction, so the assertion above can fail
  destroyed <- correct
  destroyed[["psps"]] <- scores
  expect_false(isTRUE(all.equal(predict(fit$fit.rsp, correct, combineChains = FALSE),
                                predict(fit$fit.rsp, destroyed, combineChains = FALSE))))
})

rm(testData, n.train, x, y, z, g, n.samples, n.chains, x.new, n.test)
