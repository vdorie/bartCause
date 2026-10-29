context("group.by written as an expression")

## bartCause 1.0-10 took group.by as a column name, a vector from the calling
## frame, or any expression of either, resolving it to a vector before the
## fit; a modeled grouping factor (use.ranef = TRUE) is fit by stan4bart, whose
## model frame sees only the columns of data, so the vector is handed to it as
## a column

set.seed(18)
n <- 60L
groupData <- data.frame(x1 = rnorm(n), x2 = rnorm(n), country = factor(sample(letters[1:4], n, TRUE)))
groupData$z <- rbinom(n, 1, pnorm(0.5 * groupData$x1))
groupData$y <- groupData$x1 + groupData$z + rnorm(n)

smallArgs <- list(n.samples = 10L, n.burn = 10L, n.chains = 1L, n.threads = 1L, n.trees = 10L,
                  verbose = FALSE, seed = 18L)

fitWith <- function(group.by, data = groupData, confounders = quote(x1 + x2), ...) {
  dat <- data
  eval(as.call(c(list(quote(bartCause::bartc), quote(y), quote(z), confounders, data = quote(dat),
                      group.by = group.by), smallArgs, list(...))))
}

## the grouping the sampler saw, by row
groupingOf <- function(fit, which = c("fit.rsp", "fit.trt")) {
  flist <- fit[[match.arg(which)]]$reTrms$flist
  expect_equal(length(flist), 1L)
  flist[[1L]]
}

test_that("group.by given as an expression of data reaches stan4bart as the vector 1.0-10 resolved", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")

  plain <- fitWith(quote(country))
  expect_true(inherits(plain$fit.rsp, "stan4bartFit"))
  expect_true(inherits(plain$fit.trt, "stan4bartFit"))

  same <- function(fit) {
    expect_true(inherits(fit$fit.rsp, "stan4bartFit"))
    expect_true(inherits(fit$fit.trt, "stan4bartFit"))
    expect_equal(fit$group.by, groupData$country)
    ## the same grouping reached both samplers, and so the same draws
    expect_equal(as.character(groupingOf(fit, "fit.rsp")), as.character(groupData$country))
    expect_equal(as.character(groupingOf(fit, "fit.trt")), as.character(groupData$country))
    expect_equal(fit$mu.hat.obs, plain$mu.hat.obs)
    expect_equal(fit$p.score, plain$p.score)
  }

  ## a column of data, in an expression
  same(fitWith(quote(as.factor(country))))

  ## an element of the data, named in the calling function
  elementOf <- function(dat, gname)
    bartc(response = y, treatment = z, confounders = x1 + x2, data = dat, group.by = dat[[gname]],
          n.samples = 10L, n.burn = 10L, n.chains = 1L, n.threads = 1L, n.trees = 10L, verbose = FALSE, seed = 18L)
  same(elementOf(groupData, "country"))

  ## a variable local to the calling function
  localVariable <- function(dat) {
    gid <- dat$country
    bartc(response = y, treatment = z, confounders = x1 + x2, data = dat, group.by = as.factor(gid),
          n.samples = 10L, n.burn = 10L, n.chains = 1L, n.threads = 1L, n.trees = 10L, verbose = FALSE, seed = 18L)
  }
  same(localVariable(groupData))

  ## a vector of the calling frame, with the data given by name
  fromFrame <- groupData$country
  same(bartc(y, z, x1 + x2, data = groupData, group.by = fromFrame,
             n.samples = 10L, n.burn = 10L, n.chains = 1L, n.threads = 1L, n.trees = 10L, verbose = FALSE, seed = 18L))

  ## a function of the caller's that reads a column, as the expression is read within data
  localFunction <- function(dat) {
    toFactor <- function(v) factor(v)
    bartc(response = y, treatment = z, confounders = x1 + x2, data = dat, group.by = toFactor(country),
          n.samples = 10L, n.burn = 10L, n.chains = 1L, n.threads = 1L, n.trees = 10L, verbose = FALSE, seed = 18L)
  }
  same(localFunction(groupData))
})

test_that("a group.by column of data keeps its name, and a resolved vector takes one nothing uses", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")

  plain <- fitWith(quote(country))
  expect_equal(names(plain$fit.rsp$reTrms$flist), "country")
  expect_equal(names(plain$fit.trt$reTrms$flist), "country")

  ## a column called g, a column called data, and a confounder g that the data
  ## does not hold: none is shadowed by the column the grouping factor rides in
  clash <- groupData
  clash$g <- runif(n)
  clash$data <- runif(n)
  fit <- fitWith(quote(as.factor(country)), data = clash, confounders = quote(x1 + g + data))
  expect_true(inherits(fit$fit.rsp, "stan4bartFit"))
  expect_equal(fit$fit.rsp$frame$g, clash$g)
  expect_equal(fit$fit.rsp$frame$data, clash$data)
  expect_equal(fit$fit.trt$frame$g, clash$g)
  expect_false(names(fit$fit.rsp$reTrms$flist) %in% c("g", "data"))
  expect_equal(as.character(groupingOf(fit)), as.character(clash$country))

  local <- function(dat) {
    g <- seq_len(nrow(dat))
    bartc(response = y, treatment = z, confounders = x1 + g, data = dat, group.by = as.factor(country),
          n.samples = 10L, n.burn = 10L, n.chains = 1L, n.threads = 1L, n.trees = 10L, verbose = FALSE, seed = 18L)
  }
  fit <- local(groupData)
  expect_equal(fit$fit.rsp$frame$g, seq_len(n))
  expect_equal(fit$fit.trt$frame$g, seq_len(n))

  ## a column shadowed in the calling frame by a different vector of the same name
  shadowed <- function(dat) {
    country <- rev(dat$country)
    bartc(response = y, treatment = z, confounders = x1 + x2, data = dat, group.by = country,
          n.samples = 10L, n.burn = 10L, n.chains = 1L, n.threads = 1L, n.trees = 10L, verbose = FALSE, seed = 18L)
  }
  fit <- shadowed(groupData)
  expect_equal(fit$group.by, factor(rev(as.character(groupData$country)), levels = levels(groupData$country)))
  expect_equal(as.character(groupingOf(fit)), rev(as.character(groupData$country)))
})

test_that("group.by as an expression with a glm treatment model and with a subset", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")
  skip_if_not_installed("lme4")

  fit <- fitWith(quote(as.factor(country)), method.trt = "glm")
  expect_true(inherits(fit$fit.trt, "glmerMod"))
  expect_equal(names(fit$fit.trt@flist), "g")
  expect_equal(as.character(fit$fit.trt@flist[[1L]]), as.character(groupData$country))
  
  ## the rows kept are the model frame's to choose, so the grouping column is
  ## handed over whole
  keep <- which(groupData$x1 > -1)
  fit <- fitWith(quote(as.factor(country)), method.trt = "none", subset = keep)
  expect_equal(as.character(groupingOf(fit, "fit.rsp")), as.character(groupData$country[keep]))
})

test_that("subset written as an expression of the columns of data fits with group.by, modeled or not", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")
  
  ## the expression is read within data, as the model frame reads it, and
  ## no variable of the calling frame is called x1
  keep <- groupData$x1 > -0.5
  for (use.ranef in c(TRUE, FALSE)) {
    label <- paste("use.ranef =", use.ranef)
    given <- fitWith(quote(country), subset = keep, use.ranef = use.ranef)
    fit <- fitWith(quote(country), subset = quote(x1 > -0.5), use.ranef = use.ranef)
    
    expect_equal(inherits(fit$fit.rsp, "stan4bartFit"), use.ranef, label = label)
    expect_equal(length(fit$p.score), sum(keep), label = label)
    expect_equal(fit$group.by, groupData$country[keep], label = label)
    ## the rows are the same as when the logical vector is given
    expect_equal(fit$group.by, given$group.by, label = label)
    expect_equal(fit$mu.hat.obs, given$mu.hat.obs, label = label)
    expect_equal(fit$p.score, given$p.score, label = label)
  }
  
  ## it is also read in the calling frame where data has no such column
  localCut <- function(dat) {
    threshold <- -0.5
    bartc(response = y, treatment = z, confounders = x1 + x2, data = dat, group.by = country, subset = x1 > threshold,
          n.samples = 10L, n.burn = 10L, n.chains = 1L, n.threads = 1L, n.trees = 10L, verbose = FALSE, seed = 18L)
  }
  expect_equal(localCut(groupData)$group.by, groupData$country[keep])
  
  ## and group.by as an expression, with a glm treatment model
  skip_if_not_installed("lme4")
  fit <- fitWith(quote(as.factor(country)), subset = quote(x1 > -0.5), method.trt = "glm")
  expect_equal(fit$group.by, groupData$country[keep])
})

test_that("subset written as an expression of the columns of data fits without group.by", {
  skip_on_cran()
  
  keep <- groupData$x1 > -0.5
  for (method.trt in c("bart", "glm", "none")) {
    given <- fitWith(NULL, subset = keep, method.trt = method.trt)
    fit <- fitWith(NULL, subset = quote(x1 > -0.5), method.trt = method.trt)
    expect_equal(nrow(fit$data.rsp@x), sum(keep), label = method.trt)
    expect_equal(fit$mu.hat.obs, given$mu.hat.obs, label = method.trt)
  }
})

test_that("a group.by that is not one value for each row is refused in each builder of bare vectors", {
  ## the treatment and the response model each build their data from bare
  ## vectors, with and without a parametric part, and each checks the length
  y <- groupData$y; z <- groupData$z
  x <- as.matrix(groupData["x1"]); parametric <- as.matrix(groupData["x2"])
  ## a length that does not divide the rows, which would otherwise be refused by the assignment alone
  short <- groupData$country[1:7]
  pattern <- "one value for each row of data (60), but its length is 7"
  
  fitBare <- function(...)
    bartc(y, z, x, group.by = short, n.samples = 10L, n.burn = 10L, n.chains = 1L, n.threads = 1L, n.trees = 10L,
          verbose = FALSE, ...)
  ## the propensity model first, the response model alone with method.trt = "none"
  expect_error(fitBare(), pattern, fixed = TRUE)
  expect_error(fitBare(method.trt = "none"), pattern, fixed = TRUE)
  expect_error(fitBare(parametric = parametric), pattern, fixed = TRUE)
  expect_error(fitBare(parametric = parametric, method.trt = "none"), pattern, fixed = TRUE)
  ## and the grouping entered as a fixed effect
  expect_error(fitBare(parametric = parametric, method.trt = "none", use.ranef = FALSE), pattern, fixed = TRUE)
})

test_that("predict reads the grouping column a resolved group.by was fit under", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")
  
  plain <- fitWith(quote(country), keepTrees = TRUE)
  resolved <- fitWith(quote(as.factor(country)), keepTrees = TRUE)
  expect_equal(names(plain$fit.rsp$reTrms$flist), "country")
  expect_equal(names(resolved$fit.rsp$reTrms$flist), "g")
  
  rows <- 1:12
  for (type in c("p.score", "mu.1", "mu.0", "icate"))
    expect_equal(predict(resolved, groupData[rows, c("x1", "x2")], group.by = groupData$country[rows], type = type),
                 predict(plain, groupData[rows, c("x1", "x2")], group.by = groupData$country[rows], type = type),
                 label = type)
})

test_that("use.ranef = FALSE reads group.by as a fixed effect, in whatever form it was written", {
  skip_on_cran()

  plain <- fitWith(quote(country), use.ranef = FALSE)
  expect_false(inherits(plain$fit.rsp, "stan4bartFit"))
  expect_equal(plain$group.by, groupData$country)

  elementOf <- function(dat, gname)
    bartc(response = y, treatment = z, confounders = x1 + x2, data = dat, group.by = dat[[gname]], use.ranef = FALSE,
          n.samples = 10L, n.burn = 10L, n.chains = 1L, n.threads = 1L, n.trees = 10L, verbose = FALSE, seed = 18L)
  fit <- elementOf(groupData, "country")
  expect_false(inherits(fit$fit.rsp, "stan4bartFit"))
  expect_equal(fit$group.by, groupData$country)
  expect_equal(fit$mu.hat.obs, plain$mu.hat.obs)

  fit <- fitWith(quote(as.factor(country)), use.ranef = FALSE)
  expect_equal(fit$mu.hat.obs, plain$mu.hat.obs)
})

test_that("a string given as group.by is refused with the column to pass", {
  ## 1.0-10 refused it for not being the length of data; naming the column says what to do instead
  expect_error(fitWith("country"), "'group.by' is the string \"country\", which names a column of data; pass the column itself, group.by = country, not its name", fixed = TRUE)
  expect_error(fitWith("country", use.ranef = FALSE), "pass the column itself, group.by = country", fixed = TRUE)
  expect_error(fitWith("country", method.trt = "glm"), "pass the column itself, group.by = country", fixed = TRUE)

  ## a string that names no column is a length error, as is any vector that is not one value per row
  expect_error(fitWith("elsewhere"), "'group.by' must have one value for each row of data (60), but its length is 1", fixed = TRUE)
  expect_error(fitWith(quote(country[1:5])), "one value for each row of data (60), but its length is 5", fixed = TRUE)

  ## and so is a bare vector, which has no data frame to name a column of
  x <- as.matrix(groupData[c("x1", "x2")])
  expect_error(bartc(groupData$y, groupData$z, x, group.by = "country", n.samples = 10L, n.burn = 10L, n.chains = 1L,
                     n.threads = 1L, n.trees = 10L, verbose = FALSE),
               "one value for each row of data (60), but its length is 1", fixed = TRUE)
})
