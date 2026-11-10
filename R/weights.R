## bartc takes two weight vectors. 'weights' are population weights: a row's
## share of the population the effects are averaged over, read only in those
## averages (and by a glm propensity model). 'weights.rsp' are the response
## fit's own weights, read by dbarts under its weights.type. Each is one value
## per row of data, read as a model frame reads its weights: within data first,
## then in the caller's frame. Each is checked against the rows of data and cut
## to the subset, so that it runs over every row an effect is estimated for,
## those with a missing response included. Returns NULL where the argument is
## absent or NULL, else a list of the values on the rows modeled ('value'), on
## every row of data ('full'), and the positions of the rows modeled ('rows').
resolveRowWeights <- function(matchedCall, argName, env)
{
  expr <- matchedCall[[argName]]
  if (is.null(expr)) return(NULL)

  data <- if (!is.null(matchedCall[["data"]])) eval(matchedCall[["data"]], env) else NULL
  value <- eval(expr, data, env)
  if (is.null(value)) return(NULL)

  if (!is.numeric(value) || (!is.null(dim(value)) && NCOL(value) != 1L))
    stop("'", argName, "' must be a numeric vector with one value for each row",
         if (argName == "weights") paste0("; for a survey design, pass weights(design), with the strata among ",
                                          "the confounders and the clusters in group.by") else "",
         call. = FALSE)
  value <- as.vector(value)

  n <- if (is.data.frame(data)) nrow(data) else NROW(eval(matchedCall[["response"]], data, env))
  if (length(value) != n)
    stop("'", argName, "' must have one value for each row (", n, "), but its length is ", length(value),
         call. = FALSE)

  full <- value
  rows <- seq_len(n)
  if (!is.null(matchedCall[["subset"]])) {
    rows <- normalizeBCFSubset(eval(matchedCall[["subset"]], data, env), n,
                               if (is.data.frame(data)) rownames(data) else NULL)
    value <- value[rows]
  }

  listRows <- function(bad) {
    bad <- rows[bad]
    paste0(if (length(bad) > 1L) "rows " else "row ",
           paste0(head(bad, 5L), collapse = ", "),
           if (length(bad) > 5L) paste0(" and ", length(bad) - 5L, " more") else "")
  }
  if (anyNA(value))
    stop("'", argName, "' is missing at ", listRows(is.na(value)),
         "; give every row a weight, or leave the rows out with subset", call. = FALSE)
  if (any(!is.finite(value) | value < 0))
    stop("'", argName, "' must be finite and non-negative, which ", listRows(!is.finite(value) | value < 0),
         if (sum(!is.finite(value) | value < 0) > 1L) " are not" else " is not", call. = FALSE)
  if (sum(value) <= 0)
    stop("'", argName, "' must be positive for at least one of the rows modeled", call. = FALSE)

  list(value = value, full = full, rows = rows)
}

## Whether the response, on the rows modeled, is 0/1 valued, which dbarts fits
## by probit when no family is given.
responseIsZeroOne <- function(matchedCall, env)
{
  data <- if (!is.null(matchedCall[["data"]])) eval(matchedCall[["data"]], env) else NULL
  y <- eval(matchedCall[["response"]], data, env)
  if (!is.null(matchedCall[["subset"]]))
    y <- y[normalizeBCFSubset(eval(matchedCall[["subset"]], data, env), NROW(y),
                              if (is.data.frame(data)) rownames(data) else NULL)]
  y <- as.vector(y)
  is.numeric(y) && all(y[!is.na(y)] %in% c(0, 1))
}

## A name for a variable that the fits read by name, taken by no column of data
## and by no variable the call names.
uniqueName <- function(stem, matchedCall, env)
{
  data <- if (!is.null(matchedCall[["data"]])) eval(matchedCall[["data"]], env) else NULL
  taken <- c(names(data), all.vars(matchedCall))
  name <- stem
  while (name %in% taken) name <- paste0(name, stem)
  name
}

## The population weights of a fit over every row an effect is estimated for,
## or NULL for equal weights. A fit made before bartc stored them, by bartCause
## 1.0-10, read its weights back from the response fit, which there had them
## as precisions; it is read as it was.
populationWeights <- function(object)
{
  if ("weights" %in% names(object)) return(object[["weights"]])

  weights <- if (inherits(object$fit.rsp, "stan4bartFit")) object$fit.rsp$weights else object$data.rsp@weights
  if (length(weights) == 0L) NULL else weights
}

## Whether the population weights were stored over every row, which a fit made
## before bartc stored them does not have for its missing-response rows.
populationWeightsAreStored <- function(object)
  "weights" %in% names(object)

## The response fit's own weights over every row an effect is estimated for, or
## NULL for an unweighted fit. A fit made before bartc stored them had one
## vector in both jobs, the response fit's precisions.
responseFitWeights <- function(object)
  if ("weights.rsp" %in% names(object)) object[["weights.rsp"]] else populationWeights(object)

## Whether the response fit reads its weights as counts of identical rows
## (dbarts's weights.type = "frequency"); stan4bart reads them as precisions.
responseFitWeightsAreCounts <- function(object)
{
  fit <- object$fit.rsp
  if (inherits(fit, "stan4bartFit") || is.null(fit)) return(FALSE)
  identical(fit[["weights.type"]], "frequency")
}

## The factor each row's noise standard deviation is sigma times, in the
## response fit as fit: a precision weight w gives sigma / sqrt(w), and a row
## of weight 0 has no noise scale, so its draws are NA. A count c is c
## identical rows, each drawn as one unit at sigma; where the draws are of a
## row's average over its units ('average'), that average is at
## sigma / sqrt(c). NULL where every row is at sigma.
noiseScale <- function(object, average = FALSE)
{
  weights <- responseFitWeights(object)
  if (is.null(weights)) return(NULL)

  scale <- if (responseFitWeightsAreCounts(object) && !average) rep_len(1, length(weights)) else 1 / sqrt(weights)
  scale[weights <= 0] <- NA_real_
  scale
}

## The noise variance, in units of sigma^2, of a weighted average of rows
## whose shares are 'shares' (summing to one), each row's noise in the
## response fit as fit: sum(a_i^2 / w_i). Kish's 1 / n_eff for an unweighted
## fit, it is free of the population weights' scale, and of the fit weights'
## too when they are precisions, which sigma's units follow. Rows of no share
## add nothing.
averageNoiseVariance <- function(shares, fitWeights)
{
  if (is.null(fitWeights)) return(sum(shares^2))
  rows <- shares > 0
  sum(shares[rows]^2 / fitWeights[rows])
}

## Whether a group of a grouped fit has no population weight over the rows
## its estimand averages ('weights' and 'trt' on the group's rows), so that
## its estimate is NA, the other groups estimated.
groupHasNoWeight <- function(weights, trt, estimand)
{
  if (is.null(weights)) return(FALSE)
  rows <- switch(estimand[1L], ate = rep_len(TRUE, length(trt)), att = trt > 0, atc = trt <= 0)
  !(sum(weights[rows]) > 0)
}

## The treatment on the rows modeled, read from its expression as written
## before any fit, or NULL where it cannot be read as a 0/1 vector one value for
## each row of data; the fits then find any problem with it.
treatmentBeforeFit <- function(matchedCall, env, n, rows)
{
  tryCatch({
    data <- if (!is.null(matchedCall[["data"]])) eval(matchedCall[["data"]], env) else NULL
    z <- eval(matchedCall[["treatment"]], data, env)
    if (is.logical(z)) z <- as.numeric(z)
    if (!is.numeric(z) || !is.null(dim(z)) || length(z) != n || anyNA(z)) NULL else as.vector(z)[rows]
  }, error = function(e) NULL)
}

## An estimand of no population weight over the rows it averages is an error.
stopIfNoWeight <- function(weights, trt, estimand)
{
  if (groupHasNoWeight(weights, trt, estimand))
    stop("the population weights are 0 on every row the estimate averages over", call. = FALSE)
}

## The estimates of a grouped propensity-weighted fit, NULL for the groups of
## no weight, given NA draws in the shape of the others; every group of no
## weight is an estimand of no weight, an error.
fillEmptyGroups <- function(est)
{
  empty <- vapply(est, is.null, NA)
  if (all(empty))
    stop("the population weights are 0 on every row the estimate averages over", call. = FALSE)
  template <- est[[which(!empty)[1L]]]
  template[] <- NA_real_
  est[empty] <- list(template)
  est
}

## The shares of the population weights over 'rows', a logical vector over
## every row; NULL weights are equal.
weightShares <- function(weights, rows)
{
  if (is.null(weights)) return(rep_len(1 / sum(rows), sum(rows)))
  weights <- weights[rows]
  total <- sum(weights)
  if (!(total > 0))
    stop("the population weights are 0 on every row the estimate averages over", call. = FALSE)
  weights / total
}

## The response over every row an effect is estimated for, NA where it is
## missing; the response fit holds the complete rows only.
observedResponse <- function(object)
{
  y <- if (inherits(object$fit.rsp, "stan4bartFit")) object$fit.rsp$y else object$data.rsp@y
  n <- length(object$trt)
  if (length(y) == n || !any(object$missingRows)) return(y)
  result <- rep_len(NA_real_, n)
  result[!object$missingRows] <- y
  result
}
