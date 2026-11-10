summaryATEFunctions <- list(
  getPATEEstimate.tmle.nosamp = function(samples.tmle)
{
  est <- samples.tmle[["est"]]
  sd  <- samples.tmle[["se"]]
  
  c(estimate = est, sd = sd)
},

getPATEIntervalFunction.tmle.nosamp = function(ci.style)
{
  if (ci.style %in% c("quant", "hpd"))
    stop("ci style '", ci.style, "' requires full posterior distribution to be calculated")
  
  function(ci.probs, estimate)
    estimate[["estimate"]] + estimate[["sd"]] * qnorm(ci.probs)
},

getPATEEstimate.tmle.var.exp = function(samples.tmle, samples.icate, weights)
{
  samples.pate <- if (length(dim(samples.tmle)) == 2L) samples.tmle[,"est"] else as.vector(samples.tmle[,,"est"])
  
  n.obs     <- ncol(samples.icate)
  n.samples <- nrow(samples.icate)
  
  est <- mean(samples.pate)
  var.between.samples <- var(samples.pate)
  if (!is.null(weights))
    variance.samples <- apply(weights * t(samples.icate - samples.pate)^2, 2L, sum) * n.obs / sum(weights) / (n.obs - 1)
  else
    variance.samples <- apply(t(samples.icate - samples.pate)^2, 2L, sum) / (n.obs - 1)
  var.within.samples <- mean(variance.samples)
  
  c(estimate = est, sd = sqrt(var.between.samples + var.within.samples))
},

getPATEEstimate.tmle.psi = function(samples.tmle)
{
  samples.pate <- if (length(dim(samples.tmle)) == 2L) samples.tmle[,"est"] else as.vector(samples.tmle[,,"est"])
  
  est <- mean(samples.pate)
  var.between.samples <- var(samples.pate)
  variance.samples <- if (length(dim(samples.tmle)) == 2L) samples.tmle[,"se"]^2 else samples.tmle[,,"se"]^2
  var.within.samples <- mean(variance.samples)
  
  c(estimate = est, sd = sqrt(var.between.samples + var.within.samples))
},


getPATEIntervalFunction.tmle = function(ci.style)
{
  switch(ci.style,
         norm  = function(ci.probs, estimate) estimate[["estimate"]] + estimate[["sd"]] * qnorm(ci.probs),
         quant = function(ci.probs, samples.tmle) {
           samples.pate <- if (length(dim(samples.tmle)) == 2L) samples.tmle[,"est"] else as.vector(samples.tmle[,,"est"])
           unname(quantile(samples.pate, ci.probs))
         },
         hpd   = function(ci.level, samples.tmle) {
           samples.pate <- if (length(dim(samples.tmle)) == 2L) samples.tmle[,"est"] else as.vector(samples.tmle[,,"est"])
           getHPDInterval(samples.pate, ci.level)
         })
},

getPATEEstimate.bart.var.exp = function(samples.icate, weights)
{
  n.obs <- ncol(samples.icate)
  n.samples <- nrow(samples.icate)
  
  if (!is.null(weights)) {
    samples.cate <- apply(weights * t(samples.icate), 2L, sum)
    est <- mean(samples.cate)
    
    var.between.samples <- var(samples.cate)
    variance.samples <- apply(weights * t(samples.icate - samples.cate)^2, 2L, sum) * n.obs / sum(weights) / (n.obs - 1)
  } else {
    est <- mean(samples.icate)
    
    cate.samples <- apply(samples.icate, 1L, mean)
    var.between.samples <- var(cate.samples)
    # Each should be a sample of the variance across people
    variance.samples <- apply(samples.icate, 1L, var)
  }
  var.within.samples  <- mean(variance.samples)
  
  c(estimate = est, sd = sqrt(var.between.samples + var.within.samples))
},

## 'weights' are the population weights' shares of the rows averaged over, NULL
## for equal shares, and 'noise.var' the noise variance of that average in
## units of sigma^2 (see averageNoiseVariance); the observed and the
## counterfactual draw each add it.
getPATEEstimate.bart.ppd = function(samples.icate, weights, sigma, samples.obs, samples.cf, noise.var)
{
  if (!is.null(weights)) {
    samples.cate <- as.vector(samples.icate %*% weights)
  } else {
    samples.cate <- rowMeans(samples.icate)
  }
  est <- mean(samples.cate)
  
  sd <- sqrt(var(samples.cate) +
    if (!is.null(sigma))
      2 * mean(sigma^2) * noise.var
    else {
      samples.var <- samples.obs * (1 - samples.obs) + samples.cf * (1 - samples.cf)
      if (!is.null(weights))
        mean(samples.var %*% weights^2)
      else
        mean(apply(samples.var, 2L, mean) / ncol(samples.var))
    })
  c(estimate = est, sd = sd)
},

getPATEIntervalFunction.bart = function(ci.style)
{
  switch(ci.style,
         norm  = function(ci.probs, estimate) estimate[["estimate"]] + estimate[["sd"]] * qnorm(ci.probs),
         quant = function(ci.probs, samples.pate) unname(quantile(samples.pate, ci.probs)),
         hpd   = function(ci.level, samples.pate) getHPDInterval(samples.pate, ci.level))
},

## As getPATEEstimate.bart.ppd, with one draw per row; a row with a missing
## response ('missing.obs') adds the noise of its observed draw as well, in
## units of sigma^2 'missing.noise.var'.
getSATEEstimate.bart = function(samples.iscate, weights, sigma, samples.obs, samples.cf, noise.var, missing.obs,
                                missing.noise.var)
{
  if (!is.null(weights)) {
    samples.scate <- as.vector(samples.iscate %*% weights)
  } else {
    samples.scate <- rowMeans(samples.iscate)
  }
  est <- mean(samples.scate)
  
  sd <- sqrt(var(samples.scate) +
    if (!is.null(sigma)) {
      mean(sigma^2) * (noise.var + missing.noise.var)
    } else {
      samples.var <- samples.cf * (1 - samples.cf)
      if (any(missing.obs))
        samples.var[,missing.obs] <- samples.var[,missing.obs] + samples.obs[,missing.obs] * (1 - samples.obs[,missing.obs])
      if (!is.null(weights))
        mean(samples.var %*% weights^2)
      else
        mean(apply(samples.var, 2L, mean) / ncol(samples.var))
    })
  c(estimate = est, sd = sd)
},

getSATEIntervalFunction.bart = function(ci.style)
{
  switch(ci.style,
         norm  = function(ci.probs, estimate) estimate[["estimate"]] + estimate[["sd"]] * qnorm(ci.probs),
         quant = function(ci.probs, samples.sate) unname(quantile(samples.sate, ci.probs)),
         hpd   = function(ci.level, samples.sate) getHPDInterval(samples.sate, ci.level))
},

getCATEEstimate.bart = function(samples.icate, weights)
{
  if (!is.null(weights)) {
    samples.cate <- as.vector(samples.icate %*% weights)
  } else {
    samples.cate <- rowMeans(samples.icate)
  }
  c(estimate = mean(samples.cate), sd = sd(samples.cate))
},

getCATEIntervalFunction.bart = function(ci.style)
{
  switch(ci.style,
         norm  = function(ci.probs, estimate) estimate[["estimate"]] + estimate[["sd"]] * qnorm(ci.probs),
         quant = function(ci.probs, samples.cate) unname(quantile(samples.cate, ci.probs)),
         hpd   = function(ci.level, samples.cate) getHPDInterval(samples.cate, ci.level))
})

getATEEstimates <- function(object, target, ci.style, ci.level, pate.style)
{
  x <- NULL
  ci.probs <- evalx((1 - ci.level) / 2, c(x, 1 - x))
  
  # for R CMD check
  getPATEEstimate.tmle.nosamp <- getPATEIntervalFunction.tmle.nosamp <- getPATEEstimate.tmle.var.exp <-
    getPATEEstimate.tmle.psi <- getPATEIntervalFunction.tmle <- getPATEEstimate.bart.var.exp <-
    getPATEEstimate.bart.ppd <- getPATEIntervalFunction.bart <- getSATEEstimate.bart <- 
    getSATEIntervalFunction.bart <- getCATEEstimate.bart <- getCATEIntervalFunction.bart <- NULL
  for (functionName in names(summaryATEFunctions))
    assign(functionName, summaryATEFunctions[[functionName]])
  
  # this is annoyingly complicated in order to handle grouped data
  group.effects <- if (is.null(object$group.effects)) FALSE else object$group.effects
  groupEstimates <- !is.null(object$group.by) && group.effects
  if (target == "pate") {
    if (object$method.rsp %in% c("tmle", "p.weight")) {
      if ((!groupEstimates && is.null(dim(object$est))) ||
          ( groupEstimates && is.null(dim(object$est[[1L]]))))
      {
        getATEEstimate <- getPATEEstimate.tmle.nosamp
        getATEInterval <- getPATEIntervalFunction.tmle.nosamp(ci.style)
      } else {
        getATEEstimate <- if (pate.style == "var.exp") getPATEEstimate.tmle.var.exp else getPATEEstimate.tmle.psi
        getATEInterval <- getPATEIntervalFunction.tmle(ci.style)
      }
    } else {
      getATEEstimate <- if (pate.style == "var.exp") getPATEEstimate.bart.var.exp else getPATEEstimate.bart.ppd
      getATEInterval <- getPATEIntervalFunction.bart(ci.style)
    }
  } else if (target == "sate") {
    getATEEstimate <- getSATEEstimate.bart
    getATEInterval <- getSATEIntervalFunction.bart(ci.style)
  } else if (target == "cate") {
    getATEEstimate <- getCATEEstimate.bart
    getATEInterval <- getCATEIntervalFunction.bart(ci.style)
  }
  
  estimateVariables <- names(formals(getATEEstimate))
  intervalVariables <- names(formals(getATEInterval))
  
  estimateCall <- quote(getATEEstimate())
  intervalCall <- quote(getATEInterval())
  
  pop  <- populationWeights(object)
  fitw <- responseFitWeights(object)
  n <- length(object$trt)
  y <- observedResponse(object)
  missingObs <- is.na(y)
  
  inferentialSubset <- switch(object$estimand,
                              ate = rep(TRUE, n),
                              att = object$trt >  0,
                              atc = object$trt <= 0)
  keepSubset <- object$commonSup.sub
  
  if (responseIsBinary(object))
    estimateVariables <- setdiff(estimateVariables, "sigma")
  
  # load whatever we might need
  for (i in seq_along(estimateVariables)) {
    varName <- estimateVariables[i]
    assign(varName, switch(varName,
           samples.tmle  = object$est,
           samples.icate = extract(object, "icate", "all"),
           samples.cate  = extract(object, "cate"),
           samples.sate  = extract(object, "sate"),
           samples.iscate = {
             ## a missing response is the posterior of its observed outcome,
             ## its noise counted below
             iscate <- multiplyArrayByVec(subtractArrayFromVec(y, extract(object, "mu.cf", "all")), 2 * object$trt - 1)
             if (any(missingObs)) iscate[,missingObs] <- extract(object, "icate", "all")[,missingObs]
             iscate
           },
           sigma         = extract(object, "sigma"),
           samples.obs   = extract(object, "mu.obs", "all"),
           samples.cf    = extract(object, "mu.cf", "all"),
           weights       = NULL,
           noise.var     = NULL,
           missing.obs   = NULL,
           missing.noise.var = NULL
    ))
    estimateCall[[i + 1L]] <- str2lang(varName)
  }
  
  if (responseIsBinary(object) && "sigma" %in% names(formals(getATEEstimate)))
    estimateCall["sigma"] <- list(NULL)
  
  for (i in seq_along(intervalVariables)) {
    varName <- intervalVariables[i]
    switch(varName,
           samples.pate = samples.pate <- extract(object, "pate"),
           samples.sate = samples.sate <- extract(object, "sate"),
           samples.cate = samples.cate <- extract(object, "cate"))
    intervalCall[[i + 1L]] <- str2lang(varName)
  }
  
  rowVariables <- intersect(estimateVariables, c("samples.icate", "samples.obs", "samples.cf", "samples.iscate"))
  ## the estimate over a set of rows (a logical vector over every row): the
  ## population weights' shares of those rows, and the noise of their average
  ## in the response fit as fit
  estimateOver <- function(rows) {
    shares <- if (is.null(pop)) rep_len(1 / sum(rows), sum(rows)) else weightShares(pop, rows)
    fitRows <- if (is.null(fitw)) NULL else fitw[rows]
    noise.var <- if (is.null(pop) && is.null(fitw)) 1 / sum(rows) else averageNoiseVariance(shares, fitRows)
    ## the PATE's draws of a row with a missing response are both predictive
    ## already, so only the SATE adds the noise of its observed draw
    missing.obs <- missingObs[rows]
    missing.noise.var <- if (any(missing.obs)) averageNoiseVariance(shares * missing.obs, fitRows) else 0
    weights <- if (is.null(pop)) NULL else shares
    for (varName in rowVariables)
      assign(varName, get(varName)[,rows,drop = FALSE])
    estimate <- eval(estimateCall)
    interval <- eval(intervalCall)
    names(interval) <- c("ci.lower", "ci.upper")
    c(estimate, interval)
  }
  environment(estimateOver) <- environment()
  
  if (!groupEstimates) {
    estimates <- as.data.frame(t(estimateOver(inferentialSubset & keepSubset)))
    row.names(estimates) <- object$estimand
  } else {
    estimates <- as.data.frame(t(sapply(seq_along(levels(object$group.by)), function(i) {
      level <- levels(object$group.by)[i]
      levelObs <- object$group.by == level & inferentialSubset & keepSubset
      
      # go through and subset variables that are grouped, sometimes in lists
      # sometimes in vectors
      for (varName in c("samples.tmle", "samples.cate", "samples.pate", "samples.sate", "samples.scate"))
        if (varName %in% estimateVariables || varName %in% intervalVariables)
          assign(varName, get(varName)[[i]])
      environment(estimateOver) <- environment()
      
      ## a group of no population weight is no part of the population: NA,
      ## the other groups estimated
      if (!is.null(pop) && !(sum(pop[levelObs]) > 0))
        return(c(estimate = NA_real_, sd = NA_real_, ci.lower = NA_real_, ci.upper = NA_real_, n = sum(levelObs),
                 n.eff = NA_real_))
      c(estimateOver(levelObs), n = sum(levelObs), if (!is.null(pop)) c(n.eff = kishSampleSize(pop[levelObs])))
    })))
    row.names(estimates) <- levels(object$group.by)
    
    # combine back to get whole if possible
    if (object$method.rsp %in% c("bart", "bcf")) {
      ## a group's share of the whole is its share of the population weights,
      ## over the rows each estimate averages
      allObs <- inferentialSubset & keepSubset
      levelObs <- lapply(levels(object$group.by), function(level) object$group.by == level & allObs)
      weights.g <- if (is.null(pop)) estimates$n / sum(allObs)
                   else vapply(levelObs, function(rows) sum(pop[rows]), 0) / sum(pop[allObs])
      ## a group of no population weight adds nothing, its draws being NA;
      ## the whole of no weight is an error
      if (!is.null(pop)) weightShares(pop, allObs)
      inTotal <- if (is.null(pop)) rep_len(TRUE, length(weights.g)) else weights.g > 0
      for (varName in c("samples.cate", "samples.pate", "samples.sate"))
        if (varName %in% estimateVariables || varName %in% intervalVariables)
          assign(varName, colSums(t(sapply(get(varName)[inTotal], function(x) x)) * weights.g[inTotal]))
      
      totalName <- "total"
      while (totalName %in% rownames(estimates))
        totalName <- paste0(totalName, ".")
      estimates[totalName,] <- c(estimateOver(allObs), sum(allObs), if (!is.null(pop)) kishSampleSize(pop[allObs]))
    }
  }
  
  estimates
}

## Kish's effective sample size of a set of weights.
kishSampleSize <- function(weights) sum(weights)^2 / sum(weights^2)

## Kish's effective sample size of the population weights over the rows the
## estimate averages, NULL for an unweighted fit.
getEffectiveSampleSize <- function(object)
{
  pop <- populationWeights(object)
  if (is.null(pop)) return(NULL)
  rows <- switch(object$estimand,
                 ate = rep(TRUE, length(object$trt)),
                 att = object$trt >  0,
                 atc = object$trt <= 0) & object$commonSup.sub
  kishSampleSize(pop[rows])
}

summary.bartcFit <- function(object,
                             target = c("pate", "sate", "cate"),
                             ci.style = c("norm", "quant", "hpd"), ci.level = 0.95,
                             pate.style = c("ppd", "var.exp"),
                             ...)
{
  if (!is.character(target) || target[1L] %not_in% eval(formals(summary.bartcFit)$target))
    stop("target must be in '", paste0(eval(formals(summary.bartcFit)$target), collapse = "', '"), "'")
  target <- target[1L]
  
  if (!is.character(ci.style) || ci.style[1L] %not_in% eval(formals(summary.bartcFit)$ci.style))
    stop("ci.style must be in '", paste0(eval(formals(summary.bartcFit)$ci.style), collapse = "', '"), "'")
  ci.style <- ci.style[1L]
  
  if (!is.numeric(ci.level))
    stop("ci.level must be numeric")
  ci.level <- ci.level[1L]
  if (ci.level <= 0 || ci.level >= 1)
    stop("ci.level must be in (0, 1)")
  
  if (!is.character(pate.style) || pate.style[1L] %not_in% eval(formals(summary.bartcFit)$pate.style))
    stop("pate.style must be in '", paste0(eval(formals(summary.bartcFit)$pate.style), collapse = "', '"), "'")
  pate.style <- pate.style[1L]
  
  if (target != "pate" && object$method.rsp %in% c("tmle", "p.weight"))
    stop("target '", target, "' not supported for method '", object$method.rsp)
  
  
  ## the summary reads the draws several ways; they are read out of the fits
  ## once, into this call's copy of the fit, which then holds them as stored
  object <- withDraws(object, "mu")
  draw.dims <- dim(object$mu.hat.obs)
  n.samples <- if (length(draw.dims) > 2L) draw.dims[2L] else draw.dims[1L]
  n.obs     <- draw.dims[length(draw.dims)]
  
  
  estimates <- getATEEstimates(object, target, ci.style, ci.level, pate.style)
  
  result <-
    namedList(call       = object$call,
              method.rsp = object$method.rsp,
              method.trt = object$method.trt,
              ci.info = namedList(target, ci.style, ci.level, pate.style),
              n.obs,
              n.samples,
              n.chains   = object$n.chains,
              commonSup.rule = object$commonSup.rule,
              estimates)
  
  result$n.eff <- getEffectiveSampleSize(object)
  
  groupEstimates <- !is.null(object$group.by) && !is.null(object$group.effects) && object$group.effects
  if (object$commonSup.rule != "none") {
    result$commonSup.cut <- object$commonSup.cut
    if (!groupEstimates) {
      n.cut <- as.data.frame(t(c(trt = sum(object$trt & !object$commonSup.sub), ctl = sum(!object$trt & !object$commonSup.sub))))
      row.names(n.cut) <- ""
    } else {
      n.cut <- as.data.frame(t(sapply(seq_along(levels(object$group.by)), function(i) {
        level <- levels(object$group.by)[i]
        levelObs <- object$group.by == level
        c(trt = sum( object$trt[levelObs] & !object$commonSup.sub[levelObs]),
          ctl = sum(!object$trt[levelObs] & !object$commonSup.sub[levelObs]))
      })))
      row.names(n.cut) <- levels(object$group.by)
    }

    ## the share of the population weight that the cut removes from the rows
    ## the estimand averages over, since the cut changes the population
    ## the estimate describes
    pop <- populationWeights(object)
    if (!is.null(pop)) {
      estimandRows <- switch(object$estimand,
                             ate = rep_len(TRUE, length(object$trt)),
                             att = object$trt > 0,
                             atc = object$trt <= 0)
      total <- sum(pop[estimandRows])
      if (total > 0)
        result$weight.cut <- sum(pop[estimandRows & !object$commonSup.sub]) / total
      ## and for each group of a grouped fit, in the table of rows cut
      if (groupEstimates)
        n.cut$weight.cut <- vapply(levels(object$group.by), function(level) {
          rows <- estimandRows & object$group.by == level
          if (sum(pop[rows]) > 0) sum(pop[rows & !object$commonSup.sub]) / sum(pop[rows]) else NA_real_
        }, 0)
    }
    result$n.cut <- n.cut
  }
  class(result) <- "bartcFit.summary"
  
  result
}

print.bartcFit.summary <- function(x, ...)
{
  dotsList <- list(...)
  digits <- if (!is.null(dotsList$digits)) dotsList$digits else max(3L, getOption("digits") - 3L)
  
  if (!is.null(x$call))
    cat("Call: ", paste0(deparse(x$call), collapse = "\n        "), "\n\n", sep = "", ...)
  
  cat("Causal inference model fit by:\n",
      "  model.rsp: ", x$method.rsp, "\n",
      "  model.trt: ", x$method.trt, "\n\n", sep = "", ...)
  
  target <- switch(x$ci.info$target,
                   pate = "population average",
                   sate = "sample average",
                   cate = "conditional average")
  cat("Treatment effect (", target, "):\n", ..., sep = "")
  
  print(x$estimates, digits = digits, ...)
  if (!is.null(x$estimates[["n"]]) && any(x$estimates[["n"]] <= 10L) && x$n.obs > 10L)
    cat("  if (n < 10) group-size estimates may be unstable\n", ...)
  cat("Estimates fit from ", x$n.obs - if (x$commonSup.rule != "none") sum(x$n.cut[c("trt", "ctl")]) else 0, " total observations\n", sep = "")
  if (!is.null(x$n.eff) && is.null(x$estimates[["n.eff"]]))
    cat("Effective sample size of the population weights over the rows averaged (Kish): ",
        format(x$n.eff, digits = digits), "\n", sep = "")
  
  cat(round(100 * x$ci.info$ci.level, digits), "% credible interval calculated by: ",
      switch(x$ci.info$ci.style,
             norm = "normal approximation",
             quant = "empirical quantiles",
             hpd   = "highest posterior density"),
      "\n", sep = "", ...)
  if (x$ci.info$target == "pate" && x$ci.info$ci.style == "norm")
    cat("  population TE approximated by: ", if (x$ci.info$pate.style == "var.exp") "variance expansion" else "posterior predictive distribution", "\n", sep = "")
  
  if (x$commonSup.rule != "none") {
    cat("\nCommon support enforced by cutting using '", x$commonSup.rule, "' rule, cutoff: ", x$commonSup.cut,
        "\nSuppressed observations:\n", sep = "")
    print(x$n.cut)
    if (!is.null(x$weight.cut) && x$weight.cut > 0)
      cat("The rows cut carry ", format(100 * x$weight.cut, digits = digits),
          "% of the population weight; the estimate describes the weighted population that remains\n",
          sep = "", ...)
    cat("\n")
  }
  
  cat("Result based on ", x$n.samples, " posterior samples", 
      if (!is.null(x$n.chains) && x$n.chains > 1L) paste0(" times ", x$n.chains, " chains") else "",
      "\n", sep = "", ...)
  
  invisible(x)
}

print.bartcFit <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if (!is.null(x$call))
    cat("Call:\n", paste0(deparse(x$call), collapse = "\n  "), "\n\n", sep = "")
  
  if (x$method.rsp %in% c("p.weight", "tmle")) {
    target <- "pate"
    target_str <- "population average"
  } else {
    target <- "cate"
    target_str <- "conditional average"
  }
  cat("Treatment effect (", x$estimand, ", ", target_str, "): ", format(fitted(x, target), digits = digits), "\n", sep = "", ...)
  
  invisible(x)
}

