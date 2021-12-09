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
    cate.samples <- apply(weights * t(samples.icate), 2L, sum)
    est <- mean(cate.samples)
    
    var.between.samples <- var(cate.samples)
    variance.samples <- apply(weights * t(samples.icate - cate.samples)^2, 2L, sum) * n.obs / sum(weights) / (n.obs - 1)
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

getPATEEstimate.bart.ppd = function(samples.cate, weights, sigma, samples.obs, samples.cf, n.obs)
{
  est <- mean(samples.cate)
  
  sd <- sqrt(var(samples.cate) +
    if (!is.null(sigma))
      2 * mean(sigma^2) / n.obs
    else {
      samples.var <- samples.obs * (1 - samples.obs) + samples.cf * (1 - samples.cf)
      if (!is.null(weights))
        mean(apply(weights * samples.var, 2L, sum) / n.obs)
      else
        mean(apply(samples.var, 2L, mean) / n.obs)
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

getSATEEstimate.bart = function(samples.cate, weights, sigma, samples.cf, n.obs)
{
  est <- mean(samples.cate)
  
  sd <- sqrt(var(samples.cate) +
    if (!is.null(sigma))
      mean(sigma^2) / n.obs
    else {
      samples.var <- samples.cf * (1 - samples.cf)
      if (!is.null(weights))
        mean(apply(weights * samples.var, 2L, sum) / n.obs)
      else
        mean(apply(samples.var, 2L, mean) / n.obs)
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

getCATEEstimate.bart = function(samples.cate)
{
  c(estimate = mean(samples.cate), sd = sd(samples.cate))
},

getCATEIntervalFunction.bart = function(ci.style)
{
  switch(ci.style,
         norm  = function(ci.probs, estimate) estimate[["estimate"]] + estimate[["sd"]] * qnorm(ci.probs),
         quant = function(ci.probs, samples.sate) unname(quantile(samples.cate, ci.probs)),
         hpd   = function(ci.level, samples.sate) getHPDInterval(samples.cate, ci.level))
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
  
  weights <- if (inherits(object$fit.rsp, "stan4bartFit")) object$fit.rsp$weights else object$data.rsp@weights
  if (!is.null(weights) && length(weights) == 0L) weights <- NULL
  y <- if (inherits(object$fit.rsp, "stan4bartFit")) object$fit.rsp$y else object$data.rsp@y
  
  inferentialSubset <- switch(object$estimand,
                              ate = rep(TRUE, length(y)),
                              att = object$trt >  0,
                              atc = object$trt <= 0)
  inferentialSubset <- inferentialSubset
  keepSubset <- object$commonSup.sub
  n.obs <- sum(inferentialSubset & keepSubset)
  
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
           sigma         = extract(object, "sigma"),
           samples.obs   = extract(object, "mu.obs", "all"),
           samples.cf    = extract(object, "mu.cf", "all"),
           weights       = weights,
           n.obs = n.obs
    ))
    estimateCall[[i + 1L]] <- str2lang(varName)
  }
  
  if (responseIsBinary(object))
    estimateCall["sigma"] <- list(NULL)
  
  for (i in seq_along(intervalVariables)) {
    varName <- intervalVariables[i]
    switch(varName,
           samples.pate = samples.pate <- extract(object, "pate"))
    intervalCall[[i + 1L]] <- str2lang(varName)
  }
  
  if (!groupEstimates) {
    if (!is.null(weights)) {
      weights <- weights[inferentialSubset]
      weights <- weights / sum(weights)
    }
    for (varName in intersect(estimateVariables, c("samples.icate", "samples.obs", "samples.cf")))
      assign(varName, get(varName)[,keepSubset & inferentialSubset,drop = FALSE])
    
    estimate <- eval(estimateCall)
    interval <- eval(intervalCall)
    names(interval) <- c("ci.lower", "ci.upper")
    estimates <- as.data.frame(t(c(estimate, interval)))
    row.names(estimates) <- object$estimand
  } else {
    estimates <- as.data.frame(t(sapply(seq_along(levels(object$group.by)), function(i) {
      level <- levels(object$group.by)[i]
      levelObs <- object$group.by == level & inferentialSubset & keepSubset
      
      if (!is.null(weights)) {
        weights <- weights[levelObs]
        weights <- weights / sum(weights)
      }
      
      # go through and subset variables that are grouped, sometimes in lists
      # sometimes in vectors
      for (varName in c("samples.tmle", "samples.cate", "samples.pate"))
        if (varName %in% estimateVariables || varName %in% intervalVariables)
          assign(varName, get(varName)[[i]])
      for (varName in c("samples.icate", "samples.obs", "samples.cf")) {
        if (varName %in% estimateVariables || varName %in% intervalVariables) {
          var <- get(varName)
          assign(varName, ifelse_3(is.null(dim(var)), length(dim(var)) == 2L,
                                   var[levelObs], var[,levelObs,drop=FALSE], var[,,levelObs,drop=FALSE]))
        }
      }
      n.obs <- sum(levelObs)
      
      estimate <- eval(estimateCall)
      interval <- eval(intervalCall)
      names(interval) <- c("ci.lower", "ci.upper")
      c(estimate, interval, n = n.obs)
    })))
    row.names(estimates) <- levels(object$group.by)
    
    # combine back to get whole if possible
    if (object$method.rsp %in% "bart") {
      weights.g <- estimates$n / n.obs
      for (varName in c("samples.cate", "samples.pate"))
        if (varName %in% estimateVariables || varName %in% intervalVariables)
          assign(varName, colSums(t(sapply(get(varName), function(x) x)) * weights.g))
      estimate <- eval(estimateCall)
      interval <- eval(intervalCall)
      
      totalName <- "total"
      while (totalName %in% rownames(estimates))
        totalName <- paste0(totalName, ".")
      estimates[totalName,] <- c(estimate, interval, n.obs)
    }
  }
  
  estimates
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
  
  
  n.samples <- if (length(dim(object$mu.hat.obs)) > 2L) dim(object$mu.hat.obs)[2L] else dim(object$mu.hat.obs)[1L]
  n.obs     <- dim(object$mu.hat.obs)[length(dim(object$mu.hat.obs))] 
  
  
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
      result$estimates$n <- result$estimates$n - rowSums(n.cut)
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
  cat("Estimates fit from ", x$n.obs - if (x$commonSup.rule != "none") sum(x$n.cut) else 0, " total observations\n", sep = "")
  
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

