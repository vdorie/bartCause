getTreatmentEffectEstimate <- function(samples, ci.style, ci.level, subset = NULL) {
  if (!is.null(subset)) samples <- samples[subset]
  
  est <- mean(samples)
  sd <- sd(samples)
  probs <- evalx((1 - ci.level) / 2, c(x, 1 - x))
  
  ci <- switch(ci.style,
               norm  = est + sd(samples) * qnorm(probs),
               quant = unname(quantile(samples, probs)),
               hpd   = getHPDInterval(samples, ci.level))
  c(estimate = est, sd = sd, ci.lower = ci[1L], ci.upper = ci[2L])
}

summary.cibartFit <- function(object, ci.style = c("norm", "quant", "hpd"), ci.level = 0.95, ...)
{
  if (!is.character(ci.style))
    stop("ci.style must be in '", paste0(eval(formals(summary.cibartFit)$ci.style), collapse = "', '"), "'")
  ci.style <- ci.style[1L]
  
  if (!is.numeric(ci.level))
    stop("ci.level must be numeric")
  ci.level <- ci.level[1L]
  if (ci.level <= 0 || ci.level >= 1)
    stop("ci.level must be in (0, 1)")
  
  samples <- extract(object, "est", combineChains = TRUE)
  
  numSamples <- if (is.null(object$group.by)) ncol(samples[[1L]]) else length(samples)
  numObservations <- dim(object$samples.indiv.diff)[1L]
  
  if (is.null(object$group.by)) {
    estimates <- as.data.frame(t(getTreatmentEffectEstimate(samples, ci.style, ci.level)))
    row.names(estimates) <- object$estimand
  } else {
    estimates <- as.data.frame(t(sapply(seq_along(levels(object$group.by)), function(i) {
      level <- levels(object$group.by)[i]
      levelObs <- object$group.by == level
      c(getTreatmentEffectEstimate(samples[[i]], ci.style, ci.level, levelObs), n = sum(levelObs))
    })))
    row.names(estimates) <- paste0(object$estimand, ".", levels(object$group.by))
  }
  
  result <- list(call       = object$call,
                 method.rsp = object$method.rsp,
                 method.trt = object$method.trt,
                 ci.style   = ci.style,
                 ci.level   = ci.level,
                 numObservations = numObservations,
                 numSamples = numSamples,
                 estimates  = estimates,
                 n.chains   = object$n.chains)
  class(result) <- "cibartFit.summary"
  
  result
}

print.cibartFit.summary <- function(object, digits = max(3L, getOption("digits") - 3L), ...)
{
  dotsList <- list(...)
  if (!is.null(object$call))
    cat("Call: ", paste0(deparse(object$call), collapse = "\n      "), "\n\n", sep = "", ...)
  
  cat("Causal inference model fit by:\n",
      "  model.rsp: ", object$method.rsp, "\n",
      "  model.trt: ", object$method.trt, "\n\n", sep = "", ...)
  
  cat("Treatment effect:\n", ...)
  
  print(object$estimates, digits = digits, ...)
  if (!is.null(object$estimates[["n"]]) && any(object$estimates[["n"]] <= 10L) && object$numObservations > 10L)
    cat("  if (n < 10) group-size estimates may be unstable\n", ...)
  cat("\n")
  
  cat(round(100 * object$ci.level, digits), "% credible interval calculated by: ",
      switch(object$ci.style,
             norm = "normal approximation",
             quant = "empirical quantiles",
             hpd   = "highest posterior density"),
      "\n", sep = "", ...)
  cat("Result based on ", object$numSamples, " posterior samples", 
      if (!is.null(object$n.chains)) paste0(" across ", object$n.chains, " chains") else "",
      "\n", sep = "", ...)
  
  invisible(object)
}

print.cibartFit <- function(object, ...)
{
  cat("cibart fit using methods:\n",
      "  model.rsp: ", object$method.rsp, "\n",
      "  model.trt: ", object$method.trt, "\n\n", sep = "", ...)
  
  invisible(object)
}
