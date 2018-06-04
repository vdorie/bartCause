getTreatmentEffectEstimate <- function(samples, ci.style, ci.level) {
  x <- NULL ## R CMD check
  
  est <- mean(samples)
  sd <- sd(samples)
  probs <- evalx((1 - ci.level) / 2, c(x, 1 - x))
  
  ci <- switch(ci.style,
               norm  = est + sd * qnorm(probs),
               quant = unname(quantile(samples, probs)),
               hpd   = getHPDInterval(samples, ci.level))
  c(estimate = est, sd = sd, ci.lower = ci[1L], ci.upper = ci[2L])
}

summary.bartcFit <- function(object, ci.style = c("norm", "quant", "hpd"), ci.level = 0.95, ...)
{
  if (!is.character(ci.style))
    stop("ci.style must be in '", paste0(eval(formals(summary.bartcFit)$ci.style), collapse = "', '"), "'")
  ci.style <- ci.style[1L]
  
  if (!is.numeric(ci.level))
    stop("ci.level must be numeric")
  ci.level <- ci.level[1L]
  if (ci.level <= 0 || ci.level >= 1)
    stop("ci.level must be in (0, 1)")
  
  samples <- extract(object, "est", combineChains = TRUE)
  
  numSamples <- if (!is.null(object$group.by)) length(samples[[1L]]) else length(samples)
  numObservations <- dim(object$samples.indiv.diff)[1L]
  
  if (is.null(object$group.by)) {
    estimates <- as.data.frame(t(getTreatmentEffectEstimate(samples, ci.style, ci.level)))
    row.names(estimates) <- object$estimand
  } else {
    estimates <- as.data.frame(t(sapply(seq_along(levels(object$group.by)), function(i) {
      level <- levels(object$group.by)[i]
      levelObs <- object$group.by == level
      c(getTreatmentEffectEstimate(samples[[i]], ci.style, ci.level), n = sum(levelObs))
    })))
    row.names(estimates) <- paste0(object$estimand, ".", levels(object$group.by))
  }
  
  result <-
    namedList(call       = object$call,
              method.rsp = object$method.rsp,
              method.trt = object$method.trt,
              ci.style,
              ci.level,
              numObservations,
              numSamples,
              estimates,
              n.chains   = object$n.chains,
              commonSup.rule = object$commonSup.rule)
  
  if (object$commonSup.rule != "none") {
    result$commonSup.cut <- object$commonSup.cut
    if (is.null(object$group.by)) {
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
  
  cat("Treatment effect:\n", ...)
  
  print(x$estimates, digits = digits, ...)
  if (!is.null(x$estimates[["n"]]) && any(x$estimates[["n"]] <= 10L) && x$numObservations > 10L)
    cat("  if (n < 10) group-size estimates may be unstable\n", ...)
  cat("\n")
  
  cat(round(100 * x$ci.level, digits), "% credible interval calculated by: ",
      switch(x$ci.style,
             norm = "normal approximation",
             quant = "empirical quantiles",
             hpd   = "highest posterior density"),
      "\n", sep = "", ...)
  
  if (x$commonSup.rule != "none") {
    cat("\nCommon support enforced by cutting using '", x$commonSup.rule, "' rule, cutoff: ", x$commonSup.cut,
        "\nSuppressed observations:\n", sep = "")
    print(x$n.cut)
    cat("\n")
  }
  
  cat("Result based on ", x$numSamples, " posterior samples", 
      if (!is.null(x$n.chains)) paste0(" across ", x$n.chains, " chains") else "",
      "\n", sep = "", ...)
  
  invisible(x)
}

print.bartcFit <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if (!is.null(x$call))
    cat("Call:\n", paste0(deparse(x$call), collapse = "\n  "), "\n\n", sep = "")
  
  cat("Treatment effect (", x$estimand, "): ", format(fitted(x, "est"), digits = digits), "\n", sep = "", ...)
  
  invisible(x)
}
