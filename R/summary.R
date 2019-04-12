
getPATEEstimate <- function(object, ci.style, ci.level, pate.style) {
  x <- NULL ## R CMD check
  probs <- evalx((1 - ci.level) / 2, c(x, 1 - x))
  
  weights <- object$data.rsp@weights
  weights <- if (!is.null(weights)) weights / sum(weights) else evalx(length(object$data.rsp@y), rep(1 / x, x))
  
  if (object$method.rsp == "tmle") {
    samples.tmle <- extract(object, "pate")
    
    
    if (!is.null(dim(samples.tmle))) {
      samples.pate <- if (length(dim(samples.tmle)) == 2L) samples.tmle[,"est"] else as.vector(samples.tmle[,,"est"])
      est <- mean(samples.pate)
      var.between.samples <- var(samples.pate)
      
      if (pate.style == "var.exp") {
        samples.icate <- extract(object, "icate")
        n <- nrow(samples.icate)
        L <- ncol(samples.icate)
        
        if (FALSE) { # more explicit calculations
          C_li <- samples.icate
          Cbar_l <- samples.pate
          term2 <- var(Cbar_l)
          varsL <- apply((t(C_li) - Cbar_l)^2, 1L, mean) * (n / ((n - 1) * L))
          term1 <- mean(varsL)
          sd <- sqrt(term1 + term2)
        }
        
        variance.samples <- apply(weights * t(t(samples.icate) - samples.pate)^2, 2L, sum) * (n / ((n - 1) * L))
      } else {
        variance.samples <- if (length(dim(samples.tmle)) == 2L) samples.tmle[,"se"]^2 else samples.tmle[,,"se"]^2
      }
      var.within.samples <- mean(variance.samples)
      
      sd <- sqrt(var.between.samples + var.within.samples)
    } else {
      if (ci.style %in% c("quant", "hpd"))
        stop("ci style '", ci.style, "' requires posterior distribution of TMLE transformation")
      est <- samples[["est"]]
      sd  <- samples[["se"]]
    }
    ci <- switch(ci.style,
                 norm  = est + sd * qnorm(probs),
                 quant = unname(quantile(samples, probs)),
                 hpd   = getHPDInterval(samples, ci.level))
    
    return(c(estimate = est, sd = sd, ci.lower = ci[1L], ci.upper = ci[2L]))
  }
  
  if (ci.style == "norm" && pate.style == "var.exp") {
    samples <- extract(object, "icate")
    n <- nrow(samples)
    L <- ncol(samples)
    
    est <- mean(weights %*% samples)
    
    cate.samples <- apply(weights * samples, 2L, sum)
    var.between.samples <- var(cate.samples)
    variance.samples <- apply(weights * t(t(samples) - cate.samples)^2, 2L, sum) * (n / ((n - 1) * L))
    var.within.samples  <- mean(variance.samples)
    
    if (FALSE) {
      C_li <- extract(object, "icate") # has dims n x L
      Cbar_l <- apply(C_li, 2L, mean)
      term2 <- var(Cbar_l)
      varsL <- apply(C_li, 2L, var) / ncol(C_li)
      term1 <- mean(varsL)
      sd <- sqrt(term1 + term2)
      # sd <- sqrt(var(apply(samples, 2L, mean)) + mean(apply(samples, 2L, var) / L))
    }
    
    sd <- sqrt(var.between.samples + var.within.samples)
    ci <- est + sd * qnorm(probs)
    return(c(estimate = est, sd = sd, ci.lower = ci[1L], ci.upper = ci[2L]))
  }
  
  samples <- extract(object, "cate") # already weighted
  est <- mean(samples)
  
  sd <- sqrt(var(samples) +
    if (!is.null(object$fit.rsp$sigma))
      2 * mean(object$fit.rsp$sigma^2) / length(object$data.rsp@y)
    else {
      mu.o <- flattenSamples(object$mu.hat.obs)
      mu.c <- flattenSamples(object$mu.hat.cf)
      mean(apply(weights * (mu.o * (1 - mu.o) + mu.c * (1 - mu.c)), 2L, sum) / length(object$data.rsp@y))
    })
  
  # there might be a better way of getting quantile or hpd intervals without
  # drawing from the PPD
  
  ci <- switch(ci.style,
               norm  = est + sd * qnorm(probs),
               quant = { samples <- extract(object, "pate"); unname(quantile(samples, probs)) },
               hpd   = { samples <- extract(object, "pate"); getHPDInterval(samples, ci.level) })
  
  c(estimate = est, sd = sd, ci.lower = ci[1L], ci.upper = ci[2L])
}

getSATEEstimate <- function(object, ci.style, ci.level) {
  x <- NULL
  probs <- evalx((1 - ci.level) / 2, c(x, 1 - x))
  
  weights <- object$data.rsp@weights
  weights <- if (!is.null(weights)) weights / sum(weights) else evalx(length(object$data.rsp@y), rep(1 / x, x))
  
  samples <- extract(object, "cate")
  est <- mean(samples)
  sd  <- sqrt(var(samples) +
    if (!is.null(object$fit.rsp$sigma))
      mean(object$fit.rsp$sigma^2) / length(object$data.rsp@y)
    else {
      mu.c <- flattenSamples(object$mu.hat.cf)
      mean(apply(weights * mu.c * (1 - mu.c), 2L, sum) / length(object$data.rsp@y))
    })
  
  ci <- switch(ci.style,
               norm  = est + sd * qnorm(probs),
               quant = { samples <- extract(object, "sate"); unname(quantile(samples, probs)) },
               hpd   = { samples <- extract(object, "sate"); getHPDInterval(samples, ci.level) })
  
  c(estimate = est, sd = sd, ci.lower = ci[1L], ci.upper = ci[2L])
}

getCATEEstimate <- function(object, ci.style, ci.level) {
  x <- NULL
  probs <- evalx((1 - ci.level) / 2, c(x, 1 - x))
  
  samples <- extract(object, "cate")
  est <- mean(samples)
  sd  <- sd(samples)
  
  ci <- switch(ci.style,
               norm  = est + sd * qnorm(probs),
               quant = unname(quantile(samples, probs)),
               hpd   = getHPDInterval(samples, ci.level))
  c(estimate = est, sd = sd, ci.lower = ci[1L], ci.upper = ci[2L])
}

summary.bartcFit <- function(object,
                             target = c("pate", "sate", "cate"),
                             ci.style = c("norm", "quant", "hpd"), ci.level = 0.95,
                             pate.style = c("var.exp", "ppd"),
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
  
  
  numSamples      <- prod(dim(object$mu.hat.obs)[-1L])
  numObservations <- dim(object$mu.hat.obs)[1L]
  
  getTreatmentEffectEstimate <- switch(target,
    pate = function(object, ci.style, ci.level) getPATEEstimate(object, ci.style, ci.level, pate.style),
    sate = getSATEEstimate,
    cate = getCATEEstimate)
  
  
  if (is.null(object$group.by)) {
    estimates <- as.data.frame(t(getTreatmentEffectEstimate(object, ci.style, ci.level)))
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
              ci.info = namedList(target, ci.style, ci.level, pate.style),
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
  
  cat("Treatment effect (", x$ci.info$target, "):\n", ..., sep = "")
  
  print(x$estimates, digits = digits, ...)
  if (!is.null(x$estimates[["n"]]) && any(x$estimates[["n"]] <= 10L) && x$numObservations > 10L)
    cat("  if (n < 10) group-size estimates may be unstable\n", ...)
  cat("Estimates fit from ", x$numObservations - if (x$commonSup.rule != "none") sum(x$n.cut) else 0, " total observations\n", sep = "")
  
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
  
  cat("Result based on ", x$numSamples, " posterior samples", 
      if (!is.null(x$n.chains)) paste0(" across ", x$n.chains, " chains") else "",
      "\n", sep = "", ...)
  
  invisible(x)
}

print.bartcFit <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if (!is.null(x$call))
    cat("Call:\n", paste0(deparse(x$call), collapse = "\n  "), "\n\n", sep = "")
  
  cat("Treatment effect (", x$estimand, ", cate): ", format(fitted(x, "cate"), digits = digits), "\n", sep = "", ...)
  
  invisible(x)
}
