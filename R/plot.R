plot.sigma <- function(x, main = "Traceplot sigma", xlab = "iteration", ylab = "sigma", lty = 1:x$n.chains, ...) {
  if (!is(x, "cibartFit")) stop("plot.sigma requires an object of class 'cibartFit'")
  
  first.sigma <- x$fit.rsp$first.sigma
  
  if (is.null(dim(x$fit.rsp$sigma))) {
    sigma <- c(x$fit.rsp$sigma, sigma)
    numBurnIn  <- length(first.sigma)
    numSamples <- length(sigma)
  } else {
    sigma <- cbind(x$fit.rsp$first.sigma, sigma)
    numBurnIn  <- ncol(first.sigma)
    numSamples <- ncol(sigma)
  }
  
  plot(NULL, xlim = c(1L, numSamples), ylim = range(sigma), main = main, xlab = xlab, ylab = ylab, ...)
  abline(v = numBurnIn, col = "red", lwd = 0.5)
  
  if (NCOL(sigma) > 1L) {
    for (i in seq_len(x$n.chains))
      lines(seq_len(numSamples), sigma[i,], lty = lty[i], ...)
  } else {
    lines(seq_len(numSamples), sigma, lty = lty[1L], ...)
  }
  
  invisible(NULL)
}

plot.est <- function(x, main = paste("Traceplot", x$estimand), xlab = "iteration", ylab = x$estimand, lty = 1:x$n.chains, ...) {
  if (!is(x, "cibartFit")) stop("plot.set requires an object of class 'cibartFit'")
  
  samples <- extract(x, "est", combineChains = FALSE)
  numSamples <- if (!is.null(dim(samples))) dim(samples)[2L] else length(samples)

  plot(NULL, xlim = c(1L, numSamples), ylim = range(samples), main = main, xlab = xlab, ylab = ylab, ...)
  
  if (!is.null(dim(samples))) {
    for (i in seq_len(x$n.chains))
      lines(seq_len(numSamples), samples[i,], lty = lty[i], ...)
  } else {
    lines(seq_len(numSamples), samples, lty = lty[1L], ...)
  }
  
  invisible(NULL)
}

plot.indiv <- function(x, main = "Histogram Individual Effects", xlab = "treatment effect", breaks = 20, ...)
{
  if (!is(x, "cibartFit")) stop("plot.indiv requires an object of class 'cibartFit'")
  
  samples <- extract(x, "indiv.diff", combineChains = TRUE)
  
  hist(apply(samples, 1L, mean), main = main, xlab = xlab, breaks = breaks, ...)
}

if (FALSE) plot.support <- function(x, main = "Histogram Common Support Stat", xlab = "common support statistic", breaks = 20, ...)
{
  if (!is(x, "cibartFit")) stop("plot.support requires an object of class 'cibartFit'")
  
  commonSup.stat <- getCommonSupportStatistic(x$sd.obs, x$sd.cf, x$commonSup.rule, x$commonSup.cut)
  commonSup.cut  <- getCommonSupportCutoff(x$sd.obs, x$sd.cf, x$commonSup.rule, x$commonSup.cut)
  
  hist(commonSup.stat, main = main, xlab = xlab, breaks = breaks, ...)
  abline(v = commonSup.cut, col = "red", ...)
}

