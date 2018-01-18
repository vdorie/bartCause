plot_sigma <- function(x, main = "Traceplot sigma", xlab = "iteration", ylab = "sigma", lty = 1:x$n.chains, ...) {
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

plot_est <- function(x, main = paste("Traceplot", x$estimand), xlab = "iteration", ylab = x$estimand, lty = 1:x$n.chains, ...) {
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

plot_indiv <- function(x, main = "Histogram Individual Effects", xlab = "treatment effect", breaks = 20, ...)
{
  if (!is(x, "cibartFit")) stop("plot.indiv requires an object of class 'cibartFit'")
  
  samples <- extract(x, "indiv.diff", combineChains = TRUE)
  
  hist(apply(samples, 1L, mean), main = main, xlab = xlab, breaks = breaks, ...)
}

plot_support <- function(x, main = "Common Support Scatterplot", xvar = "pca.1", yvar = "pca.2",
                         pch.trt = 21, bg.trt = "black",
                         pch.ctl = pch.trt, bg.ctl = NA,
                         pch.sup = pch.trt, bg.sup = NA, col.sup = "red", cex.sup = 1.5,
                         legend.x = "topleft", legend.y = NULL,
                         ...)
{
  if (!is(x, "cibartFit")) stop("plot.support requires an object of class 'cibartFit'")
  if (x$commonSup.rule == "none") stop("common support plot requires support rule other than 'none'")
  
  matchedCall <- match.call()
  
  x.matrix <- x$data.rsp@x
  
  getColumn <- function(var, lab) {
    callingEnv <- parent.frame(1L)
    
    val <- NULL
    if (is.character(var) && startsWith(var, "pca.")) {
      if (is.null(callingEnv$x.rot)) {
        decomp <- eigen(crossprod(x.matrix))
        callingEnv$x.rot <- x.matrix %*% decomp$vectors
      }
      
      if (is.character(var) && startsWith(var, "pca.")) {
        col <- suppressWarnings(as.integer(substr(var, 5L, nchar(var))))
        if (is.na(col)) stop("unable to parse pca column '", var, "'")
        if (col > ncol(callingEnv$x.rot)) stop("column '", var, "' greater than number of columns")
        
        val <- callingEnv$x.rot[,col]
        if (is.null(lab)) lab <- paste0("principal component ", col)
      }
    } else if (var == "css") {
      val <- with(x, getCommonSupportStatistic(sd.obs, sd.cf, commonSup.rule, commonSup.cut))
      if (is.null(lab)) lab <- paste0(x$commonSup.rule, " common support stat")
    } else if (var == "p.score") {
      if (is.null(x$p.score)) stop("p.score not present in fit")
      val <- x$p.score
      if (is.null(lab)) lab <- "propensity score"
    } else if (is.numeric(var)) {
      if (length(var) == 1L) {
        var <- as.integer(round(var, 0))
        val <- x.matrix[,var]
        if (is.null(lab)) lab <- colnames(x.matrix)[var]
      } else {
        val <- var
        if (is.null(lab)) lab <- "given vector"
      }
    } else if (is.character(var)) {
      if (var %not_in% colnames(x.matrix))
        stop("unrecognized variable '", var, "'")
      val <- x.matrix[,var]
      if (is.null(lab)) lab <- var
    } else {
      stop("unrecognized variable '", var, "'")
    }
    
    list(val, lab)
  }
  
  if (is.null(matchedCall$xlab)) xlab <- NULL
  if (is.null(matchedCall$ylab)) ylab <- NULL
  
  x.val <- y.val <- NULL
  
  massign[x.val, xlab] <- getColumn(xvar, xlab)
  massign[y.val, ylab] <- getColumn(yvar, ylab)
  
  plot(x.val, y.val, main = main,
       pch = ifelse(x$trt > 0, pch.trt, pch.ctl), bg = ifelse(x$trt > 0, bg.trt, bg.ctl), xlab = xlab, ylab = ylab, ...)
  if (any(x$commonSup.sub == FALSE))
    points(x.val[!x$commonSup.sub], y.val[!x$commonSup.sub], pch = pch.sup, col = col.sup, bg = bg.sup, cex = cex.sup, ...)
  
  if (xvar == "css") {
    cut <- with(x, getCommonSupportCutoff(sd.obs, sd.cf, commonSup.rule, commonSup.cut, trt))
    if (is.null(matchedCall$lwd)) lwd <- 0.75
    abline(v = cut, col = col.sup, lwd = lwd, ...)
    if (length(cut) > 1L) {
      mtext(names(cut)[1L], 1, 1.2 * par("mgp")[2L], at = cut[1L], cex = par("cex"), col = col.sup, las = 2, ...)
      mtext(names(cut)[2L], 1, 1.2 * par("mgp")[2L], at = cut[2L], cex = par("cex"), col = col.sup, las = 2, ...)
    } else {
      mtext("cut", 1, 1.2 * par("mgp")[2L], at = cut[1L], cex = par("cex"), col = col.sup, las = 2, ...)
    }
  }
  if (yvar == "css") {
    cut <- with(x, getCommonSupportCutoff(sd.obs, sd.cf, commonSup.rule, commonSup.cut, trt))
    if (is.null(matchedCall$lwd)) lwd <- 0.75
    abline(h = cut, col = col.sup, lwd = lwd, ...)
    if (length(cut) > 1L) {
      mtext(names(cut)[1L], 2, 1.2 * par("mgp")[2L], at = cut[1L], cex = par("cex"), col = col.sup, las = 2, ...)
      mtext(names(cut)[2L], 2, 1.2 * par("mgp")[2L], at = cut[2L], cex = par("cex"), col = col.sup, las = 2, ...)
    } else {
      mtext("cut", 2, 1.2 * par("mgp")[2L], at = cut[1L], cex = par("cex"), col = col.sup, las = 2, ...)
    }
  }
  
  if (!is.null(legend.x) && !is.na(legend.x)) {
    if (is.null(matchedCall$col)) col <- "black"
    legend(legend.x, legend.y, legend = c("trt", "ctl", "drop"),
           pch = c(pch.trt, pch.ctl, pch.sup), col = c(col, col, col.sup),
           pt.bg = c(bg.trt, bg.ctl, bg.sup), ...)
  }
  
  invisible(NULL)
}

