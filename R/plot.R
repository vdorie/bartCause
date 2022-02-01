plot_sigma <- function(x, main = "Traceplot sigma", xlab = "iteration", ylab = "sigma", lty = 1:x$n.chains, ...) {
  if (!inherits(x, "bartcFit"))
    stop("plot.sigma requires an object of class 'bartcFit'")
  
  if (responseIsBinary(x))
    stop("residual standard deviation plot requires a continuous response")
  
  if (inherits(x$fit.rsp, "stan4bartFit")) {
    warmup.sigma <- t(extract(x$fit.rsp, "sigma", combine_chains = FALSE, include_warmup = "only"))
    sample.sigma <- t(extract(x$fit.rsp, "sigma", combine_chains = FALSE, include_warmup = FALSE))
  } else {
    warmup.sigma <- x$fit.rsp$first.sigma
    sample.sigma <- x$fit.rsp$sigma
  }
  
  if (is.null(dim(sample.sigma))) {
    sigma <- c(warmup.sigma, sample.sigma)
    numBurnIn  <- length(warmup.sigma)
    numSamples <- length(sample.sigma)
  } else {
    sigma <- cbind(warmup.sigma, sample.sigma)
    numBurnIn  <- ncol(warmup.sigma)
    numSamples <- ncol(sample.sigma)
  }
  
  plot(NULL, xlim = c(1L, numSamples), ylim = range(sigma), main = main, xlab = xlab, ylab = ylab, ...)
  abline(v = numBurnIn, col = "red", lwd = 0.5)

  numTotalSamples <- numBurnIn + numSamples
  
  if (NCOL(sigma) > 1L) {
    for (i in seq_len(x$n.chains))
      lines(seq_len(numTotalSamples), sigma[i,], lty = lty[i], ...)
  } else {
    lines(seq_len(numTotalSamples), sigma, lty = lty[1L], ...)
  }
  
  invisible(NULL)
}

plot_est <- function(x, main = paste("Traceplot", x$estimand),
                     xlab = "iteration", ylab = x$estimand,
                     lty = 1:x$n.chains, col = NULL, ...) {
  if (!inherits(x, "bartcFit")) stop("plot.set requires an object of class 'bartcFit'")
  
  estType <- if (x$method.rsp %in% c("tmle", "p.weight")) "pate" else "cate"
  samples <-  extract(x, type = estType, combineChains = FALSE)
  if (!is.list(samples)) samples <- list(samples)
  
  if (is.null(col)) { col <- if (is.null(x$group.by)) 1 else seq_len(nlevels(x$group.by)) }
  col <- rep_len(col, length(samples))
  
  numSamples <- if (!is.null(dim(samples[[1L]]))) dim(samples[[1L]])[2L] else length(samples[[1L]])
  
  plot(NULL, xlim = c(1L, numSamples), ylim = range(sapply(samples, range)), main = main, xlab = xlab, ylab = ylab, ...)
  
  for (j in seq_along(samples)) {
    for (i in seq_len(x$n.chains)) {
      if (!is.null(dim(samples[[j]]))) {
        lines(seq_len(numSamples), samples[[j]][i,], lty = lty[i], col = col[j], ...)
      } else {
        lines(seq_len(numSamples), samples[[j]], lty = lty[1L], col = col[j], ...)
      }
    }
  }
  
  invisible(NULL)
}

plot_indiv <- function(x, main = "Histogram Individual Quantities",
                       type = c("icate", "mu.obs", "mu.cf", "mu.0", "mu.1", "y.cf", "y.0", "y.1", "ite"),
                       xlab = "treatment effect", breaks = 20, ...)
{
  if (!inherits(x, "bartcFit")) stop("plot.indiv requires an object of class 'bartcFit'")
  
  if (!is.character(type) || type[1L] %not_in% eval(formals(plot_indiv)$type))
    stop("type must be in '", paste0(eval(formals(plot_indiv)$type), collapse = "', '"), "'")
  type <- type[1L]
  
  icates <- fitted(x, type = type)
  
  hist(icates, main = main, xlab = xlab, breaks = breaks, ...)
}

plot_support <- function(x, main = "Common Support Scatterplot",
                         xvar = "tree.1", yvar = "tree.2",
                         sample = c("inferential", "all"),
                         xlab = NULL, ylab = NULL,
                         pch.trt = 21, bg.trt = "black",
                         pch.ctl = pch.trt, bg.ctl = NA,
                         pch.sup = pch.trt, bg.sup = NA, col.sup = "red", cex.sup = 1.5,
                         legend.x = "topleft", legend.y = NULL,
                         ...)
{
  if (!inherits(x, "bartcFit")) stop("plot.support requires an object of class 'bartcFit'")
  if (x$commonSup.rule == "none") stop("common support plot requires support rule other than 'none'")
  
  if (!is.character(sample) || sample[1L] %not_in% eval(formals(plot_support)$sample))
    stop("sample must be in '", paste0(eval(formals(plot_support)$sample), collapse = "', '"), "'")
  sample <- sample[1L]
  
  labelMap <- list(p.score = "propensity score",
                   mu.obs  = expression(hat(mu)(z)),
                   mu.cf   = expression(hat(mu)(1-z)),
                   mu.1    = expression(hat(mu)(1)),
                   mu.0    = expression(hat(mu)(0)),
                   y.cf    = expression(hat(y)(1-z)),
                   y.1     = expression(yz + hat(y)(1)*(1 - z)),
                   y.0     = expression(y(1-z) + hat(y)(0) * z),
                   ite     = expression(hat(y)(1) - hat(y)(0)),
                   icate   = expression(hat(mu)(1) - hat(mu)(0)),
                   p.weights = "p.weights")
  
  matchedCall <- match.call()
    
  x.matrix <- x$data.rsp@x
  
  subset <- rep_len(TRUE, nrow(x.matrix))
  if (sample == "inferential") {
    if (x$estimand == "att") subset <- x$trt > 0
    else if (x$estimand == "atc") subset <- x$trt <= 0
  }
  
  x.matrix <- x.matrix[subset,,drop = FALSE]
  
  getColumn <- function(var, lab) {
    callingEnv <- parent.frame(1L)
    
    val <- NULL
    if (is.character(var) && startsWith(var, "pca.")) {
      if (is.null(callingEnv$x.rot)) {
        decomp <- eigen(crossprod(x.matrix))
        callingEnv$x.rot <- x.matrix %*% decomp$vectors
      }
      
      col <- suppressWarnings(as.integer(substr(var, 5L, nchar(var))))
      if (is.na(col)) stop("unable to parse pca column '", var, "'")
      if (col > ncol(callingEnv$x.rot)) stop("column '", var, "' greater than number of columns")
      
      val <- callingEnv$x.rot[,col]
      if (is.null(lab)) lab <- paste0("principal component ", col)
    } else if (is.character(var) && startsWith(var, "tree.")) {
      if (is.null(callingEnv$treeVars)) {
        if (!requireNamespace("rpart", quietly = TRUE))
          stop("tree plots require the package 'rpart'; please install it to use this feature", call. = FALSE)
        df <- data.frame(.y.hat = fitted(x, type = "icate", sample = sample), x.matrix)
        tree <- rpart::rpart(.y.hat ~ ., df, ...)
        contCols <- apply(x.matrix, 2L, function(col) length(unique(col)) > 2L)
        callingEnv$treeVars <- tree$variable.importance[names(tree$variable.importance) %in% names(which(contCols))]
      }
      
      col <- suppressWarnings(as.integer(substr(var, 6L, nchar(var))))
      if (is.na(col)) stop("unable to parse tree column '", var, "'")
      if (col > length(callingEnv$treeVars)) stop("column '", var, "' greater than the number of continuous columns")
      
      val <- x.matrix[,names(callingEnv$treeVars)[col]]
      if (is.null(lab)) lab <- paste0(names(callingEnv$treeVars)[col], " (by tree)")
    } else if (var == "css") {
      val <- with(x, getCommonSupportStatistic(sd.obs, sd.cf, commonSup.rule, commonSup.cut))[subset]
      if (is.null(lab)) lab <- paste0(x$commonSup.rule, " common support stat")
    } else if (var %not_in% c("pate", "sate", "cate") && var %in% eval(formals(fitted.bartcFit)$type)) {
      val <- fitted(x, type = var, sample = sample)
      if (is.null(lab)) lab <- labelMap[[var]]
    } else if (var %in% c("y", "y.obs")) {
      val <- x$data.rsp@y[subset]
      if (is.null(lab)) lab <- expression(y(z))
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
  
  trt <- x$trt[subset]
  
  x.val <- y.val <- NULL
  massign[x.val, xlab] <- getColumn(xvar, xlab)
  massign[y.val, ylab] <- getColumn(yvar, ylab)
  plot(x.val, y.val, main = main,
       pch = ifelse(trt > 0, pch.trt, pch.ctl), bg = ifelse(trt > 0, bg.trt, bg.ctl), xlab = xlab, ylab = ylab, ...)
  if (any(x$commonSup.sub == FALSE)) {
    commonSup.sub <- x$commonSup.sub[subset]
    points(x.val[!commonSup.sub], y.val[!commonSup.sub], pch = pch.sup, col = col.sup, bg = bg.sup, cex = cex.sup, ...)
  }
  
  if (xvar == "css") {
    cut <- with(x, getCommonSupportCutoff(sd.obs, sd.cf, commonSup.rule, commonSup.cut, trt, missingRows))
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
    cut <- with(x, getCommonSupportCutoff(sd.obs, sd.cf, commonSup.rule, commonSup.cut, trt, missingRows))
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

