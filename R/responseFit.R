getBartResponseFit <- function(response, treatment, confounders, data, subset, weights, estimand,
                               group.by = NULL, use.ranef = TRUE,
                               commonSup.rule, commonSup.cut, p.score, crossvalidate = FALSE, calculateEstimates = TRUE, ...)
{
  treatmentIsMissing    <- missing(treatment)
  responseIsMissing     <- missing(response)
  confoundersAreMissing <- missing(confounders)
  weightsAreMissing     <- missing(weights)
  dataAreMissing        <- missing(data)
  
  matchedCall <- match.call()
  callingEnv <- parent.frame(1L)
  
  if (treatmentIsMissing)
    stop("'treatment' variable must be specified")
  if (responseIsMissing)
    stop("'response' variable must be specified")
  if (confoundersAreMissing)
    stop("'confounders' variable must be specified")
  
  if (!is.character(estimand) || estimand[1L] %not_in% c("ate", "att", "atc"))
    stop("estimand must be one of 'ate', 'att', or 'atc'")
  estimand <- estimand[1L]
  
  dbartsDataCall <- NULL; treatmentName <- NULL; missingRows <- NULL
  if (!dataAreMissing && is.data.frame(data)) {
    evalEnv <- NULL
    dataCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getResponseDataCall)), "fn", quote(dbarts::dbartsData))
    dataCall <- addCallDefaults(dataCall, eval(quoteInNamespace(getBartResponseFit)))
    massign[dbartsDataCall, evalEnv, treatmentName, missingRows] <- eval(dataCall, envir = callingEnv)
  } else {
    df <- NULL
    literalCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getResponseLiteralCall)), "fn", quote(dbarts::dbartsData))
    literalCall <- addCallDefaults(literalCall, eval(quoteInNamespace(getBartResponseFit)))
    
    dataEnv <- if (dataAreMissing) callingEnv else list2env(data, parent = callingEnv)
    
    massign[dbartsDataCall, df, treatmentName, missingRows] <- eval(literalCall, envir = dataEnv)
    
    evalEnv <- sys.frame(sys.nframe())
  }
  
  responseData <- eval(dbartsDataCall, envir = evalEnv)
  n <- length(missingRows)
  n.obs <- nrow(responseData@x)
  
  missingData <- NULL
  if (any(missingRows)) {
    ## replace response with inverted missing ness so we can get a data object for use later
    evalEnv[[deparse(dbartsDataCall$data)]][[deparse(dbartsDataCall[[2L]][[2L]])]] <- 
      ifelse(missingRows, 0, NA)
    missingData <- eval(dbartsDataCall, envir = evalEnv)
    n.mis <- nrow(missingData@x)
  }
  
  if (is.null(missingData)) {
    responseData@x.test <- responseData@x
    responseData@x.test[,treatmentName] <- 1 - responseData@x.test[,treatmentName]
  } else {
    ## structure so that first part is a counterfactual estimate ordered as are all observations
    responseData@x.test <- matrix(0, n + n.mis, ncol(responseData@x), dimnames = dimnames(responseData@x))
    responseData@x.test[which(!missingRows),] <- responseData@x
    responseData@x.test[which( missingRows),] <- missingData@x
    responseData@x.test[seq.int(n + 1L, n + n.mis),] <- missingData@x
    
    cfRows <- rep_len(TRUE, n)
    responseData@x.test[cfRows,treatmentName] <- 1 - responseData@x.test[cfRows,treatmentName]
  }
  
  ## redirect to pull in any args passed
  use.ranef <- !is.null(matchedCall[["group.by"]]) && use.ranef
  if (!use.ranef) {
    bartCall <- if (!use.ranef) redirectCall(matchedCall, dbarts::bart2)
  } else {
    group.by <- eval(redirectCall(matchedCall, quoteInNamespace(getGroupBy)), envir = callingEnv)
    if (!is.null(missingData)) {
      group.by.test <- c(group.by[!missingRows], group.by[missingRows], group.by[missingRows])
      group.by <- group.by[!missingRows]
    } else {
      group.by.test <- group.by
    }
    matchedCall$group.by      <- group.by
    matchedCall$group.by.test <- group.by.test
    
    bartCall <- redirectCall(matchedCall, dbarts::rbart_vi)
  }
  
  invalidArgs <- names(bartCall)[-1L] %not_in% names(eval(formals(eval(bartCall[[1L]])))) &
                 names(bartCall)[-1L] %not_in% names(eval(formals(dbarts::dbartsControl)))

  bartCall <- bartCall[c(1L, 1L + which(!invalidArgs))]
    
  bartCall$formula <- quote(responseData)
  bartCall$data    <- NULL
  bartCall$verbose <- FALSE
  if (is.null(bartCall[["n.chains"]])) bartCall[["n.chains"]] <- 10L
  
  evalEnv <- new.env(parent = callingEnv)
  evalEnv[["responseData"]] <- responseData
  
  
  responseIsBinary <- unique(responseData@y)
  responseIsBinary <- length(responseIsBinary) == 2L && min(responseIsBinary) == 0 && max(responseIsBinary) == 1
  
  if (crossvalidate)
    bartCall <- optimizeBARTCall(bartCall, evalEnv)
  
  bartFit <- eval(bartCall, envir = evalEnv)
  
  if (is.null(missingData)) {
    trt <- responseData@x[,treatmentName]
  } else {
    trt <- numeric(n)
    trt[!missingRows] <- responseData@x[,treatmentName]
    trt[ missingRows] <- missingData@x[,treatmentName]
  }
  treatmentRows <- trt > 0
  
  combineChains <- if (is.null(matchedCall[["combineChains"]])) FALSE else list(...)[["combineChains"]]
  mu.hat.train <- extract(bartFit, sample = "train", combineChains = combineChains)
  mu.hat.test  <- extract(bartFit, sample = "test",  combineChains = combineChains)
  
  if (is.null(missingData)) {
    mu.hat.obs <- mu.hat.train
    mu.hat.cf  <- mu.hat.test      
  } else {
    # input dims are n.chains x n.samples x n.obs
    # perm to n.obs x n.chains x n.samples and then perm back
    if (length(dim(mu.hat.train)) > 2L) {
      mu.hat.train <- aperm(mu.hat.train, c(3L, 1L, 2L))
      mu.hat.test  <- aperm(mu.hat.test, c(3L, 1L, 2L))
    } else {
      mu.hat.train <- t(mu.hat.train)
      mu.hat.test  <- t(mu.hat.test)
    }
    mu.hat.obs <- array(0, c(n, dim(mu.hat.train)[-1L]))
    addDimsToSubset(mu.hat.obs[!missingRows] <- mu.hat.train)
    addDimsToSubset(mu.hat.obs[ missingRows] <- mu.hat.test[seq.int(n + 1L, n + n.mis), drop = FALSE])
    
    addDimsToSubset(mu.hat.cf <- mu.hat.test[seq_len(n)])
    
    if (length(dim(mu.hat.train)) > 2L) {
      mu.hat.train <- aperm(mu.hat.train, c(2L, 3L, 1L))
      mu.hat.test  <- aperm(mu.hat.test, c(2L, 3L, 1L))
      mu.hat.obs   <- aperm(mu.hat.obs, c(2L, 3L, 1L))
      mu.hat.cf    <- aperm(mu.hat.cf, c(2L, 3L, 1L))
    } else {
      mu.hat.train <- t(mu.hat.train)
      mu.hat.test  <- t(mu.hat.test)
      mu.hat.obs   <- t(mu.hat.obs)
      mu.hat.cf    <- t(mu.hat.cf)
    }
  }
  
  sd.obs <- apply(mu.hat.obs, length(dim(mu.hat.obs)), sd)
  sd.cf  <- apply(mu.hat.cf,  length(dim(mu.hat.obs)), sd)
  
  commonSup.sub <- getCommonSupportSubset(sd.obs, sd.cf, commonSup.rule, commonSup.cut, trt, missingRows)
  
  if (is.null(bartFit[["y"]])) bartFit[["y"]] <- responseData@y
  
  result <- namedList(fit = bartFit, data = responseData, mu.hat.obs, mu.hat.cf, name.trt = treatmentName, trt, sd.obs, sd.cf, commonSup.sub, missingRows, est = NULL, fitPars = NULL)
  
  if (crossvalidate)
    result[["k"]] <- bartCall[["k"]]
  
  result
}

boundValues <- function(x, bounds){
  x[x > max(bounds)] <- max(bounds)
  x[x < min(bounds)] <- min(bounds)
  x
}

# expects inputs with permuted dims: n.obs x n.chains x n.samples
getPWeightEstimates <- function(y, z, weights, estimand, mu.hat.0, mu.hat.1, p.score, yBounds, p.scoreBounds)
{
  flattenSamples.perm <- function(y) {
    x <- NULL ## R CMD check
    if (!is.null(dim(y)) && length(dim(y)) > 2L) evalx(dim(y), matrix(y, nrow = x[1L], ncol = x[2L] * x[3L])) else y
  }

  if (!is.character(estimand) || estimand[1L] %not_in% c("ate", "att", "atc"))
    stop("estimand must be one of 'ate', 'att', or 'atc'")
  estimand <- estimand[1L]
  
  if (!is.null(weights)) {
    weights <- rep_len(weights, length(y))
    weights <- weights / sum(weights)
  }
  
  m <- min(y, na.rm = TRUE)
  M <- max(y, na.rm = TRUE)
  
  r <- range(y)
  r <- r + 0.1 * c(-abs(r[1L]), abs(r[2L]))
  y.st <- boundValues(y, r)
  r.st <- range(y.st)
  y.st <- (y.st - min(r.st)) / diff(r.st)
  
  origDims <- dim(mu.hat.0)
  
  ## map mu.hat to (0, 1)
  mu.hat.0 <- flattenSamples.perm(boundValues((boundValues(mu.hat.0, c(m, M)) - m) / (M - m), yBounds))
  mu.hat.1 <- flattenSamples.perm(boundValues((boundValues(mu.hat.1, c(m, M)) - m) / (M - m), yBounds))
  
  icate <- mu.hat.1 - mu.hat.0
  
  p.score <- boundValues(p.score, p.scoreBounds)
  
  if (!is.null(dim(p.score))) {
    if (!all(dim(p.score) == origDims))
      stop("dimensions of p.score samples must match that of observations")
      
    p.score <- flattenSamples.perm(p.score)
  }
    
  getPWeightEstimate <- getPWeightFunction(estimand, weights, icate, p.score)
  tmleFuncs <- getTMLEFunctions(estimand, weights)
  mu.hat.1.deriv <- tmleFuncs$mu.hat.1.deriv
  mu.hat.0.deriv <- tmleFuncs$mu.hat.0.deriv
  if (!is.null(weights)) {
    icBody <- switch(estimand,
      att = quote((length(y) * weights * a.weight * (y - mu.hat) + z * t(t(icate) - psi)) / sum(p.score * weights)),
      atc = quote((length(y) * weights * a.weight * (y - mu.hat) + (1 - z) * t(t(icate) - psi)) / sum((1 - p.score) * weights)),
      ate = quote(length(y) * weights * a.weight * (y - mu.hat) + t(t(icate) - psi)))
  } else {
    icBody <- switch(estimand,
      att = quote((a.weight * (y - mu.hat) + z * t(t(icate) - psi)) / mean(z)),
      atc = quote((a.weight * (y - mu.hat) + (1 - z) * t(t(icate) - psi)) / mean(1 - z)),
      ate = quote(a.weight * (y - mu.hat) + t(t(icate) - psi)))
  }
  getIC <- function(y, mu.hat, icate, psi, a.weight) { }
  body(getIC) <- icBody
  
  mu.hat <- mu.hat.1 * z + mu.hat.0 * (1 - z)
  
  psi <- getPWeightEstimate(z, weights, icate, p.score)
    
  a.weight <- z * mu.hat.1.deriv(z, weights, p.score) + (1 - z) * mu.hat.0.deriv(z, weights, p.score)
  ic <- getIC(y.st, mu.hat, icate, psi, a.weight)
  
  se <- apply(ic, 2L, sd) / sqrt(length(y))
  result <- c(psi * (M - m), sd(ic) / sqrt(length(y)))
  
  if (!is.null(origDims) && length(origDims) > 2L)
    array(c(psi * (M - m), se), c(origDims[2L], origDims[3L], 2L), dimnames = list(NULL, NULL, c("est", "se")))
  else 
    matrix(c(psi * (M - m), se), length(psi), 2L, dimnames = list(NULL, c("est", "se")))
}

# TODO: rewrite this so it doesn't have to permute/transpose the samples
getPWeightResponseFit <-
  function(response, treatment, confounders, data, subset, weights, estimand,
           group.by = NULL, use.ranef = TRUE, group.effects = FALSE,
           p.score, samples.p.score,
           yBounds = c(.005, .995), p.scoreBounds = c(0.025, 0.975), ...)
{
  dataAreMissing    <- missing(data)
  weightsAreMissing <- missing(weights)
  
  matchedCall <- match.call()
  callingEnv <- parent.frame(1L)
  
  if (is.null(matchedCall$p.score) || is.null(p.score))
    stop("propensity score weighting only possible if propensity score provided")
  
  if (weightsAreMissing) {
    weights <- NULL
  } else if (!dataAreMissing) {
    weights <- eval(matchedCall$weights, envir = data)
  }
  
  bartCall <- redirectCall(matchedCall, quoteInNamespace(getBartResponseFit))
  bartCall$calculateEstimates <- FALSE
  
  fit <- data <- mu.hat.obs <- mu.hat.cf <- name.trt <- trt <- sd.obs <- sd.cf <- commonSup.sub <- missingRows <- NULL
  assignAll(eval(bartCall, envir = callingEnv))
  
  treatmentRows <- trt > 0
  
  mu.hat.obs.orig <- mu.hat.obs
  mu.hat.cf.orig  <- mu.hat.cf
  # input dims are n.chains x n.samples x n.obs
  # permuate to n.obs x n.chains x n.samples
  if (length(dim(mu.hat.obs)) > 2L) {
    mu.hat.obs <- aperm(mu.hat.obs, c(3L, 1L, 2L))
    mu.hat.cf  <- aperm(mu.hat.cf,  c(3L, 1L, 2L))
  } else {
    mu.hat.obs <- t(mu.hat.obs)
    mu.hat.cf  <- t(mu.hat.cf)
  }
  
  mu.hat.1 <- mu.hat.obs * trt       + mu.hat.cf * (1 - trt)
  mu.hat.0 <- mu.hat.obs * (1 - trt) + mu.hat.cf * trt
  
  p.score <- if (!is.null(matchedCall$samples.p.score) && !is.null(samples.p.score)) samples.p.score else p.score
  if (!is.null(dim(p.score)) && length(dim(p.score)) < length(dim(mu.hat.obs.orig))) {
    # chains were collapsed
    n.chains  <- dim(mu.hat.obs.orig)[1L]
    n.samples <- dim(mu.hat.obs.orig)[2L]
    n.obs     <- dim(mu.hat.obs.orig)[3L]
    p.score   <- aperm(array(p.score, c(n.chains, n.obs, n.samples)), c(3L, 1L, 2L))
  } else {
    if (!is.null(dim(p.score)))
      p.score <- if (length(dim(p.score)) > 2L) aperm(p.score, c(3L, 1L, 2L)) else t(p.score)
  }
  
  if (is.null(matchedCall[["group.by"]]) || !group.effects) {
    if (any(commonSup.sub != TRUE)) {
      addDimsToSubset(mu.hat.0 <- mu.hat.0[commonSup.sub, drop = FALSE])
      addDimsToSubset(mu.hat.1 <- mu.hat.1[commonSup.sub, drop = FALSE])
      addDimsToSubset(p.score <- p.score[commonSup.sub, drop = FALSE])
      
      if (!is.null(weights)) weights <- weights[commonSup.sub]
    }
      
    est <- getPWeightEstimates(fit$y[commonSup.sub], trt[commonSup.sub], weights, estimand, mu.hat.0, mu.hat.1, p.score, yBounds, p.scoreBounds)
  } else {
    # we might have been given fixed effects which would live in a data frame, but we need
    # a literal to estimate within subsets
    group.by <- eval(redirectCall(matchedCall, quoteInNamespace(getGroupBy)), envir = callingEnv)
    
    est <- lapply(levels(group.by), function(level) {
      levelRows <- group.by == level & commonSup.sub
      
      addDimsToSubset(mu.hat.0 <- mu.hat.0[levelRows, drop = FALSE])
      addDimsToSubset(mu.hat.1 <- mu.hat.1[levelRows, drop = FALSE])
      addDimsToSubset(p.score <- p.score[levelRows, drop = FALSE])
      
      if (!is.null(weights)) weights <- weights[levelRows]
      
      getPWeightEstimates(fit$y[levelRows], trt[levelRows], weights, estimand, mu.hat.0, mu.hat.1, p.score, yBounds, p.scoreBounds)
    })
    names(est) <- levels(group.by)
  }
  
  namedList(fit, data, mu.hat.obs = mu.hat.obs.orig, mu.hat.cf = mu.hat.cf.orig,
            name.trt, trt, sd.obs, sd.cf, commonSup.sub, missingRows, est,
            fitPars = namedList(yBounds, p.scoreBounds))
}

getTMLEEstimates <- function(y, z, weights, estimand, mu.hat.0, mu.hat.1, p.score, yBounds, p.scoreBounds, depsilon, maxIter, n.threads)
{
  if (!is.character(estimand) || estimand[1L] %not_in% c("ate", "att", "atc"))
    stop("estimand must be one of 'ate', 'att', or 'atc'")
  estimand <- estimand[1L]
  
  flattenSamples.perm <- function(y) {
    x <- NULL ## R CMD check
    if (!is.null(dim(y)) && length(dim(y)) > 2L) evalx(dim(y), matrix(y, nrow = x[1L], ncol = x[2L] * x[3L])) else y
  }

  
  if (anyNA(y)) {
    completeRows <- !is.na(y)
    y <- y[completeRows]
    z <- z[completeRows]
    if (!is.null(weights)) weights <- weights[completeRows]
    addDimsToSubset(mu.hat.0 <- mu.hat.0[completeRows, drop = FALSE])
    addDimsToSubset(mu.hat.1 <- mu.hat.1[completeRows, drop = FALSE])
    addDimsToSubset(p.score <- p.score[completeRows, drop = FALSE])
  }
  
  tmle <- NULL
  if (is.null(weights) && is(tryCatch(tmle <- tmle::tmle, error = function(e) e), "error"))
    warning("tmle package not found; install for up-to-date results with method.rsp = 'tmle'")
  
  if (!is.null(tmle)) {
    if (is.null(dim(mu.hat.0))) { 
      result <- tmle(Y = y, A = z, W = matrix(0, length(y), 1L), Q = cbind(Q0W = mu.hat.0, Q1W = mu.hat.1), g1W = p.score)
      result <- unlist(result$estimates[[switch(estimand, ate = "ATE", att = "ATT", atc = "ATC")]][c("psi", "var.psi")])
      names(result) <- c("est", "se")
      result["se"] <- sqrt(result["se"])
      result
    } else {
      p.score <- boundValues(flattenSamples.perm(p.score), p.scoreBounds)
      W <- matrix(0, length(y), 1L)
      Q <- aperm(array(c(flattenSamples.perm(mu.hat.0), flattenSamples.perm(mu.hat.1)),
                       c(length(y), prod(dim(mu.hat.0)[-1L]), 2L), dimnames = list(NULL, NULL, c("Q0W", "Q1W"))),
                 c(1L, 3L, 2L))
      
      if (n.threads == 1L) {
        result <- t(sapply(seq_len(dim(Q)[3L]), function(i) {
          res <- tmle(Y = y, A = z, W = W, Q = Q[,,i], g1W = if (!is.null(dim(p.score))) p.score[,i] else p.score)
          unlist(res$estimates[[switch(estimand, ate = "ATE", att = "ATT", atc = "ATC")]][c("psi", "var.psi")])
        }))
      } else {
        cluster <- makeCluster(n.threads)
        
        clusterExport(cluster, c("y", "z", "W", "estimand"), sys.frame(sys.nframe()))
        
        numSamples <- dim(Q)[3L]
        numSamplesPerThread <- numSamples %/% n.threads + if (numSamples %% n.threads != 0L) 1L else 0L
        numFullThreads <- n.threads + numSamples - numSamplesPerThread * n.threads
        
        data.list <- lapply(seq_len(n.threads), function(i) {
          start <- 1L + if (i <= numFullThreads) (i - 1L) * numSamplesPerThread else numFullThreads * numSamplesPerThread + (i - numFullThreads - 1L) * (numSamplesPerThread - 1L)
          ind <- seq.int(start, length.out = numSamplesPerThread - if (i <= numFullThreads) 0L else 1L)
          list(Q = Q[,,ind,drop = FALSE], p.score = if (!is.null(dim(p.score))) p.score[,ind,drop = FALSE] else p.score)
        })
        
        tryResult <- tryCatch(results.list <- clusterApply(cluster, data.list, function(x) {
          Q <- x$Q
          p.score <- x$p.score
          sapply(seq_len(dim(Q)[3L]), function(i) {
            res <- tmle(Y = y, A = z, W = W, Q = Q[,,i], g1W = if (!is.null(dim(p.score))) p.score[,i] else p.score)
            unlist(res$estimates[[switch(estimand, ate = "ATE", att = "ATT", atc = "ATC")]][c("psi", "var.psi")])
          })
        }), error = I)
    
        stopCluster(cluster)
        
        result <- t(matrix(unlist(results.list), 2L, numSamples))
      }
      result[,2L] <- sqrt(result[,2L])
      if (length(dim(mu.hat.0)) > 2L) {
        result <- array(result, c(dim(mu.hat.0)[-1L], 2L), dimnames = list(NULL, NULL, c("est", "se")))
      } else {
        colnames(result) <- c("est", "se")
      }
    }
    return(result)
  }
  
  if (!is.null(weights)) {
    weights <- rep_len(weights, length(y))
    weights <- weights / sum(weights)
  }
  
  r <- range(y)
  r <- r + 0.1 * c(-abs(r[1L]), abs(r[2L]))
  y.st <- boundValues(y, r)
  mu.hat.0.st <- boundValues(mu.hat.0, r)
  mu.hat.1.st <- boundValues(mu.hat.1, r)
  r.st <- range(y.st)
  
  y.st <- (y.st - min(r.st)) / diff(r.st)
  mu.hat.0.st <- qlogis(boundValues((mu.hat.0.st - min(r.st)) / diff(r.st), yBounds))
  mu.hat.1.st <- qlogis(boundValues((mu.hat.1.st - min(r.st)) / diff(r.st), yBounds))
  
  
  mu.hat.0.samp <- flattenSamples.perm(mu.hat.0.st)
  mu.hat.1.samp <- flattenSamples.perm(mu.hat.1.st)
  
  p.score.samp <- boundValues(flattenSamples.perm(p.score), p.scoreBounds)
  
  origDims <- dim(mu.hat.0)
  
  getPWeightEstimate <- getPWeightFunction(estimand, weights, numeric(), numeric())
  
  mu.hat.0.deriv <- mu.hat.1.deriv <- p.score.deriv <- getIC <- calcLoss <- NULL
  assignAll(getTMLEFunctions(estimand, weights))
  
  result <- t(sapply(seq_len(ncol(mu.hat.0.samp)), function(i) {
    mu.hat.0 <- mu.hat.0.samp[,i]
    mu.hat.1 <- mu.hat.1.samp[,i]
    mu.hat <- mu.hat.1 * z + mu.hat.0 * (1 - z)
    
    p.score <- if (!is.null(dim(p.score.samp))) p.score.samp[,i] else p.score.samp
    p.score.st <- boundValues(p.score, p.scoreBounds)
    
    H1W <- z / p.score.st
    H0W <- (1 - z) / (1 - p.score.st)
    
    suppressWarnings(epsilon <- coef(glm(y.st ~ -1 + offset(mu.hat) + H0W + H1W, family = binomial)))
    epsilon[is.na(epsilon)] <- 0 
  
    mu.hat.0 <- plogis(mu.hat.0 + epsilon["H0W"] / (1 - p.score.st))
    mu.hat.1 <- plogis(mu.hat.1 + epsilon["H1W"] / p.score.st)
    
    icate <- mu.hat.1 - mu.hat.0
    mu.hat <- mu.hat.1 * z + mu.hat.0 * (1 - z)
    
    psi <- getPWeightEstimate(z, weights, icate, p.score)
    psi.prev <- psi
    
    a.weight <- z * mu.hat.1.deriv(z, weights, p.score) + (1 - z) * mu.hat.0.deriv(z, weights, p.score)
    ic.prev <- ic <- getIC(y.st, mu.hat, icate, psi, a.weight)
    
    if (mean(ic) > 0) depsilon <- -depsilon
    
    loss.prev <- Inf
    loss <- calcLoss(y.st, z, mu.hat, p.score, weights)
    if (is.nan(loss) || is.na(loss) || is.infinite(loss)) return(psi)
      
    iter <- 0L
    while (loss.prev > loss && iter < maxIter)
    {
      p.score.prev <- p.score
      p.score <- boundValues(plogis(qlogis(p.score.prev) - depsilon * p.score.deriv(z, weights, p.score.prev, icate, psi.prev)), p.scoreBounds)
      
      mu.hat.0.prev <- mu.hat.0
      mu.hat.1.prev <- mu.hat.1
      mu.hat.0 <- boundValues(plogis(qlogis(mu.hat.0.prev) - depsilon * mu.hat.0.deriv(z, weights, p.score.prev)), yBounds)
      mu.hat.1 <- boundValues(plogis(qlogis(mu.hat.1.prev) - depsilon * mu.hat.1.deriv(z, weights, p.score.prev)), yBounds)
      icate <- mu.hat.1 - mu.hat.0
      mu.hat <- mu.hat.1 * z + mu.hat.0 * (1 - z)
      
      psi.prev <- psi
      psi <- getPWeightEstimate(z, weights, icate, p.score)
      
      loss.prev <- loss
      loss <- calcLoss(y.st, z, mu.hat, p.score, weights)
      
      ic.prev <- ic
      ic <- getIC(y.st, mu.hat, icate, psi, a.weight)
      
      if (is.nan(loss) || is.infinite(loss) || is.na(loss)) loss <- Inf
    
      iter <- iter + 1L
    }
    
    if (is.infinite(loss) || loss.prev < loss)
      c(psi.prev, sd(ic.prev) / sqrt(length(y.st)))
    else
      c(psi, sd(ic) / sqrt(length(y.st)))
  }))
  
  result[,1L] <- result[,1L] * (max(r.st) - min(r.st))
  colnames(result) <- c("est", "se")

  
  if (!is.null(origDims) && length(origDims) > 2L)
    result <- array(result, c(origDims[2L], origDims[3L], 2L), dimnames = list(NULL, NULL, c("est", "se")))
  
  result
}

getTMLEResponseFit <-
  function(response, treatment, confounders, data, subset, weights, estimand,
           group.by = NULL, use.ranef = TRUE, group.effects = FALSE,
           p.score, samples.p.score,
           verbose, posteriorOfTMLE = TRUE,
           yBounds = c(.005, .995), p.scoreBounds = c(0.025, 0.975), depsilon = 0.001, maxIter = max(1000, 2 / depsilon), ...)
{
  dataAreMissing    <- missing(data)
  weightsAreMissing <- missing(weights)
  
  matchedCall <- match.call()
  callingEnv <- parent.frame(1L)
  
  if (is.null(matchedCall$p.score) || is.null(p.score))
    stop("TMLE only possible if propensity score provided")
  
  if (weightsAreMissing) {
    weights <- NULL
  } else if (!dataAreMissing) {
    weights <- eval(matchedCall$weights, envir = data)
  }
  
  bartCall <- redirectCall(matchedCall, quoteInNamespace(getBartResponseFit))
  bartCall$calculateEstimates <- FALSE
  
  fit <- data <- mu.hat.obs <- mu.hat.cf <- name.trt <- trt <- sd.obs <- sd.cf <- commonSup.sub <- missingRows <- NULL
  assignAll(eval(bartCall, envir = callingEnv))
  
  mu.hat.obs.orig <- mu.hat.obs
  mu.hat.cf.orig  <- mu.hat.cf
  # input dims are n.chains x n.samples x n.obs
  if (length(dim(mu.hat.obs)) > 2L) {
    mu.hat.obs <- aperm(mu.hat.obs, c(3L, 1L, 2L))
    mu.hat.cf  <- aperm(mu.hat.cf,  c(3L, 1L, 2L))
  } else {
    mu.hat.obs <- t(mu.hat.obs)
    mu.hat.cf  <- t(mu.hat.cf)
  }
  
  treatmentRows <- trt > 0
  
  mu.hat.1 <- mu.hat.obs * trt       + mu.hat.cf * (1 - trt)
  mu.hat.0 <- mu.hat.obs * (1 - trt) + mu.hat.cf * trt
   
  p.score <- if (!is.null(matchedCall$samples.p.score) && !is.null(samples.p.score)) samples.p.score else p.score
  if (!is.null(dim(p.score)) && length(dim(p.score)) < length(dim(mu.hat.obs.orig))) {
    n.chains  <- dim(mu.hat.obs.orig)[1L]
    n.samples <- dim(mu.hat.obs.orig)[2L]
    n.obs     <- dim(mu.hat.obs.orig)[3L]
    p.score <- aperm(array(p.score, c(n.samples, n.chains, n.obs)), c(2L, 1L, 3L))
  }
  if (!is.null(dim(p.score)))
    p.score <- if (length(dim(p.score)) > 2L) aperm(p.score, c(3L, 1L, 2L)) else t(p.score)
  
  maxIter <- round(maxIter, 0)
  
  if (verbose)
    cat("calculating TMLE adjustment\n")
  
  if (is.null(matchedCall$group.by) || !group.effects) {
    if (any(commonSup.sub != TRUE)) {
      addDimsToSubset(mu.hat.0 <- mu.hat.0[commonSup.sub, drop = FALSE])
      addDimsToSubset(mu.hat.1 <- mu.hat.1[commonSup.sub, drop = FALSE])
      
      addDimsToSubset(p.score <- p.score[commonSup.sub, drop = FALSE])
      
      if (!is.null(weights)) weights <- weights[commonSup.sub]
    }
    
    if (posteriorOfTMLE) {
      n.threads <- if ("n.threads" %in% names(list(...))) list(...)[["n.threads"]] else dbarts::guessNumCores()
      est <- getTMLEEstimates(fit$y[commonSup.sub], trt[commonSup.sub], weights, estimand,
                              mu.hat.0, mu.hat.1, p.score, yBounds, p.scoreBounds, depsilon, maxIter,
                              n.threads = n.threads)
    } else {
      est <- getTMLEEstimates(fit$y[commonSup.sub], trt[commonSup.sub], weights, estimand,
                              apply(mu.hat.0, 1L, mean),
                              apply(mu.hat.1, 1L, mean),
                              if (!is.null(dim(p.score))) apply(p.score, 1L, mean) else p.score,
                              yBounds, p.scoreBounds, depsilon, maxIter, n.threads = 1L)
    }
  } else {
    group.by <- eval(redirectCall(matchedCall, quoteInNamespace(getGroupBy)), envir = callingEnv)
    
    est <- lapply(levels(group.by), function(level) {
      levelRows <- group.by == level & commonSup.sub
      
      addDimsToSubset(mu.hat.0 <- mu.hat.0[levelRows, drop = FALSE])
      addDimsToSubset(mu.hat.1 <- mu.hat.1[levelRows, drop = FALSE])
      
      addDimsToSubset(p.score <- p.score[levelRows, drop = FALSE])
      
      if (!is.null(weights)) weights <- weights[levelRows]
      
      if (posteriorOfTMLE) {
        n.threads <- if ("n.threads" %in% names(list(...))) list(...)[["n.threads"]] else dbarts::guessNumCores()
        getTMLEEstimates(fit$y[levelRows], trt[levelRows], weights, estimand,
                         mu.hat.0, mu.hat.1, p.score, yBounds, p.scoreBounds, depsilon, maxIter,
                         n.threads = n.threads)
      } else {
        est <- getTMLEEstimates(fit$y[levelRows], trt[levelRows], weights, estimand,
                                apply(mu.hat.0, 1L, mean),
                                apply(mu.hat.1, 1L, mean),
                                if (!is.null(dim(p.score))) apply(p.score, 1L, mean) else p.score,
                                yBounds, p.scoreBounds, depsilon, maxIter, n.threads = 1L)
      }
    })
    names(est) <- levels(group.by)
  }

  namedList(fit, data , mu.hat.obs = mu.hat.obs.orig, mu.hat.cf = mu.hat.cf.orig,
            name.trt, trt, sd.obs, sd.cf, commonSup.sub,
            missingRows, est, fitPars = namedList(yBounds, p.scoreBounds, depsilon, maxIter))
}
