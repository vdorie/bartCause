getBartEstimates <- function(treatmentRows, weights, estimand, indiv.diff, commonSup.sub)
{
  x <- NULL ## R CMD check
  
  if (!is.character(estimand) || estimand[1L] %not_in% c("ate", "att", "atc"))
    stop("estimand must be one of 'ate', 'att', or 'atc'")
  estimand <- estimand[1L]
  
  origDims <- dim(indiv.diff)
    
  if (!is.null(weights)) {
    weights <- rep_len(weights, origDims[1L])
    weights <- weights / sum(weights)
  }
  
  result <- 
    if (length(origDims) > 2L) {
      apply(evalx(if (is.null(weights)) indiv.diff else indiv.diff * weights,
                  switch(estimand,
                         att = x[ treatmentRows & commonSup.sub,,],
                         atc = x[!treatmentRows & commonSup.sub,,],
                         ate = x[commonSup.sub,,])),
            c(2L, 3L), mean)
    } else {
      apply(evalx(if (is.null(weights)) indiv.diff else indiv.diff * weights,
                  switch(estimand,
                         att = x[ treatmentRows & commonSup.sub,],
                         atc = x[!treatmentRows & commonSup.sub,],
                         ate = x[commonSup.sub,])),
            2L, mean)
    }
  
  if (!is.null(origDims) && length(origDims) > 2L)
    matrix(result, origDims[2L], origDims[3L])
  else 
    result
}

getBartResponseFit <- function(response, treatment, confounders, data, subset, weights, estimand, group.by,
                               commonSup.rule, commonSup.cut, p.score, use.rbart, calculateEstimates = TRUE, ...)
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
    massign[dbartsDataCall, evalEnv, treatmentName, missingRows] <- eval(dataCall, envir = callingEnv)
  } else {
    df <- NULL
    literalCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getResponseLiteralCall)), "fn", quote(dbarts::dbartsData))
    
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
  
  ## redict to pull in any args passed
  use.rbart <- !is.null(matchedCall$group.by) && use.rbart
  bartCall <- if (!use.rbart) redirectCall(matchedCall, dbarts::bart2) else redirectCall(matchedCall, dbarts::rbart_vi)
  invalidArgs <- names(bartCall)[-1L] %not_in% names(eval(formals(eval(bartCall[[1L]])))) &
                 names(bartCall)[-1L] %not_in% names(eval(formals(dbarts::dbartsControl)))
  if (any(invalidArgs))
    bartCall <- bartCall[c(1L, 1L + which(!invalidArgs))]
    
  bartCall$formula <- quote(responseData)
  bartCall$data    <- NULL
  if (is.null(bartCall$keepTrees)) bartCall$keepTrees <- FALSE
  
  ## it is possible that some extra args were given for an xbart treatment call
  dotsList <- list(...)
  if (!is.null(matchedCall$n.burn) && is.list(dotsList[["n.burn"]]) && length(dotsList[["n.burn"]]) > 1L) {
    bartCall$n.burn  <- dotsList[["n.burn"]][[2L]]
  }
  if (!is.null(matchedCall$n.samples) && is.list(dotsList[["n.samples"]]) && length(dotsList[["n.samples"]]) > 1L) {
    bartCall$n.samples  <- dotsList[["n.samples"]][[2L]]
  }
  
  bartFit <- eval(bartCall)
  
  weights <- if (weightsAreMissing) NULL else eval(matchedCall$weights, envir = data)
  if (!is.null(weights)) weights <- weights / sum(weights)
  
  if (is.null(missingData)) {
    trt <- responseData@x[,treatmentName]
  } else {
    trt <- numeric(n)
    trt[!missingRows] <- responseData@x[,treatmentName]
    trt[ missingRows] <- missingData@x[,treatmentName]
  }
  treatmentRows <- trt > 0
    
  if (length(dim(bartFit$yhat.train)) > 2L) {
    yhat.train <- aperm(bartFit$yhat.train, c(3L, 1L, 2L))
    yhat.test  <- aperm(bartFit$yhat.test, c(3L, 1L, 2L))
  } else {
    yhat.train <- t(bartFit$yhat.train)
    yhat.test  <- t(bartFit$yhat.test)
  }
  
  if (is.null(missingData)) {
    yhat.obs <- yhat.train
    yhat.cf  <- yhat.test      
  } else {
    yhat.obs <- array(0, c(n, dim(yhat.train)[-1L]))
    addDimsToSubset(yhat.obs[!missingRows] <- yhat.train)
    addDimsToSubset(yhat.obs[ missingRows] <- yhat.test[seq.int(n + 1L, n + n.mis),drop = FALSE])
    
    addDimsToSubset(yhat.cf <- yhat.test[seq_len(n)])
  }
  
  sd.obs <- apply(yhat.obs, 1L, sd)
  sd.cf  <- apply(yhat.cf,  1L, sd)
  
  commonSup.sub <- getCommonSupportSubset(sd.obs, sd.cf, commonSup.rule, commonSup.cut, trt, missingRows)
  
  samples.est <- NULL
  if (calculateEstimates) {
    samples.indiv.diff <- (yhat.obs - yhat.cf) * ifelse(treatmentRows, 1, -1)
    
    if (is.null(matchedCall$group.by)) {
      samples.est <- getBartEstimates(treatmentRows, weights, estimand, samples.indiv.diff, commonSup.sub)
    } else {
      samples.est <- lapply(levels(group.by), function(level) {
        levelRows <- group.by == level
        if (!is.null(weights)) weights <- weights[levelRows]
        
        getBartEstimates(treatmentRows[levelRows], weights, estimand,
                         if (length(dim(bartFit$yhat.train)) > 2L) samples.indiv.diff[levelRows,,] else samples.indiv.diff[levelRows,],
                         commonSup.sub[levelRows])
      })
      names(samples.est) <- levels(group.by)
    }
  }
  
  namedList(fit = bartFit, data = responseData, yhat.obs, yhat.cf, samples.est, name.trt = treatmentName, trt, sd.obs, sd.cf, commonSup.sub, missingRows, fitPars = NULL)
}

boundValues <- function(x, bounds){
  x[x > max(bounds)] <- max(bounds)
  x[x < min(bounds)] <- min(bounds)
  x
}

getPWeightEstimates <- function(y, z, weights, estimand, yhat.0, yhat.1, p.score, yBounds, p.scoreBounds)
{
  if (!is.character(estimand) || estimand[1L] %not_in% c("ate", "att", "atc"))
    stop("estimand must be one of 'ate', 'att', or 'atc'")
  estimand <- estimand[1L]
  
  if (!is.null(weights)) {
    weights <- rep_len(weights, length(y))
    weights <- weights / sum(weights)
  }
  
  m <- min(y, na.rm = TRUE)
  M <- max(y, na.rm = TRUE)
  
  ## map yhat to (0, 1)
  yhat.0 <- boundValues((boundValues(yhat.0, c(m, M)) - m) / (M - m), yBounds)
  yhat.1 <- boundValues((boundValues(yhat.1, c(m, M)) - m) / (M - m), yBounds)
  
  origDims <- dim(yhat.0)
  
  indiv.diff <- flattenSamples(yhat.1) - flattenSamples(yhat.0)
  
  p.score <- boundValues(p.score, p.scoreBounds)
  
  if (!is.null(dim(p.score))) {
    if (!all(dim(p.score) == origDims))
      stop("dimensions of p.score samples must match that of observations")
      
    p.score <- flattenSamples(p.score)
  }
    
  getPWeightEstimate <- getPWeightFunction(estimand, weights, indiv.diff, p.score)
  
  result <- getPWeightEstimate(z, weights, indiv.diff, p.score) * (M - m)
  
  if (!is.null(origDims) && length(origDims) > 2L)
    matrix(result, origDims[2L], origDims[3L])
  else 
    result
}

getPWeightResponseFit <-
  function(response, treatment, confounders, data, subset, weights, estimand, group.by, p.score, samples.p.score, use.rbart,
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
  
  bartFit <- responseData <- yhat.obs <- yhat.cf <- name.trt <- trt <- sd.obs <- sd.cf <- commonSup.sub <- missingRows <- NULL
  massign[bartFit, responseData, yhat.obs, yhat.cf,, name.trt, trt, sd.obs, sd.cf, commonSup.sub,missingRows] <- eval(bartCall, envir = callingEnv)
  
  treatmentRows <- trt > 0
  
  yhat.1 <- yhat.obs * trt       + yhat.cf * (1 - trt)
  yhat.0 <- yhat.obs * (1 - trt) + yhat.cf * trt
  
  p.score <- if (!is.null(matchedCall$samples.p.score) && !is.null(samples.p.score)) samples.p.score else p.score
  
  if (is.null(matchedCall$group.by)) {
    if (any(commonSup.sub != TRUE)) {
      addDimsToSubset(yhat.0 <- yhat.0[commonSup.sub, drop = FALSE])
      addDimsToSubset(yhat.1 <- yhat.1[commonSup.sub, drop = FALSE])
           
      p.score <- addDimsToSubset(p.score[commonSup.sub, drop = FALSE])
      
      if (!is.null(weights)) weights <- weights[commonSup.sub]
    }
      
    samples.est <- getPWeightEstimates(bartFit$y[commonSup.sub], trt[commonSup.sub], weights, estimand, yhat.0, yhat.1, p.score, yBounds, p.scoreBounds)
  } else {
    samples.est <- lapply(levels(group.by), function(level) {
      levelRows <- group.by == level & commonSup.sub
      
      addDimsToSubset(yhat.0 <- yhat.0[levelRows, drop = FALSE])
      addDimsToSubset(yhat.1 <- yhat.1[levelRows, drop = FALSE])
      addDimsToSubset(p.score <- p.score[levelRows, drop = FALSE])
      
      if (!is.null(weights)) weights <- weights[levelRows]
      
      getPWeightEstimates(bartFit$y[levelRows], trt[levelRows], weights, estimand, yhat.0, yhat.1, p.score, yBounds, p.scoreBounds)
    })
    names(samples.est) <- levels(group.by)
  }
  
  namedList(fit = bartFit, data = responseData, yhat.obs, yhat.cf, samples.est, name.trt, trt, sd.obs, sd.cf, commonSup.sub, missingRows, fitPars = namedList(yBounds, p.scoreBounds))
}

getTMLEEstimates <- function(y, z, weights, estimand, yhat.0, yhat.1, p.score, yBounds, p.scoreBounds, depsilon, maxIter)
{
  if (!is.character(estimand) || estimand[1L] %not_in% c("ate", "att", "atc"))
    stop("estimand must be one of 'ate', 'att', or 'atc'")
  estimand <- estimand[1L]
  
  if (anyNA(y)) {
    completeRows <- !is.na(y)
    y <- y[completeRows]
    z <- z[completeRows]
    if (!is.null(weights)) weights <- weights[completeRows]
    addDimsToSubset(yhat.0 <- yhat.0[completeRows, drop = FALSE])
    addDimsToSubset(yhat.1 <- yhat.1[completeRows, drop = FALSE])
    addDimsToSubset(p.score <- p.score[completeRows, drop = FALSE])
  }
  
  if (!is.null(weights)) {
    weights <- rep_len(weights, length(y))
    weights <- weights / sum(weights)
  } else {
    yhat.0 <- flattenSamples(yhat.0)
    yhat.1 <- flattenSamples(yhat.1)
    p.score <- boundValues(flattenSamples(p.score), p.scoreBounds)
    return(sapply(seq_len(ncol(yhat.0)), function(i) {
      res <- tmle::tmle(Y = y, A = z, W = matrix(0, length(y), 1L), Q = cbind(Q0W = yhat.0[,i], Q1W = yhat.1[,i]), g1W = if (!is.null(dim(p.score))) p.score[,i] else p.score)
      res$estimates[[switch(estimand, ate = "ATE", att = "ATT", atc = "ATC")]]$psi
    }))
  }
  
  r <- range(y)
  r <- r + 0.1 * c(-abs(r[1L]), abs(r[2L]))
  y.st <- boundValues(y, r)
  yhat.0.st <- boundValues(yhat.0, r)
  yhat.1.st <- boundValues(yhat.1, r)
  r.st <- range(y.st)
  
  y.st <- (y.st - min(r.st)) / diff(r.st)
  yhat.0.st <- qlogis(boundValues((yhat.0.st - min(r.st)) / diff(r.st), yBounds))
  yhat.1.st <- qlogis(boundValues((yhat.1.st - min(r.st)) / diff(r.st), yBounds))
  
  
  yhat.0.samp <- flattenSamples(yhat.0.st)
  yhat.1.samp <- flattenSamples(yhat.1.st)
  
  p.score.samp <- boundValues(flattenSamples(p.score), p.scoreBounds)
  
  origDims <- dim(yhat.0)
  
  
  getPWeightEstimate <- getPWeightFunction(estimand, weights, numeric(), numeric())
  
  getYhat0Deriv <- getYhat1Deriv <- getPScoreDeriv <- getIC <- calcLoss <- NULL
  massign[getYhat0Deriv, getYhat1Deriv, getPScoreDeriv, getIC, calcLoss] <-
    getTMLEFunctions(estimand, weights)
  
  result <- sapply(seq_len(ncol(yhat.0.samp)), function(i) {
    yhat.0 <- yhat.0.samp[,i]
    yhat.1 <- yhat.1.samp[,i]
    yhat <- yhat.1 * z + yhat.0 * (1 - z)
    
    p.score <- if (!is.null(dim(p.score.samp))) p.score.samp[,i] else p.score.samp
    p.score.st <- boundValues(p.score, p.scoreBounds)
    
    H1W <- z / p.score.st
    H0W <- (1 - z) / (1 - p.score.st)
    
    suppressWarnings(epsilon <- coef(glm(y.st ~ -1 + offset(yhat) + H0W + H1W, family = binomial)))
    epsilon[is.na(epsilon)] <- 0 
  
    yhat.0 <- plogis(yhat.0 + epsilon["H0W"] / (1 - p.score.st))
    yhat.1 <- plogis(yhat.1 + epsilon["H1W"] / p.score.st)
    
    indiv.diff <- yhat.1 - yhat.0
    yhat <- yhat.1 * z + yhat.0 * (1 - z)
    
    psi <- getPWeightEstimate(z, weights, indiv.diff, p.score)
    psi.prev <- psi
    
    a.weight <- z * getYhat1Deriv(z, weights, p.score) + (1 - z) * getYhat0Deriv(z, weights, p.score)
    ic <- getIC(y.st, yhat, indiv.diff, psi, a.weight)
    
    if (mean(ic) > 0) depsilon <- -depsilon
    
    loss.prev <- Inf
    loss <- calcLoss(y.st, z, yhat, p.score, weights)
    if (is.nan(loss) || is.na(loss) || is.infinite(loss)) return(psi)
      
    iter <- 0L
    while (loss.prev > loss && iter < maxIter)
    {
      p.score.prev <- p.score
      p.score <- boundValues(plogis(qlogis(p.score.prev) - depsilon * getPScoreDeriv(z, weights, p.score.prev, indiv.diff, psi.prev)), p.scoreBounds)
      
      yhat.0.prev <- yhat.0
      yhat.1.prev <- yhat.1
      yhat.0 <- boundValues(plogis(qlogis(yhat.0.prev) - depsilon * getYhat0Deriv(z, weights, p.score.prev)), yBounds)
      yhat.1 <- boundValues(plogis(qlogis(yhat.1.prev) - depsilon * getYhat1Deriv(z, weights, p.score.prev)), yBounds)
      indiv.diff <- yhat.1 - yhat.0
      yhat <- yhat.1 * z + yhat.0 * (1 - z)
      
      psi.prev <- psi
      psi <- getPWeightEstimate(z, weights, indiv.diff, p.score)
      
      loss.prev <- loss
      loss <- calcLoss(y.st, z, yhat, p.score, weights)
      
      ## ic <- getIC(y, yhat, indiv.diff, psi, a.weight)
      
      if (is.nan(loss) || is.infinite(loss) || is.na(loss)) loss <- Inf
    
      iter <- iter + 1L
    }
    
    if (is.infinite(loss) || loss.prev < loss) psi.prev else psi
  })
  
  if (!is.null(origDims) && length(origDims) > 2L)
    matrix(result, origDims[2L], origDims[3L]) * (max(r.st) - min(r.st))
  else
    result * (max(r.st) - min(r.st))
}

getTMLEResponseFit <-
  function(response, treatment, confounders, data, subset, weights, estimand, group.by, p.score, samples.p.score, use.rbart,
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
  
  bartFit <- responseData <- yhat.obs <- yhat.cf <- name.trt <- trt <- sd.obs <- sd.cf <- commonSup.sub <- missingRows <- NULL
  massign[bartFit, responseData, yhat.obs, yhat.cf,, name.trt, trt, sd.obs, sd.cf, commonSup.sub, missingRows] <- eval(bartCall, envir = callingEnv)
  
  treatmentRows <- trt > 0
  
  yhat.1 <- yhat.obs * trt       + yhat.cf * (1 - trt)
  yhat.0 <- yhat.obs * (1 - trt) + yhat.cf * trt
   
  p.score <- if (!is.null(matchedCall$samples.p.score) && !is.null(samples.p.score)) samples.p.score else p.score
  
  maxIter <- round(maxIter, 0)
  
  if (is.null(matchedCall$group.by)) {
    if (any(commonSup.sub != TRUE)) {
      addDimsToSubset(yhat.0 <- yhat.0[commonSup.sub, drop = FALSE])
      addDimsToSubset(yhat.1 <- yhat.1[commonSup.sub, drop = FALSE])
      
      addDimsToSubset(p.score <- p.score[commonSup.sub, drop = FALSE])
      
      if (!is.null(weights)) weights <- weights[commonSup.sub]
    }
    
    samples.est <- getTMLEEstimates(bartFit$y[commonSup.sub], trt[commonSup.sub], weights, estimand, yhat.0, yhat.1, p.score, yBounds, p.scoreBounds, depsilon, maxIter)
  } else {
    samples.est <- lapply(levels(group.by), function(level) {
      levelRows <- group.by == level & commonSup.sub
      
      addDimsToSubset(yhat.0 <- yhat.0[levelRows, drop = FALSE])
      addDimsToSubset(yhat.1 <- yhat.1[levelRows, drop = FALSE])
      
      addDimsToSubset(p.score <- p.score[levelRows, drop = FALSE])
      
      if (!is.null(weights)) weights <- weights[levelRows]
      
      getTMLEEstimates(bartFit$y[levelRows], trt[levelRows], weights, estimand, yhat.0, yhat.1, p.score, yBounds, p.scoreBounds, depsilon, maxIter)
    })
    names(samples.est) <- levels(group.by)
  }

  namedList(fit = bartFit, data = responseData, yhat.obs, yhat.cf, samples.est, name.trt, trt, sd.obs, sd.cf, commonSup.sub, missingRows, fitPars = namedList(yBounds, p.scoreBounds, depsilon, maxIter))
}
