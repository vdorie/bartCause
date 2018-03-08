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
                               commonSup.rule, commonSup.cut, p.score, calculateEstimates = TRUE, ...)
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
  
  dbartsDataCall <- NULL; treatmentName <- NULL
  if (!dataAreMissing && is.data.frame(data)) {
    evalEnv <- NULL
    dataCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getResponseDataCall)), "fn", quote(dbarts::dbartsData))
    
    massign[dbartsDataCall, evalEnv, treatmentName] <- eval(dataCall, envir = callingEnv)
  } else {
    df <- NULL
    literalCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getResponseLiteralCall)), "fn", quote(dbarts::dbartsData))
    
    dataEnv <- if (dataAreMissing) callingEnv else list2env(data, parent = callingEnv)
    
    massign[dbartsDataCall, df, treatmentName] <- eval(literalCall, envir = dataEnv)
    
    evalEnv <- sys.frame(sys.nframe())
  }
  
  responseData <- eval(dbartsDataCall, envir = evalEnv)
  
  trt <- responseData@x[,treatmentName]
  treatmentRows <- trt > 0
  responseData@x.test <- responseData@x
  responseData@x.test[,treatmentName] <- 1 - responseData@x.test[,treatmentName]
  
  ## redict to pull in any args passed
  bartCall <- redirectCall(matchedCall, dbarts::bart2)
  bartCall$formula <- quote(responseData)
  bartCall$data    <- NULL
  bartCall$keepTrees <- FALSE
  
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
  
  if (length(dim(bartFit$yhat.train)) > 2L) {
    yhat.train <- aperm(bartFit$yhat.train, c(3L, 1L, 2L))
    yhat.test  <- aperm(bartFit$yhat.test, c(3L, 1L, 2L))
    
    samples.indiv.diff <- yhat.train - yhat.test
    samples.indiv.diff[!treatmentRows,,] <- -samples.indiv.diff[!treatmentRows,,]
  } else {
    yhat.train <- t(bartFit$yhat.train)
    yhat.test  <- t(bartFit$yhat.test)
    
    samples.indiv.diff <- yhat.train - yhat.test
    samples.indiv.diff[!treatmentRows,] <- -samples.indiv.diff[!treatmentRows,]
  }
  
  sd.obs <- apply(yhat.train, 1L, sd)
  sd.cf  <- apply(yhat.test,  1L, sd)
  
  commonSup.sub <- getCommonSupportSubset(sd.obs, sd.cf, commonSup.rule, commonSup.cut, trt)
  
  samples.est <- NULL
  if (calculateEstimates) {
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
  
  namedList(fit = bartFit, data = responseData, samples.est, samples.indiv.diff, name.trt = treatmentName, trt, sd.obs, sd.cf, commonSup.sub)
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
  
  m <- min(y)
  M <- max(y)
  
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
  function(response, treatment, confounders, data, subset, weights, estimand, group.by, p.score, samples.p.score,
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
  
  bartFit <- responseData <- samples.indiv.diff <- name.trt <- trt <- sd.obs <- sd.cf <- commonSup.sub <- NULL
  massign[bartFit, responseData,, samples.indiv.diff, name.trt, trt, sd.obs, sd.cf, commonSup.sub] <- eval(bartCall, envir = callingEnv)
  
  treatmentRows <- trt > 0
  
  if (length(dim(bartFit$yhat.train)) > 2L) {
    trainingSamples <- aperm(bartFit$yhat.train, c(3L, 1L, 2L))
    testSamples     <- aperm(bartFit$yhat.test, c(3L, 1L, 2L))
    
    yhat.0 <- yhat.1 <- trainingSamples
    
    yhat.0[ treatmentRows,,] <- testSamples[ treatmentRows,,]
    yhat.1[!treatmentRows,,] <- testSamples[!treatmentRows,,]
  } else {
    trainingSamples <- t(bartFit$yhat.train)
    testSamples     <- t(bartFit$yhat.test)
    
    yhat.0 <- yhat.1 <- trainingSamples
    
    yhat.0[ treatmentRows,] <- testSamples[ treatmentRows,]
    yhat.1[!treatmentRows,] <- testSamples[!treatmentRows,]
  }
  
  p.score <- if (!is.null(matchedCall$samples.p.score) && !is.null(samples.p.score)) samples.p.score else p.score
  
  if (is.null(matchedCall$group.by)) {
    if (any(commonSup.sub != TRUE)) {
      if (length(dim(yhat.0)) > 2L) {
        yhat.0 <- yhat.0[commonSup.sub,,]; yhat.1 <- yhat.0[commonSup.sub,,]
      } else {
        yhat.0 <- yhat.0[commonSup.sub,];  yhat.1 <- yhat.1[commonSup.sub,]
      }
      
      p.score <- if (is.null(dim(p.score)))
        p.score[commonSup.sub]
      else {
        if (length(dim(p.score)) > 2L) p.score[commonSup.sub,,] else p.score[commonSup.sub,]
      }
      
      if (!is.null(weights)) weights <- weights[commonSup.sub]
    }
      
    samples.est <- getPWeightEstimates(bartFit$y[commonSup.sub], trt[commonSup.sub], weights, estimand, yhat.0, yhat.1, p.score, yBounds, p.scoreBounds)
  } else {
    samples.est <- lapply(levels(group.by), function(level) {
      levelRows <- group.by == level & commonSup.sub
      
      yhat.0 <- if (length(dim(yhat.0)) > 2L) yhat.0[levelRows,,] else yhat.0[levelRows,]
      yhat.1 <- if (length(dim(yhat.0)) > 2L) yhat.1[levelRows,,] else yhat.1[levelRows,]
      
      p.score <- if (is.null(dim(p.score))) {
        p.score[levelRows]
      } else {
        if (length(dim(p.score)) > 2L) p.score[levelRows,,] else p.score[levelRows,]
      }
      if (!is.null(weights)) weights <- weights[levelRows]
      
      getPWeightEstimates(bartFit$y[levelRows], trt[levelRows], weights, estimand, yhat.0, yhat.1, p.score, yBounds, p.scoreBounds)
    })
    names(samples.est) <- levels(group.by)
  }
  
  namedList(fit = bartFit, data = responseData, samples.est, samples.indiv.diff, name.trt, trt, sd.obs, sd.cf, commonSup.sub)
}

getTMLEEstimates <- function(y, z, weights, estimand, yhat.0, yhat.1, p.score, yBounds, p.scoreBounds, depsilon, maxIter)
{
  if (!is.character(estimand) || estimand[1L] %not_in% c("ate", "att", "atc"))
    stop("estimand must be one of 'ate', 'att', or 'atc'")
  estimand <- estimand[1L]
  
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
  function(response, treatment, confounders, data, subset, weights, estimand, group.by, p.score, samples.p.score,
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
  
  bartFit <- responseData <- samples.indiv.diff <- name.trt <- trt <- sd.obs <- sd.cf <- commonSup.sub <- NULL
  massign[bartFit, responseData,, samples.indiv.diff, name.trt, trt, sd.obs, sd.cf, commonSup.sub] <- eval(bartCall, envir = callingEnv)
  
  treatmentRows <- trt > 0
  
  if (length(dim(bartFit$yhat.train)) > 2L) {
    trainingSamples <- aperm(bartFit$yhat.train, c(3L, 1L, 2L))
    testSamples     <- aperm(bartFit$yhat.test, c(3L, 1L, 2L))
    
    yhat.0 <- yhat.1 <- trainingSamples
    
    yhat.0[ treatmentRows,,] <- testSamples[ treatmentRows,,]
    yhat.1[!treatmentRows,,] <- testSamples[!treatmentRows,,]
  } else {
    trainingSamples <- t(bartFit$yhat.train)
    testSamples     <- t(bartFit$yhat.test)
    
    yhat.0 <- yhat.1 <- trainingSamples
    
    yhat.0[ treatmentRows,] <- testSamples[ treatmentRows,]
    yhat.1[!treatmentRows,] <- testSamples[!treatmentRows,]
  }
  
  p.score <- if (!is.null(matchedCall$samples.p.score) && !is.null(samples.p.score)) samples.p.score else p.score
  
  maxIter <- round(maxIter, 0)
  
  if (is.null(matchedCall$group.by)) {
    if (any(commonSup.sub != TRUE)) {
      if (length(dim(yhat.0)) > 2L) {
        yhat.0 <- yhat.0[commonSup.sub,,]; yhat.1 <- yhat.0[commonSup.sub,,]
      } else {
        yhat.0 <- yhat.0[commonSup.sub,];  yhat.1 <- yhat.1[commonSup.sub,]
      }
      
      p.score <- if (is.null(dim(p.score)))
        p.score[commonSup.sub]
      else {
        if (length(dim(p.score)) > 2L) p.score[commonSup.sub,,] else p.score[commonSup.sub,]
      }
      
      if (!is.null(weights)) weights <- weights[commonSup.sub]
    }
    
    samples.est <- getTMLEEstimates(bartFit$y[commonSup.sub], trt[commonSup.sub], weights, estimand, yhat.0, yhat.1, p.score, yBounds, p.scoreBounds, depsilon, maxIter)
  } else {
    samples.est <- lapply(levels(group.by), function(level) {
      levelRows <- group.by == level & commonSup.sub
      yhat.0 <- if (length(dim(yhat.0)) > 2L) yhat.0[levelRows,,] else yhat.0[levelRows,]
      yhat.1 <- if (length(dim(yhat.0)) > 2L) yhat.1[levelRows,,] else yhat.1[levelRows,]
      
      p.score <- if (is.null(dim(p.score))) {
        p.score[levelRows]
      } else {
        if (length(dim(p.score)) > 2L) p.score[levelRows,,] else p.score[levelRows,]
      }
      if (!is.null(weights)) weights <- weights[levelRows]
      
      getTMLEEstimates(bartFit$y[levelRows], trt[levelRows], weights, estimand, yhat.0, yhat.1, p.score, yBounds, p.scoreBounds, depsilon, maxIter)
    })
    names(samples.est) <- levels(group.by)
  }

  namedList(fit = bartFit, data = responseData, samples.est, samples.indiv.diff, name.trt, trt, sd.obs, sd.cf, commonSup.sub)
}
