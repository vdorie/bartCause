getGLMTreatmentFit <- function(treatment, confounders, data, subset, weights, ...)
{
  treatmentIsMissing    <- missing(treatment)
  confoundersAreMissing <- missing(confounders)
  dataAreMissing        <- missing(data)
  
  matchedCall <- match.call()
  callingEnv <- parent.frame(1L)
  
  if (treatmentIsMissing)
    stop("'treatment' variable must be specified")
  if (confoundersAreMissing)
    stop("'confounders' variable must be specified")
  
  glmCall <- evalEnv <- NULL
  if (!dataAreMissing && is.data.frame(data)) {
    evalEnv <- NULL
    dataCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentDataCall)), "fn", quote(stats::glm))
    
    massign[glmCall, evalEnv] <- eval(dataCall, envir = callingEnv)
  } else {
    df <- NULL
    literalCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentLiteralCall)), "fn", quote(stats::glm))
    
    dataEnv <- if (dataAreMissing) callingEnv else list2env(data, parent = callingEnv)
    
    massign[glmCall, df] <- eval(literalCall, envir = dataEnv)
    
    evalEnv <- sys.frame(sys.nframe())
  }
  
  glmCall <- addCallArgument(glmCall, 2L, quote(stats::binomial))
  glmCall <- addCallArguments(glmCall, list(...))
  
  glmFit <- eval(glmCall, envir = evalEnv)
  
  list(fit = glmFit, p.score = fitted(glmFit), samples = NULL)
}

getBartTreatmentFit <- function(treatment, confounders, data, subset, weights, ...)
{
  treatmentIsMissing    <- missing(treatment)
  confoundersAreMissing <- missing(confounders)
  dataAreMissing        <- missing(data)
  
  matchedCall <- match.call()
  callingEnv <- parent.frame(1L)
  
  if (treatmentIsMissing)
    stop("'treatment' variable must be specified")
  if (confoundersAreMissing)
    stop("'confounders' variable must be specified")
  
  bartCall <- NULL
  if (!dataAreMissing && is.data.frame(data)) {
    evalEnv <- NULL
    dataCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentDataCall)), "fn", quote(dbarts::bart2))
    
    massign[bartCall, evalEnv] <- eval(dataCall, envir = callingEnv)
  } else {
    df <- NULL
    literalCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentLiteralCall)), "fn", quote(dbarts::bart2))
    
    dataEnv <- if (dataAreMissing) callingEnv else list2env(data, parent = callingEnv)
    
    massign[bartCall, df] <- eval(literalCall, envir = dataEnv)
    
    evalEnv <- sys.frame(sys.nframe())
  }
  
  bartCall$verbose <- FALSE
  bartCall <- addCallArguments(bartCall, list(...))
  
  bartFit <- eval(bartCall, envir = evalEnv)
  samples <- evalx(pnorm(bartFit$yhat.train), if (length(dim(x)) > 2L) aperm(x, c(3L, 1L, 2L)) else t(x))
  
  list(fit = bartFit,
       p.score = apply(samples, 1L, mean),
       samples = samples)
}

TEST_K <- c(0.5, 1, 2, 4, 8)

getBartXValTreatmentFit <- function(treatment, confounders, data, subset, weights, ...)
{
  treatmentIsMissing    <- missing(treatment)
  confoundersAreMissing <- missing(confounders)
  dataAreMissing        <- missing(data)
  
  matchedCall <- match.call()
  callingEnv <- parent.frame(1L)
  
  if (treatmentIsMissing)
    stop("'treatment' variable must be specified")
  if (confoundersAreMissing)
    stop("'confounders' variable must be specified")
  
  xbartCall <- NULL
  if (!dataAreMissing && is.data.frame(data)) {
    evalEnv <- NULL
    dataCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentDataCall)), "fn", quote(dbarts::xbart))
    
    massign[xbartCall, evalEnv] <- eval(dataCall, envir = callingEnv)
  } else {
    df <- NULL
    literalCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentLiteralCall)), "fn", quote(dbarts::xbart))
    
    dataEnv <- if (dataAreMissing) callingEnv else list2env(data, parent = callingEnv)
    
    massign[xbartCall, df] <- eval(literalCall, envir = dataEnv)
    
    evalEnv <- sys.frame(sys.nframe())
  }
  xbartCall$verbose <- FALSE
  xbartCall$k <- TEST_K
  bartCall <- xbartCall
  bartCall[[1L]] <- quote(dbarts::bart2)
  
  dotsList <- list(...)
  
  xbartCall <- addCallArguments(xbartCall, dotsList)
  bartCall  <- addCallArguments(bartCall, dotsList)
  
  if (!is.null(matchedCall$n.burn) && is.list(dotsList[["n.burn"]]) && length(dotsList[["n.burn"]]) > 1L) {
    xbartCall$n.burn <- dotsList[["n.burn"]][[1L]]
    bartCall$n.burn  <- dotsList[["n.burn"]][[2L]]
  }
  if (!is.null(matchedCall$n.samples) && is.list(dotsList[["n.samples"]]) && length(dotsList[["n.samples"]]) > 1L) {
    xbartCall$n.samples <- dotsList[["n.samples"]][[1L]]
    bartCall$n.samples  <- dotsList[["n.samples"]][[2L]]
  }
  
  xbartFit <- eval(xbartCall, envir = evalEnv)
  
  bartCall$k <- TEST_K[which.min(apply(xbartFit, 2L, mean))]
  
  bartFit <- eval(bartCall, envir = evalEnv)
  samples <- evalx(pnorm(bartFit$yhat.train), if (length(dim(x)) > 2L) aperm(x, c(3L, 1L, 2L)) else t(x))
  
  list(fit = bartFit,
       p.score = apply(samples, 1L, mean),
       samples = samples)
}

