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
    
    evalEnv <- new.env(parent = callingEnv)
    evalEnv[["df"]] <- df
  }
  extraArgs <- matchedCall[names(matchedCall) %not_in% names(glmCall) | names(matchedCall) == ""]
  
  glmCall <- addCallArgument(glmCall, 2L, quote(stats::binomial))
  glmCall <- addCallArguments(glmCall, extraArgs)
  
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
    
    evalEnv <- new.env(parent = callingEnv)
    evalEnv[["df"]] <- df
  }
  extraArgs <- matchedCall[names(matchedCall) %not_in% names(bartCall) | names(matchedCall) == ""]
  
  bartCall$verbose <- FALSE
  bartCall <- addCallArguments(bartCall, extraArgs)
  if (!is.null(bartCall[["n.chains"]])) bartCall[["n.chains"]] <- 10L
  
  bartFit <- eval(bartCall, envir = evalEnv)
  x <- NULL ## R CMD check
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
    
    evalEnv <- new.env(parent = callingEnv)
    evalEnv[["df"]] <- df
  }
  extraArgs <- matchedCall[names(matchedCall) %not_in% names(xbartCall) | names(matchedCall) == ""]
  dotsList <- list(...)
  
  xbartCall[["verbose"]] <- FALSE
  xbartCall[["k"]] <- TEST_K
  bartCall <- xbartCall
  bartCall[[1L]] <- quote(dbarts::bart2)
  
  xbartCall <- addCallArguments(xbartCall, extraArgs)
  bartCall  <- addCallArguments(bartCall, extraArgs)
  
  # these are common to both and could have been specified twice in a list format
  for (argName in c("n.samples", "n.burn", "n.threads", "n.trees", "k", "power", "base")) {
    if (!is.null(matchedCall[[argName]]) && is.list(dotsList[[argName]]) && length(dotsList[[argName]]) > 1L) {
      xbartCall[[argName]] <- dotsList[[argName]][[1L]]
      bartCall[[argName]]  <- dotsList[[argName]][[2L]]
    }
  }
  
  xbartFit <- eval(xbartCall, envir = evalEnv)
  xval <- apply(xbartFit, seq.int(2L, length(dim(xbartFit))), mean)
  if (is.null(dim(xval))) {
    dim(xval) <- length(xval)
    dimnames(xval) <- dimnames(xbartFit)[-1L]
  }
  xvalOpt <- which.min(xval)
  
  xvalInd <- getArrayIndicesForOffset(xvalOpt, dim(xval))
  for (i.dim in seq_len(length(dim(xval)))) {
    parName <- names(dimnames(xval))[i.dim]
    
    bartCall[[parName]] <-
      if (is.numeric(xbartCall[[parName]])) xbartCall[[parName]][xvalInd[i.dim]]
      else                                  dotsList[[parName]][xvalInd[i.dim]]
  }
  
  if (!is.null(bartCall[["n.chains"]])) bartCall[["n.chains"]] <- 10L
  
  bartFit <- eval(bartCall, envir = evalEnv)
  x <- NULL ## R CMD check
  samples <- evalx(pnorm(bartFit$yhat.train), if (length(dim(x)) > 2L) aperm(x, c(3L, 1L, 2L)) else t(x))
  
  list(fit = bartFit,
       p.score = apply(samples, 1L, mean),
       samples = samples)
}

