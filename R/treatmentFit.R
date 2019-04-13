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

getBartTreatmentFit <- function(treatment, confounders, data, subset, weights, crossvalidateBinary = FALSE, ...)
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
  if (is.null(bartCall[["n.chains"]])) bartCall[["n.chains"]] <- 10L
  # if (is.null(bartCall[["k"]])) bartCall[["k"]] <- "chi(1.25, Inf)"
  if (crossvalidateBinary)
    bartCall <- optimizeBARTCall(bartCall, evalEnv)
  
  bartFit <- eval(bartCall, envir = evalEnv)
  x <- NULL ## R CMD check
  samples <- evalx(pnorm(bartFit$yhat.train), if (length(dim(x)) > 2L) aperm(x, c(3L, 1L, 2L)) else t(x))
  
  list(fit = bartFit,
       p.score = apply(samples, 1L, mean),
       samples = samples)
}

