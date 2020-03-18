getGLMTreatmentFit <- function(treatment, confounders, data, subset, weights, group.by = NULL, use.ranef = TRUE, ...)
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
  
  if (is.null(matchedCall[["group.by"]]) || !use.ranef) {
    fn <- quote(stats::glm)
  } else {
    if (requireNamespace("lme4", quietly = TRUE) == FALSE)
      stop("random effect model for glm treatment requires lme4 package to be available")
    fn <- quote(lme4::glmer)
  }
  
  glmCall <- evalEnv <- NULL
  if (!dataAreMissing && is.data.frame(data)) {
    evalEnv <- NULL
    
    dataCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentDataCall)), "fn", fn)
    dataCall <- addCallDefaults(dataCall, eval(quoteInNamespace(getGLMTreatmentFit)))
    dataCall[["use.lmer"]] <- !is.null(matchedCall[["group.by"]]) && use.ranef
    
    massign[glmCall, evalEnv] <- eval(dataCall, envir = callingEnv)
  } else {
    df <- NULL
    literalCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentLiteralCall)), "fn", fn)
    literalCall <- addCallDefaults(literalCall, eval(quoteInNamespace(getGLMTreatmentFit)))
    
    dataEnv <- if (dataAreMissing) callingEnv else list2env(data, parent = callingEnv)
    literalCall[["use.lmer"]] <- !is.null(matchedCall$group.by) && use.ranef
    
    massign[glmCall, df] <- eval(literalCall, envir = dataEnv)
    
    evalEnv <- new.env(parent = callingEnv)
    evalEnv[["df"]] <- df
  }
  extraArgs <- matchedCall[names(matchedCall) %not_in% names(glmCall) | names(matchedCall) == ""]
  
  glmCall <- addCallArguments(glmCall, extraArgs)
  if (is.null(glmCall[["family"]])) glmCall[["family"]] <- quote(stats::binomial)
  
  glmFit <- eval(glmCall, envir = evalEnv)
  
  list(fit = glmFit, p.score = fitted(glmFit), samples = NULL)
}

getBartTreatmentFit <- function(treatment, confounders, data, subset, weights, group.by = NULL, use.ranef = TRUE,
                                crossvalidateBinary = FALSE, ...)
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
  
  fn <- if (is.null(matchedCall[["group.by"]]) || !use.ranef)
      quote(dbarts::bart2)
    else 
      quote(dbarts::rbart_vi)
  
  bartCall <- NULL
  if (!dataAreMissing && is.data.frame(data)) {
    evalEnv <- NULL
    dataCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentDataCall)), "fn", fn)
    dataCall <- addCallDefaults(dataCall, eval(quoteInNamespace(getBartTreatmentFit)))
    dataCall[["use.lmer"]] <- FALSE
    
    massign[bartCall, evalEnv] <- eval(dataCall, envir = callingEnv)
  } else {
    df <- NULL
    literalCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentLiteralCall)), "fn", fn)
    literalCall <- addCallDefaults(literalCall, eval(quoteInNamespace(getBartTreatmentFit)))
    literalCall[["use.lmer"]] <- FALSE
    
    dataEnv <- if (dataAreMissing) callingEnv else list2env(data, parent = callingEnv)
    
    massign[bartCall, df] <- eval(literalCall, envir = dataEnv)
    
    evalEnv <- new.env(parent = callingEnv)
    evalEnv[["df"]] <- df
  }
  extraArgs <- matchedCall[names(matchedCall) %not_in% names(bartCall) | names(matchedCall) == ""]
  
  bartCall$verbose <- FALSE
  bartCall <- addCallArguments(bartCall, extraArgs)
  if (is.null(bartCall[["n.chains"]])) bartCall[["n.chains"]] <- 10L
  
  if (crossvalidateBinary)
    bartCall <- optimizeBARTCall(bartCall, evalEnv)
  
  bartFit <- eval(bartCall, envir = evalEnv)
  combineChains <- if (is.null(matchedCall[["combineChains"]])) FALSE else list(...)[["combineChains"]]
  samples <- extract(bartFit, combineChains = combineChains)
  
  list(fit = bartFit,
       p.score = apply(samples, length(dim(samples)), mean),
       samples = samples)
}

