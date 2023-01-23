getGLMTreatmentFit <- function(response, treatment, confounders, parametric, data, subset, weights, group.by = NULL, use.ranef = TRUE, ...)
{
  treatmentIsMissing    <- missing(treatment)
  confoundersAreMissing <- missing(confounders)
  parametricIsMissing   <- missing(parametric)
  dataAreMissing        <- missing(data)
  
  matchedCall <- match.call()
  callingEnv <- parent.frame(1L)
  
  if (treatmentIsMissing)
    stop("'treatment' variable must be specified")
  if (confoundersAreMissing)
    stop("'confounders' variable must be specified")
  
  useLmer <-
    (!is.null(matchedCall[["group.by"]]) && use.ranef) ||
    (!is.null(matchedCall[["parametric"]]) && anyBars(matchedCall[["parametric"]]))
  
  if (!useLmer) {
    fn <- quote(stats::glm)
  } else {
    if (requireNamespace("lme4", quietly = TRUE) == FALSE)
      stop("random effect model for glm treatment requires lme4 package to be available")
    fn <- quote(lme4::glmer)
  }

  glmCall <- evalEnv <- NULL
  if (!dataAreMissing && is.data.frame(data)) {
    dataEnv <- NULL
    if (grepl("(?<![\\w\\.])\\.(?![\\w\\.])", deparse(matchedCall$confounders), perl = TRUE)) {
      mfCall <- quote(stats::model.frame(a ~ b, data = data))
      mfCall[[2L]][[2L]] <- matchedCall$response
      mfCall[[2L]][[3L]] <- matchedCall$confounders
      mf <- eval(mfCall)
      data <- data[,setdiff(colnames(data), deparse(attr(attr(mf, "terms"), "variables")[[2L]]))]
      dataEnv <- new.env(parent = callingEnv)
      dataEnv[[deparse(matchedCall$data)]] <- data
    }
    
    evalEnv <- NULL
    dataCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentDataCall)), "fn", fn)
    dataCall <- addCallDefaults(dataCall, eval(quoteInNamespace(getGLMTreatmentFit)))
    dataCall[["use.lmer"]] <- useLmer
    
    massign[glmCall, evalEnv] <- eval(dataCall, envir = callingEnv)
    
    if (!is.null(dataEnv)) evalEnv <- dataEnv

    df <- evalEnv[[as.character(dataCall$data)]]
    treatment <- df[[matchedCall$treatment]]
    if (!all(treatment %in% c(0, 1)))
      stop("response must be in {0, 1}")
  } else {
    df <- NULL
    literalCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentLiteralCall)), "fn", fn)
    literalCall <- addCallDefaults(literalCall, eval(quoteInNamespace(getGLMTreatmentFit)))
    
    dataEnv <- if (dataAreMissing) callingEnv else list2env(data, parent = callingEnv)
    literalCall[["use.lmer"]] <- useLmer
    
    massign[glmCall, df] <- eval(literalCall, envir = dataEnv)

    treatment <- df[[deparse(glmCall[[2L]][[2L]])]]
    if (!all(treatment %in% c(0, 1)))
      stop("response must be in {0, 1}")
    
    evalEnv <- new.env(parent = callingEnv)
    evalEnv[["df"]] <- df
  }
  
  extraArgs <- matchedCall[names(matchedCall) %not_in% names(glmCall) | names(matchedCall) == ""]
  
  glmCall <- addCallArguments(glmCall, extraArgs)
  if (is.null(glmCall[["family"]])) glmCall[["family"]] <- quote(stats::binomial)
  
  glmFit <- eval(glmCall, envir = evalEnv)
  
  list(fit = glmFit, p.score = fitted(glmFit), samples = NULL)
}

getBartTreatmentFit <- function(response, treatment, confounders, parametric, data, subset, weights, group.by = NULL, use.ranef = TRUE,
                                crossvalidate = FALSE, ...)
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
  
  bartMethod <- "bart"
  fn <- quote(dbarts::bart2)
  if (!is.null(matchedCall[["parametric"]])) {
    if (!is.null(matchedCall[["group.by"]]))
      stop("`group.by` must be missing or NULL if `parametric` is supplied; for varying intercepts, add (1 | group) to parametric equation")
    if (requireNamespace("stan4bart", quietly = TRUE) == FALSE)
      stop("semiparametric BART treatment model requires stan4bart package to be available")
    fn <- quote(stan4bart::stan4bart)
    bartMethod <- "stan4bart"
  } else if (!is.null(matchedCall[["group.by"]]) && use.ranef) {
    fn <- quote(dbarts::rbart_vi)
    bartMethod <- "rbart"
  }
  
  if (crossvalidate && bartMethod %not_in% "bart")
    stop("crossvalidation not yet supported for varying intercept or semiparametric BART models")
  
  bartCall <- NULL
  if (!dataAreMissing && is.data.frame(data)) {
    # if the confounders contain a '.', they can end up including the response variable.
    # This regex looks for dots that are not preceeded or followed by valid variable name
    # name characters. It can miss some cases, e.g. '_.' doesn't parse, but it should work
    # most of the time.
    dataEnv <- NULL
    if (grepl("(?<![\\w\\.])\\.(?![\\w\\.])", deparse(matchedCall$confounders), perl = TRUE)) {
      mfCall <- quote(stats::model.frame(a ~ b, data = data))
      mfCall[[2L]][[2L]] <- matchedCall$response
      mfCall[[2L]][[3L]] <- matchedCall$confounders
      mf <- eval(mfCall)
      data <- data[,setdiff(colnames(data), deparse(attr(attr(mf, "terms"), "variables")[[2L]]))]
      dataEnv <- new.env(parent = callingEnv)
      dataEnv[[deparse(matchedCall$data)]] <- data
    }
    dataCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentDataCall)), "fn", fn)
    dataCall <- addCallDefaults(dataCall, eval(quoteInNamespace(getBartTreatmentFit)))
    dataCall[["use.lmer"]] <- FALSE
    
    massign[bartCall, evalEnv] <- eval(dataCall, envir = callingEnv)
    
    if (!is.null(dataEnv)) evalEnv <- dataEnv

    df <- evalEnv[[as.character(dataCall$data)]]
    treatment <- df[[matchedCall$treatment]]
    if (!all(treatment %in% c(0, 1)))
      stop("response must be in {0, 1}")
  } else {
    df <- NULL
    literalCall <- addCallArgument(redirectCall(matchedCall, quoteInNamespace(getTreatmentLiteralCall)), "fn", fn)
    literalCall <- addCallDefaults(literalCall, eval(quoteInNamespace(getBartTreatmentFit)))
    literalCall[["use.lmer"]] <- FALSE
    
    dataEnv <- if (dataAreMissing) callingEnv else list2env(data, parent = callingEnv)
    
    massign[bartCall, df] <- eval(literalCall, envir = dataEnv)

    treatment <- df[[deparse(bartCall[[2L]][[2L]])]]
    if (!all(treatment %in% c(0, 1)))
      stop("response must be in {0, 1}")
    
    evalEnv <- new.env(parent = callingEnv)
    evalEnv[["df"]] <- df
  }
  extraArgs <- matchedCall[names(matchedCall) %not_in% names(bartCall) | names(matchedCall) == ""]
  
  bartCall$verbose <- if (bartMethod %in% "stan4bart") -1L else FALSE
  bartCall <- addCallArguments(bartCall, extraArgs)
  
  chainsArgument <- if (bartMethod %in% "stan4bart") "chains" else "n.chains"
  if (is.null(bartCall[[chainsArgument]])) bartCall[[chainsArgument]] <- 10L
  
  if (crossvalidate)
    bartCall <- optimizeBARTCall(bartCall, evalEnv)
  
  bartFit <- eval(bartCall, envir = evalEnv)
  combineChains <- if (is.null(matchedCall[["combineChains"]])) FALSE else list(...)[["combineChains"]]
  
  if (bartMethod %in% "stan4bart")
    samples <- extract(bartFit, combine_chains = combineChains)
  else
    samples <- extract(bartFit, combineChains = combineChains)
  
  result <- 
    list(fit = bartFit,
         p.score = apply(samples, length(dim(samples)), mean),
         samples = samples)
  
  if (crossvalidate)
    result[["k"]] <- bartCall[["k"]]
  
  result
}

