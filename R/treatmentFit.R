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
  
  ## a parametric equation or a modeled group intercept both make the propensity
  ## model semiparametric, and stan4bart is the only sampler that fits one; its
  ## binary family is probit, matching the dbarts binary route
  bartMethod <- "bart"
  fn <- quote(dbarts::bart2)
  if (!is.null(matchedCall[["parametric"]]) || (!is.null(matchedCall[["group.by"]]) && use.ranef)) {
    if (requireNamespace("stan4bart", quietly = TRUE) == FALSE)
      stop("semiparametric BART treatment model, including a varying intercept from 'group.by' with use.ranef = TRUE, requires stan4bart package to be available; pass use.ranef = FALSE to enter the grouping factor as a fixed effect instead")
    fn <- quote(stan4bart::stan4bart)
    bartMethod <- "stan4bart"
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
  
  if (bartMethod %in% "stan4bart") {
    bartCall <- addStan4BartSamplingArguments(bartCall, matchedCall, callingEnv)
  } else if (is.null(bartCall[["n.chains"]])) {
    bartCall[["n.chains"]] <- 10L
  }

  ## a propensity model has no counterfactual surface, so the treatment name the
  ## redirect carries over is dropped rather than fit as a test set
  if (bartMethod %in% "stan4bart") bartCall[["treatment"]] <- NULL

  ## the treatment is binary, so every BART-backed propensity model is a probit,
  ## and dbarts refuses a weighted probit (no tractable latent form) whether it
  ## is reached directly or through stan4bart; the weights carry the design
  ## information in the treatment-effect estimators (p.weights, tmle).
  if (bartMethod %in% c("bart", "stan4bart") && !is.null(bartCall[["weights"]])) {
    bartCall[["weights"]] <- NULL
    message("propensity score model is fit unweighted; weights enter the treatment-effect estimators")
  }

  if (crossvalidate)
    bartCall <- optimizeBARTCall(bartCall, evalEnv)
  
  bartFit <- eval(bartCall, envir = evalEnv)
  combineChains <- if (is.null(matchedCall[["combineChains"]])) FALSE else list(...)[["combineChains"]]
  
  if (bartMethod %in% "stan4bart") {
    samples <- extract(bartFit, combine_chains = combineChains)
    if (length(dim(samples)) == 3L) {
      samples <- aperm(samples, c(3L, 2L, 1L))
    } else {
      samples <- t(samples)
    }
  } else {
    samples <- extract(bartFit, combineChains = combineChains)
  }
  
  result <- 
    list(fit = bartFit,
         p.score = apply(samples, length(dim(samples)), mean),
         samples = samples)
  
  if (crossvalidate)
    result[["k"]] <- bartCall[["k"]]
  
  result
}

