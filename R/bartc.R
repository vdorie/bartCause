bartc <- function(
  response, treatment, confounders, parametric, data, subset, weights,
  method.rsp = c("bart", "bcf", "tmle", "p.weight"),
  method.trt = c("bart", "glm", "none"),
  estimand   = c("ate", "att", "atc"),
  group.by = NULL,
  commonSup.rule = c("none", "sd", "chisq"),
  commonSup.cut  = c(NA_real_, 1, 0.05),
  args.rsp = list(), args.trt = list(),
  p.scoreAsCovariate = TRUE, use.ranef = TRUE, group.effects = FALSE,
  crossvalidate = FALSE,
  keepCall = TRUE, verbose = TRUE,
  seed = NA_integer_, ...
)
{
  matchedCall <- match.call()
  sysCall     <- sys.call()
  callingEnv  <- parent.frame(1L)
  
  # some dots arg can get eaten by R's argument matching algorithm, like 'k' for keepCall
  mismatchedArgs.sys <- names(sysCall) %not_in% names(matchedCall) & names(sysCall) != "" &
                        names(sysCall) %in% names(formals(dbarts::bart2))
  if (any(mismatchedArgs.sys)) {
    mismatchedArgs.mc  <- sapply(matchedCall, function(x)
      any(sapply(which(mismatchedArgs.sys), function(j) x == sysCall[[j]])))
    
    oldNames <- names(matchedCall)[mismatchedArgs.mc]
    newNames <- names(sysCall)[mismatchedArgs.sys]
    names(matchedCall)[mismatchedArgs.mc] <- newNames
    for (i in seq_along(oldNames)) {
      oldValue <- get(oldNames[i])
      assign(oldNames[i], eval(formals(bartc)[[oldNames[i]]]))
    }
  }
  
  givenCall <- if (keepCall) matchedCall else call("NULL")
  matchedCall$verbose <- NULL  
  
  group.byIsLiteral <- FALSE
  
  ## check validity of character vector arguments by comparing to function prototype
  for (argName in c("method.rsp", "estimand", "commonSup.rule")) {
    arg <- get(argName)
    if (!is.character(arg) || arg[1L] %not_in% eval(formals(bartCause::bartc)[[argName]]))
      stop(argName, " must be in '", paste0(eval(formals(bartCause::bartc)[[argName]]), collapse = "', '"), "'")
    assign(argName, arg[1L])
  }
  
  if (length(crossvalidate) != 1L ||
      (!is.logical(crossvalidate) && !is.character(crossvalidate)) ||
      (is.logical(crossvalidate) && is.na(crossvalidate)) ||
      (is.character(crossvalidate) && crossvalidate %not_in% c("rsp", "trt")))
    stop("crossvalidate must be one of TRUE, FALSE, 'rsp', or 'trt'")
  
  if (!is.na(seed)) {
    oldSeed <- .GlobalEnv[[".Random.seed"]]
    set.seed(seed)
    matchedCall[["seed"]] <- NULL
  }
  
  fit.trt <- p.score <- samples.p.score <- NULL
  if (is.numeric(method.trt)) {
    if (!is.null(dim(method.trt))) {
      samples.p.score <- method.trt
      p.score <- apply(samples.p.score, 1L, mean)
    } else  {
      p.score <- method.trt
    }
    method.trt <- "given"
  } else if (is.character(method.trt)) {
    method.trt <- method.trt[1L]
    if (method.trt %not_in% eval(formals(bartCause::bartc)$method.trt))
      stop("method.trt must be in '", paste0(eval(formals(bartCause::bartc)$method.trt), collapse = "', '"), "'")
    
    if (method.trt %not_in% c("none", "glm") && !is.null(matchedCall[["group.by"]]) && use.ranef) {
      group.by <- eval(redirectCall(matchedCall, quoteInNamespace(getGroupBy)), envir = callingEnv)
      group.byIsLiteral <- TRUE
  
      matchedCall[["group.by"]] <- group.by
    }
    
    treatmentCall <- switch(method.trt,
      glm       = redirectCall(matchedCall, quoteInNamespace(getGLMTreatmentFit)),
      bart      = redirectCall(matchedCall, quoteInNamespace(getBartTreatmentFit)),
      none      = NULL)
    
    if (!is.null(treatmentCall)) {
      if (!is.null(args.trt) && length(args.trt) > 0L)
        treatmentCall[names(matchedCall[["args.trt"]])[-1L]] <- matchedCall[["args.trt"]][-1L]
    
      if (!is.null(treatmentCall[["crossvalidate"]]))
        treatmentCall[["crossvalidate"]] <- if (is.logical(crossvalidate)) crossvalidate else crossvalidate == "trt"
      
      if (!is.na(seed) && method.trt == "bart")
        treatmentCall$seed <- sample.int(.Machine$integer.max, 1L)
      
      if (verbose) cat("fitting treatment model via method '", method.trt, "'\n", sep = "")
      
      massign[fit.trt, p.score, samples.p.score] <- eval(treatmentCall, envir = callingEnv)
    }
  } else {
    stop("method.trt must be in '", paste0(eval(formals(bartCause::bartc)$method.trt), collapse = "', '"), "' or a fixed vector")
  }
    
  if (!is.logical(p.scoreAsCovariate) || length(p.scoreAsCovariate) != 1L || is.na(p.scoreAsCovariate)) {
    ## a numeric vector here is almost always p.score = <vector> partial-matching
    ## this formal, since bartc() has no p.score of its own
    hint <- if (is.numeric(p.scoreAsCovariate))
      "; to supply propensity scores directly, use method.trt = <vector>" else ""
    stop("p.scoreAsCovariate must be a single TRUE or FALSE", hint)
  }
  if (method.rsp %in% c("p.weight", "tmle") && method.trt == "none")
    stop("response method '", method.rsp, "' requires propensity score estimation")
  if (method.rsp %in% c("bart", "bcf") && p.scoreAsCovariate == FALSE && method.trt != "none")
    warning("for response method '", method.rsp, "', propensity score not used unless included as covariate")
  if (!is.null(matchedCall$p.scoreAsCovariate) && p.scoreAsCovariate == TRUE && method.trt == "none")
    warning("p.scoreAsCovariate == TRUE requires method.trt != 'none'")
  
  
  if (!group.byIsLiteral && !is.null(matchedCall[["group.by"]]) && use.ranef) {
    group.by <- eval(redirectCall(matchedCall, quoteInNamespace(getGroupBy)), envir = callingEnv)
    group.byIsLiteral <- TRUE
    
    matchedCall$group.by <- group.by
  }
  
  responseCall <- switch(method.rsp,
    bcf      = redirectCall(matchedCall, quoteInNamespace(getBCFResponseFit)),
    bart     = redirectCall(matchedCall, quoteInNamespace(getBartResponseFit)),
    p.weight = redirectCall(matchedCall, quoteInNamespace(getPWeightResponseFit)),
    tmle     = redirectCall(matchedCall, quoteInNamespace(getTMLEResponseFit)))
  
  argsToAdd <- names(matchedCall) %in% names(formals(getBartResponseFit)) & names(matchedCall) %not_in% names(responseCall)
  if (any(argsToAdd)) for (argName in names(matchedCall)[argsToAdd])
    responseCall[[argName]] <- matchedCall[[argName]]
  
  if (!is.null(matchedCall$commonSup.rule)) {
    if (is.null(matchedCall$commonSup.cut))
      commonSup.cut <- eval(formals(bartCause::bartc)$commonSup.cut)[match(commonSup.rule, eval(formals(bartCause::bartc)$commonSup.rule))]
    responseCall$commonSup.rule <- commonSup.rule[1L]
    responseCall$commonSup.cut <- commonSup.cut[1L]
  } else {
    responseCall$commonSup.rule <- "none"
    responseCall$commonSup.cut  <- NA_real_
  }
  
  responseCall <- addCallDefaults(responseCall, bartCause::bartc)
  if ("verbose" %in% names(responseCall)) responseCall[[which(names(responseCall) == "verbose")]] <- verbose
  
  evalEnv <- callingEnv
  if ((p.scoreAsCovariate || method.rsp %in% c("tmle", "p.weight")) && !is.null(p.score)) {
    evalEnv <- new.env(parent = callingEnv)
    pScoreArgName <- "ps"
    if (!is.null(matchedCall$data))
      while (pScoreArgName %in% names(data)) pScoreArgName <- paste0(pScoreArgName, "ps")
    evalEnv[[pScoreArgName]] <- p.score
    
    responseCall$p.score <- as.symbol(pScoreArgName)
    
    if ("samples.p.score" %in% names(formals(eval(responseCall[[1L]])))) {
      evalEnv$samples.p.score <- samples.p.score
      responseCall$samples.p.score <- quote(samples.p.score)
    }
  }
  
  if (!is.null(args.rsp) && length(args.rsp) > 0L)
    responseCall[names(matchedCall[["args.rsp"]])[-1L]] <- matchedCall[["args.rsp"]][-1L]
  
  ## crossvalidation tunes the tree prior with xbart, which has no multi-forest
  ## form; refuse it here rather than let the argument reach the fitter as a
  ## silent no-op. 'trt' still crossvalidates the treatment model only.
  if (method.rsp == "bcf" && (isTRUE(crossvalidate) || identical(crossvalidate, "rsp")))
    stop("crossvalidate is not supported for response method 'bcf'; the crossvalidation engine has no multi-forest form, so the tree prior's 'k' cannot be tuned. Use method.rsp = 'bart' to crossvalidate the response model")

  if (!is.null(responseCall[["crossvalidate"]]))
    responseCall[["crossvalidate"]] <- if (is.logical(crossvalidate)) crossvalidate else crossvalidate == "rsp"
  
  if (!is.na(seed))
    responseCall$seed <- sample.int(.Machine$integer.max, 1L)
  
  if (verbose) cat("fitting response model via method '", method.rsp, "'\n", sep = "")
  
  fit <- data <- mu.hat.obs <- mu.hat.cf <- name.trt <- name.p.score <- trt <- sd.obs <-
    sd.cf <- commonSup.sub <- missingRows <- est <- fitPars <- NULL
  assignAll(eval(responseCall, envir = evalEnv))

  ## name.p.score is the design column the response fitter actually put the
  ## propensity score in, carried out the same way name.trt is: predict must
  ## not re-derive it, since no rule over the column names can tell the score
  ## from a confounder whose own name starts with the score's stem
  result <- namedList(fit.rsp = fit, data.rsp = data, fit.trt, mu.hat.obs, mu.hat.cf, p.score, samples.p.score,
                      method.rsp, method.trt, estimand,
                      commonSup.rule, commonSup.cut,
                      name.trt, name.p.score, trt,
                      sd.obs, sd.cf, commonSup.sub, missingRows, est, fitPars,
                      call = givenCall)
  if (!is.null(matchedCall[["group.by"]])) {
    result[["group.by"]] <-
      if (group.byIsLiteral)
        group.by
      else
        eval(redirectCall(matchedCall, quoteInNamespace(getGroupBy)), envir = callingEnv)
    if (use.ranef) result[["use.ranef"]] <- use.ranef
    if (group.effects) result[["group.effects"]] <- group.effects
  }
  ## dbarts >= 1.0-0 combines chains in fit$yhat.train (now 2-D), so derive the
  ## chain count from bartCause's own 3-D mu.hat.obs [n.chains, n.samples, n.obs]
  result$n.chains <- if (length(dim(mu.hat.obs)) > 2L) dim(mu.hat.obs)[1L] else 1L
  
  if (!is.na(seed)) {
    result$seed <- .GlobalEnv$.Random.seed
    
    if (!is.null(oldSeed))
      .GlobalEnv[[".Random.seed"]] <- oldSeed
    else
      rm(list = ".Random.seed", envir = .GlobalEnv) # unset the seed
  } else {
    if (!exists(".Random.seed", .GlobalEnv)) runif(1L)
    result$seed <- .GlobalEnv$.Random.seed
  }
   
  class(result) <- "bartcFit"
  result
}

responseIsBinary <- function(object) {
  if (inherits(object$fit.rsp, "stan4bartFit")) {
    object$fit.rsp$family$family != "gaussian"
  } else {
    is.null(object$fit.rsp[["sigma"]])
  }
}

