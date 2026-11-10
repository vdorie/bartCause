bartc <- function(
  response, treatment, confounders, parametric, data, subset, weights, weights.rsp,
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
                        names(sysCall) %in% names(formals(dbarts::bart))
  if (any(mismatchedArgs.sys)) {
    ## the formal that took each name is the one before '...' it is a prefix of;
    ## the arguments are not told apart by their values, which two can share
    newNames <- names(sysCall)[mismatchedArgs.sys]
    ownFormals <- names(formals(bartc))
    ownFormals <- ownFormals[seq_len(match("...", ownFormals) - 1L)]
    ownFormals <- ownFormals[ownFormals %in% names(matchedCall)]
    oldNames <- vapply(newNames, function(name) ownFormals[startsWith(ownFormals, name)][1L], "")
    
    names(matchedCall)[match(oldNames[!is.na(oldNames)], names(matchedCall))] <- newNames[!is.na(oldNames)]
    ## the value stays an unevaluated expression in matchedCall for the fits to
    ## read, which is where the prior vocabulary (chi, normal, ...) is defined;
    ## reading the formal here would force it in the caller's frame
    for (oldName in oldNames[!is.na(oldNames)])
      assign(oldName, eval(formals(bartc)[[oldName]]))
  }
  
  ## the fourth argument is 'parametric', so data given by position lands
  ## there, and the formulas it was meant for then fail on their variables
  refuseDataGivenByPosition(matchedCall, callingEnv)
  if ("tree.prior" %in% names(matchedCall)) refuseRetiredTreePrior(matchedCall[["tree.prior"]])
  
  givenCall <- if (keepCall) matchedCall else call("NULL")
  matchedCall$verbose <- NULL  
  
  ## check validity of character vector arguments by comparing to function prototype
  for (argName in c("method.rsp", "estimand", "commonSup.rule")) {
    arg <- get(argName)
    if (!is.character(arg) || arg[1L] %not_in% eval(formals(bartCause::bartc)[[argName]]))
      stop(argName, " must be in '", paste0(eval(formals(bartCause::bartc)[[argName]]), collapse = "', '"), "'")
    assign(argName, arg[1L])
  }
  
  ## weights with tmle need the package: say so before fitting anything
  if (method.rsp == "tmle" && !is.null(matchedCall[["weights"]])) getTMLEFunction(TRUE)
  
  if (length(crossvalidate) != 1L ||
      (!is.logical(crossvalidate) && !is.character(crossvalidate)) ||
      (is.logical(crossvalidate) && is.na(crossvalidate)) ||
      (is.character(crossvalidate) && crossvalidate %not_in% c("rsp", "trt")))
    stop("crossvalidate must be one of TRUE, FALSE, 'rsp', or 'trt'")
  
  # the fits below get their own seeds, and dbarts reads seed = NA as a
  # retired spelling, so the forwarded call never carries this one
  if (!is.null(matchedCall[["seed"]])) matchedCall[["seed"]] <- NULL
  if (!is.na(seed)) {
    oldSeed <- .GlobalEnv[[".Random.seed"]]
    set.seed(seed)
  }
  
  ## an argument taken out of a fit's call because another governs is named
  ## once, with the reason, in a warning before the response model is fit
  overridden <- character()
  takeOverridden <- function(call) {
    overridden <<- c(overridden, attr(call, "overridden"))
    attr(call, "overridden") <- NULL
    call
  }
  ## each fit is given a seed drawn from this one, which governs where an
  ## rngSeed is also given, as in 1.0-10; the fits would refuse the pair
  dropRngSeed <- function(call) {
    if (is.na(seed) || is.null(call[["rngSeed"]])) return(call)
    overridden <<- c(overridden, c(rngSeed = "seed gives each fit its own seed, which governs"))
    call[["rngSeed"]] <- NULL
    call
  }
  
  ## the population weights and the response fit's own weights are resolved
  ## once, over the rows modeled; the fits read the response fit's weights
  ## from its expression as written, and the averages read the population's
  popWeights <- resolveRowWeights(matchedCall, "weights", callingEnv)
  rspWeights <- resolveRowWeights(matchedCall, "weights.rsp", callingEnv)
  ## NULL weights are absent, and the fits' calls do not carry them
  if (is.null(popWeights) && "weights" %in% names(matchedCall)) matchedCall[["weights"]] <- NULL
  if (is.null(rspWeights) && "weights.rsp" %in% names(matchedCall)) matchedCall[["weights.rsp"]] <- NULL
  ## an estimand that the weights and the treatment leave no weight is refused
  ## before anything is fit; a common support rule that leaves none, after
  if (!is.null(popWeights)) {
    z <- treatmentBeforeFit(matchedCall, callingEnv, length(popWeights$full), popWeights$rows)
    if (!is.null(z)) stopIfNoWeight(popWeights$value, z, estimand)
  }
  argsRspNames <- if (hasArgsList(matchedCall[["args.rsp"]], args.rsp)) argsListNames(matchedCall[["args.rsp"]], args.rsp) else character()
  familyGiven <- !is.null(matchedCall[["family"]]) || "family" %in% argsRspNames
  if (!is.null(rspWeights) && length(rspWeights$value) > 0L && !all(rspWeights$value %in% c(0, 1)) &&
      !familyGiven && responseIsZeroOne(matchedCall, callingEnv))
    stop("a 0/1 response is fit by probit BART, which takes weights.rsp of 0 or 1 only, having no weighted ",
         "form; counts of identical rows are written out as repeated rows, or fit with family = \"logistic\" ",
         "and weights.type = \"frequency\"; population weights go in 'weights', which weight the average of ",
         "the unit effects, with the variables they were built from among the confounders", call. = FALSE)
  
  ## a row out of a continuous fit is refused before any fit where the response
  ## plainly is one: without a family, a response fit is binary exactly when its
  ## response is 0/1 valued, and a family that says gaussian is continuous; any
  ## other family is read after the fit
  if (!is.null(rspWeights) &&
      (if (!familyGiven) !responseIsZeroOne(matchedCall, callingEnv)
       else familyIsPlainlyGaussian(givenResponseFamily(matchedCall, args.rsp), callingEnv)))
    refuseRowsOutOfResponseFit(rspWeights, popWeights)
  
  ## the response fit's weights come by weights.rsp alone, and their reading
  ## says how it reads them; both refused before anything is fit
  if ("weights" %in% argsRspNames)
    stop("args.rsp gives 'weights': the response fit's own weights are bartc's weights.rsp, and the ",
         "population weights that every average reads are 'weights'", call. = FALSE)
  if (is.null(rspWeights) && (!is.null(matchedCall[["weights.type"]]) || "weights.type" %in% argsRspNames))
    stop("weights.type says how the response fit reads weights.rsp, and without weights.rsp the response ",
         "fit is unweighted; drop weights.type, or give weights.rsp", call. = FALSE)
  
  ## the grouping factor is resolved once, so the treatment fit, the response
  ## fit and the group-level estimands all see the same levels; the call keeps
  ## the user's expression, which is what the model formulas name
  group.by <- if (!is.null(matchedCall[["group.by"]]))
    eval(redirectCall(matchedCall, quoteInNamespace(getGroupBy)), envir = callingEnv)
  else
    NULL
  
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
    
    treatmentCall <- switch(method.trt,
      glm       = redirectCall(matchedCall, quoteInNamespace(getGLMTreatmentFit)),
      bart      = redirectCall(matchedCall, quoteInNamespace(getBartTreatmentFit)),
      none      = NULL)
    
    trtEnv <- callingEnv
    if (!is.null(treatmentCall)) {
      ## weights.type says how the response fit reads its weights; the
      ## propensity fit is not given them
      if (!is.null(treatmentCall[["weights.type"]])) treatmentCall[["weights.type"]] <- NULL
      if (method.trt == "glm" && !is.null(popWeights)) {
        ## glm's prior weights, at mean 1 over the rows modeled: its starting
        ## values are taken from the weights, and at a survey's scale it
        ## diverges; the estimate does not depend on their scale
        trtEnv <- new.env(parent = callingEnv)
        weightsName <- uniqueName("w", matchedCall, callingEnv)
        trtEnv[[weightsName]] <- popWeights$full / mean(popWeights$value)
        treatmentCall[["weights"]] <- as.symbol(weightsName)
      } else if (!is.null(treatmentCall[["weights"]])) {
        ## a BART propensity model is a probit, which has no weighted form
        treatmentCall[["weights"]] <- NULL
        if (method.rsp %in% c("p.weight", "tmle"))
          message("the propensity score model is fit unweighted by BART, so its scores are the sample's; ",
                  "if sampling depended on treatment, method.trt = \"glm\" fits it with the weights")
      }
      if (hasArgsList(matchedCall[["args.trt"]], args.trt))
        treatmentCall <- takeOverridden(addArgsList(treatmentCall, matchedCall[["args.trt"]], args.trt, "args.trt"))
      treatmentCall <- dropRngSeed(treatmentCall)
    
      if (!is.null(treatmentCall[["crossvalidate"]]))
        treatmentCall[["crossvalidate"]] <- if (is.logical(crossvalidate)) crossvalidate else crossvalidate == "trt"
      
      if (!is.na(seed) && method.trt == "bart")
        treatmentCall$seed <- sample.int(.Machine$integer.max, 1L)
      
      if (verbose) cat("fitting treatment model via method '", method.trt, "'\n", sep = "")
      
      massign[fit.trt, p.score, samples.p.score] <- eval(treatmentCall, envir = trtEnv)
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
  
  ## the response fit takes weights.rsp, and is unweighted without it
  if (!is.null(matchedCall[["weights.rsp"]]))
    responseCall[["weights"]] <- matchedCall[["weights.rsp"]]
  else if (!is.null(responseCall[["weights"]]))
    responseCall[["weights"]] <- NULL
  
  evalEnv <- new.env(parent = callingEnv)
  ## the propensity-weighted estimators average with the population weights,
  ## handed over on the rows modeled
  if (method.rsp %in% c("tmle", "p.weight") && !is.null(popWeights)) {
    evalEnv$pop.weights <- popWeights$value
    responseCall$pop.weights <- quote(pop.weights)
  }
  if ((p.scoreAsCovariate || method.rsp %in% c("tmle", "p.weight")) && !is.null(p.score)) {
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
  
  if (hasArgsList(matchedCall[["args.rsp"]], args.rsp))
    responseCall <- takeOverridden(addArgsList(responseCall, matchedCall[["args.rsp"]], args.rsp, "args.rsp"))
  responseCall <- dropRngSeed(responseCall)
  ## each argument is named once for each reason, though both fits dropped it
  warnArgumentsNotUsed(overridden[!duplicated(paste(names(overridden), overridden))])
  
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
  ## both weight vectors are kept over every row an effect is estimated for,
  ## the missing-response rows included, NULL where not given
  ## the posterior draws are kept once, in fit.rsp and fit.trt: mu.hat.obs,
  ## mu.hat.cf and samples.p.score are read out of them when asked for (see
  ## R/draws.R). The exception is a matrix of propensity draws that was given,
  ## which no fit holds. combineChains records how the fits' draws were laid out.
  if (!is.null(popWeights))
    stopIfNoWeight(popWeights$value[commonSup.sub], trt[commonSup.sub], estimand)
  combineChains <- combineChainsRecord(mu.hat.obs, samples.p.score)
  result <- namedList(fit.rsp = fit, data.rsp = data, fit.trt, p.score,
                      method.rsp, method.trt, estimand,
                      commonSup.rule, commonSup.cut,
                      name.trt, name.p.score, trt,
                      sd.obs, sd.cf, commonSup.sub, missingRows, est, fitPars,
                      weights = if (!is.null(popWeights)) popWeights$value else NULL,
                      weights.rsp = if (!is.null(rspWeights)) rspWeights$value else NULL,
                      combineChains,
                      call = givenCall)
  if (method.trt == "given" && !is.null(samples.p.score))
    result[["samples.p.score"]] <- samples.p.score
  
  ## a row out of a continuous response fit has no noise scale for its
  ## counterfactual, so it cannot be in an average; with a family given, the
  ## fit says whether it is binary
  if (!is.null(rspWeights) && !responseIsBinary(result))
    refuseRowsOutOfResponseFit(rspWeights, popWeights)
  if (!is.null(matchedCall[["group.by"]])) {
    result[["group.by"]] <- group.by
    if (use.ranef) result[["use.ranef"]] <- use.ranef
    if (group.effects) result[["group.effects"]] <- group.effects
  }
  ## combineChains = TRUE collapses mu.hat.obs to 2-D, so its own dims cannot
  ## be trusted for the chain count; every response fitter's own object
  ## (dbarts's bart/bartBCF, both of which record n.chains directly) survives
  ## that collapse, so read it from there and fall back to the old dim-based
  ## guess only when a fit object does not carry one
  result$n.chains <- if (!is.null(fit[["n.chains"]])) fit[["n.chains"]]
                      else if (length(dim(mu.hat.obs)) > 2L) dim(mu.hat.obs)[1L] else 1L
  
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
    ## a held sigma lives only in the fit's $fixed, so it is not binary
    is.null(object$fit.rsp[["sigma"]]) && is.null(heldSigma(object))
  }
}

## The residual standard deviation a dbarts response fit holds fixed
## (family = gaussian(sigma = fixed(v)), or aft), a single number, else NULL.
heldSigma <- function(object) {
  fit <- object$fit.rsp
  if (inherits(fit, "stan4bartFit") || is.null(fit)) return(NULL)
  fit[["fixed"]][["sigma"]]
}


## A row with weights.rsp 0 is out of a continuous response fit, so it has no
## noise scale for its counterfactual and cannot be in an average its
## population weight is positive in.
refuseRowsOutOfResponseFit <- function(rspWeights, popWeights)
{
  outOfFit <- rspWeights$value <= 0 & (if (is.null(popWeights)) TRUE else popWeights$value > 0)
  if (!any(outOfFit)) return(invisible(NULL))
  bad <- rspWeights$rows[outOfFit]
  stop(if (length(bad) > 1L) "rows " else "row ", paste0(head(bad, 5L), collapse = ", "),
       if (length(bad) > 5L) paste0(" and ", length(bad) - 5L, " more") else "",
       if (length(bad) > 1L) " have" else " has",
       " weights.rsp 0 but a positive population weight: a row out of the response fit has no ",
       "noise scale for its counterfactual; give it a positive weights.rsp, or a population weight of 0; a row ",
       "that is to be in the treatment model only is given an NA response and a population weight of 0",
       call. = FALSE)
}
