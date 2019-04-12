flattenSamples <- function(y) {
  x <- NULL ## R CMD check
  if (!is.null(dim(y)) && length(dim(y)) > 2L) evalx(dim(y), matrix(y, nrow = x[1L], ncol = x[2L] * x[3L])) else y
}

bartc <- function(
  response, treatment, confounders, data, subset, weights,
  method.rsp = c("bart", "tmle", "bcf", "p.weight"),
  method.trt = c("bart", "glm", "none", "bart.xval"),
  estimand   = c("ate", "att", "atc"),
  group.by = NULL,
  commonSup.rule = c("none", "sd", "chisq"),
  commonSup.cut  = c(NA_real_, 1, 0.05),
  args.rsp = list(), args.trt = list(),
  p.scoreAsCovariate = TRUE, use.rbart = FALSE,
  keepCall = TRUE, verbose = TRUE, ...
)
{
  matchedCall    <- match.call()
  sysCall        <- sys.call()
  callingEnv <- parent.frame(1L)
  
  # some dots arg can get eaten by R's argument matching algorithm, like 'k' for keepCall
  mismatchedArgs.sys <- names(sysCall) %not_in% names(matchedCall) & names(matchedCall) != "" &
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
  
  group.by <- eval(redirectCall(matchedCall, quoteInNamespace(getGroupBy)), envir = callingEnv)
  
  matchedCall$group.by <- group.by
  matchedCall$verbose <- NULL
  
  ## check validity of character vector arguments by comparing to function prototype
  for (argName in c("method.rsp", "estimand", "commonSup.rule")) {
    arg <- get(argName)
    if (!is.character(arg) || arg[1L] %not_in% eval(formals(bartCause::bartc)[[argName]]))
      stop(argName, " must be in '", paste0(eval(formals(bartCause::bartc)[[argName]]), collapse = "', '"), "'")
    assign(argName, arg[1L])
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
    
    
    treatmentCall <- switch(method.trt,
      glm       = redirectCall(matchedCall, quoteInNamespace(getGLMTreatmentFit)),
      bart      = redirectCall(matchedCall, quoteInNamespace(getBartTreatmentFit)),
      bart.xval = redirectCall(matchedCall, quoteInNamespace(getBartXValTreatmentFit)),
      none      = NULL)
    if (!is.null(args.trt) && length(args.trt) > 0L)
      treatmentCall[names(matchedCall[["args.trt"]])[-1L]] <- matchedCall[["args.trt"]][-1L]
    
    if (!is.null(treatmentCall)) {
      
      if (verbose) cat("fitting treatment model via method '", method.trt, "'\n", sep = "")
      
      massign[fit.trt, p.score, samples.p.score] <- eval(treatmentCall, envir = callingEnv)
    }
  } else {
    stop("method.trt must be in '", paste0(eval(formals(bartCause::bartc)$method.trt), collapse = "', '"), "' or a fixed vector")
  }
    
  if (is.na(p.scoreAsCovariate))
    stop("p.scoreAsCovariate must be TRUE or FALSE")
  if (method.rsp %in% c("p.weight", "tmle") && method.trt == "none")
    stop("response method '", method.rsp, "' requires propensity score estimation")
  if (method.rsp == "bart" && p.scoreAsCovariate == FALSE && method.trt != "none")
    warning("for response method 'bart', propensity score not used unless included as covariate")
  if (!is.null(matchedCall$p.scoreAsCovariate) && p.scoreAsCovariate == TRUE && method.trt == "none")
    warning("p.scoreAsCovariate == TRUE requires method.trt != 'none'")
  
  responseCall <- switch(method.rsp,
    bcf      = redirectCall(matchedCall, quoteInNamespace(bcf)),
    bart     = redirectCall(matchedCall, quoteInNamespace(getBartResponseFit)),
    p.weight = redirectCall(matchedCall, quoteInNamespace(getPWeightResponseFit)),
    tmle     = redirectCall(matchedCall, quoteInNamespace(getTMLEResponseFit)))
  
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
  
  evalEnv <- callingEnv
  if (p.scoreAsCovariate && !is.null(p.score)) {
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
  
  if (verbose) cat("fitting response model via method '", method.rsp, "'\n", sep = "")
  
  fit <- data <- mu.hat.obs <- mu.hat.cf <- name.trt <- trt <- sd.obs <-
    sd.cf <- commonSup.sub <- missingRows <- est <- fitPars <- NULL
  assignAll(eval(responseCall, envir = evalEnv))
  
  result <- namedList(fit.rsp = fit, data.rsp = data, fit.trt, mu.hat.obs, mu.hat.cf, p.score, samples.p.score,
                      method.rsp, method.trt, estimand, group.by,
                      commonSup.rule, commonSup.cut,
                      name.trt, trt,
                      sd.obs, sd.cf, commonSup.sub, missingRows, est, fitPars,
                      call = givenCall)
  result$n.chains <- if (length(dim(fit$yhat.train) > 2L)) dim(fit$yhat.train)[1L] else 1L
  
  if (!exists(".Random.seed", .GlobalEnv)) runif(1)
  result$seed <- .Random.seed
  
  class(result) <- "bartcFit"
  result
}

