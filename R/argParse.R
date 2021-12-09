getGroupBy <- function(data, subset, group.by)
{
  dataAreMissing <- missing(data)
  subsetIsMissing <- missing(subset)
  
  matchedCall <- match.call()
  
  if (is.null(matchedCall$group.by)) return(NULL)
  
  tryResult <- tryCatch(group.by.literal <- group.by, error = function(e) e)
  
  if (!dataAreMissing && inherits(tryResult, "error"))
    group.by <- eval(matchedCall$group.by, envir = data)
  
  if (!subsetIsMissing) group.by <- group.by[subset]
  
  as.factor(group.by)
}

# set up call to look inside 'data'
getTreatmentDataCall <- function(fn, treatment, confounders, parametric, data, subset, weights, group.by, use.ranef, use.lmer)
{
  matchedCall <- match.call()
  tryResult <- tryCatch(confounders.literal <- confounders, error = function(e) e, warning = function(w) w)
  if (!inherits(tryResult, "error") && !inherits(tryResult, "warning")) {
    if (is.language(confounders.literal))
      matchedCall$confounders <- confounders.literal
    else if (is.character(confounders.literal))
      matchedCall$confounders <- str2lang(confounders.literal)
  }
  if (!is.null(matchedCall[["parametric"]])) {
    tryResult <- tryCatch(parametric.literal <- parametric, error = function(e) e, warning = function(w) w)
    if (!inherits(tryResult, "error") && !inherits(tryResult, "warning")) {
      if (is.language(parametric.literal))
        matchedCall$parametric <- parametric.literal
      else if (is.character(parametric.literal))
        matchedCall$parametric <- str2lang(parametric.literal)
    }
  }
  
  if (is.null(matchedCall[["parametric"]])) {
    if (is.null(matchedCall[["group.by"]]) || (use.ranef && !use.lmer)) {
      formula <- a ~ b
      formula[[2L]] <- matchedCall$treatment
      formula[[3L]] <- matchedCall$confounders
    } else if (!use.ranef) {
      # add as fixed effects
      formula <- a ~ b + c
      formula[[2L]] <- matchedCall$treatment
      formula[[3L]][[2L]] <- matchedCall$confounders
      formula[[3L]][[3L]] <- matchedCall$group.by
    } else {
      formula <- a ~ b + (1 | c)
      formula[[2L]] <- matchedCall$treatment
      formula[[3L]][[2L]] <- matchedCall$confounders
      formula[[3L]][[3L]][[2L]][[3L]] <- matchedCall$group.by
    }
  } else {
    if (!is.null(matchedCall[["group.by"]]))
      stop("`group.by` must be missing or NULL if `parametric` is supplied; for varying intercepts, add (1 | group) to parametric equation")
    if (!use.lmer) {
      formula <- treatment ~ parametrics + bart(nonParametrics)
      formula[[2L]] <- matchedCall$treatment
      formula[[3L]][[2L]] <- matchedCall$parametric
      formula[[3L]][[3L]][[2L]] <- matchedCall$confounders
    } else {
      formula <- treatment ~ parametrics + nonParametrics
      formula[[2L]] <- matchedCall$treatment
      formula[[3L]][[2L]] <- matchedCall$parametric
      formula[[3L]][[3L]] <- matchedCall$confounders
    }
  } 
  
  environment(formula) <- parent.frame(1L)
  
  fn <- matchedCall$fn; matchedCall$fn <- NULL
  result <- redirectCall(matchedCall, fn)
  result <- addCallArgument(result, 1L, formula)
  
  list(call = result, env = parent.frame(1L))
}

getResponseDataCall <- function(fn, response, treatment, confounders, parametric, data, subset, weights, p.score, group.by, use.ranef)
{
  matchedCall <- match.call()
  tryResult <- tryCatch(confounders.literal <- confounders, error = function(e) e, warning = function(w) w)
  if (!inherits(tryResult, "error") && !inherits(tryResult, "warning")) {
    if (is.language(confounders.literal))
      matchedCall$confounders <- confounders.literal
    else if (is.character(confounders.literal))
      matchedCall$confounders <- str2lang(confounders.literal)
  }
  if (!is.null(matchedCall[["parametric"]])) {
    tryResult <- tryCatch(parametric.literal <- parametric, error = function(e) e, warning = function(w) w)
    if (!inherits(tryResult, "error") && !inherits(tryResult, "warning")) {
      if (is.language(parametric.literal))
        matchedCall$parametric <- parametric.literal
      else if (is.character(parametric.literal))
        matchedCall$parametric <- str2lang(parametric.literal)
    }
  }
  
  if (is.null(matchedCall$p.score)) {
    evalEnv <- parent.frame(1L)
    
    if (is.null(matchedCall[["parametric"]])) {
      if (is.null(matchedCall[["group.by"]]) || use.ranef) {
        formula <- a ~ b
        formula[[2L]] <- matchedCall$response
        formula[[3L]] <- quote(a + b)
        formula[[3L]][[2L]] <- matchedCall$confounders
        formula[[3L]][[3L]] <- matchedCall$treatment
      } else {
        formula <- a ~ b + c + d
        formula[[2L]] <- matchedCall$response
        formula[[3L]][[2L]][[2L]] <- matchedCall$confounders
        formula[[3L]][[2L]][[3L]] <- matchedCall$treatment
        formula[[3L]][[3L]] <- matchedCall$group.by
      }
    } else {
      if (!is.null(matchedCall[["group.by"]]))
        stop("group.by must be missing or NULL if parametric is supplied; for varying intercepts, add (1 | group) to parametric equation")
      formula <- response ~ treatment + bart(confounders) + parametric
      # ~(response, RHS)
      formula[[2L]] <- matchedCall$response
      # formula[[3L]] := +(treatment + bart(confounders), parametric)
      formula[[3L]][[3L]] <- matchedCall$parametric
      formula[[3L]][[2L]][[3L]][[2L]] <- matchedCall$confounders
      formula[[3L]][[2L]][[2L]] <- matchedCall$treatment
    }
  } else {
    # if the p.score is present it was likely estimated (or just given) and thus not
    # present in 'data' or data's environment
    
    evalEnv <- parent.frame(1L)
    # check to see if p.score is in the calling environment
    p.scoreEval <- tryCatch(p.score, error = function(e) e)
    if (!inherits(p.scoreEval, "error")) {
      # add it to data, copy data into a new environment
      pScoreName <- "ps"
      while (pScoreName %in% names(data))
        pScoreName <- paste0(pScoreName, "ps")
      
      evalEnv <- new.env(parent = parent.frame(1L))
      data[[pScoreName]] <- p.scoreEval
      
      evalEnv[["data"]] <- data
      
      matchedCall$data <- quote(data) # going to redirect to a different data object
    } else {
      pScoreName <- deparse(matchedCall$p.score)
    }
    
    if (is.null(matchedCall[["parametric"]])) {
      if (is.null(matchedCall[["group.by"]]) || use.ranef) {
        formula <- a ~ b
        formula[[2L]] <- matchedCall$response
        formula[[3L]] <- quote(a + b)
        formula[[3L]][[2L]] <- quote(a + b)
        formula[[3L]][[2L]][[2L]] <- matchedCall$confounders
        formula[[3L]][[2L]][[3L]] <- str2lang(pScoreName)
        formula[[3L]][[3L]] <- matchedCall$treatment
      } else {
        formula <- a ~ b + c + d + e
        formula[[2L]] <- matchedCall$response
        formula[[3L]][[2L]][[2L]][[2L]] <- matchedCall$confounders
        formula[[3L]][[2L]][[2L]][[3L]] <- str2lang(pScoreName)
        formula[[3L]][[2L]][[3L]]       <- matchedCall$treatment
        formula[[3L]][[3L]]             <- matchedCall$group.by
      }
    } else {
      if (!is.null(matchedCall[["group.by"]]))
        stop("group.by must be missing or NULL if parametric is supplied; for varying intercepts, add (1 | group) to parametric equation")
      
      formula <- response ~ treatment + p.score + bart(confounders) + parametric
      # ~(response, RHS)
      formula[[2L]] <- matchedCall$response
      formula[[3L]][[3L]] <- matchedCall$parametric
      formula[[3L]][[2L]][[3L]][[2L]] <- matchedCall$confounders
      formula[[3L]][[2L]][[2L]][[3L]] <- str2lang(pScoreName)
      formula[[3L]][[2L]][[2L]][[2L]] <- matchedCall$treatment
    }
  }
  
  environment(formula) <- evalEnv
  
  fn <- matchedCall$fn; matchedCall$fn <- NULL
  result <- redirectCall(matchedCall, fn)
  result <- addCallArgument(result, 1L, formula)
  
  #responseVar <- as.vector(evalEnv[[deparse(result$data)]][[result[[2L]][[2L]]]])
  responseVar <- as.vector(get(deparse(result$data), envir = evalEnv)[[result[[2L]][[2L]]]])
  if (length(responseVar) == 0L) browser()
  list(call = result, env = evalEnv, trt = deparse(matchedCall$treatment), missingRows = is.na(responseVar))
}

# treat args as literals
getTreatmentLiteralCall <- function(fn, treatment, confounders, parametric, subset, weights, group.by, use.ranef, use.lmer)
{
  matchedCall <- match.call()
  if (is.null(matchedCall[["group.by"]])) group.by <- NULL
  
  x <- NULL # R CMD check
  treatmentName <- "z"
  
  if (is.null(matchedCall[["parametric"]])) {
    confounderNames <- colnames(confounders)
    
    if (is.null(confounderNames))
      confounderNames <- paste0("V", seq_len(NCOL(confounders)))
    
    while (treatmentName %in% confounderNames)
      treatmentName <- paste0(treatmentName, "z")
    
    df <- as.data.frame(cbind(treatment, confounders))
    colnames(df) <- c(treatmentName, confounderNames)
    
    if (!is.null(group.by)) {
      group.byName <- "g"
      while (group.byName %in% colnames(df))
        group.byName <- paste0(group.byName, "g")
      df[[group.byName]] <- group.by
      
      if (!use.ranef) confounderNames <- c(confounderNames, group.byName)
    }
    
    if (is.null(group.by) || !use.ranef || !use.lmer) {
      formula <- a ~ b
      formula[[2L]] <- str2lang(treatmentName)
      formula[[3L]] <- str2lang(paste0(confounderNames, collapse = " + "))
    } else {   
      formula <- a ~ b + (1 | c)
      formula[[2L]] <- str2lang(treatmentName)
      formula[[3L]][[2L]] <- str2lang(paste0(setdiff(colnames(df), c(treatmentName, group.byName)), collapse = " + "))
      formula[[3L]][[3L]][[2L]][[3L]] <- str2lang(group.byName)
    }
  } else {
    if (!is.null(group.by))
      stop("group.by must be missing or NULL if parametric is supplied; for varying intercepts, add (1 | group) to parametric equation")
    
    confounderNames <- colnames(confounders)
    parametricNames <- colnames(parametric)
    
    if (is.null(confounderNames))
      confounderNames <- paste0("V", seq_len(NCOL(confounders)), "_bart")
    if (is.null(parametricNames))
      parametricNames <- paste0("V", seq_len(NCOL(parametric)))
    
    nameCollidedConfoundersExpression <- evalx(confounderNames, quote(x[x %in% parametricNames]))
    evalx(nameCollidedConfoundersExpression, x <- paste0(x, "_bart"), forceX = TRUE)
    
    while (treatmentName %in% confounderNames || treatmentName %in% parametricNames)
      treatmentName <- paste0(treatmentName, "z")
     
    df <- as.data.frame(cbind(treatment, confounders, parametric))
    colnames(df) <- c(treatmentName, confounderNames, parametricNames)
    
    if (!use.lmer) {
      formula <- treatment ~ parametrics + bart(nonParametrics)
      formula[[2L]] <- str2lang(treatmentName)
      formula[[3L]][[2L]] <- str2lang(paste0(parametricNames, collapse = " + "))
      formula[[3L]][[3L]][[2L]] <- str2lang(paste0(confounderNames, collapse = " + "))
    } else {
      formula <- treatment ~ allTerms
      formula[[2L]] <- str2lang(treatmentName)
      formula[[3L]] <- str2lang(paste0(c(parametricNames, confounderNames), collapse = " + "))
    }
  }
    
  result <- quote(functionName(formula, data = df))
  result[[1L]] <- matchedCall$fn
  result[[2L]] <- formula
  
  if (!is.null(matchedCall$subset)) result$subset <- subset
  if (!is.null(matchedCall$weights)) result$weights <- weights
  
  list(call = result, df = df)
}

getResponseLiteralCall <- function(fn, response, treatment, confounders, parametric, subset, weights, p.score, group.by, use.ranef)
{
  matchedCall <- match.call()
  if (is.null(matchedCall[["group.by"]])) group.by <- NULL
  
  x <- NULL # R CMD check
  
  responseName <- "y"
  treatmentName <- "z"
  
  if (is.null(matchedCall[["parametric"]])) {
    confounderNames <- colnames(confounders)
    
    if (is.null(confounderNames))
      confounderNames <- paste0("V", seq_len(NCOL(confounders)))
    
    while (responseName %in% confounderNames)
      responseName <- paste0(responseName, "y")
    while (treatmentName %in% confounderNames)
      treatmentName <- paste0(treatmentName, "z")
    
    df <- as.data.frame(cbind(response, treatment, confounders))
    colnames(df) <- c(responseName, treatmentName, confounderNames)
    
    if (!is.null(matchedCall$p.score)) {
      pScoreName <- "ps"
      while (pScoreName %in% names(df))
        pScoreName <- paste0(pScoreName, "ps")
      
      if (!is.null(matchedCall$subset)) {
        df[[pScoreName]] <- numeric(nrow(df))
        df[[pScoreName]][subset] <- p.score
      } else {
        df[[pScoreName]] <- p.score
      }
    }
    
    if (!is.null(group.by) && !use.ranef) {
      group.byName <- "g"
      while (group.byName %in% colnames(df))
        group.byName <- paste0(group.byName, "g")
      df[[group.byName]] <- group.by
    }
    
    formula <- a ~ b
    formula[[2L]] <- str2lang(responseName)
    formula[[3L]] <- str2lang(paste0(setdiff(colnames(df), responseName), collapse = " + "))
  } else {
    if (!is.null(group.by))
      stop("group.by must be missing or NULL if parametric is supplied; for varying intercepts, add (1 | group) to parametric equation")
    
    confounderNames <- colnames(confounders)
    parametricNames <- colnames(parametric)
    
    if (is.null(confounderNames))
      confounderNames <- paste0("V", seq_len(ncol(confounders)), "_bart")
    if (is.null(parametricNames))
      parametricNames <- paste0("V", seq_len(ncol(parametric)))
    
    nameCollidedConfoundersExpression <- evalx(confounderNames, quote(x[x %in% parametricNames]))
    evalx(nameCollidedConfoundersExpression, x <- paste0(x, "_bart"), forceX = TRUE)
    
    while (responseName %in% confounderNames || responseName %in% parametricNames)
      responseName <- paste0(responseName, "y")
    
    while (treatmentName %in% confounderNames || treatmentName %in% parametricNames)
      treatmentName <- paste0(treatmentName, "z")
    
    df <- as.data.frame(cbind(response, treatment, confounders, parametric))
    colnames(df) <- c(responseName, treatmentName, confounderNames, parametricNames)
    
    if (!is.null(matchedCall$p.score)) {
      pScoreName <- "ps"
      while (pScoreName %in% colnames(df))
        pScoreName <- paste0(pScoreName, "ps")
      
      if (!is.null(matchedCall$subset)) {
        df <- cbind(df, numeric())
        colnames(df)[ncol(df)] <- pScoreName
        df[subset,pScoreName] <- p.score
      } else {
        df <- cbind(df, p.score)
        colnames(df)[ncol(df)] <- pScoreName
      }
    }
    
    formula <- response ~ parametrics + bart(nonParametrics)
    formula[[2L]] <- str2lang(responseName)
    formula[[3L]][[2L]] <- str2lang(paste0(setdiff(colnames(df), c(responseName, confounderNames)), collapse = " + "))
    formula[[3L]][[3L]][[2L]] <- str2lang(paste0(confounderNames, collapse = " + "))
  }
  
  result <- quote(functionName(formula, data = df))
  result[[1L]] <- matchedCall$fn
  result[[2L]] <- formula
  
  if (!is.null(matchedCall$subset))  result$subset <- subset
  if (!is.null(matchedCall$weights)) result$weights <- weights
   
  list(call = result, df = df, trt = treatmentName, missingRows = is.na(as.vector(df[,responseName])))
}

