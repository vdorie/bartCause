## group.by rides in the formula rather than as a separate argument: a modeled
## intercept becomes an lmer-style (1 | g) term that stan4bart or lme4 reads,
## an unmodeled one just another right-hand-side term
addGroupTerm <- function(formula, group.by, use.ranef)
{
  if (is.null(group.by)) return(formula)
  
  term <- group.by
  if (use.ranef) {
    term <- quote((1 | g))
    term[[2L]][[3L]] <- group.by
  }
  formula[[3L]] <- call("+", formula[[3L]], term)
  
  formula
}

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
  
  groupIsPresent <- !is.null(matchedCall[["group.by"]])
  
  if (is.null(matchedCall[["parametric"]])) {
    if (!groupIsPresent || !use.ranef || use.lmer) {
      formula <- a ~ b
      formula[[2L]] <- matchedCall$treatment
      formula[[3L]] <- matchedCall$confounders
    } else {
      ## a modeled intercept with no parametric equation is fit by stan4bart,
      ## which takes the nonparametric part as a bart() term
      formula <- a ~ bart(b)
      formula[[2L]] <- matchedCall$treatment
      formula[[3L]][[2L]] <- matchedCall$confounders
    }
  } else {
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
  
  formula <- addGroupTerm(formula, if (groupIsPresent) matchedCall$group.by else NULL, use.ranef)
  
  environment(formula) <- parent.frame(1L)
  
  fn <- matchedCall$fn; matchedCall$fn <- NULL
  result <- redirectCall(matchedCall, fn)
  result <- addCallArgument(result, 1L, formula)
  
  list(call = result, env = parent.frame(1L))
}

getResponseDataCall <- function(fn, response, treatment, confounders, parametric, data, subset, weights, p.score, group.by, use.ranef)
{
  matchedCall <- match.call()
  groupIsPresent <- !is.null(matchedCall[["group.by"]])
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
      if (!groupIsPresent || !use.ranef) {
        formula <- a ~ b
        formula[[2L]] <- matchedCall$response
        formula[[3L]] <- quote(a + b)
        formula[[3L]][[2L]] <- matchedCall$confounders
        formula[[3L]][[3L]] <- matchedCall$treatment
      } else {
        formula <- a ~ bart(b + c)
        formula[[2L]] <- matchedCall$response
        formula[[3L]][[2L]][[2L]] <- matchedCall$confounders
        formula[[3L]][[2L]][[3L]] <- matchedCall$treatment
      }
    } else {
      formula <- response ~ treatment + bart(confounders + treatment) + parametric
      # ~(response, RHS)
      formula[[2L]] <- matchedCall$response
      # formula[[3L]] := +(treatment + bart(confounders), parametric)
      formula[[3L]][[3L]] <- matchedCall$parametric
      formula[[3L]][[2L]][[3L]][[2L]][[3L]] <- matchedCall$treatment
      formula[[3L]][[2L]][[3L]][[2L]][[2L]] <- matchedCall$confounders
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

      ## p.scoreEval is fit on the subset only (its length is length(subset), not
      ## nrow(data)); assigning it straight into the full-length column errors when
      ## the lengths don't divide and silently recycles/misaligns when they do
      ## (subset can refer to 'data', so resolve it the same way bcf.R does).
      ## Place it at the subset positions instead, mirroring getResponseLiteralCall.
      if (!is.null(matchedCall$subset)) {
        subsetValue <- eval(matchedCall$subset, data, evalEnv)
        data[[pScoreName]] <- numeric(nrow(data))
        data[[pScoreName]][subsetValue] <- p.scoreEval
      } else {
        data[[pScoreName]] <- p.scoreEval
      }

      evalEnv <- new.env(parent = parent.frame(1L))
      evalEnv[["data"]] <- data
      
      matchedCall$data <- quote(data) # going to redirect to a different data object
    } else {
      pScoreName <- deparse(matchedCall$p.score)
    }
    
    if (is.null(matchedCall[["parametric"]])) {
      if (!groupIsPresent || !use.ranef) {
        formula <- a ~ b
        formula[[2L]] <- matchedCall$response
        formula[[3L]] <- quote(a + b)
        formula[[3L]][[2L]] <- quote(a + b)
        formula[[3L]][[2L]][[2L]] <- matchedCall$confounders
        formula[[3L]][[2L]][[3L]] <- str2lang(pScoreName)
        formula[[3L]][[3L]] <- matchedCall$treatment
      } else {
        formula <- a ~ bart(b + c + d)
        formula[[2L]] <- matchedCall$response
        formula[[3L]][[2L]][[2L]][[2L]] <- matchedCall$confounders
        formula[[3L]][[2L]][[2L]][[3L]] <- str2lang(pScoreName)
        formula[[3L]][[2L]][[3L]]       <- matchedCall$treatment
      }
    } else {
      if (exists("pScoreName")) {
        formula <- response ~ treatment + p.score + bart(confounders + treatment + p.score) + parametric
        formula[[2L]] <- matchedCall$response
        # formula[[3L]] is all of RHS
        # modify parse tree from end of RHS back, since the tails of binary ops are scalars
        formula[[3L]][[3L]] <- matchedCall$parametric
        
        # formula[[3L]][[2L]][[3L]] - all of bart(); formula[[3L]][[2L]][[3L]][[2L]] - what's inside
        formula[[3L]][[2L]][[3L]][[2L]][[3L]] <- str2lang(pScoreName)
        formula[[3L]][[2L]][[3L]][[2L]][[2L]][[3L]] <- matchedCall$treatment
        formula[[3L]][[2L]][[3L]][[2L]][[2L]][[2L]] <- matchedCall$confounders

        # linear model part
        formula[[3L]][[2L]][[2L]][[3L]] <- str2lang(pScoreName)
        formula[[3L]][[2L]][[2L]][[2L]] <- matchedCall$treatment
      } else {
        formula <- response ~ treatment + bart(confounders + treatment) + parametric
        formula[[2L]] <- matchedCall$response
        # formula[[3L]] is all of RHS
        # modify parse tree from end of RHS back, since the tails of binary ops are scalars
        formula[[3L]][[3L]] <- matchedCall$parametric
        
        # formula[[3L]][[2L]][[3L]] - all of bart(); formula[[3L]][[2L]][[3L]][[2L]] - what's inside
        formula[[3L]][[2L]][[3L]][[2L]][[3L]] <- matchedCall$treatment
        formula[[3L]][[2L]][[3L]][[2L]][[2L]] <- matchedCall$confounders

        # linear model part
        formula[[3L]][[2L]][[2L]] <- matchedCall$treatment
      }
    }
  }
  
  formula <- addGroupTerm(formula, if (groupIsPresent) matchedCall$group.by else NULL, use.ranef)
  
  environment(formula) <- evalEnv
  
  fn <- matchedCall$fn; matchedCall$fn <- NULL
  result <- redirectCall(matchedCall, fn)
  result <- addCallArgument(result, 1L, formula)
  
  #responseVar <- as.vector(evalEnv[[deparse(result$data)]][[result[[2L]][[2L]]]])
  responseVar <- as.vector(get(deparse(result$data), envir = evalEnv)[[result[[2L]][[2L]]]])
  ## the propensity score's column name is RETURNED rather than re-derived by a
  ## fitter: with a confounder named "psps" no ladder over colnames(data@x) can
  ## tell the score from the confounder, and the bcf fitter has to know which
  ## column to keep out of the treatment forest. Appended, so the four-target
  ## massign calls that consume this list positionally never see it.
  list(call = result, env = evalEnv, trt = deparse(matchedCall$treatment), missingRows = is.na(responseVar),
       p.score = if (is.null(matchedCall$p.score)) NULL else pScoreName)
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
    }
    
    if (is.null(group.by) || !use.ranef || use.lmer) {
      formula <- a ~ b
      formula[[2L]] <- str2lang(treatmentName)
      formula[[3L]] <- str2lang(paste0(confounderNames, collapse = " + "))
    } else {   
      formula <- a ~ bart(b)
      formula[[2L]] <- str2lang(treatmentName)
      formula[[3L]][[2L]] <- str2lang(paste0(confounderNames, collapse = " + "))
    }
  } else {
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
    
    if (!is.null(group.by)) {
      group.byName <- "g"
      while (group.byName %in% colnames(df))
        group.byName <- paste0(group.byName, "g")
      df[[group.byName]] <- group.by
    }
    
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
  
  formula <- addGroupTerm(formula, if (!is.null(group.by)) str2lang(group.byName) else NULL, use.ranef)
    
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
    
    if (!is.null(group.by)) {
      group.byName <- "g"
      while (group.byName %in% colnames(df))
        group.byName <- paste0(group.byName, "g")
      df[[group.byName]] <- group.by
    }
    
    modelNames <- setdiff(colnames(df), if (is.null(group.by)) responseName else c(responseName, group.byName))
    if (is.null(group.by) || !use.ranef) {
      formula <- a ~ b
      formula[[2L]] <- str2lang(responseName)
      formula[[3L]] <- str2lang(paste0(modelNames, collapse = " + "))
    } else {
      formula <- a ~ bart(b)
      formula[[2L]] <- str2lang(responseName)
      formula[[3L]][[2L]] <- str2lang(paste0(modelNames, collapse = " + "))
    }
  } else {
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
    
    if (!is.null(group.by)) {
      group.byName <- "g"
      while (group.byName %in% colnames(df))
        group.byName <- paste0(group.byName, "g")
      df[[group.byName]] <- group.by
    }
    
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

      allParametricNames <- c(treatmentName, pScoreName, parametricNames)
      nonParametricNames <- c(confounderNames, treatmentName, pScoreName)
    } else {
      allParametricNames <- c(treatmentName, parametricNames)
      nonParametricNames <- c(confounderNames, treatmentName)
    }


    formula <- response ~ parametrics + bart(nonParametrics)
    formula[[2L]] <- str2lang(responseName)
    formula[[3L]][[2L]] <- str2lang(paste0(allParametricNames, collapse = " + "))
    formula[[3L]][[3L]][[2L]] <- str2lang(paste0(nonParametricNames, collapse = " + "))
  }
  
  formula <- addGroupTerm(formula, if (!is.null(group.by)) str2lang(group.byName) else NULL, use.ranef)
  
  result <- quote(functionName(formula, data = df))
  result[[1L]] <- matchedCall$fn
  result[[2L]] <- formula
  
  if (!is.null(matchedCall$subset))  result$subset <- subset
  if (!is.null(matchedCall$weights)) result$weights <- weights

  ## as in getResponseDataCall: the resolved propensity-score column name is
  ## appended for the bcf fitter's moderator exclusion, and positional consumers
  ## of the first four elements are unaffected
  list(call = result, df = df, trt = treatmentName, missingRows = is.na(as.vector(df[,responseName])),
       p.score = if (is.null(matchedCall$p.score)) NULL else pScoreName)
}

