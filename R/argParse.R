getGroupBy <- function(data, subset, group.by)
{
  dataAreMissing <- missing(data)
  subsetIsMissing <- missing(subset)
  
  matchedCall <- match.call()
  
  if (is.null(matchedCall$group.by)) return(NULL)
  
  tryResult <- tryCatch(group.by.literal <- group.by, error = function(e) e)
  
  if (!dataAreMissing && is(tryResult, "error"))
    group.by <- eval(matchedCall$group.by, envir = data)
  
  if (!subsetIsMissing) group.by <- group.by[subset]
  
  as.factor(group.by)
}

## set up call to look inside 'data'
getTreatmentDataCall <- function(fn, treatment, confounders, data, subset, weights, group.by, use.ranef, use.lmer)
{
  matchedCall <- match.call()
  tryResult <- tryCatch(confounders.literal <- confounders, error = function(e) e)
  if (!is(tryResult, "error")) {
    if (is.language(confounders.literal))
      matchedCall$confounders <- confounders.literal
    else if (is.character(confounders.literal))
      matchedCall$confounders <- str2lang(confounders.literal)
  }
  
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
  environment(formula) <- parent.frame(1L)
  
  fn <- matchedCall$fn; matchedCall$fn <- NULL
  result <- redirectCall(matchedCall, fn)
  result <- addCallArgument(result, 1L, formula)
  
  list(call = result, env = parent.frame(1L))
}

getResponseDataCall <- function(fn, response, treatment, confounders, data, subset, weights, p.score, group.by, use.ranef)
{
  matchedCall <- match.call()
  tryResult <- tryCatch(confounders.literal <- confounders, error = function(e) e)
  if (!is(tryResult, "error")) {
    if (is.language(confounders.literal))
      matchedCall$confounders <- confounders.literal
    else if (is.character(confounders.literal))
      matchedCall$confounders <- str2lang(confounders.literal)
  }
  
  if (is.null(matchedCall$p.score)) {
    evalEnv <- parent.frame(1L)
    
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
    ## if the p.score is present it was likely estimated (or just given) and thus not
    ## present in 'data' or data's environment
    
    evalEnv <- parent.frame(1L)
    ## check to see if p.score is in the calling environment
    p.scoreEval <- tryCatch(p.score, error = function(e) e)
    if (!is(p.scoreEval, "error")) {
      ## add it to data, copy data into a new environment
      pScoreName <- "ps"
      while (pScoreName %in% names(data))
        pScoreName <- paste0(pScoreName, "ps")
      
      evalEnv <- new.env(parent = parent.frame(1L))
      data[[pScoreName]] <- p.scoreEval
      
      evalEnv[["data"]] <- data
      
      matchedCall$data <- quote(data) ## going to redirect to a different data object
    } else {
      pScoreName <- deparse(matchedCall$p.score)
    }
    
    if (is.null(matchedCall[["group.by"]]) || use.ranef) {
      formula <- a ~ b
      formula[[2L]] <- matchedCall$response
      formula[[3L]] <- quote(a + b)
      formula[[3L]][[2L]] <- quote(a + b)
      formula[[3L]][[2L]][[2L]] <- matchedCall$confounders
      formula[[3L]][[2L]][[3L]] <- parse(text = pScoreName)[[1L]]
      formula[[3L]][[3L]] <- matchedCall$treatment
    } else {
      formula <- a ~ b + c + d + e
      formula[[2L]] <- matchedCall$response
      formula[[3L]][[2L]][[2L]][[2L]] <- matchedCall$confounders
      formula[[3L]][[2L]][[2L]][[3L]] <- parse(text = pScoreName)[[1L]]
      formula[[3L]][[2L]][[3L]]       <- matchedCall$treatment
      formula[[3L]][[3L]]             <- matchedCall$group.by
    }
  }
  
  environment(formula) <- evalEnv
  
  fn <- matchedCall$fn; matchedCall$fn <- NULL
  result <- redirectCall(matchedCall, fn)
  result <- addCallArgument(result, 1L, formula)
  
  responseVar <- evalEnv[[deparse(result$data)]][[result[[2L]][[2L]]]]
  list(call = result, env = evalEnv, trt = deparse(matchedCall$treatment), missingRows = is.na(responseVar))
}

## treat args as literals
getTreatmentLiteralCall <- function(fn, treatment, confounders, subset, weights, group.by, use.ranef, use.lmer)
{
  matchedCall <- match.call()
  if (is.null(matchedCall[["group.by"]])) group.by <- NULL
  
  x <- NULL ## R CMD check
  
  treatmentName <- "z"
  while (treatmentName %in% colnames(confounders))
    treatmentName <- paste0(treatmentName, "z")
  df <- as.data.frame(cbind(confounders, treatment))
  names(df)[ncol(df)] <- treatmentName
  
  formulaTerms <- evalx(colnames(df), x[x != treatmentName])
  if (!is.null(group.by)) {
    group.byName <- "g"
    while (group.byName %in% colnames(df))
      group.byName <- paste0(group.byName, "g")
    df[[group.byName]] <- group.by
    
    if (!use.ranef) formulaTerms <- c(formulaTerms, group.byName)
  }
  
  if (is.null(group.by) || !use.ranef || !use.lmer) {
    formula <- a ~ b
    formula[[2L]] <- parse(text = treatmentName)[[1L]]
    formula[[3L]] <- parse(text = paste0(formulaTerms, collapse = " + "))[[1L]]
  } else {   
    formula <- a ~ b + (1 | c)
    formula[[2L]] <- parse(text = treatmentName)[[1L]]
    formula[[3L]][[2L]] <- parse(text = paste0(evalx(colnames(df), x[x %not_in% c(treatmentName, group.byName)]), collapse = " + "))[[1L]]
    formula[[3L]][[3L]][[2L]][[3L]] <- parse(text = group.byName)[[1L]]
  }
    
  ## ls is temp
  result <- quote(ls(formula, data = df))
  result[[1L]] <- matchedCall$fn
  result[[2L]] <- formula
  if (!is.null(matchedCall$subset)) result$subset <- subset
  if (!is.null(matchedCall$weights)) result$weights <- weights
  
  list(call = result, df = df)
}

getResponseLiteralCall <- function(fn, response, treatment, confounders, subset, weights, p.score, group.by, use.ranef)
{
  matchedCall <- match.call()
  if (is.null(matchedCall[["group.by"]])) group.by <- NULL
  
  x <- NULL ## R CMD check
  
  df <- as.data.frame(cbind(confounders, response, treatment))
  responseName <- "y"
  while (responseName %in% colnames(df))
    responseName <- paste0(responseName, "y")
  treatmentName <- "z"
  while (treatmentName %in% colnames(df))
    treatmentName <- paste0(treatmentName, "z")
  names(df)[length(df) - 1L] <- responseName
  names(df)[length(df)]      <- treatmentName
  
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
  
  formulaTerms <- evalx(colnames(df), x[x != responseName])
  if (!is.null(group.by) && !use.ranef) {
    group.byName <- "g"
    while (group.byName %in% colnames(df))
      group.byName <- paste0(group.byName, "g")
    df[[group.byName]] <- group.by
    
    formulaTerms <- c(formulaTerms, group.byName)
  }
  
  formula <- a ~ b
  formula[[2L]] <- parse(text = responseName)[[1L]]
  formula[[3L]] <- parse(text = paste0(formulaTerms, collapse = " + "))[[1L]]
  
  ## ls is temp
  result <- quote(ls(formula, data = df))
  result[[1L]] <- matchedCall$fn
  result[[2L]] <- formula
  if (!is.null(matchedCall$subset))  result$subset <- subset
  if (!is.null(matchedCall$weights)) result$weights <- weights
  
  list(call = result, df = df, trt = treatmentName, missingRows = is.na(df[[responseName]]))
}

