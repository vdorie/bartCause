getGroupBy <- function(confounders, data, subset, group.by)
{
  dataAreMissing <- missing(data)
  subsetIsMissing <- missing(subset)
  
  matchedCall <- match.call()
  
  if (is.null(matchedCall$group.by)) return(NULL)
  
  if (!dataAreMissing) {
    group.by <- data[[matchedCall$group.by]]
  }
  
  if (!subsetIsMissing) group.by <- group.by[subset]
  
  as.factor(group.by)
}

## set up call to look inside 'data'
getTreatmentDataCall <- function(fn, treatment, confounders, data, subset, weights)
{
  matchedCall <- match.call()
  
  formula <- a ~ b
  formula[[2L]] <- matchedCall$treatment
  formula[[3L]] <- matchedCall$confounders
  environment(formula) <- parent.frame(1L)
  
  fn <- matchedCall$fn; matchedCall$fn <- NULL
  result <- redirectCall(matchedCall, fn)
  result <- addCallArgument(result, 1L, formula)
  
  list(call = result, env = parent.frame(1L))
}

getResponseDataCall <- function(fn, response, treatment, confounders, data, subset, weights, p.score)
{
  matchedCall <- match.call()
  
  if (is.null(matchedCall$p.score)) {
    evalEnv <- parent.frame(1L)
    
    formula <- a ~ b
    formula[[2L]] <- matchedCall$response
    formula[[3L]] <- quote(a + b)
    formula[[3L]][[2L]] <- matchedCall$confounders
    formula[[3L]][[3L]] <- matchedCall$treatment
  } else {
    ## if the p.score is preset it was likely estimated (or just given) and thus not
    ## preset in 'data' or data's environment
    
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
    
    formula[[3L]][[2L]] <- quote(a + b)
    formula[[3L]][[2L]][[2L]] <- matchedCall$confounders
    formula[[3L]][[2L]][[3L]] <- parse(text = pScoreName)[[1L]]
    formula[[3L]][[3L]] <- matchedCall$treatment
  }
  
  environment(formula) <- evalEnv
  
  fn <- matchedCall$fn; matchedCall$fn <- NULL
  result <- redirectCall(matchedCall, fn)
  result <- addCallArgument(result, 1L, formula)
  
  list(call = result, env = evalEnv, trt = deparse(matchedCall$treatment))
}

## treat args as literals
getTreatmentLiteralCall <- function(fn, treatment, confounders, subset, weights)
{
  matchedCall <- match.call()
  
  treatmentName <- "z"
  while (treatmentName %in% colnames(confounders))
    treatmentName <- paste0(treatmentName, "z")
  df <- as.data.frame(cbind(confounders, treatment))
  names(df)[ncol(df)] <- treatmentName
  
  formula <- a ~ b
  formula[[2L]] <- parse(text = treatmentName)[[1L]]
  formula[[3L]] <- parse(text = paste0(evalx(colnames(df), x[x != treatmentName]), collapse = " + "))[[1L]]
  
  ## ls is temp
  result <- quote(ls(formula, data = df))
  result[[1L]] <- matchedCall$fn
  result[[2L]] <- formula
  if (!is.null(matchedCall$subset)) result$subset <- subset
  if (!is.null(matchedCall$weights)) result$weights <- weights
  
  list(call = result, df = df)
}

getResponseLiteralCall <- function(fn, response, treatment, confounders, subset, weights, p.score)
{
  matchedCall <- match.call()
  
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
    
    df[[pScoreName]] <- p.score
  }
  
  formula <- a ~ b
  formula[[2L]] <- parse(text = responseName)[[1L]]
  formula[[3L]] <- parse(text = paste0(evalx(colnames(df), x[x != responseName]), collapse = " + "))[[1L]]
  
  ## ls is temp
  result <- quote(ls(formula, data = df))
  result[[1L]] <- matchedCall$fn
  result[[2L]] <- formula
  if (!is.null(matchedCall$subset))  result$subset <- subset
  if (!is.null(matchedCall$weights)) result$weights <- weights
  
  list(call = result, df = df, trt = treatmentName)
}

