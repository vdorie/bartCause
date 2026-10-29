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

## What group.by names, by the rule bartCause 1.0-10 used: its value where the
## caller can see it, else its value within data, which is how a bare column
## name or a function of columns is read. Whatever the expression needs that is
## not a column is found from the caller's frame. A grouping factor is one value
## per row of data, so the value is checked against it before any subset is
## applied.
getGroupBy <- function(data, subset, group.by)
{
  dataAreMissing <- missing(data)
  subsetIsMissing <- missing(subset)
  
  matchedCall <- match.call()
  
  if (is.null(matchedCall$group.by)) return(NULL)
  
  tryResult <- tryCatch(group.by.literal <- group.by, error = function(e) e)
  
  if (!dataAreMissing && inherits(tryResult, "error"))
    group.by <- eval(matchedCall$group.by, envir = data, enclos = parent.frame())
  
  if (!dataAreMissing && is.data.frame(data))
    checkGroupByLength(group.by, nrow(data), data)
  
  ## subset is read as the model frame reads it, within data first, so an
  ## expression of columns does not have to be found in the caller's frame
  if (!subsetIsMissing) {
    subset <- if (dataAreMissing) subset else eval(matchedCall$subset, envir = data, enclos = parent.frame())
    group.by <- group.by[subset]
  }
  
  as.factor(group.by)
}

## A grouping factor has a value for each of the 'n' rows. A single string is
## almost always the name of a column of 'data' given where the column itself
## was wanted, which would otherwise be recycled to one group.
checkGroupByLength <- function(group.by, n, data = NULL)
{
  if (length(group.by) == n) return(invisible(NULL))
  
  if (is.character(group.by) && length(group.by) == 1L && is.data.frame(data) && group.by %in% names(data))
    stop("'group.by' is the string \"", group.by, "\", which names a column of data; pass the column itself, ",
         "group.by = ", deparse(as.name(group.by), backtick = TRUE), ", not its name", call. = FALSE)
  
  stop("'group.by' must have one value for each row of data (", n, "), but its length is ", length(group.by), call. = FALSE)
}

## A modeled grouping factor (use.ranef = TRUE) rides the formula as
## (1 | <group.by as written>), and the mixed-model fitters look that grouping
## factor up in the model frame, where a function of columns, an element of the
## caller's list or a variable of the caller's is not to be found. As 1.0-10
## did, group.by is resolved to a vector first (getGroupBy); it is handed to
## the fit as a column of a copy of 'data' under a name nothing in the call
## uses, the copy bound under another such name in a child of 'env', and
## 'dataCall', the call to a data-call builder, is pointed at both. A group.by
## that is a bare column of data already carrying that value is left alone, so
## that the fit names its grouping factor as the user did. Returns the call and
## the environment to evaluate it in.
resolveGroupColumn <- function(matchedCall, dataCall, data, env, use.ranef)
{
  unchanged <- list(call = dataCall, env = env)
  groupExpr <- matchedCall[["group.by"]]
  if (is.null(groupExpr) || !isTRUE(use.ranef)) return(unchanged)
  
  ## the rows are those of data: a subset is the model frame's to apply
  groupCall <- redirectCall(matchedCall, quoteInNamespace(getGroupBy))
  if (!is.null(groupCall[["subset"]])) groupCall[["subset"]] <- NULL
  value <- eval(groupCall, envir = env)
  checkGroupByLength(value, nrow(data), data)
  
  if (is.symbol(groupExpr) && as.character(groupExpr) %in% names(data) &&
      identical(as.factor(data[[as.character(groupExpr)]]), value))
    return(unchanged)
  
  taken <- c(names(data), all.vars(matchedCall))
  columnName <- "g"
  while (columnName %in% taken) columnName <- paste0(columnName, "g")
  dataName <- "data"
  while (dataName %in% c(taken, columnName)) dataName <- paste0(dataName, "data")
  
  data[[columnName]] <- value
  dataEnv <- new.env(parent = env)
  dataEnv[[dataName]] <- data
  
  dataCall[["data"]] <- as.symbol(dataName)
  dataCall[["group.by"]] <- as.symbol(columnName)
  list(call = dataCall, env = dataEnv)
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
      checkGroupByLength(group.by, nrow(df))
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
      checkGroupByLength(group.by, nrow(df))
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
      checkGroupByLength(group.by, nrow(df))
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
      checkGroupByLength(group.by, nrow(df))
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

## dbarts and stan4bart name their sampling controls differently, and stan4bart
## takes no '...', so an argument left untranslated here is dropped by the call
## builders and the fit quietly runs stan4bart's own defaults - on the bartCause
## defaults, about five times the work that was asked for. 'args' is the
## fitter's matched call, read for the dbarts spellings; anything already
## written in stan4bart's own vocabulary came from the caller and wins.
addStan4BartSamplingArguments <- function(call, args, env, defaultChains = 10L, ownFormals = NULL,
                                             residualScale = TRUE)
{
  argNames <- names(args)
  
  ## the arguments that go on to be translated below or read by the caller are
  ## accounted for; any other name stan4bart has no formal for is dropped by the
  ## call builders, which 'ownFormals' (the caller's own arguments) lets this say.
  ## sigdf and sigquant are dbarts' prior on the residual scale, which the
  ## parametric component draws here, so they have nowhere to go.
  if (!is.null(ownFormals)) {
    accountedFor <- c("n.samples", "n.burn", "n.chains", "n.thin", "n.threads", "n.trees",
                      "k", "power", "base", "keepTrees", "combineChains",
                      names(formals(stan4bart::stan4bart)), ownFormals)
    givenNames <- argNames[argNames != ""]
    notUsed <- character()
    resid <- intersect(c("sigdf", "sigquant"), givenNames)
    ## a propensity model is binary, so only the response fit reports them
    if (residualScale)
      notUsed[resid] <- "the parametric component draws the residual standard deviation, so it has no prior to set"
    dropped <- setdiff(givenNames, c(accountedFor, resid))
    notUsed[dropped] <- "not an argument of stan4bart"
    warnArgumentsNotUsed(notUsed)
  }
  getArg <- function(name)
    if (name %in% argNames) eval(args[[name]], env) else NULL
  
  n.samples <- getArg("n.samples")
  n.burn    <- getArg("n.burn")
  n.chains  <- getArg("n.chains")
  
  if (is.null(call[["chains"]]))
    call[["chains"]] <- if (!is.null(n.chains)) as.integer(n.chains) else defaultChains
  
  ## stan4bart's 'iter' counts the warmup draws that dbarts' 'n.samples' excludes
  if (is.null(call[["iter"]])) {
    if (is.null(n.samples)) n.samples <- eval(formals(dbarts::bart)$n.samples)
    if (is.null(n.burn))    n.burn    <- eval(formals(dbarts::bart)$n.burn)
    
    ## dbarts' n.samples is a sweep budget of which n.samples %/% n.thin draws
    ## are kept, whereas every one of stan4bart's post-warmup iterations is
    ## kept, so the budget is divided here
    n.thin <- getArg("n.thin")
    n.kept <- as.integer(n.samples)
    ## an unusable n.thin is left for stan4bart's check of 'skip' to refuse
    if (is.null(call[["skip"]]) && is.numeric(n.thin) && length(n.thin) == 1L && !is.na(n.thin) && n.thin >= 1)
      n.kept <- n.kept %/% as.integer(n.thin)
    
    call[["iter"]] <- n.kept + as.integer(n.burn)
    if (is.null(call[["warmup"]])) call[["warmup"]] <- as.integer(n.burn)
  } else if (is.null(call[["warmup"]]) && !is.null(n.burn)) {
    call[["warmup"]] <- as.integer(n.burn)
  }
  
  ## 'skip' is dbarts' 'n.thin' for both blocks - that many transitions per kept
  ## draw
  n.thin <- getArg("n.thin")
  if (is.null(call[["skip"]]) && !is.null(n.thin))
    call[["skip"]] <- as.integer(n.thin)
  
  ## dbarts threads chains inside its own engine at no startup cost; stan4bart
  ## builds a cluster, about 1.25 seconds flat, and runs the chains serially
  ## unless told otherwise. It repays that only with more than one chain and a
  ## long enough fit, which is where the bartCause defaults sit.
  if (is.null(call[["cores"]])) {
    chains <- as.integer(eval(call[["chains"]], env))
    iter   <- as.integer(eval(call[["iter"]], env))
    n.threads <- getArg("n.threads")
    if (is.null(n.threads)) n.threads <- dbarts::guessNumCores()
    n.threads <- as.integer(n.threads)
    
    if (!is.na(n.threads) && n.threads > 1L && chains > 1L && chains * iter >= 4000L) {
      call[["cores"]] <- min(n.threads, chains)
      ## a clustered stan4bart fit seeds its chains off the clock unless it is
      ## given a seed, so hand it one drawn from R's generator and keep the fit
      ## reproducible under set.seed; the draws depend on that seed and the
      ## chain count, not on how many cores ran them
      if (is.null(call[["seed"]])) call[["seed"]] <- sample.int(.Machine$integer.max, 1L)
    }
  }
  
  ## the tree prior reaches stan4bart's forest through 'bart_args'; a bart_args
  ## that is not a literal list is the caller's own object and is left as given
  bartArgNames <- c("n.trees", "k", "power", "base", "keepTrees")
  bartArgNames <- bartArgNames[bartArgNames %in% argNames]
  if (length(bartArgNames) > 0L) {
    bartArgs <- call[["bart_args"]]
    if (is.null(bartArgs)) bartArgs <- quote(list())
    if (is.call(bartArgs) && bartArgs[[1L]] == quote(list)) {
      for (name in bartArgNames) {
        if (name %in% names(bartArgs)) next
        ## a prior written out in the call, k = chi(1.25, Inf), is stan4bart's to
        ## read with its own vocabulary, as it reads one written in bart_args;
        ## evaluated here it would find no chi() to call
        if (name == "k" && isPriorCall(args[[name]])) {
          bartArgs[[name]] <- args[[name]]
          next
        }
        value <- getArg(name)
        ## a call that is the value, as k = quote(chi(1.25, Inf)) gives, is put
        ## in bart_args as written and read as a prior there; a bare name would
        ## be looked up by stan4bart, which has no object of that name
        if (name == "k" && is.symbol(value))
          stop("'k' is a bare name, which is not taken as a prior; write k = chi(1.25, Inf) directly in the call, ",
               "or give dbartsPriors$chi(1.25, Inf) or the string \"chi(1.25, Inf)\" instead", call. = FALSE)
        bartArgs[[name]] <- value
      }
      call[["bart_args"]] <- bartArgs
    }
  }
  
  call
}

## A call to a function of dbarts' prior vocabulary, chi(1.25, Inf) or fixed(2)
## as written in the call, which only the fits know how to evaluate
isPriorCall <- function(expr)
  is.call(expr) && is.symbol(expr[[1L]]) && as.character(expr[[1L]]) %in% names(dbarts::dbartsPriors)

## bartCause 1.0-10 passed these names to dbarts, whose 0.9-x bart2 took them as
## arguments; dbarts 1.0-0 retired them, so a call that still names one is
## rewritten into the form that carries it: power, base and split.probs ride
## tree.prior = cgm(), sigdf and sigquant ride family = gaussian(sigma =
## chisq(df, quant)), rngSeed is seed. A name that has no 1.0-0 form, or that
## cannot apply, is returned in 'notUsed' (name -> reason) rather than dropped
## unseen. 'call' is a call or anything as.list() takes apart whose first
## element is the function. 'binary' is a response with no residual scale to
## give a prior to. Writing the successor as well as the retired name is
## refused, since the two would have to agree and nothing says which wins.
translateRetiredBartArguments <- function(call, binary)
{
  elements <- as.list(call)
  argNames <- names(elements)
  if (is.null(argNames)) argNames <- rep.int("", length(elements))
  notUsed <- character()
  
  refuseBoth <- function(retired, successor)
    stop(paste0("'", retired, "'", collapse = ", "), " and '", successor, "' cannot both be given; ",
         "dbarts 1.0-0 takes the one in '", successor, "'", call. = FALSE)
  
  treeNames <- intersect(c("power", "base", "split.probs"), argNames)
  if (length(treeNames) > 0L) {
    if ("tree.prior" %in% argNames) refuseBoth(treeNames, "tree.prior")
    treeCall <- as.call(c(list(quote(cgm)), elements[match(treeNames, argNames)]))
    elements[match(treeNames, argNames)] <- NULL
    elements[["tree.prior"]] <- treeCall
    argNames <- names(elements)
  }
  
  residNames <- intersect(c("sigdf", "sigquant"), argNames)
  if (length(residNames) > 0L) {
    if (binary) {
      notUsed[residNames] <- "a binary response has no residual scale to give a prior"
    } else {
      if ("family" %in% argNames) refuseBoth(residNames, "family")
      ## the defaults of 0.9-x's bart2
      sigdf    <- if ("sigdf" %in% residNames) elements[[match("sigdf", argNames)]] else 3
      sigquant <- if ("sigquant" %in% residNames) elements[[match("sigquant", argNames)]] else 0.9
      elements[match(residNames, argNames)] <- NULL
      elements[["family"]] <- call("gaussian", sigma = call("chisq", sigdf, sigquant))
    }
    if (binary) elements[match(residNames, argNames)] <- NULL
    argNames <- names(elements)
  }
  
  if ("rngSeed" %in% argNames) {
    if ("seed" %in% argNames) refuseBoth("rngSeed", "seed")
    names(elements)[match("rngSeed", argNames)] <- "seed"
    argNames <- names(elements)
  }
  
  noForm <- intersect(c("rngKind", "rngNormalKind"), argNames)
  if (length(noForm) > 0L) {
    notUsed[noForm] <- "dbarts 1.0-0 has no choice of random number generator"
    elements[match(noForm, argNames)] <- NULL
  }
  
  list(call = as.call(elements), notUsed = notUsed)
}

## One warning for a call, naming each argument that did not reach the sampler
## and why, in the order given. 'notUsed' is a character vector of reasons
## named by argument.
warnArgumentsNotUsed <- function(notUsed)
{
  if (length(notUsed) == 0L) return(invisible(NULL))
  
  reasons <- unique(notUsed)
  parts <- vapply(reasons, function(reason)
    paste0(paste0("'", names(notUsed)[notUsed == reason], "'", collapse = ", "), " (", reason, ")"), "")
  warning("arguments not used: ", paste0(parts, collapse = "; "), call. = FALSE)
}

## What bartc's '...' carries for the response method it was asked for -
## maxIter and posteriorOfTMLE for "tmle", a forest's moderators for "bcf" -
## travels to the BART fits as well, which have no use for it and are not to
## report it as dropped: "p.weight" and "tmle" fit their response through
## getBartResponseFit, and the propensity fit sees whatever the response does.
responseMethodArguments <- function()
  unique(unlist(lapply(list(getBartResponseFit, getStan4BartResponseFit, getPWeightResponseFit,
                            getTMLEResponseFit, getBCFResponseFit, bcf),
                       function(fn) names(formals(fn)))))

## args.rsp and args.trt hand arguments to a fit. Written as a list() call, each
## element is an expression that the fit evaluates as though it were written in
## its own call, as 1.0-10 took it. Any other expression - a variable, a
## function's result - is evaluated here to a list of values, which go into the
## call as they are, a symbol or call among them quoted so that the fit
## receives it rather than evaluating it. 'args' is that value, 'argsExpr' the
## expression it came from.
addArgsList <- function(call, argsExpr, args, argName)
{
  if (is.call(argsExpr) && identical(argsExpr[[1L]], quote(list))) {
    call[names(argsExpr)[-1L]] <- argsExpr[-1L]
    return(call)
  }
  
  if (!is.list(args))
    stop("'", argName, "' must be a list of arguments, each element named for the argument it gives", call. = FALSE)
  argNames <- names(args)
  if (is.null(argNames) || anyNA(argNames) || any(argNames == ""))
    stop("every element of '", argName, "' must be named for the argument it gives", call. = FALSE)
  
  for (i in seq_along(args)) {
    value <- args[i]
    if (is.language(value[[1L]])) value[[1L]] <- call("quote", value[[1L]])
    call[argNames[i]] <- value
  }
  call
}

## bartc(y, z, x1 + x2, d) binds d to 'parametric', the fourth argument, and
## fails on x1 with nothing to say that d was taken for something else. A data
## frame is a legitimate 'parametric' only for the vector form of the call,
## where the confounders are values too, so the call is refused when data is
## missing, the fourth argument evaluates to a data frame and the confounders
## are not values the caller can see. Each is evaluated in the caller's frame
## where an error or a warning means "not a value there", and nothing else is
## forced.
refuseDataGivenByPosition <- function(matchedCall, env)
{
  if (is.null(matchedCall[["parametric"]]) || !is.null(matchedCall[["data"]])) return(invisible(NULL))
  
  valueOf <- function(expr) tryCatch(eval(expr, env), error = function(e) NULL, warning = function(w) NULL)
  
  fourth <- valueOf(matchedCall[["parametric"]])
  if (!is.data.frame(fourth)) return(invisible(NULL))
  confounders <- valueOf(matchedCall[["confounders"]])
  if (is.numeric(confounders) || is.data.frame(confounders)) return(invisible(NULL))
  
  stop("the fourth argument of bartc is 'parametric', which was given a data frame; ",
       "'data' must be passed by name, data = ", paste0(deparse(matchedCall[["parametric"]]), collapse = " "), call. = FALSE)
}
