## allows multiple assignment
massign <- structure(NA, class = "lval")

"[<-.lval" <- function(x, ..., value) {
  callingEnv <- parent.frame(1L)
  
  args <- as.list(match.call())
  args <- args[-c(1L, 2L, length(args))]
  argNames <- names(args)
  
  valueNames <- names(value)
  # length(value) <- length(args)
  if (length(value) < length(args)) {
    value <- rep_len(value, length(args))
    if (!is.null(valueNames)) names(value) <- rep_len(valueNames, length(args))
  }
  
  varNames <- as.character(args)
  if (any(argNames != ""))
    varNames[argNames != ""] <- argNames[argNames != ""]
  
  varNamesNoSkip <- varNames[varNames != ""]
  
  duplicateVarNames <- duplicated(varNamesNoSkip)
  if (any(duplicateVarNames))
    warning("names on left-hand-side of assignment appear more than once: ",
            paste0(varNamesNoSkip[duplicateVarNames], collapse = ", "),
            "; result undefined", sep = "")
  
  
  ## for unnamed rhs, we go entirely by position
  if (is.null(valueNames)) {
    if (any(argNames != ""))
      warning("right-hand-side of assignment is unnamed; using position only")
    
    for (i in seq_along(varNames)) {
      var <- args[[i]]
      
      if (!missing(var)) {
        #cat("assigning '", varNames[i], " value (", value[[1L]], ") in env '", format(callingEnv), "'\n")
        assign(varNames[i], value[[1L]], envir = callingEnv)
        #eval.parent(substitute(var <- val, list(var = varNames[i], val = value[[1L]])))
      }
      
      value <- value[-1L]
    }
    return(x)
  }
  
  
  ## go through named args first
  for (i in seq_along(varNames)) {
    if (argNames[i] == "") next
    
    varName <- varNames[i]
    valueName <- as.character(args[[i]])
    
    sel <- valueNames == valueName
    numMatches <- sum(sel)
    
    if (numMatches == 0L)
      stop("'", valueName, "' not present in right-hand-side of assignment", sep = "")
    
    if (numMatches > 1L) {
      warning("'", valueName, "' present multiple times in right-hand-side of assignment; only first will be used", sep = "")
      selectionIndex <- which.max(sel)
      sel <- logical(length(value))
      sel[selectionIndex] <- TRUE
    }
    
    
    #cat("assigning '", varName, " value (", value[sel][[1L]], ") in env '", format(callingEnv), "'\n")
    assign(varName, value[sel][[1L]], envir = callingEnv)
    #eval.parent(substitute(varName <- val, list(varName = varName, val = value[sel][[1L]])))
    ## check to see if the value is named later, if not pop it off
    if (i == length(argNames) || !any(valueName %in% args[seq.int(i + 1L, length(args))])) {
      ## tail(args, -i))) {
      value <- value[!sel]
      valueNames <- valueNames[!sel]
    }
  }
  
  ## now for unnamed args
  for (i in seq_along(args)) {
    if (argNames[i] != "") next
    
    var <- args[[i]]
    
    if (!missing(var)) {
      #cat("assigning '", as.character(var), " value (", value[[1L]], ") in env '", format(callingEnv), "'\n")
      assign(as.character(var), value[[1L]], envir = callingEnv)
      #eval.parent(substitute(var <- val, list(var = var, val = value[[1L]])))
    }
    
    ## pop selected values 
    value <- value[-1L]
    valueNames <- valueNames[-1L]
  }
  x
}

unpack <- structure(NA, class = "named_lval")

"[<-.named_lval" <- function(x, ..., value) {
  mc <- match.call()
  
  mc[[1L]] <- as.symbol("[<-.lval")
  sel <- seq.int(3L, length(mc) - 1L)
  varNames <- as.character(mc[sel])
  names(mc)[sel] <- as.character(mc[sel])
  
  if (!all(varNames %in% names(value)))
    mc <- mc[-sel[!(varNames %in% names(value))]]
  
  #cat("calling ", as.character(mc), "\n")
  eval(mc, parent.frame(1L))
}
