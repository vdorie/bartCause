coerceOrError <- function(x, type)
{
  mc <- match.call()
  
  if (is.null(x)) stop("'", mc[[2L]], "' cannot be NULL")
  
  func <- switch(type, logical = as.logical, integer = as.integer, numeric = as.numeric)
  result <- tryCatch(func(x), warning = function(e) e)
  if (inherits(result, "warning")) stop("'", mc[[2L]], "' must be coercible to type: ", type)
  
  result
}

"%not_in%" <- function(x, table) match(x, table, nomatch = 0L) <= 0L

evalx.recurse <- function(x, e) {
  if (length(e) == 0L || typeof(e) == "symbol") return(e)
  
  for (i in seq_along(e)) {
    if (!is.language(e[[i]])) next
    
    e[[i]] <- if (e[[i]] == "x") {
        # prematurely evaluate recursive calls to evalx
        if (is.language(x) && !is.symbol(x) && length(x) > 0L && x[[1L]] == quote(evalx)) eval(x, parent.frame()) else x
      } else {
        evalx.recurse(x, e[[i]])
      }
  }
  
  e
}

# evaluates the expression 'e' by after first replacing all instances of 'x' with the expression x
# if force is TRUE, x will be evaluated and not passed on as an expression. This is useful when
# x itself evaluates to the expession desired
#
# for example:
# complicated_x <- quote(big_variable_name[big_variable_name > 0])
# evalx(complicated_x, x <- x + 5, forceX = TRUE)

# resolves to:
# big_variable_name[big_variable_name > 0] <-
#    big_variable_name[big_variable_name > 0] + 5
#
# in contrast to:
evalx <- function(x, e, forceX = FALSE) {
  mc <- match.call()
  callingEnv <- parent.frame()
  
  if (!forceX) {
    e <- evalx.recurse(mc$x, mc$e)
    eval(e, callingEnv)
  } else {
    x <- force(x)
    e <- evalx.recurse(x, mc$e)
    eval(e, callingEnv)
  }
}

"%as_x_in%" <- function(a, expr) {
  mc <- match.call()
  calling_env <- parent.frame()
  
  expr <- evalx.recurse(mc$a, mc$expr)
  eval(expr, calling_env)
}


ifelse_3 <- function(a, b, c, d, e) {
  mc <- match.call(); env <- parent.frame()
  if (eval(mc[["a"]], env)) {
    c
  } else if (eval(mc[["b"]], env)) {
    d
  } else {
    e
  }
}

ifelse_n <- function(n, ...) {
  mc <- match.call(); env <- parent.frame()
  
  for (i in seq_len(n - 1L))
    if (eval(mc[[i + 2L]], env)) return(eval(mc[[n + 1L + i]], env))
  eval(mc[[2L * n - 1L]], env)
}

redirectCall <- function(call, fn, ...)
{
  matchedCall <- match.call()
  extraArgs <- if (length(matchedCall) > 3L) as.character(matchedCall[-c(1L, 2L, 3L)]) else character()
  
  originalFn <- eval(call[[1L]])
  call[[1L]] <- if (is.function(fn)) matchedCall[[3L]] else fn
  if (length(extraArgs) == 0L) {
    fn <- if (is.function(fn)) fn else eval(fn)
    
    argsToKeep <- names(call)[-1L] %in% names(formals(fn))
    if (any(names(formals(originalFn)) == "...") && any(names(formals(fn)) == "..."))
      argsToKeep <- argsToKeep | names(call)[-1L] %not_in% names(formals(originalFn))
    
    call <- call[c(TRUE, argsToKeep)]
  } else {
    matchIndices <- match(extraArgs, names(call), nomatch = 0L)
    
    call <- call[c(1L, matchIndices)]
  }
  
  call
}

addCallDefaults <- function(call, fn)
{
  possibleArgs <- names(formals(eval(call[[1L]])))
  evalx(possibleArgs, if (any(x == "...")) x <- x[x != "..."])
  
  currentArgs <- names(call)[-1L]
  fnFormals   <- formals(fn)
  
  ## prune down to just those with defaults
  fnFormals <- fnFormals[sapply(fnFormals, function(x) !is.symbol(x))]
  
  formalsToAdd <- names(fnFormals) %in% possibleArgs & names(fnFormals) %not_in% currentArgs
  if (any(formalsToAdd)) {
    fnFormals <- fnFormals[formalsToAdd]
    for (i in seq_along(fnFormals)) {
      if (!is.null(fnFormals[[i]])) call[[names(fnFormals)[i]]] <- fnFormals[[i]]
    }
  }
  
  call
}

addCallArgument <- function(call, position, argument)
{
  if (is.character(position)) {
    name <- position
    position <- length(call) + 1L
  } else {
    position <- as.integer(position) + 1L
    if (position <= length(call)) for (i in seq.int(length(call), position)) {
      call[[i + 1L]] <- call[[i]]
      names(call)[i + 1L] <- names(call)[i]
    }
    name <- ""
  }
  call[[position]] <- argument
  names(call)[position] <- name
  call
}

addCallArguments <- function(call, args, replace = TRUE)
{
  fnFormals <- formals(eval(call[[1L]]))
  
  oldArgs <- names(args) != "" & names(args) %in% names(call)
  if (replace == TRUE & any(oldArgs))
    call[names(call) != "" & names(call) %in% names(args)] <- args[oldArgs]
  args <- args[!oldArgs]
  
  if (inherits(args, "call"))
    args <- args[-1L]

  for (i in seq_along(args)) {
    if (!is.null(names(args)) && names(args)[i] != "" && names(args)[i] %not_in% names(fnFormals))
      next
    
    call[[length(call) + 1L]] <- args[[i]]
    if (!is.null(names(args)) && names(args)[i] != "")
      names(call)[length(call)] <- names(args)[i]
  }
  call
}

pruneCallArguments <- function(call, ignoreDots = FALSE)
{
  fnFormals <- formals(eval(call[[1L]]))
  if (!ignoreDots && any(names(fnFormals) == "...")) return(call)
  
  for (i in seq.int(length(call), 2L)) {
    if (names(call)[i] == "") next
    if (names(call)[i] %not_in% names(fnFormals)) {
      if (i < length(call)) for (j in seq.int(i, length(call))) {
        call[[j]] <- call[[j + 1]]
        names(call)[j] <- names(call)[j + 1]
      }
      names(call)[length(call)] <- ""
      call[[length(call)]] <- NULL
    }
  }
  
  call
}

subTermInLanguage <- function(lang, oldTerm, newTerm)
{
  if (length(lang) == 1L && is.symbol(lang))
    return(if (lang == oldTerm) newTerm else lang)

  for (i in seq_along(lang)) {
    if (is.symbol(lang[[i]])) {
      if (lang[[i]] == oldTerm) lang[[i]] <- newTerm
    } else if (is.language(lang[[i]])) {
      lang[[i]] <- subTermInLanguage(lang[[i]], oldTerm, newTerm)
    }
  }
  lang
}

setDefaultsFromFormals <- function(call, formals, ...)
{
  argsToReplace <- list(...)
  matchIndices <- match(argsToReplace, names(call), nomatch = 0L)
  missingFormals <- match(argsToReplace[matchIndices == 0L], names(formals))

  if (length(missingFormals) == 0L) return(call)
  
  call[seq.int(length(missingFormals)) + length(call)] <- formals[missingFormals]
  call
}

is.formula <- function(x) is.language(x) && x[[1L]] == '~'

## from lme4
namedList <- function(...) {
  result <- list(...)
  substituteNames <- sapply(substitute(list(...)), deparse)[-1L]
  if (is.null(resultNames <- names(result))) resultNames <- substituteNames
  if (any(noNames <- resultNames == "")) resultNames[noNames] <- substituteNames[noNames]
  setNames(result, resultNames)
}

## use this to produce calls of the form
##  dbarts:::functionName
## so that we can evaluate non-exported functions in
## the user's environment
quoteInNamespace <- function(name, character.only = FALSE) {
  result <- quote(a + b)
  result[[1L]] <- as.symbol(":::")
  result[[2L]] <- as.symbol("bartCause")
  
  result[[3L]] <- if (character.only) name else match.call()[[2]]
  result
}

## silly function to handle subsetting when there are (possibly) multiple
## chains - goes through the parse tree and adds the correct number of commas
addDimsToSubset <- function(e) {
  subDims <- function(e, env) {
    if (is.call(e) && e[[1L]] == "[") {
      temp <- quote(dim(a))
      temp[[2L]] <- e[[2L]]
     
      dims <- eval(temp, env)
      if (is.null(dims)) return(e)
      
      temp <- if (length(dims) > 2L) quote(a[b,,]) else quote(a[b,])
      
      temp[[2L]] <- e[[2L]]
      temp[[3L]] <- e[[3L]]
      if (any(names(e) %in% "drop")) temp[["drop"]] <- e[["drop"]]
      return(temp)
    }
    
    if (!is.symbol(e)) for (i in seq_along(e)) e[[i]] <- subDims(e[[i]], env)
    
    e
  }
  
  e <- match.call()$e
  env <- parent.frame()
  
  tryResult <- tryCatch(result <- eval(subDims(e, env), env), error = function(e) e)
  if (inherits(tryResult, "error")) stop(tryResult)
  result
}

getArrayIndicesForOffset <- function(i, d)
{
  res <- rep(NA, length(d))
  stride <- prod(d[-length(d)])
  if (length(d) > 1L) for (j in seq.int(length(d), 2L)) {
    res[j] <- (i - 1L) %/% stride + 1L
    i <- (i - 1L) %% stride + 1L
    stride <- stride %/% d[j - 1L]
  }
  res[1L] <- i
  res
}

anyBars <- function(expr)
  any(c("|", "||") %in% all.names(expr))

issueWarningForUnknownArguments <- function() {
  matchedCall <- match.call(
    sys.function(sys.parent(1L)),
    sys.call(sys.parent(1L))
  )
  matchedFunction <- sys.function(sys.parent(1L))

  x <- NULL # for R CMD check
  if (
    length(
      unknownArgs <-
        names(matchedCall)[-1L]
        %as_x_in%
        x[x %not_in% names(formals(matchedFunction))]
    ) > 0L
  )
  {
    warning(
      'In ', deparse(matchedCall), ' :\n  ',
     'called with unknown argument(s): "',
      paste0(unknownArgs, collapse = '", "'),
      '"',
      call. = FALSE
    )
  }
}
