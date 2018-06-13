getBARTFitForSubset <- function(object, observedSubset) {
  z <- ifelse(observedSubset, 1, 0)
  object$yhat.obs * z + object$yhat.cf * (1 - z)
}

fitted.bartcFit <-
  function(object,
           value = c("est", "y", "y0", "y1", "indiv.diff", "p.score"),
           sample = c("inferential", "all"),
           ...)
{
  if (!is.character(value) || value[1L] %not_in% eval(formals(fitted.bartcFit)$value))
    stop("value must be in '", paste0(eval(formals(fitted.bartcFit)$value), collapse = "', '"), "'")
  value <- value[1L]
  
  if (!is.character(sample) || sample[1L] %not_in% eval(formals(fitted.bartcFit)$sample))
    stop("sample must be in '", paste0(eval(formals(fitted.bartcFit)$sample), collapse = "', '"), "'")
  sample <- sample[1L]
  
  if (value == "est")
    return(if (!is.null(object$group.by)) sapply(object$samples.est, mean) else mean(object$samples.est))
  
  result <-
    switch(value,
           y           = apply(flattenSamples(object$yhat.obs), 1L, mean),
           y0          = apply(flattenSamples(getBARTFitForSubset(object, !object$trt)), 1L, mean),
           y1          = apply(flattenSamples(getBARTFitForSubset(object,  object$trt)), 1L, mean),
           indiv.diff  = apply(flattenSamples((object$yhat.obs - object$yhat.cf) * ifelse(object$trt, 1, -1)), 1L, mean),
           p.score     = object$p.score) 
  
  if (is.null(result)) return(NULL)
  
  subset <- rep_len(TRUE, length(result))
  if (sample == "inferential") {
    if (object$estimand == "att") subset <- object$trt
    else if (object$estimand == "atc") subset <- !object$trt
  }
  
  result[subset]
}

extract.bartcFit <-
  function(object,
           value = c("est", "y", "y0", "y1", "indiv.diff", "p.score"),
           sample = c("inferential", "all"),
           combineChains = TRUE,
           ...)
{
  value <- value[1L]
  if (value %not_in% eval(formals(extract.bartcFit)$value))
    stop("value must be in '", paste0(eval(formals(extract.bartcFit)$value), collapse = "', '"), "'")
  
  sample <- sample[1L]
  if (sample %not_in% eval(formals(extract.bartcFit)$sample))
    stop("sample must be in '", paste0(eval(formals(extract.bartcFit)$sample), collapse = "', '"), "'")
  
  if (value == "est") {
    if (!is.null(object$group.by))
      return(if (combineChains) lapply(object$samples.est, as.vector) else object$samples.est)
    else
      return(if (combineChains) as.vector(object$samples.est) else object$samples.est)
  }
  
  result <-
    switch(value,
           y           = object$yhat.obs,
           y0          = getBARTFitForSubset(object, !object$trt),
           y1          = getBARTFitForSubset(object,  object$trt),
           indiv.diff  = (object$yhat.obs - object$yhat.cf) * ifelse(object$trt, 1, -1),
           p.score     = object$samples.p.score)
  
  if (is.null(result)) return(NULL)
  
  if (combineChains) result <- flattenSamples(result)
  
  subset <- rep_len(TRUE, dim(result)[1L])
  if (sample == "inferential") {
    if (object$estimand == "att") subset <- object$trt
    else if (object$estimand == "atc") subset <- !object$trt
  }
  
  if (length(dim(result)) > 2L)
    result[subset,,]
  else
    result[subset,]
}

extract <- function(object, ...) UseMethod("extract")

