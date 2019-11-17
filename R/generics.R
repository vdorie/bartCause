averageDifferences <- function(samples.indiv.diff, treatmentRows, weights, estimand, commonSup.sub)
{
  x <- NULL ## R CMD check
  
  if (!is.character(estimand) || estimand[1L] %not_in% c("ate", "att", "atc"))
    stop("estimand must be one of 'ate', 'att', or 'atc'")
  estimand <- estimand[1L]
  
  origDims <- dim(samples.indiv.diff)
    
  if (!is.null(weights)) {
    weights <- rep_len(weights, origDims[1L])
    weights <- weights / sum(weights)
  }
  
  result <- 
    if (length(origDims) > 2L) {
      apply(evalx(if (is.null(weights)) samples.indiv.diff else samples.indiv.diff * weights,
                  switch(estimand,
                         att = x[ treatmentRows & commonSup.sub,,],
                         atc = x[!treatmentRows & commonSup.sub,,],
                         ate = x[commonSup.sub,,])),
            c(2L, 3L), mean)
    } else {
      apply(evalx(if (is.null(weights)) samples.indiv.diff else samples.indiv.diff * weights,
                  switch(estimand,
                         att = x[ treatmentRows & commonSup.sub,],
                         atc = x[!treatmentRows & commonSup.sub,],
                         ate = x[commonSup.sub,])),
            2L, mean)
    }
  
  if (!is.null(origDims) && length(origDims) > 2L)
    matrix(result, origDims[2L], origDims[3L])
  else 
    result
}

getEstimateSamples <- function(samples.indiv.diff, treatmentRows, weights, estimand, group.by, commonSup.sub) {
  if (is.null(group.by)) {
    samples.est <- averageDifferences(samples.indiv.diff, treatmentRows, weights, estimand, commonSup.sub)
  } else {
    samples.est <- lapply(levels(group.by), function(level) {
      levelRows <- group.by == level
      if (!is.null(weights)) weights <- weights[levelRows]
      
      averageDifferences(if (length(dim(samples.indiv.diff)) > 2L) samples.indiv.diff[levelRows,,] else samples.indiv.diff[levelRows,],
                         treatmentRows[levelRows], weights, estimand, commonSup.sub[levelRows])
    })
    names(samples.est) <- levels(group.by)
  }
  samples.est
}

predict.bartcFit <-
  function(object,
           newdata,
           value = c("mu.1", "mu.0", "icate", "p.score"),
           combineChains = TRUE, ...)
{
  value <- value[1L]
  if (value %not_in% eval(formals(predict.bartcFit)$value))
    stop("value must be in '", paste0(eval(formals(predict.bartcFit)$value), collapse = "', '"), "'")
  
  x.new <- as.data.frame(if (is.null(dim(newdata)) && length(newdata) > 0L)
                           matrix(newdata, ncol = length(newdata))
                         else
                           newdata)
  
  if (value == "p.score") {
    if (object$method.trt %in% c("given", "none"))
      stop("predict(value = 'p.score', ...) requires method.trt to specify a model")
    
    if (object$method.trt == "glm")
      p.score <- predict(object$fit.trt, x.new, type = "response", ...)
    else {
      if (is.null(object$fit.trt$fit))
        stop("predict with method.trt = '", object$method.trt, "' requires treatment model to be fit with keepTrees == TRUE")
      p.score <- pnorm(predict(object$fit.trt, x.new, ...))
      if (length(dim(p.score)) > 2L) p.score <- aperm(p.score, c(3L, 1L, 2L))
      if (combineChains) p.score <- flattenSamples(p.score)
    }
    return(p.score)
  }
  
  if (object$method.rsp != "bart")
    stop("predict(value = '", value, "', ...) requires method.rsp == 'bart'; other methods not designed to make individual predictions")
  
  if (is.null(object$fit.rsp$fit))
    stop("predict with method.rsp = 'bart' requires response model to be fit with keepTrees == TRUE")
  
  
  p.scoreName <- "ps"
  while (paste0(p.scoreName, "ps") %in% colnames(object$data.rsp@x)) p.scoreName <- paste0(p.scoreName, "ps")
  
  p.scoreAsCovariate <- !is.null(object$p.score) && p.scoreName %in% colnames(object$data.rsp@x)
  if (p.scoreAsCovariate && object$method.trt == "given")
    stop("predict requires fitting propensity scores to use in response model, however no treatment model exists");
    
  if (p.scoreAsCovariate) {
    if (object$method.trt == "glm")
      p.score <- predict(object$fit.trt, x.new, type = "response", ...)
    else {
      if (is.null(object$fit.trt$fit))
        stop("predict with method.trt = '", object$method.trt, "' and propensity scores as a covariate requires treatment model to be fit with keepTrees == TRUE")
      p.score <- pnorm(predict(object$fit.trt, x.new, ...))
    }
    if (!is.null(dim(p.score))) p.score <- apply(p.score, length(dim(p.score)), mean)
    x.new[[p.scoreName]] <- p.score
  }
  
  responseIsBinary <- is.null(object$fit.rsp[["sigma"]])
  T <- if (responseIsBinary) pnorm else function(x) x
  
  result <- switch(value,
    mu.1 = { x.new[[object$name.trt]] <- 1; T(predict(object$fit.rsp, x.new, ...)) },
    mu.0 = { x.new[[object$name.trt]] <- 0; T(predict(object$fit.rsp, x.new, ...)) },
    icate = {
      x.new[[object$name.trt]] <- 1
      mu.hat.1 <- predict(object$fit.rsp, x.new, ...)
      x.new[[object$name.trt]] <- 0
      T(mu.hat.1) - T(predict(object$fit.rsp, x.new, ...))
    })
  
  if (length(dim(result)) > 2L) result <- aperm(result, c(3L, 1L, 2L))
  
  if (combineChains) flattenSamples(result) else result
}

fitted.bartcFit <-
  function(object, 
           value = c("pate", "sate", "cate", "mu", "mu.0", "mu.1", "y.0", "y.1", "icate", "ite", "p.score", "p.weights"),
           sample = c("inferential", "all"),
           ...)
{
  if (!is.character(value) || value[1L] %not_in% eval(formals(fitted.bartcFit)$value))
    stop("value must be in '", paste0(eval(formals(fitted.bartcFit)$value), collapse = "', '"), "'")
  value <- value[1L]
  
  if (!is.character(sample) || sample[1L] %not_in% eval(formals(fitted.bartcFit)$sample))
    stop("sample must be in '", paste0(eval(formals(fitted.bartcFit)$sample), collapse = "', '"), "'")
  sample <- sample[1L]
  
  if (value == "p.weights" && is.null(object$p.score))
    stop("p.score cannot be NULL to obtain fitted p.weights")
  
  if (value == "p.score") return(object$p.score)
  
  result <- extract(object, value, sample, ...)
  
  if (!is.null(object$group.by) && value %in% c("pate", "sate", "cate"))
    return(sapply(result, function(result.i) {
      if (object$method.rsp %in% c("tmle", "p.weight") && value == "pate")
        return(mean(result.i))
      
      ifelse_3(!is.null(dim(result.i)), value != "p.score",
               apply(result.i, 1L, mean), mean(result.i), result.i)
    }))
  
  if (object$method.rsp %in% c("tmle", "p.weight") && value == "pate")
    return(mean(result))
  
  ifelse_3(!is.null(dim(result)), value != "p.score",
           apply(result, 1L, mean), mean(result), result)
}

extract.bartcFit <-
  function(object,
           value = c("pate", "sate", "cate", "mu", "mu.0", "mu.1", "y.0", "y.1", "icate", "ite", "p.score", "p.weights"),
           sample = c("inferential", "all"),
           combineChains = TRUE,
           ...)
{
  if (!is.character(value) || value[1L] %not_in% eval(formals(extract.bartcFit)$value))
    stop("value must be in '", paste0(eval(formals(extract.bartcFit)$value), collapse = "', '"), "'")
  value <- value[1L]
  
  if (!is.character(sample) || sample[1L] %not_in% eval(formals(extract.bartcFit)$sample))
    stop("sample must be in '", paste0(eval(formals(extract.bartcFit)$sample), collapse = "', '"), "'")
  sample <- sample[1L]
  
  if (value == "p.weights" && is.null(object$p.score))
    stop("p.score cannot be NULL to extract p.weights")
  
  if (object$method.rsp %in% c("p.weight", "tmle")) {
    if (value %in% c("sate", "cate"))
      stop("method '", object$method.rsp, "' does not produce estimates of ", value)
    if (value %in% c("mu", "mu.0", "mu.1", "y.0", "y.1"))
      warning("for method '", object$method.rsp, "' value '", value, "' does not have a clear interpretation")
    
    if (value == "pate") {
      result <- object$est
      return(
        if (is.null(object$group.by))
          ifelse_3(is.null(dim(result)), length(dim(result)) == 2L, result["est"], result[,"est"], result[,,"est"])
        else lapply(result, function(result.i)
          ifelse_3(is.null(dim(result.i)), length(dim(result.i)) == 2L, result.i["est"], result.i[,"est"], result.i[,,"est"]))
      )
    }
  }
  
  weights <- object$data.rsp@weights
  if (!is.null(weights)) weights <- weights / sum(weights)
  
  oldSeed <- .GlobalEnv$.Random.seed
  .GlobalEnv$.Random.seed <- object$seed
  
  n.obs     <- length(object$data.rsp@y)
  n.samples <- tail(dim(object$mu.hat.cf), 1L)
  n.chains  <- object$n.chains
  trtSign <- ifelse(object$trt == 1, 1, -1)
  
  responseIsBinary <- is.null(object$fit.rsp[["sigma"]])
  
  if (value %in% c("pate", "sate", "y.0", "y.1", "ite")) {
    y.cf <- with(object,
      if (responseIsBinary)
        rbinom(length(mu.hat.cf), 1L, mu.hat.cf)
      else
        t(t(flattenSamples(mu.hat.cf)) + rnorm(n.obs * length(fit.rsp$sigma), 0, fit.rsp$sigma))
    )
    if (object$n.chains > 1L) y.cf <- array(y.cf, c(n.obs, n.chains, n.samples))
  }
  if (value == "pate") {
    y.obs <- with(object,
      if (responseIsBinary)
        rbinom(length(mu.hat.obs), 1L, mu.hat.obs)
      else
        t(t(flattenSamples(mu.hat.obs)) + rnorm(n.obs * length(fit.rsp$sigma), 0, fit.rsp$sigma))
    )
    if (object$n.chains > 1L) y.obs <- array(y.obs, c(n.obs, n.chains, n.samples))
  }
  .GlobalEnv$.Random.seed <- oldSeed
  
  if (value %in% c("pate", "sate", "cate")) {
    
    samples.indiv.diff <-
      with(object, switch(value,
                          pate = y.obs - y.cf,
                          sate = mu.hat.obs - y.cf,
                          cate = mu.hat.obs - mu.hat.cf)) * trtSign
    
    result <- with(object,
      getEstimateSamples(samples.indiv.diff, trt > 0, weights, estimand, group.by, commonSup.sub))
    
    if (!is.null(object$group.by))
      return(if (combineChains) lapply(result, as.vector) else result)
    else
      return(if (combineChains) as.vector(result) else result)
  }
  
  result <-
    with(object, switch(value,
           mu          = mu.hat.obs,
           mu.1        = mu.hat.obs * trt       + mu.hat.cf * (1 - trt),
           mu.0        = mu.hat.obs * (1 - trt) + mu.hat.cf * trt,
           y.1         = mu.hat.obs * trt       +      y.cf * (1 - trt),
           y.0         = mu.hat.obs * (1 - trt) +      y.cf * trt,
           ite         = (mu.hat.obs -      y.cf) * trtSign,
           icate       = (mu.hat.obs - mu.hat.cf) * trtSign,
           p.score     = object$samples.p.score,
           p.weights   = getPWeights(estimand, trt, weights, if (!is.null(samples.p.score)) samples.p.score else p.score, fitPars$p.scoreBounds)))
  
  if (is.null(result)) return(NULL)
  
  if (combineChains) result <- flattenSamples(result)
  
  subset <- rep_len(TRUE, dim(result)[1L])
  if (sample == "inferential") {
    if (object$estimand == "att") subset <- object$trt > 0
    else if (object$estimand == "atc") subset <- object$trt <= 0
  }
  
  if (length(dim(result)) > 2L)
    result[subset,,]
  else
    result[subset,]
}

extract <- function(object, ...) UseMethod("extract")


refit.bartcFit <- function(object, newresp = NULL,
                           commonSup.rule = c("none", "sd", "chisq"),
                           commonSup.cut  = c(NA_real_, 1, 0.05), ...)
{
  matchedCall <- match.call()
  if (!is.null(newresp)) warning("'newresp' argument ignored, provided only for generic signature compatibility")
  
  if (!is.null(matchedCall$commonSup.rule)) {
     if (is.null(matchedCall$commonSup.cut))
       commonSup.cut <- eval(formals(refit.bartcFit)$commonSup.cut)[match(commonSup.rule, eval(formals(refit.bartcFit)$commonSup.rule))]
    commonSup.rule <- commonSup.rule[1L]
    commonSup.cut <- commonSup.cut[1L]
  } else {
    commonSup.rule <- "none"
    commonSup.cut  <- NA_real_
  }
  
  object$commonSup.rule <- commonSup.rule
  object$commonSup.cut  <- commonSup.cut
  
  object$commonSup.sub <- with(object, getCommonSupportSubset(sd.obs, sd.cf, commonSup.rule, commonSup.cut, trt, missingRows))
  commonSup.sub <- object$commonSup.sub
  
  
  treatmentRows <- object$trt > 0
  weights <- object$data.rsp@weights
  if (!is.null(weights)) weights <- weights / sum(weights)
  
  responseIsBinary <- is.null(object$fit.rsp[["sigma"]])
  T <- if (responseIsBinary) pnorm else function(x) x
  
  if (object$method.rsp == "bart") {
    samples.indiv.diff <- (T(object$yhat.obs) - T(object$yhat.cf)) * ifelse(treatmentRows, 1, -1)
    
    if (is.null(object$group.by)) {
      object$samples.est <- with(object, getBartEstimates(treatmentRows, weights, estimand, samples.indiv.diff, commonSup.sub))
    } else {
      object$samples.est <- lapply(levels(object$group.by), function(level) {
        levelRows <- object$group.by == level
        if (!is.null(weights)) weights <- weights[levelRows]
        
        with(object, getBartEstimates(treatmentRows[levelRows], weights, estimand,
                                      addDimsToSubset(samples.indiv.diff[levelRows, drop = FALSE]), commonSup.sub[levelRows]))
      })
      names(object$samples.est) <- levels(object$group.by)
    }
  } else if (object$method.rsp == "p.weight") {
    yhat.1 <- with(object, yhat.obs * trt       + yhat.cf * (1 - trt))
    yhat.0 <- with(object, yhat.obs * (1 - trt) + yhat.cf * trt)
    p.score <- object$p.score
    
    if (is.null(object$group.by)) {
      if (any(object$commonSup.sub != TRUE)) {
        addDimsToSubset(yhat.0 <- yhat.0[commonSup.sub, drop = FALSE])
        addDimsToSubset(yhat.1 <- yhat.1[commonSup.sub, drop = FALSE])
           
        p.score <- addDimsToSubset(p.score[commonSup.sub, drop = FALSE])
      
        if (!is.null(weights)) weights <- weights[commonSup.sub]
      }
      
      object$samples.est <- with(object, getPWeightEstimates(data.rsp@y[commonSup.sub], trt[commonSup.sub], weights, estimand, T(yhat.0), T(yhat.1), p.score, fitPars$yBounds, fitPars$p.scoreBounds))
    } else {
      object$samples.est <- lapply(levels(object$group.by), function(level) {
        levelRows <- object$group.by == level & object$commonSup.sub
        
        addDimsToSubset(yhat.0 <- yhat.0[levelRows, drop = FALSE])
        addDimsToSubset(yhat.1 <- yhat.1[levelRows, drop = FALSE])
        addDimsToSubset(p.score <- p.score[levelRows, drop = FALSE])
      
        if (!is.null(weights)) weights <- weights[levelRows]
      
        with(object, getPWeightEstimates(data.rsp@y[levelRows], trt[levelRows], weights, estimand, T(yhat.0), T(yhat.1), p.score,
                                         fitPars$yBounds, fitPars$p.scoreBounds))
      })
      names(object$samples.est) <- levels(object$group.by)
    }
  } else if (object$method.rsp == "tmle") {
    yhat.1 <- with(object, yhat.obs * trt       + yhat.cf * (1 - trt))
    yhat.0 <- with(object, yhat.obs * (1 - trt) + yhat.cf * trt)
    p.score <- object$p.score
    
    if (is.null(object$group.by)) {
      if (any(object$commonSup.sub != TRUE)) {
        addDimsToSubset(yhat.0 <- yhat.0[commonSup.sub, drop = FALSE])
        addDimsToSubset(yhat.1 <- yhat.1[commonSup.sub, drop = FALSE])
           
        p.score <- addDimsToSubset(p.score[commonSup.sub, drop = FALSE])
      
        if (!is.null(weights)) weights <- weights[commonSup.sub]
      }
      
      object$samples.est <- with(object, getTMLEEstimates(data.rsp@y[commonSup.sub], trt[commonSup.sub], weights, estimand, T(yhat.0), T(yhat.1), p.score, fitPars$yBounds, fitPars$p.scoreBounds, fitPars$depsilon, fitPars))
    } else {
      object$samples.est <- lapply(levels(object$group.by), function(level) {
        levelRows <- object$group.by == level & object$commonSup.sub
        
        addDimsToSubset(yhat.0 <- yhat.0[levelRows, drop = FALSE])
        addDimsToSubset(yhat.1 <- yhat.1[levelRows, drop = FALSE])
        addDimsToSubset(p.score <- p.score[levelRows, drop = FALSE])
      
        if (!is.null(weights)) weights <- weights[levelRows]
      
        with(object, getTMLEEstimates(data.rsp@y[levelRows], trt[levelRows], weights, estimand, T(yhat.0), T(yhat.1), p.score,
                                      fitPars$yBounds, fitPars$p.scoreBounds, fitPars$depsilon, fitPars$maxIter))
      })
      names(object$samples.est) <- levels(object$group.by)
    }

  }
  
  invisible(object)
}

refit <- function(object, newresp, ...) UseMethod("refit")

