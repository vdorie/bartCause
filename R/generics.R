combineChains <- function(samples) {
  if (is.null(dim(samples)) || length(dim(samples)) == 2L) return(samples)
  
  x <- NULL
  res <- evalx(dim(samples), matrix(aperm(samples, c(2L, 1L, 3L)), prod(x[1L:2L]), x[3L]))
  if (!is.null(dimnames(samples)))
  dimnames(res) <- evalx(dimnames(samples), list(NULL, x[[length(x)]]))
  res
}

uncombineChains <- function(samples, n.chains) {
  if (is.null(dim(samples))) {
    if (n.chains == 1L) samples else matrix(samples, n.chains, length(samples) %/% n.chains)
  } else {
    res <- if (n.chains == 1L) samples else aperm(array(samples, c(dim(samples)[1L] %/% n.chains, n.chains, ncol(samples))), c(2L, 1L, 3L))
    if (!is.null(dimnames(samples))) dimnames(res) <- list(NULL, NULL, dimnames(samples)[[2L]])
    res
  }
}

multiplyArrayByVec <- function(x, v) {
  if (length(dim(x)) > 2L) {
    aperm(aperm(x, c(3L, 1L, 2L)) * v, c(2L, 3L, 1L))
  } else {
    t(t(x) * v)
  }
}

averageDifferences <- function(samples.indiv.diff, treatmentRows, weights, estimand, commonSup.sub)
{
  x <- NULL ## R CMD check
  
  if (!is.character(estimand) || estimand[1L] %not_in% c("ate", "att", "atc"))
    stop("estimand must be one of 'ate', 'att', or 'atc'")
  estimand <- estimand[1L]
  
  origDims <- dim(samples.indiv.diff)
  
  if (!is.null(weights)) {
    weights <- rep_len(weights, origDims[3L])
    weights <- weights / sum(weights)
    samples.indiv.diff <- multiplyArrayByVec(samples.indiv.diff, weights)
  }
  
  result <- 
    if (length(origDims) > 2L) {
      apply(switch(estimand,
                   att = samples.indiv.diff[,, treatmentRows & commonSup.sub],
                   atc = samples.indiv.diff[,,!treatmentRows & commonSup.sub],
                   ate = samples.indiv.diff[,,commonSup.sub]),
            c(1L, 2L), mean)
    } else {
      apply(switch(estimand,
                   att = samples.indiv.diff[, treatmentRows & commonSup.sub],
                   atc = samples.indiv.diff[,!treatmentRows & commonSup.sub],
                   ate = samples.indiv.diff[,commonSup.sub]),
            1L, mean)
    }
  
  if (!is.null(origDims) && length(origDims) > 2L)
    matrix(result, origDims[1L], origDims[2L])
  else 
    result
}

getEstimateSamples <- function(samples.indiv.diff, treatmentRows, weights, estimand, group.by, group.effects, commonSup.sub) {
  if (is.null(group.by) || !group.effects) {
    samples.est <- averageDifferences(samples.indiv.diff, treatmentRows, weights, estimand, commonSup.sub)
  } else {
    samples.est <- lapply(levels(group.by), function(level) {
      levelRows <- group.by == level
      if (!is.null(weights)) weights <- weights[levelRows]
      
      averageDifferences(if (length(dim(samples.indiv.diff)) > 2L) samples.indiv.diff[,,levelRows] else samples.indiv.diff[,levelRows],
                         treatmentRows[levelRows], weights, estimand, commonSup.sub[levelRows])
    })
    names(samples.est) <- levels(group.by)
  }
  samples.est
}

predict.bartcFit <-
  function(object, newdata,
           group.by,
           type = c("mu.0", "mu.1", "y.0", "y.1", "icate", "ite", "p.score"),
           combineChains = TRUE, ...)
{
  matchedCall <- match.call()
  
  if (!is.character(type) || type[1L] %not_in% eval(formals(predict.bartcFit)$type))
    stop("type must be in '", paste0(eval(formals(predict.bartcFit)$type), collapse = "', '"), "'")
  type <- type[1L]
  
  if (type != "p.score") {
    if (object$method.rsp != "bart")
      stop("predict(type = '", type, "', ...) requires method.rsp == 'bart'; other methods not designed to make individual predictions")
    
    if (is.null(object$fit.rsp$fit))
      stop("predict with method.rsp = 'bart' requires response model to be fit with keepTrees == TRUE")
    
    p.scoreName <- "ps"
    while (paste0(p.scoreName, "ps") %in% colnames(object$data.rsp@x)) p.scoreName <- paste0(p.scoreName, "ps")
    
    p.scoreAsCovariate <- !is.null(object$p.score) && p.scoreName %in% colnames(object$data.rsp@x)
    if (p.scoreAsCovariate && object$method.trt == "given")
      stop("predict requires fitting propensity scores to use in response model, however no treatment model exists");
  }
  
  x.new <- as.data.frame(if (is.null(dim(newdata)) && length(newdata) > 0L)
                           matrix(newdata, ncol = length(newdata))
                         else
                           newdata)
  
  if (!is.null(object[["group.by"]]))
    use.ranef <- !is.null(object[["use.ranef"]]) && object[["use.ranef"]]
  
  if (type == "p.score" || p.scoreAsCovariate) {
    if (object$method.trt %in% c("given", "none"))
      stop("predict(type = 'p.score', ...) requires method.trt to specify a model")
    
    if (object$method.trt == "glm") {
      if (!is.null(object[["group.by"]])) {
        # have to put group.by in correct place
        x.new.g <- x.new
        if (use.ranef) {
          # uses lmer
          x.new.g[[names(object$fit.trt@flist)]] <- group.by
          p.score <- predict(object$fit.trt, x.new.g, type = "response", ...)
        } else {
          # uses lm
          varNames <- attr(object$fit.trt$terms, "term.labels")
          x.new.g[[varNames[varNames %not_in% names(x.new)][1L]]] <- group.by
          p.score <- predict(object$fit.trt, x.new.g, type = "response", ...)
        }
      } else {
        p.score <- predict(object$fit.trt, x.new, ...)
      }
    } else {
      if (is.null(object$fit.trt$fit))
        stop("predict with method.trt = '", object$method.trt, "' requires treatment model to be fit with keepTrees == TRUE")
      
      if (!is.null(object[["group.by"]])) {
        if (use.ranef) {
          # uses rbart
          p.score <- predict(object$fit.trt, x.new, group.by, combineChains = FALSE, ...)
        } else {
          # uses base bart
          x.new.g <- x.new
          varNames <- attr(object$fit.trt$fit$data@x, "term.labels")
          x.new.g[[varNames[varNames %not_in% names(x.new)][1L]]] <- group.by
          
          p.score <- predict(object$fit.trt, x.new.g, combineChains = FALSE, ...)
        }
      } else {
        p.score <- predict(object$fit.trt, x.new, ...)
      }
    }
  }
  
  if (type == "p.score")
    return(if (combineChains) combineChains(p.score) else p.score)
  
  if (p.scoreAsCovariate) {
    if (!is.null(dim(p.score))) p.score <- apply(p.score, length(dim(p.score)), mean)
    x.new[[p.scoreName]] <- p.score
  }
  
  responseIsBinary <- is.null(object$fit.rsp[["sigma"]])
  
  if (!is.null(object$group.by)) {
    if (use.ranef) {
      predictArgs <- list(object$fit.rsp, x.new, group.by, combineChains = FALSE, ...)
    } else {
      x.new.g <- x.new
      varNames <- attr(object$fit.rsp$fit$data@x, "term.labels")
      x.new.g[[varNames[varNames %not_in% c(names(x.new), object$name.trt)][1L]]] <- group.by
      
      predictArgs <- list(object$fit.rsp, x.new.g, combineChains = FALSE, ...)
    }
  } else {
    predictArgs <- list(object$fit.rsp, x.new, combineChains = FALSE, ...)
  }
  
  if (type %in% c("mu.0", "y.0", "icate", "ite")) {
    predictArgs[[2L]][[object$name.trt]] <- 0
    mu.0 <- do.call("predict", predictArgs)
  }
  if (type %in% c("mu.1", "y.1", "icate", "ite")) {
    predictArgs[[2L]][[object$name.trt]] <- 1
    mu.1 <- do.call("predict", predictArgs)
  }
  
  if (type %in% c("y.0", "ite"))
    y.0 <- sampleFromPPD(object, mu.0)
  if (type %in% c("y.1", "ite"))
    y.1 <- sampleFromPPD(object, mu.1)
  
  result <-
    switch(type,
           mu.0  = mu.0,
           mu.1  = mu.1,
           icate = mu.1 - mu.0,
           y.0   = y.0,
           y.1   = y.1,
           ite   = y.1 - y.0)
    
  if (combineChains) combineChains(result) else result
}

fitted.bartcFit <-
  function(object, 
           type = c("pate", "sate", "cate", "mu.obs", "mu.cf", "mu.0", "mu.1",
                    "y.cf", "y.0", "y.1", "icate", "ite", "p.score", "p.weights"),
           sample = c("inferential", "all"),
           ...)
{
  if (!is.character(type) || type[1L] %not_in% eval(formals(fitted.bartcFit)$type))
    stop("type must be in '", paste0(eval(formals(fitted.bartcFit)$type), collapse = "', '"), "'")
  type <- type[1L]
  
  if (!is.character(sample) || sample[1L] %not_in% eval(formals(fitted.bartcFit)$sample))
    stop("sample must be in '", paste0(eval(formals(fitted.bartcFit)$sample), collapse = "', '"), "'")
  sample <- sample[1L]
  
  if (type == "p.weights" && is.null(object$p.score))
    stop("p.score cannot be NULL to obtain fitted p.weights")
  
  if (type == "p.score") {
    subset <- rep_len(TRUE, length(object$trt))
    if (sample == "inferential") {
      if (object$estimand == "att") subset <- object$trt > 0
      else if (object$estimand == "atc") subset <- !object$trt <= 0
    }
    return(object$p.score[subset])
  }
  
  result <- extract(object, type = type, sample = sample, ...)
  
  group.effects <- if (!is.null(object[["group.effects"]])) object[["group.effects"]] else FALSE
  if (!is.null(object$group.by) && group.effects && type %in% c("pate", "sate", "cate")) {
    return(sapply(result, function(result.i) {
      if (object$method.rsp %in% c("tmle", "p.weight") && type == "pate")
        return(mean(result.i))
      
      ifelse_3(!is.null(dim(result.i)), type != "p.score",
               apply(result.i, length(dim(result.i)), mean), mean(result.i), result.i)
    }))
  }
  
  if (object$method.rsp %in% c("tmle", "p.weight") && type == "pate")
    return(mean(result))
  
  if (!is.null(dim(result)))
    apply(result, length(dim(result)), mean)
  else
    mean(result)
}

extract.bartcFit <-
  function(object,
           type = c("pate", "sate", "cate", "mu.obs", "mu.cf", "mu.0", "mu.1", "y.cf",
                    "y.0", "y.1", "icate", "ite", "p.score", "p.weights"),
           sample = c("inferential", "all"),
           combineChains = TRUE,
           ...)
{
  if (!is.character(type) || type[1L] %not_in% eval(formals(extract.bartcFit)$type))
    stop("type must be in '", paste0(eval(formals(extract.bartcFit)$type), collapse = "', '"), "'")
  type <- type[1L]
  
  if (!is.character(sample) || sample[1L] %not_in% eval(formals(extract.bartcFit)$sample))
    stop("sample must be in '", paste0(eval(formals(extract.bartcFit)$sample), collapse = "', '"), "'")
  sample <- sample[1L]
  
  if (type == "p.weights" && is.null(object$p.score))
    stop("p.score cannot be NULL to extract p.weights")
  
  group.effects <- if (!is.null(object[["group.effects"]])) object[["group.effects"]] else FALSE
  
  if (object$method.rsp %in% c("p.weight", "tmle")) {
    if (type %in% c("sate", "cate"))
      stop("method '", object$method.rsp, "' does not produce estimates of ", type)
    if (type %in% c("mu", "mu.0", "mu.1", "y.0", "y.1"))
      warning("for method '", object$method.rsp, "' type '", type, "' does not have a clear interpretation")
    
    if (type == "pate") {
      result <- object$est
      return(
        if (is.null(object$group.by) || !group.effects)
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
  
  if (type %in% c("pate", "sate", "y.cf", "y.0", "y.1", "ite")) {
    y.obs <- object$fit.rsp$y
    y.cf <- sampleFromPPD(object, object$mu.hat.cf)
  }
  
  if (type == "pate")
    y.obs.ppd <- sampleFromPPD(object, object$mu.hat.obs)
  
  .GlobalEnv$.Random.seed <- oldSeed
  
  obsCfToTrtCtl <- function(obs, cf, trt) {
    if (length(dim(obs)) > 2L) {
      aperm(aperm(obs, c(3L, 1L, 2L)) * trt + aperm(cf, c(3L, 1L, 2L)) * (1 - trt), c(2L, 3L, 1L))
    } else {
      t(t(obs) * trt + t(cf) * (1 - trt))
    }
  }
  
  if (type %in% c("pate", "sate", "cate")) {
    
    samples.indiv.diff <- multiplyArrayByVec(with(object,
      switch(type,
             pate = y.obs.ppd - y.cf,
             sate = y.obs     - y.cf,
             cate = mu.hat.obs - mu.hat.cf)),
      trtSign)
    
    if (is.null(object$group.by)) group.by <- NULL
    result <- with(object,
      getEstimateSamples(samples.indiv.diff, trt > 0, weights, estimand, group.by, group.effects, commonSup.sub))
    
    if (!is.null(object$group.by) && group.effects)
      return(if (combineChains) lapply(result, as.vector) else result)
    else
      return(if (combineChains) as.vector(result) else result)
  }
  
  result <-
    with(object, switch(type,
           mu.obs      = mu.hat.obs,
           mu.cf       = mu.hat.cf,
           mu.1        = obsCfToTrtCtl(mu.hat.obs, mu.hat.cf, trt),
           mu.0        = obsCfToTrtCtl(mu.hat.obs, mu.hat.cf, 1 - trt),
           y.cf        = y.cf,
           y.1         = obsCfToTrtCtl(y.obs, y.cf, trt),
           y.0         = obsCfToTrtCtl(y.obs, y.cf, 1 - trt),
           ite         = multiplyArrayByVec(     y.obs -      y.cf, trtSign),
           icate       = multiplyArrayByVec(mu.hat.obs - mu.hat.cf, trtSign),
           p.score     = object$samples.p.score,
           p.weights   = getPWeights(estimand, trt, weights, if (!is.null(samples.p.score)) samples.p.score else p.score, fitPars$p.scoreBounds)))
  
  if (is.null(result)) return(NULL)
  
  if (combineChains) result <- combineChains(result)
  
  subset <- rep_len(TRUE, dim(result)[length(dim(result))])
  if (sample == "inferential") {
    if (object$estimand == "att") subset <- object$trt > 0
    else if (object$estimand == "atc") subset <- object$trt <= 0
  }
  
  if (length(dim(result)) > 2L)
    result[,,subset]
  else
    result[,subset]
}

sampleFromPPD <- function(object, ev)
{
  responseIsBinary <- is.null(object$fit.rsp[["sigma"]])
  
  if (responseIsBinary) {
    if (length(dim(ev)) > 2L) {
      result <- array(rbinom(length(ev), 1L, ev), dim(ev), dimnames = dimnames(ev))
    } else {
      result <- matrix(rbinom(length(ev), 1L, ev), nrow(ev), ncol(ev), dimnames = list(rownames(ev), colnames(ev)))
    }
  } else {
    n.obs <- dim(ev)[length(dim(ev))]
    result <- ev + rnorm(n.obs * length(object$fit.rsp$sigma), 0, rep_len(object$fit.rsp$sigma, n.obs * length(object$fit.rsp$sigma)))
  }
  
  result
}

# extract <- function(object, ...) UseMethod("extract")


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
  group.effects <- if (!is.null(object[["group.effects"]])) object[["group.effects"]] else FALSE
  group.by <- if (!is.null(object[["group.by"]])) object[["group.by"]] else NULL
  
  if (object$method.rsp == "bart") {
    samples.indiv.diff <- extract(object, value = "icate", combineChains = FALSE)
    
    object$est <- with(object,
      getEstimateSamples(samples.indiv.diff, treatmentRows, weights, estimand, group.by, group.effects, commonSup.sub))
   
  
  } else if (object$method.rsp == "p.weight") {
    mu.hat.0 <- extract(object, "mu.hat.0", combineChains = FALSE)
    mu.hat.1 <- extract(object, "mu.hat.1", combineChains = FALSE)
    if (length(dim(mu.hat.0)) > 2L) {
      mu.hat.0 <- aperm(mu.hat.0, c(3L, 1L, 2L))
      mu.hat.1 <- aperm(mu.hat.1, c(3L, 1L, 2L))
    } else {
      mu.hat.0 <- t(mu.hat.0)
      mu.hat.1 <- t(mu.hat.1)
    }
    
    p.score <- if (!is.null(object$samples.p.score)) object$samples.p.score else object$p.score
    if (!is.null(dim(p.score)) && length(dim(p.score)) < length(dim(mu.hat.0))) {
      # chains were collapsed
      n.chains  <- dim(mu.hat.0)[2L]
      n.samples <- dim(mu.hat.0)[3L]
      n.obs     <- dim(mu.hat.0)[1L]
      p.score   <- aperm(array(p.score, c(n.chains, n.obs, n.samples)), c(3L, 1L, 2L))
    } else {
      if (!is.null(dim(p.score)))
        p.score <- if (length(dim(p.score)) > 2L) aperm(p.score, c(3L, 1L, 2L)) else t(p.score)
    }
    
    if (is.null(object$group.by) || !group.effects) {
      if (any(object$commonSup.sub != TRUE)) {
        addDimsToSubset(mu.hat.0 <- mu.hat.0[commonSup.sub, drop = FALSE])
        addDimsToSubset(mu.hat.1 <- mu.hat.1[commonSup.sub, drop = FALSE])
           
        p.score <- addDimsToSubset(p.score[commonSup.sub, drop = FALSE])
      
        if (!is.null(weights)) weights <- weights[commonSup.sub]
      }
      
      object$est <- with(object, getPWeightEstimates(data.rsp@y[commonSup.sub], trt[commonSup.sub], weights, estimand, mu.hat.0, mu.hat.1, p.score, fitPars$yBounds, fitPars$p.scoreBounds))
    } else {
      object$est <- lapply(levels(object$group.by), function(level) {
        levelRows <- object$group.by == level & object$commonSup.sub
        
        addDimsToSubset(mu.hat.0 <- mu.hat.0[levelRows, drop = FALSE])
        addDimsToSubset(mu.hat.1 <- mu.hat.1[levelRows, drop = FALSE])
        addDimsToSubset(p.score  <- p.score[levelRows, drop = FALSE])
      
        if (!is.null(weights)) weights <- weights[levelRows]
      
        with(object, getPWeightEstimates(data.rsp@y[levelRows], trt[levelRows], weights, estimand, mu.hat.0, mu.hat.1, p.score,
                                         fitPars$yBounds, fitPars$p.scoreBounds))
      })
      names(object$est) <- levels(object$group.by)
    }
  } else if (object$method.rsp == "tmle") {
    mu.hat.0 <- extract(object, "mu.hat.0", combineChains = FALSE)
    mu.hat.1 <- extract(object, "mu.hat.1", combineChains = FALSE)
    if (length(dim(mu.hat.0)) > 2L) {
      mu.hat.0 <- aperm(mu.hat.0, c(3L, 1L, 2L))
      mu.hat.1 <- aperm(mu.hat.1, c(3L, 1L, 2L))
    } else {
      mu.hat.0 <- t(mu.hat.0)
      mu.hat.1 <- t(mu.hat.1)
    }
    
    p.score <- if (!is.null(object$samples.p.score)) object$samples.p.score else object$p.score
    if (!is.null(dim(p.score)) && length(dim(p.score)) < length(dim(mu.hat.0))) {
      # chains were collapsed
      n.chains  <- dim(mu.hat.0)[2L]
      n.samples <- dim(mu.hat.0)[3L]
      n.obs     <- dim(mu.hat.0)[1L]
      p.score   <- aperm(array(p.score, c(n.chains, n.obs, n.samples)), c(3L, 1L, 2L))
    } else {
      if (!is.null(dim(p.score)))
        p.score <- if (length(dim(p.score)) > 2L) aperm(p.score, c(3L, 1L, 2L)) else t(p.score)
    }
    
    if (is.null(object$group.by) || !group.effects) {
      if (any(object$commonSup.sub != TRUE)) {
        addDimsToSubset(mu.hat.0 <- mu.hat.0[commonSup.sub, drop = FALSE])
        addDimsToSubset(mu.hat.1 <- mu.hat.1[commonSup.sub, drop = FALSE])
           
        addDimsToSubset(p.score <- p.score[commonSup.sub, drop = FALSE])
      
        if (!is.null(weights)) weights <- weights[commonSup.sub]
      }
      
      object$est <- with(object, getTMLEEstimates(data.rsp@y[commonSup.sub], trt[commonSup.sub], weights, estimand, mu.hat.0, mu.hat.1, p.score, fitPars$yBounds, fitPars$p.scoreBounds, fitPars$depsilon, fitPars))
    } else {
      object$est <- lapply(levels(object$group.by), function(level) {
        levelRows <- object$group.by == level & object$commonSup.sub
        
        addDimsToSubset(yhat.0 <- yhat.0[levelRows, drop = FALSE])
        addDimsToSubset(yhat.1 <- yhat.1[levelRows, drop = FALSE])
        addDimsToSubset(p.score <- p.score[levelRows, drop = FALSE])
      
        if (!is.null(weights)) weights <- weights[levelRows]
      
        with(object, getTMLEEstimates(data.rsp@y[levelRows], trt[levelRows], weights, estimand, mu.hat.0, mu.hat.1, p.score,
                                      fitPars$yBounds, fitPars$p.scoreBounds, fitPars$depsilon, fitPars$maxIter))
      })
      names(object$est) <- levels(object$group.by)
    }

  }
  
  invisible(object)
}

refit <- function(object, newresp, ...) UseMethod("refit")

