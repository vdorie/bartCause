combineChains <- function(samples, n.chains) {
  if (is.null(dim(samples)) || length(dim(samples)) <= 1L) return(samples)
  if (length(dim(samples)) == 2L) {
    if (n.chains == 1L) return(samples)
    else return(as.vector(t(samples)))
  }
  
  x <- NULL
  res <- evalx(dim(samples), matrix(aperm(samples, c(2L, 1L, 3L)), prod(x[1L:2L]), x[3L]))
  if (!is.null(dimnames(samples)))
  dimnames(res) <- evalx(dimnames(samples), list(NULL, x[[length(x)]]))
  res
}

multiplyArrayByVec <- function(x, v) {
  if (length(dim(x)) > 2L) {
    aperm(aperm(x, c(3L, 1L, 2L)) * v, c(2L, 3L, 1L))
  } else {
    t(t(x) * v)
  }
}

# v - x
subtractArrayFromVec <- function(v, x) {
  if (length(dim(x)) > 2L) {
    aperm(v - aperm(x, c(3L, 1L, 2L)), c(2L, 3L, 1L))
  } else {
    t(v - t(x))
  }
}

## The average of the unit differences over the estimand's rows within common
## support, per draw. Population weights are shares of those rows alone, so a
## weighted ATT is the weighted mean over the treated; rows of no weight are
## left out, their draws being NA where the response fit has no noise scale.
## Rows of no weight in all are an error, or NA draws where 'emptyIsNA' (a
## group of a grouped fit, which leaves the other groups to be estimated).
averageDifferences <- function(samples.indiv.diff, treatmentRows, weights, estimand, commonSup.sub,
                               emptyIsNA = FALSE)
{
  if (!is.character(estimand) || estimand[1L] %not_in% c("ate", "att", "atc"))
    stop("estimand must be one of 'ate', 'att', or 'atc'")
  estimand <- estimand[1L]
  
  origDims <- dim(samples.indiv.diff)
  rows <- switch(estimand,
                 att = treatmentRows & commonSup.sub,
                 atc = !treatmentRows & commonSup.sub,
                 ate = commonSup.sub)
  
  if (is.null(weights)) {
    return(
      if (length(origDims) > 2L) {
        ## drop = FALSE: a one-chain fit's leading chain margin has extent 1,
        ## which R's default drop = TRUE would silently remove along with the
        ## subsetted margin, leaving apply()'s MARGIN = c(1L, 2L) result
        ## transposed relative to what matrix() below expects
        matrix(apply(samples.indiv.diff[,,rows, drop = FALSE], c(1L, 2L), mean), origDims[1L], origDims[2L])
      } else {
        apply(samples.indiv.diff[,rows, drop = FALSE], 1L, mean)
      })
  }
  
  if (emptyIsNA && !(sum(weights[rows]) > 0))
    return(if (length(origDims) > 2L) matrix(NA_real_, origDims[1L], origDims[2L]) else rep(NA_real_, origDims[1L]))
  shares <- weightShares(weights, rows)
  rows <- which(rows)[shares > 0]
  shares <- shares[shares > 0]
  if (length(origDims) > 2L) {
    result <- matrix(0, origDims[1L], origDims[2L])
    for (i in seq_len(origDims[1L]))
      result[i,] <- matrix(samples.indiv.diff[i,,rows], origDims[2L], length(rows)) %*% shares
    result
  } else {
    as.vector(samples.indiv.diff[,rows, drop = FALSE] %*% shares)
  }
}

getEstimateSamples <- function(samples.indiv.diff, treatmentRows, weights, estimand, group.by, group.effects, commonSup.sub) {
  if (is.null(group.by) || !group.effects) {
    samples.est <- averageDifferences(samples.indiv.diff, treatmentRows, weights, estimand, commonSup.sub)
  } else {
    samples.est <- lapply(levels(group.by), function(level) {
      levelRows <- group.by == level
      if (!is.null(weights)) weights <- weights[levelRows]
      
      averageDifferences(if (length(dim(samples.indiv.diff)) > 2L) samples.indiv.diff[,,levelRows, drop = FALSE] else samples.indiv.diff[,levelRows],
                         treatmentRows[levelRows], weights, estimand, commonSup.sub[levelRows], emptyIsNA = TRUE)
    })
    names(samples.est) <- levels(group.by)
  }
  samples.est
}

## A fit saved by bartCause 1.0-10 or earlier holds samplers from an earlier
## dbarts that this one cannot replay, and it predates the recorded name of the
## propensity score column; every fit made since records that name, even when
## it is NULL, so its absence marks the old ones, and the failure that follows
## without this check (a missing 'ps' column, a missing slot in a dbarts
## object) does not say what to do.
stopIfFitIsFromEarlierVersion <- function(object)
{
  if ("name.p.score" %not_in% names(object))
    stop("this fit was made by an earlier version of bartCause, with an earlier version of dbarts, ",
         "and predict cannot use it; refit it with the current versions. ",
         "print, summary, extract, fitted and refit still read it")
}

## Saved trees are what predict replays, whichever sampler fit the model, and
## bartc's 'keepTrees' reaches both the dbarts and the stan4bart fits, so it is
## the one argument to name.
noTreesMessage <- function(model, method)
  paste0("predict with ", method, " requires the ", model, " model to be fit with keepTrees = TRUE; ",
         "pass keepTrees = TRUE to bartc")

predict.bartcFit <-
  function(object, newdata,
           group.by,
           type = c("mu", "y", "mu.0", "mu.1", "y.0", "y.1", "icate", "ite", "p.score"),
           combineChains = TRUE, ...)
{
  matchedCall <- match.call()
  
  stopIfFitIsFromEarlierVersion(object)
  
  if (!is.character(type) || type[1L] %not_in% eval(formals(predict.bartcFit)$type))
    stop("type must be in '", paste0(eval(formals(predict.bartcFit)$type), collapse = "', '"), "'")
  type <- type[1L]

  if (type != "p.score") {
    ## bcf is designed to make individual predictions; what it cannot do is make
    ## them out of sample, so it gets its own refusal rather than the one below
    if (object$method.rsp == "bcf")
      stop("predict(type = '", type, "', ...) is not available for method.rsp == 'bcf'; out-of-sample mu(x) and tau(x) need per-forest saved-tree replay, which dbarts does not expose")
    if (object$method.rsp != "bart")
      stop("predict(type = '", type, "', ...) requires method.rsp == 'bart'; other methods not designed to make individual predictions")
    
    if (is.null(if (inherits(object$fit.rsp, "stan4bartFit")) object$fit.rsp$sampler.bart else object$fit.rsp$fit))
      stop(noTreesMessage("response", "method.rsp = 'bart'"))
    
    ## the column the response fitter put the score in, taken from the fit
    ## rather than re-derived: a name ladder over the design columns cannot
    ## tell the score from a confounder sharing its stem (a confounder "psps"
    ## against a score "ps" resolved to "psps", and the score was then written
    ## over the confounder). It is NULL exactly when the response model has no
    ## score column, which is what the membership test used to stand in for.
    p.scoreName <- object$name.p.score
    p.scoreAsCovariate <- !is.null(object$p.score) && !is.null(p.scoreName)
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
        p.score <- predict(object$fit.trt, x.new, type = "response", ...)
      }
    } else {
      if (is.null(if (inherits(object$fit.trt, "stan4bartFit")) object$fit.trt$sampler.bart else object$fit.trt$fit))
        stop(noTreesMessage("treatment", paste0("method.trt = '", object$method.trt, "'")))
      
      x.new.trt <- x.new
      if (!is.null(object[["group.by"]])) {
        # the grouping factor is a column of the design either way: named by the
        # random-effect term when it is modeled, and by the one design column
        # the new data does not carry when it is not
        if (use.ranef) {
          x.new.trt[[names(object$fit.trt$reTrms$flist)[1L]]] <- group.by
        } else {
          varNames <- attr(object$fit.trt$fit$data@x, "term.labels")
          x.new.trt[[varNames[varNames %not_in% names(x.new)][1L]]] <- group.by
        }
      }
      if (inherits(object$fit.trt, "stan4bartFit")) {
        # stan4bart returns n.obs x n.samples x n.chains
        p.score <- predict(object$fit.trt, x.new.trt, combine_chains = FALSE, ...)
        p.score <- aperm(p.score, c(3L, 2L, 1L))
      } else {
        p.score <- predict(object$fit.trt, x.new.trt, combineChains = FALSE, ...)
        ## dbarts >= 1.0-0 keeps the one-chain margin predict() used to drop
        p.score <- dropSingleChainDim(p.score, object$n.chains, 2L)
      }
    }
  }
  
  n.chains <- object$n.chains
  if (type == "p.score")
    return(
      if (object$method.trt %in% "glm") p.score
      else if (combineChains) combineChains(p.score, n.chains)
      else if (n.chains == 1L) addSingleChainDim(p.score, 2L)
      else p.score)
  
  if (p.scoreAsCovariate) {
    if (!is.null(dim(p.score))) p.score <- apply(p.score, length(dim(p.score)), mean)
    x.new[[p.scoreName]] <- p.score
  }
    
  if (!is.null(object$group.by)) {
    x.new.g <- x.new
    if (use.ranef) {
      x.new.g[[names(object$fit.rsp$reTrms$flist)[1L]]] <- group.by
    } else {
      varNames <- attr(object$fit.rsp$fit$data@x, "term.labels")
      x.new.g[[varNames[varNames %not_in% c(names(x.new), object$name.trt)][1L]]] <- group.by
    }
    predictArgs <- list(object$fit.rsp, x.new.g, combineChains = FALSE, ...)
  } else {
    predictArgs <- list(object$fit.rsp, x.new, combineChains = FALSE, ...)
  }
  
  if (inherits(predictArgs[[1L]], "stan4bartFit") && any(names(predictArgs) == "combineChains"))
      evalx(names(predictArgs), x[x == "combineChains"] <- "combine_chains")
  
  if (type %in% c("mu", "y")) {
    if (is.null(x.new[[object$name.trt]]) || anyNA(x.new[[object$name.trt]]))
      stop("for predict type '", type, "', newdata must have '", object$name.trt, "' column filled")
    
    mu <- do.call("predict", predictArgs)

    if (inherits(predictArgs[[1L]], "stan4bartFit")) {
      mu <- aperm(mu, c(3L, 2L, 1L))
    } else {
      ## dbarts >= 1.0-0 keeps the one-chain margin predict() used to drop
      mu <- dropSingleChainDim(mu, object$n.chains, 2L)
    }

    if (type == "y")
      y <- sampleFromPPD(object, mu)
  }

  if (type %in% c("mu.0", "y.0", "icate", "ite")) {
    predictArgs[[2L]][[object$name.trt]] <- 0
    mu.0 <- do.call("predict", predictArgs)

    if (inherits(predictArgs[[1L]], "stan4bartFit")) {
      mu.0 <- aperm(mu.0, c(3L, 2L, 1L))
    } else {
      mu.0 <- dropSingleChainDim(mu.0, object$n.chains, 2L)
    }
  }
  if (type %in% c("mu.1", "y.1", "icate", "ite")) {
    predictArgs[[2L]][[object$name.trt]] <- 1
    mu.1 <- do.call("predict", predictArgs)

    if (inherits(predictArgs[[1L]], "stan4bartFit")) {
      mu.1 <- aperm(mu.1, c(3L, 2L, 1L))
    } else {
      mu.1 <- dropSingleChainDim(mu.1, object$n.chains, 2L)
    }
  }
  
  if (type %in% c("y.0", "ite"))
    y.0 <- sampleFromPPD(object, mu.0)
  if (type %in% c("y.1", "ite"))
    y.1 <- sampleFromPPD(object, mu.1)
  
  result <-
    switch(type,
           mu    = mu,
           y     = y,
           mu.0  = mu.0,
           mu.1  = mu.1,
           icate = mu.1 - mu.0,
           y.0   = y.0,
           y.1   = y.1,
           ite   = y.1 - y.0)
    
  if (combineChains) combineChains(result, n.chains)
  else if (n.chains == 1L) addSingleChainDim(result, 2L)
  else result
}

fitted.bartcFit <-
  function(object, 
           type = c("pate", "sate", "cate", "mu.obs", "mu.cf", "mu.0", "mu.1",
                    "y.cf", "y.0", "y.1", "icate", "ite", "p.score", "p.weights"),
           sample = c("inferential", "all"),
           ...)
{
  issueWarningForUnknownArguments()

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
  
  result <- extractBartcFit(object, type, sample, combineChains = TRUE)
  
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
  else if (type == "p.weights")
    result ## no draws: one weight per observation already
  else
    mean(result)
}

extract.bartcFit <-
  function(object,
           type = c("pate", "sate", "cate", "mu.obs", "mu.cf", "mu.0", "mu.1", "y.cf",
                    "y.0", "y.1", "icate", "ite", "p.score", "p.weights", "sigma"),
           sample = c("inferential", "all"),
           combineChains = TRUE,
           ...)
{
  issueWarningForUnknownArguments()

  if (!is.character(type) || type[1L] %not_in% eval(formals(extract.bartcFit)$type))
    stop("type must be in '", paste0(eval(formals(extract.bartcFit)$type), collapse = "', '"), "'")
  type <- type[1L]

  if (!is.character(sample) || sample[1L] %not_in% eval(formals(extract.bartcFit)$sample))
    stop("sample must be in '", paste0(eval(formals(extract.bartcFit)$sample), collapse = "', '"), "'")
  sample <- sample[1L]

  result <- extractBartcFit(object, type, sample, combineChains)
  if (combineChains || object$n.chains != 1L || is.null(result)) return(result)

  ## a one-chain result keeps its chain margin, as dbarts and stan4bart do:
  ## 1 x n.samples for a scalar draw, 1 x n.samples x n.obs for a per-observation one
  oldRank <- if (type %in% c("pate", "sate", "cate", "sigma")) 1L else 2L
  if (is.list(result)) lapply(result, addSingleChainDim, oldRank)
  else addSingleChainDim(result, oldRank)
}

## Works in the stored layout, where a one-chain result has no chain margin;
## extract.bartcFit() adds the margin on the way out. Internal callers use this.
extractBartcFit <-
  function(object, type, sample = "inferential", combineChains = TRUE)
{
  
  if (type == "p.weights" && is.null(object$p.score))
    stop("p.score cannot be NULL to extract p.weights")

  n.chains <- object$n.chains
  if (type == "sigma") {
    if (responseIsBinary(object))
      stop("binary response model does not have a residual standard deviation parameter (sigma)")
    held <- heldSigma(object)
    sigma <-
      if (!is.null(held)) {
        ## a held sigma is one number; repeat it in the drawn layout, which is
        ## [n.chains, n.samples] with the chain margin dropped at one chain
        n.samples <- dim(object$mu.hat.obs)[length(dim(object$mu.hat.obs)) - 1L]
        if (n.chains > 1L) matrix(held, n.chains, n.samples) else rep(held, n.samples)
      } else if (inherits(object$fit.rsp, "stan4bartFit")) {
        t(extract(object$fit.rsp, "sigma", combine_chains = FALSE))
      } else {
        ## dbarts returns sigma pre-combined as a length n.chains*n.samples,
        ## chain-major vector (chain 1's whole run, then chain 2's, ...);
        ## reshape to [n.chains, n.samples] by filling column-major with
        ## n.samples rows (each column is one chain's contiguous block) and
        ## transposing, so combineChains = FALSE yields per-chain draws and
        ## = TRUE recombines to the original chain-major order
        s <- object$fit.rsp$sigma
        if (is.null(dim(s)) && n.chains > 1L) t(matrix(s, ncol = n.chains)) else s
      }
    ## combined is a plain vector at every chain count
    return(if (combineChains) as.vector(t(sigma)) else sigma)
  }
  
  group.effects <- if (!is.null(object[["group.effects"]])) object[["group.effects"]] else FALSE
  
  if (object$method.rsp %in% c("p.weight", "tmle")) {
    if (type %in% c("sate", "cate"))
      stop("method '", object$method.rsp, "' does not produce estimates of ", type)
    if (type %in% c("mu", "mu.0", "mu.1", "y.0", "y.1"))
      warning("for method '", object$method.rsp, "' type '", type, "' does not have a clear interpretation")
    
    if (type == "pate") {
      pateOf <- function(result) {
        result <- ifelse_3(is.null(dim(result)), length(dim(result)) == 2L, result["est"], result[,"est"], result[,,"est"])
        if (combineChains && length(dim(result)) == 2L) as.vector(t(result)) else result
      }
      result <- object$est
      return(
        if (is.null(object$group.by) || !group.effects) pateOf(result)
        else lapply(result, pateOf)
      )
    }
  }
  
  ## not 'weights', which with(object, ...) below would shadow
  pop <- populationWeights(object)
  
  oldSeed <- .GlobalEnv[[".Random.seed"]]
  .GlobalEnv$.Random.seed <- object$seed
  
  n.chains  <- object$n.chains
  trtSign <- ifelse(object$trt == 1, 1, -1)
  
  ## a unit's noise is drawn as the response fit has it; inside an average a
  ## row of counts is the average of its units
  noise <- noiseScale(object, average = type %in% c("pate", "sate"))
  if (type %in% c("pate", "sate", "y.cf", "y.0", "y.1", "ite")) {
    y.obs <- observedResponse(object)
    y.cf <- sampleFromPPD(object, object$mu.hat.cf, noise)
  }
  
  if (type == "pate" || (type %in% c("sate", "y.0", "y.1", "ite") && anyNA(y.obs)))
    y.obs.ppd <- sampleFromPPD(object, object$mu.hat.obs, noise)
  ## a missing response is drawn from its posterior predictive, as a
  ## counterfactual is
  if (type %in% c("sate", "y.0", "y.1", "ite") && anyNA(y.obs)) {
    observed <- !is.na(y.obs)
    draws <- y.obs.ppd
    if (length(dim(draws)) > 2L) draws[,,observed] <- rep(y.obs[observed], each = prod(dim(draws)[1L:2L]))
    else draws[,observed] <- rep(y.obs[observed], each = nrow(draws))
    y.obs <- draws
  }
  
  if (!is.null(oldSeed))
    .GlobalEnv$.Random.seed <- oldSeed
  else
    rm(list = ".Random.seed", envir = .GlobalEnv)
  
  obsCfToTrtCtl <- function(obs, cf, trt) {
    if (is.null(dim(obs))) {
      if (length(dim(cf)) > 2L) {
        aperm(obs * trt + aperm(cf, c(3L, 1L, 2L)) * (1 - trt), c(2L, 3L, 1L))
      } else {
        t(obs * trt + t(cf) * (1 - trt))
      }
    } else if (length(dim(obs)) > 2L) {
      aperm(aperm(obs, c(3L, 1L, 2L)) * trt + aperm(cf, c(3L, 1L, 2L)) * (1 - trt), c(2L, 3L, 1L))
    } else {
      t(t(obs) * trt + t(cf) * (1 - trt))
    }
  }
  
  if (type %in% c("pate", "sate", "cate")) {
    
    samples.indiv.diff <- multiplyArrayByVec(with(object,
      switch(type,
             pate = y.obs.ppd - y.cf,
             sate = if (is.null(dim(y.obs))) subtractArrayFromVec(y.obs, y.cf) else y.obs - y.cf,
             cate = mu.hat.obs - mu.hat.cf)),
      trtSign)
    
    if (is.null(object$group.by)) group.by <- NULL
    result <- with(object,
      getEstimateSamples(samples.indiv.diff, trt > 0, pop, estimand, group.by, group.effects, commonSup.sub))
    
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
           ite         = multiplyArrayByVec(if (is.null(dim(y.obs))) subtractArrayFromVec(y.obs, y.cf) else y.obs - y.cf, trtSign),
           icate       = multiplyArrayByVec(mu.hat.obs - mu.hat.cf, trtSign),
           p.score     = object$samples.p.score,
           p.weights   = getPWeights(estimand, trt, if (!is.null(pop)) pop / sum(pop),
                                     if (!is.null(samples.p.score)) samples.p.score else p.score, fitPars$p.scoreBounds)))
  
  if (is.null(result)) return(NULL)
  
  if (combineChains) result <- combineChains(result, n.chains)
  
  ## a score without samples (method.trt = "glm") gives weights with no draws
  subset <- rep_len(TRUE, if (is.null(dim(result))) length(result) else dim(result)[length(dim(result))])
  if (sample == "inferential") {
    if (object$estimand == "att") subset <- object$trt > 0
    else if (object$estimand == "atc") subset <- object$trt <= 0
  }
  
  if (is.null(dim(result)))
    result[subset]
  else if (length(dim(result)) > 2L)
    result[,,subset]
  else
    result[,subset]
}

## Posterior predictive draws at the expected values 'ev', the observation
## margin last. A continuous draw's noise is sigma times 'scale', one factor
## per observation (see noiseScale), NA giving NA draws; NULL is sigma for
## every observation, as at new rows.
sampleFromPPD <- function(object, ev, scale = NULL)
{
  if (responseIsBinary(object)) {
    if (length(dim(ev)) > 2L) {
      result <- array(rbinom(length(ev), 1L, ev), dim(ev), dimnames = dimnames(ev))
    } else {
      result <- matrix(rbinom(length(ev), 1L, ev), nrow(ev), ncol(ev), dimnames = list(rownames(ev), colnames(ev)))
    }
  } else {
    n.obs <- dim(ev)[length(dim(ev))]
    sigma <- extractBartcFit(object, "sigma", combineChains = FALSE)
    n.draws <- length(sigma)
    sigma <- rep_len(sigma, n.obs * n.draws)
    ## rnorm(n, 0, s) is s times standard normal draws, so an unscaled draw is
    ## the one drawn before the scale was read
    epsilon <- rnorm(length(sigma)) * sigma
    if (!is.null(scale)) epsilon <- epsilon * rep(scale, each = n.draws)
    dim(epsilon) <- dim(ev)
    result <- ev + epsilon 
  }
  
  result
}

# extract <- function(object, ...) UseMethod("extract")


refit.bartcFit <- function(object, newresp = NULL,
                           commonSup.rule = c("none", "sd", "chisq"),
                           commonSup.cut  = c(NA_real_, 1, 0.05), ...)
{
  issueWarningForUnknownArguments()

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
  responseIsStan4Bart <- inherits(object$fit.rsp, "stan4bartFit")
  y.rsp   <- if (responseIsStan4Bart) object$fit.rsp$y       else object$data.rsp@y
  ## not 'weights', which with(object, ...) below would shadow
  pop <- populationWeights(object)
  
  ## the fit keeps the response of the complete rows only, while the
  ## estimators below run over every row, as in the response fits; a fit made
  ## before the population weights were stored kept theirs for the complete
  ## rows only too
  if (object$method.rsp %in% c("p.weight", "tmle") && any(object$missingRows) && length(y.rsp) < length(object$trt)) {
    if (!is.null(pop) && !populationWeightsAreStored(object))
      stop("refit of a weighted '", object$method.rsp, "' fit with missing responses is not supported, ",
           "as the fit does not keep the weights of the rows with a missing response; ",
           "call bartc again with the new common support rule instead")
    y.all <- rep_len(NA_real_, length(object$trt))
    y.all[!object$missingRows] <- y.rsp
    y.rsp <- y.all
  }
  
  group.effects <- if (!is.null(object[["group.effects"]])) object[["group.effects"]] else FALSE
  group.by <- if (!is.null(object[["group.by"]])) object[["group.by"]] else NULL
  
  if (object$method.rsp %in% c("bart", "bcf")) {
    samples.indiv.diff <- extractBartcFit(object, "icate", sample = "all", combineChains = FALSE)
    
    object$est <- with(object,
      getEstimateSamples(samples.indiv.diff, treatmentRows, pop, estimand, group.by, group.effects, commonSup.sub))
   
  
  } else if (object$method.rsp == "p.weight") {
    mu.hat.0 <- extractBartcFit(object, "mu.0", sample = "all", combineChains = FALSE)
    mu.hat.1 <- extractBartcFit(object, "mu.1", sample = "all", combineChains = FALSE)
    if (length(dim(mu.hat.0)) > 2L) {
      mu.hat.0 <- aperm(mu.hat.0, c(3L, 1L, 2L))
      mu.hat.1 <- aperm(mu.hat.1, c(3L, 1L, 2L))
    } else {
      mu.hat.0 <- t(mu.hat.0)
      mu.hat.1 <- t(mu.hat.1)
    }
    
    ## NB: the local variable cannot be named 'p.score' -- object$p.score also
    ## exists (the flat, posterior-averaged score), and with(object, ...) below
    ## would have this local shadowed by that field, silently reusing the
    ## averaged score for every posterior draw instead of the per-draw samples
    p.score.samples <- if (!is.null(object$samples.p.score)) object$samples.p.score else object$p.score
    if (!is.null(dim(p.score.samples)) && length(dim(p.score.samples)) < length(dim(mu.hat.0))) {
      # chains were collapsed. The chain count is the fit's own record, as in
      # bartc(); mu.hat.0's second extent agrees with it only while the
      # response fit kept its chains apart, which a combineChains = TRUE fit
      # does not
      n.chains  <- if (!is.null(object$n.chains)) object$n.chains else dim(mu.hat.0)[2L]
      n.samples <- dim(mu.hat.0)[3L]
      n.obs     <- dim(mu.hat.0)[1L]
      p.score.samples <- aperm(array(p.score.samples, c(n.chains, n.obs, n.samples)), c(3L, 1L, 2L))
    } else {
      if (!is.null(dim(p.score.samples)))
        p.score.samples <- if (length(dim(p.score.samples)) > 2L) aperm(p.score.samples, c(3L, 1L, 2L)) else t(p.score.samples)
    }

    if (is.null(object$group.by) || !group.effects) {
      if (any(object$commonSup.sub != TRUE)) {
        addDimsToSubset(mu.hat.0 <- mu.hat.0[commonSup.sub, drop = FALSE])
        addDimsToSubset(mu.hat.1 <- mu.hat.1[commonSup.sub, drop = FALSE])

        p.score.samples <- addDimsToSubset(p.score.samples[commonSup.sub, drop = FALSE])

        if (!is.null(pop)) pop <- pop[commonSup.sub]
      }
      stopIfNoWeight(pop, object$trt[object$commonSup.sub], object$estimand)

      object$est <- with(object, getPWeightEstimates(y.rsp[commonSup.sub], trt[commonSup.sub], pop, estimand, mu.hat.0, mu.hat.1, p.score.samples, fitPars$yBounds, fitPars$p.scoreBounds))
    } else {
      object$est <- lapply(levels(object$group.by), function(level) {
        levelRows <- object$group.by == level & object$commonSup.sub

        addDimsToSubset(mu.hat.0 <- mu.hat.0[levelRows, drop = FALSE])
        addDimsToSubset(mu.hat.1 <- mu.hat.1[levelRows, drop = FALSE])
        addDimsToSubset(p.score.samples  <- p.score.samples[levelRows, drop = FALSE])

        if (!is.null(pop)) pop <- pop[levelRows]
        if (groupHasNoWeight(pop, object$trt[levelRows], object$estimand)) return(NULL)

        with(object, getPWeightEstimates(y.rsp[levelRows], trt[levelRows], pop, estimand, mu.hat.0, mu.hat.1, p.score.samples,
                                         fitPars$yBounds, fitPars$p.scoreBounds))
      })
      object$est <- fillEmptyGroups(object$est)
      names(object$est) <- levels(object$group.by)
    }
  } else if (object$method.rsp == "tmle") {
    mu.hat.0 <- extractBartcFit(object, "mu.0", sample = "all", combineChains = FALSE)
    mu.hat.1 <- extractBartcFit(object, "mu.1", sample = "all", combineChains = FALSE)
    if (length(dim(mu.hat.0)) > 2L) {
      mu.hat.0 <- aperm(mu.hat.0, c(3L, 1L, 2L))
      mu.hat.1 <- aperm(mu.hat.1, c(3L, 1L, 2L))
    } else {
      mu.hat.0 <- t(mu.hat.0)
      mu.hat.1 <- t(mu.hat.1)
    }
    
    ## see the p.weight branch above: keep this off the name 'p.score' so it
    ## isn't shadowed by object$p.score inside with(object, ...) below
    p.score.samples <- if (!is.null(object$samples.p.score)) object$samples.p.score else object$p.score
    if (!is.null(dim(p.score.samples)) && length(dim(p.score.samples)) < length(dim(mu.hat.0))) {
      # chains were collapsed. The chain count is the fit's own record, as in
      # bartc(); mu.hat.0's second extent agrees with it only while the
      # response fit kept its chains apart, which a combineChains = TRUE fit
      # does not
      n.chains  <- if (!is.null(object$n.chains)) object$n.chains else dim(mu.hat.0)[2L]
      n.samples <- dim(mu.hat.0)[3L]
      n.obs     <- dim(mu.hat.0)[1L]
      p.score.samples <- aperm(array(p.score.samples, c(n.chains, n.obs, n.samples)), c(3L, 1L, 2L))
    } else {
      if (!is.null(dim(p.score.samples)))
        p.score.samples <- if (length(dim(p.score.samples)) > 2L) aperm(p.score.samples, c(3L, 1L, 2L)) else t(p.score.samples)
    }

    if (is.null(object$group.by) || !group.effects) {
      if (any(object$commonSup.sub != TRUE)) {
        addDimsToSubset(mu.hat.0 <- mu.hat.0[commonSup.sub, drop = FALSE])
        addDimsToSubset(mu.hat.1 <- mu.hat.1[commonSup.sub, drop = FALSE])

        addDimsToSubset(p.score.samples <- p.score.samples[commonSup.sub, drop = FALSE])

        if (!is.null(pop)) pop <- pop[commonSup.sub]
      }
      stopIfNoWeight(pop, object$trt[object$commonSup.sub], object$estimand)

      ## fitPars doesn't carry the n.threads originally used to fit; re-derive
      ## single-threaded rather than spawn a fresh parallel cluster as a side
      ## effect of what should be a lightweight commonSup.rule recompute
      object$est <- with(object, getTMLEEstimates(y.rsp[commonSup.sub], trt[commonSup.sub], pop, estimand, mu.hat.0, mu.hat.1, p.score.samples,
                                                  fitPars$yBounds, fitPars$p.scoreBounds, fitPars$depsilon, fitPars$maxIter,
                                                  n.threads = 1L))
    } else {
      object$est <- lapply(levels(object$group.by), function(level) {
        levelRows <- object$group.by == level & object$commonSup.sub

        addDimsToSubset(mu.hat.0 <- mu.hat.0[levelRows, drop = FALSE])
        addDimsToSubset(mu.hat.1 <- mu.hat.1[levelRows, drop = FALSE])
        addDimsToSubset(p.score.samples <- p.score.samples[levelRows, drop = FALSE])

        if (!is.null(pop)) pop <- pop[levelRows]
        if (groupHasNoWeight(pop, object$trt[levelRows], object$estimand)) return(NULL)

        with(object, getTMLEEstimates(y.rsp[levelRows], trt[levelRows], pop, estimand, mu.hat.0, mu.hat.1, p.score.samples,
                                      fitPars$yBounds, fitPars$p.scoreBounds, fitPars$depsilon, fitPars$maxIter,
                                      n.threads = 1L))
      })
      object$est <- fillEmptyGroups(object$est)
      names(object$est) <- levels(object$group.by)
    }

  }
  
  invisible(object)
}

refit <- function(object, newresp, ...) UseMethod("refit")

