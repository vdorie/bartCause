## A bartc fit keeps its posterior draws once, in fit.rsp and fit.trt. The
## response fit's expected outcomes under the treatment a row got and under the
## other one (mu.hat.obs, mu.hat.cf) and the propensity draws (samples.p.score)
## are read out of those fits whenever they are asked for. The functions here
## are the one place they are read out, for bartc() when it fits and for the
## reader of a fit that does not hold them.

## The draws of a dbarts response fit, as bartc lays them out: chains by
## samples by observations, or samples by observations when combined. Rows with
## a missing response were left out of the training data and ride the test
## set after the counterfactual rows, which is where their expected outcome
## comes from; missingRows says which rows they were.
bartResponseDraws <- function(bartFit, combineChains, missingRows)
{
  mu.hat.train <- extract(bartFit, sample = "train", combineChains = combineChains)
  mu.hat.test  <- extract(bartFit, sample = "test",  combineChains = combineChains)

  if (!combineChains) {
    ## dbarts >= 1.0-0 keeps the one-chain margin extract() used to drop
    mu.hat.train <- dropSingleChainDim(mu.hat.train, bartFit$n.chains, 2L)
    mu.hat.test  <- dropSingleChainDim(mu.hat.test,  bartFit$n.chains, 2L)
  }

  if (!any(missingRows)) {
    mu.hat.obs <- mu.hat.train
    mu.hat.cf  <- mu.hat.test
    ## the counterfactual rows are the observed ones, so they share their names
    dimnames(mu.hat.cf) <- dimnames(mu.hat.obs)
  } else {
    n     <- length(missingRows)
    n.mis <- sum(missingRows)
    # input dims are n.chains x n.samples x n.obs
    # perm to n.obs x n.chains x n.samples and then perm back
    if (length(dim(mu.hat.train)) > 2L) {
      mu.hat.train <- aperm(mu.hat.train, c(3L, 1L, 2L))
      mu.hat.test  <- aperm(mu.hat.test, c(3L, 1L, 2L))
    } else {
      mu.hat.train <- t(mu.hat.train)
      mu.hat.test  <- t(mu.hat.test)
    }
    mu.hat.obs <- array(0, c(n, dim(mu.hat.train)[-1L]))
    addDimsToSubset(mu.hat.obs[!missingRows] <- mu.hat.train)
    addDimsToSubset(mu.hat.obs[ missingRows] <- mu.hat.test[seq.int(n + 1L, n + n.mis), drop = FALSE])
    
    addDimsToSubset(mu.hat.cf <- mu.hat.test[seq_len(n)])
    
    if (length(dim(mu.hat.train)) > 2L) {
      mu.hat.obs <- aperm(mu.hat.obs, c(2L, 3L, 1L))
      mu.hat.cf  <- aperm(mu.hat.cf, c(2L, 3L, 1L))
    } else {
      mu.hat.obs <- t(mu.hat.obs)
      mu.hat.cf  <- t(mu.hat.cf)
    }
  }

  list(mu.hat.obs = mu.hat.obs, mu.hat.cf = mu.hat.cf)
}

## stan4bart has n.obs x n.samples x n.chains or n.obs x (n.samples * n.chains);
## bart has n.chains x n.samples x n.obs or (n.chains * n.samples) x n.obs
stan4bartResponseDraws <- function(bartFit, combineChains)
{
  mu.hat.train <- extract(bartFit, sample = "train", combine_chains = combineChains)
  mu.hat.test  <- extract(bartFit, sample = "test",  combine_chains = combineChains)
  
  if (combineChains) {
    mu.hat.train <- t(mu.hat.train)
    mu.hat.test  <- t(mu.hat.test)
  } else {
    mu.hat.train <- aperm(mu.hat.train, c(3L, 2L, 1L))
    mu.hat.test  <- aperm(mu.hat.test,  c(3L, 2L, 1L))
  }
  
  list(mu.hat.obs = mu.hat.train, mu.hat.cf = mu.hat.test)
}

## A propensity fit's draws, observations last
treatmentDraws <- function(bartFit, combineChains)
{
  if (inherits(bartFit, "stan4bartFit")) {
    samples <- extract(bartFit, combine_chains = combineChains)
    if (length(dim(samples)) == 3L) {
      samples <- aperm(samples, c(3L, 2L, 1L))
    } else {
      samples <- t(samples)
    }
  } else {
    samples <- extract(bartFit, combineChains = combineChains)
    ## dbarts >= 1.0-0 keeps the one-chain margin extract() used to drop
    if (!combineChains) samples <- dropSingleChainDim(samples, bartFit$n.chains, 2L)
  }
  samples
}

## What bartc() records so that a reader can lay the draws out as the fit did:
## combined chains have no chain margin, and a one-chain fit left apart has none
## either except under stan4bart, so the shape of the draws is the record.
combineChainsRecord <- function(mu.hat.obs, samples.p.score)
  c(rsp = length(dim(mu.hat.obs)) == 2L,
    trt = if (is.null(samples.p.score)) NA else length(dim(samples.p.score)) == 2L)

## The elements a bartc fit does not hold: read out of the fits when asked for
computedElements <- c("mu.hat.obs", "mu.hat.cf", "samples.p.score")

## The draws, read out of the fits once. 'which' is "mu" for mu.hat.obs and
## mu.hat.cf together (they share the extract calls) and "p.score" for
## samples.p.score. An element the fit holds (one saved by an earlier version,
## or a propensity matrix that was given) is returned as held.
bartcDraws <- function(object, which = c("mu", "p.score"))
{
  which <- match.arg(which, several.ok = TRUE)
  result <- list()
  if ("mu" %in% which) {
    if (all(c("mu.hat.obs", "mu.hat.cf") %in% names(unclass(object)))) {
      result$mu.hat.obs <- .subset2(object, "mu.hat.obs")
      result$mu.hat.cf  <- .subset2(object, "mu.hat.cf")
    } else {
      fit <- .subset2(object, "fit.rsp")
      combine <- .subset2(object, "combineChains")[["rsp"]]
      draws <-
        if (inherits(fit, "bartBCF")) list(mu.hat.obs = fit[["mu.hat.obs"]], mu.hat.cf = fit[["mu.hat.cf"]])
        else if (inherits(fit, "stan4bartFit")) stan4bartResponseDraws(fit, combine)
        else bartResponseDraws(fit, combine, .subset2(object, "missingRows"))
      result$mu.hat.obs <- draws$mu.hat.obs
      result$mu.hat.cf  <- draws$mu.hat.cf
    }
  }
  if ("p.score" %in% which) {
    result["samples.p.score"] <-
      if ("samples.p.score" %in% names(unclass(object))) list(.subset2(object, "samples.p.score"))
      else if (is.null(.subset2(object, "fit.trt")) || !inherits(.subset2(object, "fit.trt"), c("bart", "stan4bartFit"))) list(NULL)
      else list(treatmentDraws(.subset2(object, "fit.trt"), .subset2(object, "combineChains")[["trt"]]))
  }
  result
}

## One of the three elements, as a plain list would answer for a name it holds
bartcElement <- function(x, name)
{
  if (name %in% computedElements && !(name %in% names(unclass(x)))) {
    if (name == "samples.p.score") bartcDraws(x, "p.score")[[name]]
    else bartcDraws(x, "mu")[[name]]
  } else {
    .subset2(x, name)
  }
}

`$.bartcFit` <- function(x, name)
{
  if (name %in% computedElements) bartcElement(x, name)
  else .subset2(x, name, exact = FALSE)
}

`[[.bartcFit` <- function(x, i, ...)
{
  if (is.character(i) && length(i) == 1L && i %in% computedElements) bartcElement(x, i)
  else NextMethod()
}
