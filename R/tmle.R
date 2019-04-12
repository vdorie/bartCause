getPWeights <- function(estimand, z, weights, p.score, p.scoreBounds)
{
  p.score <- boundValues(p.score, p.scoreBounds)
  if (!is.null(weights)) {
    switch(estimand,
           att = p.score * weights / apply(p.score * weights, 1L, sum),
           atc = (1 - p.score) * weights / apply((1 - p.score) * weights, 2L, sum),
           ate = weights)
  } else {
    switch(estimand,
           att = p.score / mean(z),
           atc = (1 - p.score) / mean(1 - z),
           ate = array(1 / length(z), dim(p.score)))
  }
}

getPWeightFunction <- function(estimand, weights, cite, p.score)
{
  fnBody <- if (!is.null(weights)) {
    if (!is.null(dim(p.score))) {
      switch(estimand,
             att = quote(apply(cite * p.score * weights, 2L, sum) / apply(p.score * weights, 2L, sum)),
             atc = quote(apply(cite * (1 - p.score) * weights, 2L, sum) / apply((1 - p.score) * weights, 2L, sum)),
             ate = quote(apply(cite * weights, 2L, sum)))
    } else {
      if (!is.null(dim(cite))) {
        switch(estimand,
               att = quote(apply(cite * p.score * weights, 2L, sum) / sum(p.score * weights)),
               atc = quote(apply(cite * (1 - p.score) * weights, 2L, sum) / sum((1 - p.score) * weights)),
               ate = quote(apply(cite * weights, 2L, sum)))
      } else {
        switch(estimand,
               att = quote(sum(cite * p.score * weights) / sum(p.score * weights)),
               atc = quote(sum(cite * (1 - p.score) * weights) / sum((1 - p.score) * weights)),
               ate = quote(sum(cite * weights)))
      }
    }
  } else {
    if (!is.null(dim(cite))) {
      switch(estimand,
             att = quote(apply(cite * p.score, 2L, mean) / mean(z)),
             atc = quote(apply(cite * (1 - p.score), 2L, mean)  / mean(1 - z)),
             ate = quote(apply(cite, 2L, mean)))
    } else {
      switch(estimand,
             att = quote(mean(cite * p.score) / mean(z)),
             atc = quote(mean(cite * (1 - p.score)) / mean(1 - z)),
             ate = quote(mean(cite)))
    }
  }
  
  result <- function(z, weights, cite, p.score) NULL
  body(result) <- fnBody
  environment(result) <- parent.frame(1L)
  
  result
}



getTMLEFunctions <- function(estimand, weights) {
  createFunctionWithBody <- function(body, ...)
  {
    matchedCall <- match.call()
    result <- function() NULL
    base::body(result) <- body
    
    args     <- sapply(seq.int(3L, length(matchedCall)), function(i) matchedCall[[i]])
    argNames <- names(matchedCall)[seq.int(3L, length(matchedCall))]
    
    swapRows <- argNames == ""
    if (any(swapRows)) {
      argNames[swapRows]  <- sapply(args[swapRows], deparse)
      for (i in which(swapRows))
        args[[i]] <- formals(createFunctionWithBody)$body
    }
    names(args) <- argNames
    
    formals(result) <- args
    environment(result) <- parent.frame(2L)
    
    result
  }
  ## for R CMD check
  a.weight <- cite <- p.score <- psi <- x <- y <- mu.hat <- z <- NULL
  if (!is.null(weights)) {
    if (estimand == "att") {
      mu.hat.0Body <- quote(-p.score / (1 - p.score))
      mu.hat.1Body <- quote(1)
      p.scoreBody <- quote(cite - psi)
      icBody <- quote((length(y) * weights * a.weight * (y - mu.hat) + z * (cite - psi)) / sum(p.score * weights))
    } else if (estimand == "atc") {
      mu.hat.0Body <- quote(1)
      mu.hat.1Body <- quote(-(1 - p.score) / p.score)
      p.scoreBody <- quote(cite - psi)
      icBody <- quote((length(y) * weights * a.weight * (y - mu.hat) + (1 - z) * (cite - psi)) / sum((1 - p.score) * weights))
    } else if (estimand == "ate") {
      mu.hat.0Body <- quote(1 - p.score / (1 - p.score))
      mu.hat.1Body <- quote(1 - (1 - p.score) / p.score)
      p.scoreBody <- quote(cite - psi)
      icBody <- quote(length(y) * weights * a.weight * (y - mu.hat) + (cite - psi))
    }
    calcLossBody <- quote(-mean(weights * (y * log(mu.hat) + (1 - y) * log(1 - mu.hat) + z * log(p.score) + (1 - z) * log(1 - p.score))))
  } else {
    if (estimand == "att") {
      mu.hat.0Body <- quote(-p.score / (1 - p.score))
      mu.hat.1Body <- quote(1)
      p.scoreBody <- quote(cite - psi)
      icBody <- quote((a.weight * (y - mu.hat) + z * (cite - psi)) / mean(z))
    } else if (estimand == "atc") {
      mu.hat.0Body <- quote(1)
      mu.hat.1Body <- quote(-(1 - p.score) / p.score)
      p.scoreBody <- quote(cite - psi)
      icBody <- quote((a.weight * (y - mu.hat) + (1 - z) * (cite - psi)) / mean(1 - z))
    } else if (estimand == "ate") {
      mu.hat.0Body <- quote(1 - p.score / (1 - p.score))
      mu.hat.1Body <- quote(1 - (1 - p.score) / p.score)
      p.scoreBody <- quote(cite - psi)
      icBody <- quote(a.weight * (y - mu.hat) + (cite - psi))
    }
    calcLossBody <- quote(-mean(y * log(mu.hat) + (1 - y) * log(1 - mu.hat) + z * log(p.score) + (1 - z) * log(1 - p.score)))
  }
  
  mu.hat.0.deriv <- createFunctionWithBody(mu.hat.0Body, z, weights, p.score)
  mu.hat.1.deriv <- createFunctionWithBody(mu.hat.1Body, z, weights, p.score)
  p.score.deriv  <- createFunctionWithBody(p.scoreBody, z, weights, p.score, cite, psi)
  getIC          <- createFunctionWithBody(icBody, y, mu.hat, cite, psi, a.weight)
  calcLoss       <- createFunctionWithBody(calcLossBody, y, z, mu.hat, p.score, weights)
  
  namedList(mu.hat.0.deriv, mu.hat.1.deriv, p.score.deriv, getIC, calcLoss)
}
