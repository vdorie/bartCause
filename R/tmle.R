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

getPWeightFunction <- function(estimand, weights, icate, p.score)
{
  fnBody <- if (!is.null(weights)) {
    if (!is.null(dim(p.score))) {
      switch(estimand,
             att = quote(apply(icate * p.score * weights, 2L, sum) / apply(p.score * weights, 2L, sum)),
             atc = quote(apply(icate * (1 - p.score) * weights, 2L, sum) / apply((1 - p.score) * weights, 2L, sum)),
             ate = quote(apply(icate * weights, 2L, sum)))
    } else {
      if (!is.null(dim(icate))) {
        switch(estimand,
               att = quote(apply(icate * p.score * weights, 2L, sum) / sum(p.score * weights)),
               atc = quote(apply(icate * (1 - p.score) * weights, 2L, sum) / sum((1 - p.score) * weights)),
               ate = quote(apply(icate * weights, 2L, sum)))
      } else {
        switch(estimand,
               att = quote(sum(icate * p.score * weights) / sum(p.score * weights)),
               atc = quote(sum(icate * (1 - p.score) * weights) / sum((1 - p.score) * weights)),
               ate = quote(sum(icate * weights)))
      }
    }
  } else {
    if (!is.null(dim(icate))) {
      switch(estimand,
             att = quote(apply(icate * p.score, 2L, mean) / mean(z)),
             atc = quote(apply(icate * (1 - p.score), 2L, mean)  / mean(1 - z)),
             ate = quote(apply(icate, 2L, mean)))
    } else {
      switch(estimand,
             att = quote(mean(icate * p.score) / mean(z)),
             atc = quote(mean(icate * (1 - p.score)) / mean(1 - z)),
             ate = quote(mean(icate)))
    }
  }
  
  result <- function(z, weights, icate, p.score) NULL
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
  a.weight <- icate <- p.score <- psi <- x <- y <- mu.hat <- z <- NULL
  if (!is.null(weights)) {
    if (estimand == "att") {
      mu.hat.0Body <- quote(-p.score / (1 - p.score))
      mu.hat.1Body <- quote(1)
      p.scoreBody <- quote(icate - psi)
      icBody <- quote((length(y) * weights * a.weight * (y - mu.hat) + z * (icate - psi)) / sum(p.score * weights))
    } else if (estimand == "atc") {
      mu.hat.0Body <- quote(1)
      mu.hat.1Body <- quote(-(1 - p.score) / p.score)
      p.scoreBody <- quote(icate - psi)
      icBody <- quote((length(y) * weights * a.weight * (y - mu.hat) + (1 - z) * (icate - psi)) / sum((1 - p.score) * weights))
    } else if (estimand == "ate") {
      mu.hat.0Body <- quote(1 - p.score / (1 - p.score))
      mu.hat.1Body <- quote(1 - (1 - p.score) / p.score)
      p.scoreBody <- quote(icate - psi)
      icBody <- quote(length(y) * weights * a.weight * (y - mu.hat) + (icate - psi))
    }
    calcLossBody <- quote(-mean(weights * (y * log(mu.hat) + (1 - y) * log(1 - mu.hat) + z * log(p.score) + (1 - z) * log(1 - p.score))))
  } else {
    if (estimand == "att") {
      mu.hat.0Body <- quote(-p.score / (1 - p.score))
      mu.hat.1Body <- quote(1)
      p.scoreBody <- quote(icate - psi)
      icBody <- quote((a.weight * (y - mu.hat) + z * (icate - psi)) / mean(z))
    } else if (estimand == "atc") {
      mu.hat.0Body <- quote(1)
      mu.hat.1Body <- quote(-(1 - p.score) / p.score)
      p.scoreBody <- quote(icate - psi)
      icBody <- quote((a.weight * (y - mu.hat) + (1 - z) * (icate - psi)) / mean(1 - z))
    } else if (estimand == "ate") {
      mu.hat.0Body <- quote(1 - p.score / (1 - p.score))
      mu.hat.1Body <- quote(1 - (1 - p.score) / p.score)
      p.scoreBody <- quote(icate - psi)
      icBody <- quote(a.weight * (y - mu.hat) + (icate - psi))
    }
    calcLossBody <- quote(-mean(y * log(mu.hat) + (1 - y) * log(1 - mu.hat) + z * log(p.score) + (1 - z) * log(1 - p.score)))
  }
  
  mu.hat.0.deriv <- createFunctionWithBody(mu.hat.0Body, z, weights, p.score)
  mu.hat.1.deriv <- createFunctionWithBody(mu.hat.1Body, z, weights, p.score)
  p.score.deriv  <- createFunctionWithBody(p.scoreBody, z, weights, p.score, icate, psi)
  getIC          <- createFunctionWithBody(icBody, y, mu.hat, icate, psi, a.weight)
  calcLoss       <- createFunctionWithBody(calcLossBody, y, z, mu.hat, p.score, weights)
  
  namedList(mu.hat.0.deriv, mu.hat.1.deriv, p.score.deriv, getIC, calcLoss)
}
