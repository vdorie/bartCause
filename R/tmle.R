getPWeightFunction <- function(estimand, weights, indiv.diff, p.score)
{
  fnBody <- if (!is.null(weights)) {
    if (!is.null(dim(p.score))) {
      switch(estimand,
             att = quote(apply(indiv.diff * p.score * weights, 2L, sum) / apply(p.score * weights, 2L, sum)),
             atc = quote(apply(indiv.diff * (1 - p.score) * weights, 2L, sum) / apply((1 - p.score) * weights, 2L, sum)),
             ate = quote(apply(indiv.diff * weights, 2L, sum)))
    } else {
      if (!is.null(dim(indiv.diff))) {
        switch(estimand,
               att = quote(apply(indiv.diff * p.score * weights, 2L, sum) / sum(p.score * weights)),
               atc = quote(apply(indiv.diff * (1 - p.score) * weights, 2L, sum) / sum((1 - p.score) * weights)),
               ate = quote(apply(indiv.diff * weights, 2L, sum)))
      } else {
        switch(estimand,
               att = quote(sum(indiv.diff * p.score * weights) / sum(p.score * weights)),
               atc = quote(sum(indiv.diff * (1 - p.score) * weights) / sum((1 - p.score) * weights)),
               ate = quote(sum(indiv.diff * weights)))
      }
    }
  } else {
    if (!is.null(dim(indiv.diff))) {
      switch(estimand,
             att = quote(apply(indiv.diff * p.score, 2L, mean) / mean(z)),
             atc = quote(apply(indiv.diff * (1 - p.score), 2L, mean)  / mean(1 - z)),
             ate = quote(apply(indiv.diff, 2L, mean)))
    } else {
      switch(estimand,
             att = quote(mean(indiv.diff * p.score) / mean(z)),
             atc = quote(mean(indiv.diff * (1 - p.score)) / mean(1 - z)),
             ate = quote(mean(indiv.diff)))
    }
  }
  
  result <- function(z, weights, indiv.diff, p.score) NULL
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
  
  if (!is.null(weights)) {
    if (estimand == "att") {
      yhat0Body <- quote(-(p.score / (1 - p.score)) / sum(p.score * weights))
      yhat1Body <- quote(1 / sum(p.score * weights))
      p.scoreBody <- quote((indiv.diff - psi) / sum(p.score * weights))
      icBody <- quote(length(y) * weights * a.weight * (y - yhat) + z * (indiv.diff - psi) / sum(p.score * weights))
    } else if (estimand == "atc") {
      yhat0Body <- quote(1 / sum((1 - p.score) * weights))
      yhat1Body <- quote(-((1 - p.score) / p.score) / sum((1 - p.score) * weights))
      p.scoreBody <- quote((indiv.diff - psi) / sum((1 - p.score) * weights))
      icBody <- quote(length(y) * weights * a.weight * (y - yhat) + (1 - z) * (indiv.diff - psi) / sum((1 - p.score) * weights))
    } else if (estimand == "ate") {
      yhat0Body <- quote(1 - p.score / (1 - p.score))
      yhat1Body <- quote(1 - (1 - p.score) / p.score)
      p.scoreBody <- quote(indiv.diff - psi)
      icBody <- quote(length(y) * weights * a.weight * (y - yhat) + (indiv.diff - psi))
    }
    calcLossBody <- quote(-mean(weights * (y * log(yhat) + (1 - y) * log(1 - yhat) + z * log(p.score) + (1 - z) * log(1 - p.score))))
  } else {
    if (estimand == "att") {
      yhat0Body <- quote(-(p.score / (1 - p.score)) / mean(z))
      yhat1Body <- quote(1 / mean(z))
      p.scoreBody <- quote((indiv.diff - psi) / mean(z))
      icBody <- quote(a.weight * (y - yhat) + z * (indiv.diff - psi) / mean(z))
    } else if (estimand == "atc") {
      yhat0Body <- quote(1 / mean(1 - z))
      yhat1Body <- quote(-((1 - p.score) / p.score) / mean(1 - z))
      p.scoreBody <- quote((indiv.diff - psi) / mean(1 - z))
      icBody <- quote(a.weight * (y - yhat) + (1 - z) * (indiv.diff - psi) / mean(1 - z))
    } else if (estimand == "ate") {
      yhat0Body <- quote(1 - p.score / (1 - p.score))
      yhat1Body <- quote(1 - (1 - p.score) / p.score)
      p.scoreBody <- quote(indiv.diff - psi)
      icBody <- quote(a.weight * (y - yhat) + (indiv.diff - psi))
    }
    calcLossBody <- quote(-mean(y * log(yhat) + (1 - y) * log(1 - yhat) + z * log(p.score) + (1 - z) * log(1 - p.score)))
  }
  
  getYhat0Deriv  <- createFunctionWithBody(yhat0Body, z, weights, p.score)
  getYhat1Deriv  <- createFunctionWithBody(yhat1Body, z, weights, p.score)
  getPScoreDeriv <- createFunctionWithBody(p.scoreBody, z, weights, p.score, indiv.diff, psi)
  getIC          <- createFunctionWithBody(icBody, y, yhat, indiv.diff, psi, a.weight)
  calcLoss       <- createFunctionWithBody(calcLossBody, y, z, yhat, p.score, weights)
  
  namedList(getYhat0Deriv, getYhat1Deriv, getPScoreDeriv, getIC, calcLoss)
}
