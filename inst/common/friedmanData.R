generateFriedmanData <- function(n, ranef = FALSE, causal = FALSE, binary = FALSE) {
  f <- function(x)
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 + 10 * x[,4] + 5 * x[,5]
  
  set.seed(99)
  sigma <- 1.0
  
  x <- matrix(runif(n * 10), n, 10)
  mu <- f(x)
  
  result <- list(x = x, sigma = sigma)
  
  if (ranef) {
    n.g.1 <- 5L
    n.g.2 <- 8L
    
    result <- within(result, {
      g.1 <- sample(n.g.1, n, replace = TRUE)
      
      Sigma.b.1 <- matrix(c(1.5^2, .2, .2, 1^2), 2)
      R.b <- chol(Sigma.b.1)
      b.1 <- matrix(rnorm(2 * n.g.1), n.g.1) %*% R.b
      rm(R.b)
      
      g.2 <- sample(n.g.2, n, replace = TRUE)
      
      Sigma.b.2 <- as.matrix(1.2)
      b.2 <- rnorm(n.g.2, 0, sqrt(Sigma.b.2))
      
      mu.fixef <- x[,4] * 10
      mu.bart <- mu - mu.fixef
      mu.ranef <- b.1[g.1,1] + x[,4] * b.1[g.1,2] + b.2[g.2]
      mu <- mu + mu.ranef
    })
  } else {
    mu.fixef <- x[,4] * 10
    mu.bart <- mu - mu.fixef 
  }
  
  if (causal) {
    result <- within(result, {
      tau <- 5
      z <- rbinom(n, 1, 0.2)
    })
    
    if (ranef) {
      result <- within(result, {
        mu.fixef.0 <- mu.fixef
        mu.fixef.1 <- mu.fixef.0 + tau
        mu.bart.0  <- mu.bart.1  <- mu.bart
        mu.ranef.0 <- mu.ranef.1 <- mu.ranef
        
        mu.0 <- mu.bart.0 + mu.fixef.0 + mu.ranef.0
        mu.1 <- mu.bart.1 + mu.fixef.1 + mu.ranef.1
      })
    } else {
      result <- within(result, {
        mu.fixef.0 <- mu.fixef
        mu.fixef.1 <- mu.fixef.0 + tau
        mu.bart.0  <- mu.bart.1  <- mu.bart
        
        mu.0 <- mu.bart.0 + mu.fixef.0
        mu.1 <- mu.bart.1 + mu.fixef.1
      })
    }
    
    if (binary) {
      result <- within(result, {
        loc   <- mean(c(mu.0, mu.1))
        scale <- sd(c(mu.0, mu.1)) / qnorm(0.15)
        mu.0 <- (mu.0 - loc) / scale
        mu.1 <- (mu.1 - loc) / scale
        
        mu.fixef.0 <- (mu.fixef.0 - loc) / scale
        mu.fixef.1 <- (mu.fixef.1 - loc) / scale
        mu.bart.0 <- mu.bart.0 / scale
        mu.bart.1 <- mu.bart.1 / scale
        if (ranef) {
          mu.ranef.0 <- mu.ranef.0 / scale
          mu.ranef.1 <- mu.ranef.1 / scale
        }
        
        rm(loc, scale)
        
        y.0 <- rbinom(n, 1L, pnorm(mu.0))
        y.1 <- rbinom(n, 1L, pnorm(mu.1))
        y <- y.1 * z + y.0 * (1 - z)
      })
    } else {
      result <- within(result, {
        y.0 <- mu.0 + rnorm(n, 0, sigma)
        y.1 <- mu.1 + rnorm(n, 0, sigma)
        y <- y.1 * z + y.0 * (1 - z)
      })
    }
    
    
    result$mu <- NULL
    result$mu.fixef <- NULL
    result$mu.ranef <- NULL
  } else {
    if (binary) {
      result <- within(result, {
        loc <- mean(mu)
        scale <- sd(mu) / qnorm(0.15)
        mu <- (mu - loc) / scale
        
        mu.fixef <- (mu.fixef - loc) / scale
        mu.bart <- mu.bart / scale
        if (ranef)
          mu.ranef <- mu.ranef / scale
        
        rm(loc, scale)
        
        y <- rbinom(n, 1L, pnorm(mu))
      })
    } else {
      result <- within(result, {
        y <- mu + rnorm(n, 0, sigma)
      })
    }
  }
  
  result
}

