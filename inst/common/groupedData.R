generateGroupedData <- function() {
  n <- 100L
  beta.z <- c(.75, -0.5,  0.25)
  beta.y <- c(.5,   1.0, -1.5)
  
  set.seed(725)
  x <- matrix(stats::rnorm(3 * n), n, 3)
  tau <- stats::rgamma(1L, 0.25 * 16 * stats::rgamma(1L, 1 * 32, 32), 16)
  
  g <- sample(3L, n, replace = TRUE)
  u.z <- rnorm(3L, 0, 0.3)
  
  p.score <- stats::pnorm(x %*% beta.z + u.z[g])
  z <- stats::rbinom(n, 1, p.score)
  
  mu.0 <- x %*% beta.y
  mu.1 <- x %*% beta.y + tau
  
  u.y <- rnorm(3L, 0, 0.5)
  
  y <- mu.0 * (1 - z) + mu.1 * z + u.y[g] + rnorm(n, 0, 2)
  
  list(x = x, z = z, y = y, p.score = p.score, g = g)
}
testData <- generateGroupedData()
rm(generateGroupedData)
