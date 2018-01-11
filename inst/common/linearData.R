generateLinearData <- function() {
  n <- 100L
  beta.z <- c(.75, -0.5,  0.25)
  beta.y <- c(.5,   1.0, -1.5)
  
  set.seed(725)
  x <- matrix(stats::rnorm(3 * n), n, 3)
  tau <- stats::rgamma(1L, 0.25 * 16 * stats::rgamma(1L, 1 * 32, 32), 16)
  
  p.score <- stats::pnorm(x %*% beta.z)
  z <- stats::rbinom(n, 1, p.score)
  
  mu.0 <- x %*% beta.y
  mu.1 <- x %*% beta.y + tau
  
  y <- mu.0 * (1 - z) + mu.1 * z + rnorm(n, 0, 2)
  
  list(x = x, z = z, y = y, p.score = p.score)
}
testData <- generateLinearData()
rm(generateLinearData)
