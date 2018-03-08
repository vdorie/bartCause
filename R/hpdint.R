getHPDInterval <- function(x, alpha)
{
  x <- sort(x)
  lh <- 1L
  rh <- ceiling(length(x) * alpha)
  
  actualAlpha <- (rh - lh + 1L) / length(x)
  
  t <- (actualAlpha - alpha) / 2
  
  x.lh <- (1 - t) * x[lh] + t * x[lh + 1L]
  x.rh <- (1 - t) * x[rh] + t * x[rh - 1L]
  dist <- x.rh - x.lh 
  
  x.lh.min <- x.lh
  x.rh.min <- x.rh
  dist.min <- dist
  
  lh <- lh + 1L
  rh <- rh + 1L
  
  while (rh <= length(x)) {
    x.lh <- (1 - t) * x[lh] + t * x[lh + 1L]
    x.rh <- (1 - t) * x[rh] + t * x[rh - 1L]
    dist <- x.rh - x.lh
    
    if (dist < dist.min) {
      x.lh.min <- x.lh
      x.rh.min <- x.rh
      dist.min <- dist
    }
    
    lh <- lh + 1L
    rh <- rh + 1L
  }
  
  c(x.lh.min, x.rh.min)
}

