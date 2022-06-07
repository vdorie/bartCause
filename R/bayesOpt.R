optimizeBARTCall <- function(bartCall, env, kRange = NULL)
{
  initCall <- bartCall
    
  initCall[["n.chains"]] <- 1L
  initCall[["n.samples"]] <-
    if (!is.null(initCall[["n.samples"]]))
     eval(initCall[["n.samples"]], env) %/% 5L
    else
      formals(dbarts::bart2)[["n.samples"]] %/% 5L
  initCall[["n.burn"]] <-
    if (!is.null(initCall[["n.burn"]]))
      eval(initCall[["n.burn"]], env) %/% 5L
    else
      formals(dbarts::bart2)[["n.burn"]] %/% 5L
  
  # flat prior
  initCall[["k"]] <- "chi(1, Inf)"
  
  fit.init <- eval(initCall, env)
  
  xbartCall <- initCall
  xbartCall[["k"]] <- NULL
  xbartCall[["n.reps"]] <- 1L
  xbartCall[["n.samples"]] <- NULL
  xbartCall[["n.burn"]] <- NULL
  xbartCall[["n.chains"]] <- NULL
  xbartCall[[1L]] <- quote(dbarts::xbart)
  
  argsToMove <- names(xbartCall) != "" & names(xbartCall) %not_in% names(formals(dbarts::xbart)) & names(xbartCall) %in% names(formals(dbarts::dbartsControl))
  if (any(argsToMove)) {
    control <- if (!is.null(xbartCall[["control"]])) xbartCall[["control"]] else quote(dbarts::dbartsControl())
    if (is.call(control)) {
      for (argName in names(xbartCall)[argsToMove]) {
        control[[argName]] <- xbartCall[[argName]]
        xbartCall[[argName]] <- NULL
      }
    } else if (inherits(control, "dbartsControl")) {
      for (argName in names(xbartCall)[argsToMove]) {
        slot(control, argName) <- eval(xbartCall[[argName]], env)
        xbartCall[[argName]] <- NULL
      }
    } else {
      stop("unrecognized control object supplied")
    }
    xbartCall[["control"]] <- control
  }
  
  wrapper <- function(k) {
    xbartCall[["k"]] <- diff(kRange) * k + min(kRange)
    res <- eval(xbartCall)
    if (!is.null(dim(res))) res[1L,] else res
  }
  evalEnv <- new.env(parent = env)
  environment(wrapper) <- evalEnv
  kRange <- range(fit.init$k)
  kRange <- 0.5 * 2.5 * diff(kRange) * c(-1, 1) + mean(kRange)
  kRange[1L] <- max(kRange[1L], 0.1)
  if (kRange[2L] < kRange[1L]) kRange[2L] <- kRange[1L] * 2
  evalEnv$kRange <- kRange
  evalEnv$xbartCall <- xbartCall
  
  bestK <- bayesOptimize(wrapper, seq(0, 1, length.out = 4))
  bestK <- diff(kRange) * bestK + min(kRange)
  bartCall[["k"]] <- bestK
  bartCall
}

getMinInRegion <- function(lh, rh, coef)
{
  # ax^3 + bx^2 + cx + d
  if (is.na(lh) || is.na(rh)) {
    if (!is.na(rh)) return(lh)
    if (!is.na(lh)) return(rh)
    return(0.5)
  }
  if (anyNA(coef)) return(if (runif(1L) < 0.5) lh else rh)
  
  f.l <- sum(lh^(0:3) * coef)
  f.r <- sum(rh^(0:3) * coef)
  a <- coef[4L]; b <- coef[3L]; c <- coef[2L]; d <- coef[1L]
  
  if (!is.na(a) != 0) {
    disc <- 4 * (b^2 - 3 * a * c)
    if (disc < 0) return(0.5 * (lh + rh))
    roots <- (-2 * b + c(-1, 1) * sqrt(disc)) / (6 * a)
    validRoots <- roots > lh & roots < rh & (6 * a * roots + 2 * b > 0)
    if (any(validRoots))
      mean(roots[validRoots])
    else
      { if (f.l < f.r) lh else rh }
  } else if (b != 0) {
    root <- -coef[1] / (2 * coef[3])
    if (root > lh && root < rh) root else { if (f.l < f.r) lh else rh }
  } else {
    { if (f.l < f.r) lh else rh }
  }
}

bayesOptimize <- function(f, x.0, n.iter = 50L, tau = 10, theta = 1, sigma.sq = .01, plot = FALSE)
{
  covFunc <- function(x, y = x) {
    outer(x, y, function(x, y) {
      m <- pmin(x + tau, y + tau)
      theta^2 * (m^3 / 3 + 0.5 * abs(x - y) * m^2)
    })
  }
  getDerivatives <- function(t, GP) {
    x <- GP$x
    result <- theta^2 * cbind(ifelse(t < x, (t + tau) * (x + tau) - 0.5 * (t + tau)^2, 0.5 * (x + tau)^2),
                              ifelse(t < x, x - t, 0),
                              ifelse(t < x, -1, 0))
    as.vector(crossprod(result, GP$v))
  }
  post.mean <- function(x.new, GP)
  #  as.vector(crossprod(solve(GP$K.l, covFunc(GP$x, x.new)), GP$K.ly))
    as.vector(crossprod(covFunc(GP$x, x.new), GP$v))
  post.var <- function(x.new, GP)
    covFunc(x.new) - crossprod(solve(GP$K.l, covFunc(GP$x, x.new)))
  post.mean.var <- function(x.new, GP) {
    L.invKxt <- solve(GP$K.l, covFunc(GP$x, x.new))
    list(mu    = as.vector(crossprod(L.invKxt, GP$K.ly)),
         Sigma = covFunc(x.new) - crossprod(L.invKxt))
  }
  getExpectedImprovement <- function(x.new, GP, eta) {
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    
    delta <- eta - post.mean(x.new, GP)
    v <- post.var(x.new, GP)[1L]
    0.5 * delta * (1 + erf(delta / sqrt(2 * v))) + sqrt(v / (2 * pi)) * exp(-0.5 * delta^2 / v)
  }
  
  d <- length(x.0)
  
  x <- x.0
  f.x <- f(x.0)
  valid_points <- !is.infinite(f.x) & !is.na(f.x)
  if (any(!valid_points)) {
    x   <-   x[valid_points]
    f.x <- f.x[valid_points]
  }
  if (length(x) <= 2L) {
    x.new <- seq(min(x.0), max(x.0), length.out = 5L) 
    x <- c(x, x.new)
    f.x <- c(f.x, f(x.new))
    valid_points <- !is.infinite(f.x) & !is.na(f.x)
    x   <-   x[valid_points]
    f.x <- f.x[valid_points]
    
    sort_order <- order(x)
    f.x <- f.x[sort_order]
    x   <- x[sort_order]
    if (length(x) <= 2L) stop("cannot find valid starting points for optimizer")
  }
  mu    <- mean(f.x)
  sigma <- sd(f.x)
  
  for (i in seq_len(n.iter)) {
    f.x.p <- (f.x - mu) / sigma
    GP <- within(list(), {
      x <- x
      K <- covFunc(x) + diag(sigma.sq, length(x))
      K.l <- t(chol(K))
      K.ly <- solve(K.l, f.x.p)
      v <- solve(t(K.l), K.ly)
    })
    
    x.uni <- unique(x)
    n.prop  <- length(x.uni) - 1L
    x.prop  <- rep(NA, n.prop)
    
    mu.prop <- rep(NA, n.prop)
    curv.prop <- rep(NA, n.prop) # curvature
    sigma.prop <- rep(NA, n.prop)
    
    for (j in seq_len(n.prop)) {
      x.j <- x.uni[j]
      f.p <- getDerivatives(x.j, GP)
      coef <- c(0.5 * (f.p[2L] - f.p[3L] * x.j), f.p[3L] / 6)
      coef <- c(f.p[1L] - 3 * coef[2L] * x.j^2 - 2 * coef[1L] * x.j, coef)
      coef <- c(post.mean(x.j, GP) - sum(x.j^(1:3) * coef), coef)
      
      x.prop[j] <- getMinInRegion(x.j, x.uni[j + 1], coef)
      
      muSigma <- post.mean.var(x.prop[j], GP)
      mu.prop[j] <- muSigma$mu
      sigma.prop[j] <- sqrt(muSigma$Sigma[1L])
      
      f.p.prop <- getDerivatives(x.prop[j], GP)
      curv.prop[j] <- f.p.prop[2L] / (1 + f.p.prop[1L]^2)^1.5
    }
    
    # make sure end points can be considered
    if (x.prop[1L] != x.uni[1L]) {
      x.prop <- c(x.uni[1L], x.prop)
      n.prop <- n.prop + 1L
      
      muSigma <- post.mean.var(x.prop[1L], GP)
      mu.prop <- c(muSigma$mu, mu.prop)
      sigma.prop <- c(sqrt(muSigma$Sigma[1L]), sigma.prop)
      f.p <- getDerivatives(x.prop[1L], GP)
      curv.prop <- c(f.p[2L] / (1 + f.p[1L]^2)^1.5, curv.prop)
    }
    if (x.prop[length(x.prop)] != x.uni[length(x.uni)]) {
      x.prop <- c(x.prop, x.uni[length(x.uni)])
      n.prop <- n.prop + 1L
      
      muSigma <- post.mean.var(x.prop[n.prop], GP)
      mu.prop <- c(mu.prop, muSigma$mu)
      sigma.prop <- c(sigma.prop, sqrt(muSigma$Sigma[1L]))
      f.p <- getDerivatives(x.prop[n.prop], GP)
      curv.prop <- c(curv.prop, f.p[2L] / (1 + f.p[1L]^2)^1.5)
    }
       
    eta <- min(mu.prop)
    ei <- rep(NA, n.prop) # expected improvement
    for (j in seq_along(x.prop))
      ei[j] <- getExpectedImprovement(x.prop[j], GP, eta)
    
    # pick half the points from those with the highest expected improvement
    # 1/4 that are the current minima
    # 1/4 that are the most uncertainty
    n.prop <- min(d, n.prop)
    n.ei <- ceiling(0.25 * n.prop)
    #n.mu <- ceiling((n.prop - n.ei) / 3)
    n.mu <- 0L
    # n.curv <- ceiling(0.5 * (n.prop - n.ei - n.mu))
    n.curv <- ceiling((n.prop - n.ei - n.mu) / 3)
    n.sigma <- n.prop - n.ei - n.mu - n.curv
    
    i.new <- order(ei, decreasing = TRUE)[seq_len(n.ei)]
    i.ei <- i.new
    
    sort_order <- order(mu.prop)
    sort_order <- sort_order[!(sort_order %in% i.new)]
    i.mu <- sort_order[seq_len(n.mu)]
    i.new <- c(i.new, i.mu)
    
    sort_order <- order(abs(curv.prop), decreasing = TRUE)
    sort_order <- sort_order[!(sort_order %in% i.new)]
    i.curv <- sort_order[seq_len(n.curv)]
    i.new <- c(i.new, i.curv)
    
    sort_order <- order(sigma.prop, decreasing = TRUE)
    sort_order <- sort_order[!(sort_order %in% i.new)]
    i.sigma <- sort_order[seq_len(n.sigma)]
    i.new <- c(i.new, i.sigma)
    
    x.new <- x.prop[i.new]
    
    x <- c(x, x.new)
    f.x <- c(f.x, f(x.new))
    
    sort_order <- order(x)
    x <- x[sort_order]
    f.x <- f.x[sort_order]
    
    valid_points <- !is.infinite(f.x) & !is.na(f.x)
    if (any(!valid_points)) {
      x   <-   x[valid_points]
      f.x <- f.x[valid_points]
    }
    
    if (i <= 5L) {
      mu <- mean(f.x)
      sigma <- sd(f.x)
    }
  }
  
  f.x.p <- (f.x - mu) / sigma
  GP <- within(list(), {
    x <- x
    K <- covFunc(x) + diag(sigma.sq, length(x))
    K.l <- t(chol(K))
    K.ly <- solve(K.l, f.x.p)
    v <- solve(t(K.l), K.ly)
  })
  
  if (plot) {    
    x.plot <- seq(min(x.0), max(x.0), length.out = 101)
    temp <- post.mean.var(x.plot, GP)
    mu.hat <- temp$mu
    cov.hat <- temp$Sigma
    
    upper <- mu.hat + 1.96 * sqrt(diag(cov.hat))
    lower <- mu.hat - 1.96 * sqrt(diag(cov.hat))
    plot(NULL, type = "n", xlim = range(x), ylim = range(c(upper, lower, f.x.p)),
         xlab = "x", ylab = "f(x)")
    polygon(c(x.plot, rev(x.plot)), c(lower, rev(upper)), col = rgb(0, 0, 0, 0.2), border = NA)
    points(x, f.x.p, pch = 20)
    lines(x.plot, mu.hat, lwd = 1.5)
    
    points(x.prop, mu.prop, pch = 1)
    
    points(x.new, mu.prop[i.new], pch = 20, col = c(rep(2, length(i.ei)), rep(3, length(i.mu)), rep(4, length(i.curv)), rep(5, length(i.sigma))))
    legend("topright", c("ei", "mu", "curv", "sig"), pch = 20, col = 2:5, bty = "n")
  }
  
  x.uni <- unique(x)
  res <- x.uni[which.min(post.mean(x.uni, GP))[1L]]
  if (length(res) != 1L || anyNA(res))
    stop("error in bayesOptimize: result is unexpectedly NA or not of length 1")
  res
}
