## Whether the installed stan4bart checks mixing: a fit carries the R-hats it
## checked. Probed once, on a small fit.
stan4bartChecksMixing <- local({
  answer <- NULL
  function() {
    if (is.null(answer)) {
      answer <<- FALSE
      if (requireNamespace("stan4bart", quietly = TRUE)) {
        set.seed(1)
        d <- data.frame(x = rnorm(40), g = factor(rep(1:4, 10)))
        d$y <- d$x + rnorm(40)
        fit <- stan4bart::stan4bart(y ~ bart(x) + (1 | g), data = d, iter = 10, warmup = 3,
                                    chains = 2, verbose = -1,
                                    bart_args = list(n.trees = 5, n.threads = 1))
        answer <<- !is.null(fit$mixing)
      }
    }
    answer
  }
})
