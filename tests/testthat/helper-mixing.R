## stan4bart warns once per fit when its chains disagree (split R-hat over 1.1).
## The fits these tests run are a handful of draws long, so it often does, and
## what they check is not the mixing. This holds back that one warning, by its
## text, and passes every other warning on; where stan4bart has no such check it
## changes nothing.
mixingWarning <- "disagree on .*; the fit is unreliable"

quietMixing <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl(mixingWarning, conditionMessage(w)))
      invokeRestart("muffleWarning")
  })
}

## test_that() with the mixing warning held back for the whole test body.
testQuietMixing <- function(desc, code) {
  code <- substitute(code)
  eval(bquote(testthat::test_that(.(desc), quietMixing(.(code)))), parent.frame())
}

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
        fit <- quietMixing(stan4bart::stan4bart(y ~ bart(x) + (1 | g), data = d, iter = 10, warmup = 3,
                                                chains = 2, verbose = -1,
                                                bart_args = list(n.trees = 5, n.threads = 1)))
        answer <<- !is.null(fit$mixing)
      }
    }
    answer
  }
})
