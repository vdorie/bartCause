## A small data frame and a bartc call on it for the tests of the forms an
## argument is written in. The arguments are given as the expressions a user
## would write, quoted, and evaluated in the caller's frame.
set.seed(19)
n <- 60L
formData <- data.frame(x1 = rnorm(n), x2 = rnorm(n), country = factor(sample(letters[1:4], n, TRUE)))
formData$x3 <- rnorm(n)
formData$z <- rbinom(n, 1, pnorm(0.5 * formData$x1))
formData$y <- formData$x1 + formData$z + rnorm(n)

small <- list(n.samples = 10L, n.burn = 10L, n.chains = 1L, n.threads = 1L, n.trees = 10L, keepTrees = TRUE,
              verbose = FALSE)
fit <- function(...) eval(as.call(c(list(quote(bartCause::bartc), quote(y), quote(z), quote(x1 + x2), data = quote(formData)),
                                    small, list(...))), parent.frame())
