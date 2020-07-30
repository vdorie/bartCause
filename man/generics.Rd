\name{bartc-generics}
\alias{bartc-generics}
\alias{extract}
\alias{fitted.bartcFit}
\alias{extract.bartcFit}
\alias{refit}
\alias{refit.bartcFit}
\alias{predict.bartcFit}
\title{Generic Methods for \code{bartcFit} Objects}
\description{
  Visual exploratory data analysis and model fitting diagnostics for causal inference models fit
  using the \code{\link{bartc}} function.
}
\usage{
\method{fitted}{bartcFit}(object,
       type = c("pate", "sate", "cate", "mu.obs", "mu.cf", "mu.0",
                "mu.1", "y.cf", "y.0", "y.1", "icate", "ite",
                "p.score", "p.weights"),
       sample = c("inferential", "all"),
       \dots)

extract(object, \dots)

\method{extract}{bartcFit}(object,
        type = c("pate", "sate", "cate", "mu.obs", "mu.cf", "mu.0",
                 "mu.1", "y.cf", "y.0", "y.1", "icate", "ite",
                 "p.score", "p.weights"),
        sample = c("inferential", "all"),
        combineChains = TRUE,
        \dots)

\method{predict}{bartcFit}(object, newdata,
        group.by,
        type = c("mu", "y", "mu.0", "mu.1", "y.0", "y.1", "icate", "ite",
                 "p.score"),
        combineChains = TRUE,
        \dots)

refit(object, newresp, \dots)

\method{refit}{bartcFit}(object,
      newresp = NULL,
      commonSup.rule = c("none", "sd", "chisq"),
      commonSup.cut  = c(NA_real_, 1, 0.05),
      \dots)

}
\arguments{
  \item{object}{
    Object of class \code{bartcFit}.
  }
  \item{type}{
    Which quantity to return. See details for a description of possible values.
  }
  \item{sample}{
    Return information for either the \code{"inferential"} (e.g. treated observations when the estimand is
    \code{att}) or \code{"all"} observations.
  }
  \item{combineChains}{
    If the models were fit with more than one chain, results retain the chain structure unless
    \code{combineChains} is \code{TRUE}.
  }
  \item{newresp}{
    Not presently used, but provided for compatibility with other definitions of the \code{refit} generic.
  }
  \item{newdata}{
    Data corresponding to the confounders in a \code{\link{bartc}} fit.
  }
  \item{group.by}{
    Optional grouping variable. See definition of \code{group.by} in \code{\link{bartc}}.
  }
  \item{commonSup.rule, commonSup.cut}{As in \code{\link{bartc}}.}
  \item{\dots}{
    Additional parameters passed up the generic method chain.
  }
}
\details{
  \code{fitted} returns the values that would serve as predictions for an object returned by the
  \code{\link{bartc}} function, while \code{extract} instead returns the full matrix or array of posterior
  samples. The possible options are:
  
  \itemize{
    \item \code{"pate", "sate", "cate"} - various target quantities; see \code{\link[=summary.bartcFit]{summary}}
    \item \code{"mu"} - predict only: expected value; requires user-supplied treatment variable in
                        \code{newdata}
    \item \code{"y"} - predict only: sample of the response; requires user-supplied treatment variable in
                        \code{newdata}
    \item \code{"mu.obs"} - (samples from the posterior of) the expected value under the observed
                            treatment condition, i.e.
                           \eqn{\hat{mu}_i(1) * z_i + \hat{mu}_i(0) * (1 - z_i)}
    \item \code{"mu.cf"} - the expected value under the counterfactual treatment condition, i.e.
                           \eqn{\hat{mu}_i(1) * (1 - z_i) + \hat{mu}_i(0) * z_i)}
    \item \code{"mu.0"} - the expected value under the control condition
    \item \code{"mu.1"} - the expected value under the treated condition
    \item \code{"y.cf"} - samples of the response under the the counterfactual treatment condition, i.e.
                          \eqn{\hat{y}_i(1 - z_i))}; values are obtained by adding noise to \code{mu.cf}
                          using the posterior predictive distribution
    \item \code{"y.0"} - observed responses under the control together with predicted under the treated, i.e.
                         \eqn{\hat{y}_i(1) * z_i + y(0) * (1 - z_i)}
    \item \code{"y.1"} - observed responses under the treatment together with predicted under the control, i.e.
                         \eqn{y_i(1) * z_i + \hat{y}(0) * (1 - z_i)}
    \item \code{"ite"} - (sample) individual treatment effect estimates, i.e.
                         \eqn{(y_i(z_i) - y_i(1 - z_i)) * (2z_i - 1)}; uses observed responses and posterior
                         predicted counterfactuals
    \item \code{"icate"} - individual conditional average treatment effect estimates, i.e.
                         \eqn{\hat{mu}_i(1) - \hat{mu}_i(0)}
    \item \code{"p.score"} - probability that each observation is assigned to the treatment group
    \item \code{"p.weights"} - weights assigned to each individual difference if the response method
                               is \code{"p.weight"}
  }
  
  \code{refit} exists to allow the same regressions to be used to calculate estimates under different
  common support rules. To refit those models on a subset, see the examples in \code{\link{bartc}}.
  
  \code{predict} allows the fitted model to be used to make predictions on an out-of-sample set.
  Requires model to be fit with \code{keepTrees} equal to \code{TRUE}. As \sQuote{\code{y}} values are
  all considered out of sample, the posterior predictive distribution is always used when relevant.
}
\value{
  For \code{fitted}, \code{extract}, and \code{predict}, a matrix, array, or vector depending on the
  dimensions of the result and the number of chains. For the following, when \code{n.chains} is one
  the dimension is dropped.
  
  \itemize{
    \item \code{"pate"}, \code{"sate"}, or \code{"cate"} - with \code{fitted}, a scalar; with
          \code{extract}, \code{n.chains x n.samples}
    \item \code{"p.score"} - depending on the fitting method, samples may or not be present; when samples
          are absent, a vector is returned for both functions; when present, the same as \code{"y"}.
    \item all other types - with \code{fitted}, a vector of length equal to the number of
          observations (\code{n.obs}); with \code{extract} or \code{predict}, a matrix or array of
          dimensions \code{n.chains x n.samples x n.obs}.
  }
  
  For \code{refit}, an object of class \code{bartcFit}.
}
\author{
Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
\code{\link{bartc}}
}
\examples{
## fit a simple linear model
n <- 100L
beta.z <- c(.75, -0.5,  0.25)
beta.y <- c(.5,   1.0, -1.5)
sigma <- 2

set.seed(725)
x <- matrix(rnorm(3 * n), n, 3)
tau <- rgamma(1L, 0.25 * 16 * rgamma(1L, 1 * 32, 32), 16)

p.score <- pnorm(x \%*\% beta.z)
z <- rbinom(n, 1, p.score)

mu.0 <- x \%*\% beta.y
mu.1 <- x \%*\% beta.y + tau

y <- mu.0 * (1 - z) + mu.1 * z + rnorm(n, 0, sigma)

# low parameters only for example
fit <- bartc(y, z, x, n.samples = 100L, n.burn = 15L, n.chains = 2L)

# compare fit to linear model
lm.fit <- lm(y ~ z + x)

plot(fitted(fit, type = "mu.obs"), fitted(lm.fit))

# rank order sample individual treatment effect estimates and plot
ites   <- extract(fit, type = "ite")
ite.m  <- apply(ites, 2, mean)
ite.sd <- apply(ites, 2, sd)
ite.lb <- ite.m - 2 * ite.sd
ite.ub <- ite.m + 2 * ite.sd

ite.o <- order(ite.m)

plot(NULL, type = "n",
     xlim = c(1, length(ite.m)), ylim = range(ite.lb, ite.ub),
     xlab = "effect order", ylab = "individual treatment effect")
lines(rbind(seq_along(ite.m), seq_along(ite.m), NA),
      rbind(ite.lb[ite.o], ite.ub[ite.o], NA), lwd = 0.5)
points(seq_along(ite.m), ite.m[ite.o], pch = 20)

}
\keyword{causal}
\keyword{bart}
\keyword{tmle}

