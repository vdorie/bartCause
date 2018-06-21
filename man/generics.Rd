\name{bartc-generics}
\alias{bartc-generics}
\alias{extract}
\alias{fitted.bartcFit}
\alias{extract.bartcFit}
\alias{refit}
\alias{refit.bartcFit}
\title{Generic Methods for \code{bartcFit} Objects}
\description{
  Visual exploratory data analysis and model fitting diagnostics for causal inference models fit
  using the \code{\link{bartc}} function.
}
\usage{
\method{fitted}{bartcFit}(object,
                          value = c("est", "y", "y0", "y1", "indiv.diff", "p.score", "p.weights"),
                          sample = c("inferential", "all"),
                           \dots)

extract(object, \dots)

\method{extract}{bartcFit}(object,
                           value = c("est", "y", "y0", "y1", "indiv.diff", "p.score", "p.weights"),
                           sample = c("inferential", "all"),
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
  \item{value}{
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
  \item{commonSup.rule, commonSup.cut}{As in \code{\link{bartc}}}
  \item{\dots}{
    Additional parameters passed up the generic method chain.
  }
}
\details{
  \code{fitted} returns the values that would serve as predictions for an object returned by the
  \code{\link{bartc}} function, while \code{extract} instead returns the full matrix or array of posterior
  samples. The possible options are:
  
  \itemize{
    \item \code{"est"} - the estimate itself, e.g. ATE
    \item \code{"y"} - predictions under the observed treatment condition, i.e.
                       \eqn{\hat{y}_i(1) * z_i + \hat{y}_i(0) * (1 - z_i)}.
    \item \code{"y0"} - predictions for all observations under the control
    \item \code{"y1"} - predictions for all observations under the treatment
    \item \code{"indiv.diff"} - for all observations the individual treatment effect estimates, i.e.
                                \eqn{\hat{y}_i(1) \hat{y}_i(0)}.
    \item \code{"p.score"} - probability that each observation is assigned to the treatment group
    \item \code{"p.weights"} - weights assigned to each individual difference if the response method
                               is \code{"p.weight"}
  }
  
  \code{refit} exists to allow the same regressions to be used to calculate estimates under different
  common support rules. To refit those models on a subset, see the examples in \code{\link{bartc}}.
}
\value{
  For \code{fitted} and \code{extract}, a matrix, array, or vector depending on the dimensions of the
  result and the number of chains. For the following, when \code{n.chains} is one the dimension is dropped.
  
  \itemize{
    \item \code{"est"} - when \code{fitted}, a scalar; when \code{extract}, \code{n.samples x n.chains}
    \item \code{"y"}, \code{"y0"}, \code{"y1"}, \code{"indiv.diff"}, \code{"p.weights"} - when
          \code{fitted}, a vector of length equal to the number of observations (\code{n.obs}); when
          \code{extract}, a matrix or array of dimensions \code{n.obs x n.samples x n.chains}.
    \item \code{"p.score"} - depending on the fitting method, samples may or not be present. When samples
          are absent, a vector is returned for both functions. When present, the same as \code{"y"}.
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
}
\keyword{causal}
\keyword{bart}
\keyword{tmle}

