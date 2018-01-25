\name{cibart-generics}
\alias{cibart-generics}
\alias{extract}
\alias{fitted.cibartFit}
\alias{extract.cibartFit}
\title{Generic Methods for \code{cibartFit} Objects}
\description{
  Visual exploratory data analysis and model fitting diagnostics for causal inference models fit
  using the \code{\link{cibart}} function.
}
\usage{
\method{fitted}{cibartFit}(object,
                           value = c("est", "y", "y0", "y1", "indiv.diff", "p.score"),
                           sample = c("inferential", "all"),
                           \dots)

extract(object, \dots)

\method{extract}{cibartFit}(object,
                            value = c("est", "y", "y0", "y1", "indiv.diff", "p.score"),
                            sample = c("inferential", "all"),
                            combineChains = TRUE,
                            \dots)
}
\arguments{
  \item{object}{
    Object of class \code{cibartFit}.
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
  \item{...}{
    Additional parameters passed up the generic method chain.
  }
}
\details{
  \code{fitted} returns the values that would serve as predictions for an object returned by the
  \code{\link{cibart}} function, while \code{extract} instead returns the full matrix or array of posterior
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
  }
}
\value{
  A matrix, array, or vector depending on the dimensions of the result and the number of chains. For the
  following, when \code{n.chains} is one the dimension is dropped.
  
  \itemize{
    \item \code{"est"} - when \code{fitted}, a scalar; when \code{extract}, \code{n.samples x n.chains}
    \item \code{"y"}, \code{"y0"}, \code{"y1"}, \code{"indiv.diff"} - when \code{fitted}, a vector of length
          equal to the number of observations (\code{n.obs}); when \code{extract}, a matrix or array of
          dimensions \code{n.obs x n.samples x n.chains}.
    \item \code{"p.score"} - depending on the fitting method, samples may or not be present. When samples
          are absent, a vector is returned for both functions. When present, the same as \code{"y"}.
  }
}
\author{
Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
\code{\link{cibart}}
}
\examples{
## fit a simple linear model
}
\keyword{causal}
\keyword{bart}
\keyword{tmle}

