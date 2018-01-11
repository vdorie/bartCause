\name{cibart}
\alias{cibart}
\title{Causal Inference using Bayesian Additive Regression Trees}
\description{
  Fits a collection of treatment and response models using the Bayesian Additive Regresssion Trees
  (BART) algorithm, producing estimates of treatment effects.
}
\usage{
cibart(response, treatment, confounders, data, subset, weights,
       method.rsp = c("bart", "pweight", "tmle"),
       method.trt = c("none", "glm", "bart", "bart.xval"),
       estimand   = c("ate", "att", "atc"),
       group.by = NULL,
       propensityScoreAsCovariate = TRUE,
       keepCall = TRUE, verbose = TRUE,
       ...)
}
\arguments{
  \item{response}{
    A vector of the continuous outcome variable, or a reference to such in the \code{data} argument.
  }
  \item{treatment}{
    A vector of the binary treatment variable, or a reference to \code{data}.
  }
  \item{confounders}{
    A matrix or data frame of covariates to be used in estimating the treatment and response model.
    Can also the right-hand-side of a formula (e.g. \code{x1 + x2 + ...}). The \code{data} argument
    will be searched when supplied
  }
  \item{data}{
    An optional data frame or named list containing the \code{response}, \code{treatment}, and
    \code{confounders}.
  }
  \item{subset}{
    An optional vector using to subset the data. Can refer to \code{data} if provided.
  }
  \item{weights}{
    An optional vector of population weights used in model fitting and estimating the treatment
    effect. Can refer to \code{data} if provided.
  }
  \item{method.rsp}{
    A character string specifying which method to use when fitting the response surface and estimating
    the treatment effect. Options are: \code{"bart"} - fit the response surface with BART and take
    the average of the individual treatment effect estimates, \code{"pweight"} - fit the resposne surface
    with BART but compute the treatment effect estimate by using a propensity score weighted sum of
    individual effects, and \code{"tmle"} - as above, but further adjust the individual estimates using
    the Targetted Minimum Loss based Estimation (TMLE) adjustment.
  }
  \item{method.trt}{
    A character string specifying which method to use when fitting the treatment assingment mechanism,
    or a vector/matrix of propensity scores. Character string options are: \code{"none"} - do no
    propensity score estimation, \code{"glm"} - fit a generalized linear model with a binomial response
    and all confounders added linearly, \code{"bart"} - fit BART directly to the treatment variable,
    and \code{"bart.xval"} - crossvalidate the treatment fit with different values of the end-node prior
    variance. Cannot be \code{"none"} if the response model requires propensity scores. When supplied as
    a matrix, it should be of dimensions equal to the number of observations times the number of samples
    used in any response model.
  }
  \item{estimand}{
    A character string specifying which causal effect to target. Options are \code{"ate"} - average
    treatment effect, \code{"att"} - average treatment effect on the treated, and \code{"atc"} - average
    treatment effect on the controls.
  }
  \item{group.by}{
    An optional factor that, when present, causes the treatment effect estimate to be calculated within
    each group.
  }
  \item{propensityScoreAsCovariate}{
    A logical such that when \code{TRUE}, the propensity score is added to the response model as a
    covariate.
  }
  \item{keepCall}{
    A logical such that when \code{FALSE}, the call to \code{cibart} is not kept. This can reduce the 
    amount of information printed by \code{\link{summary}} when passing in data as literals.
  }
  \item{verbose}{
    A logical that when \code{TRUE} prints information as the model is fit.
  }
  \item{\dots}{
    Further arguments to the treatment and response model fitting algorithms. See \code{\link{glm}}
    and \code{\link{bart2}} for reference. For \code{method.trt == "bart.xval"}, \code{n.burn} and
    \code{n.samples} can be lists supplying two sets of arguments, the first for the crossvalidation
    method (\code{\link{xbart}}) and the second for the fit itself.
  }
}
\details{
  cibart words words
}
\value{
  \code{cibart} returns an object of class \code{cibartFit}.
}
\references{
Chipman, H., George, E., and McCulloch, R. (2009)
   BART: Bayesian Additive Regression Trees.

Chipman, H., George, E., and McCulloch R. (2006)
   Bayesian Ensemble Learning. 
   Advances in Neural Information Processing Systems 19,
   Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

both of the above at:
\url{http://www.rob-mcculloch.org}

JLH Stat Sci
}
\author{
Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
\code{\link{pdbart}}
}
\examples{
## fit a simple linear model
n <- 100L
beta.z <- c(.75, -0.5,  0.25)
beta.y <- c(.5,   1.0, -1.5)
  
set.seed(725)
x <- matrix(rnorm(3 * n), n, 3)
tau <- rgamma(1L, 0.25 * 16 * rgamma(1L, 1 * 32, 32), 16)

p.score <- pnorm(x %*% beta.z)
z <- rbinom(n, 1, p.score)

mu.0 <- x %*% beta.y
mu.1 <- x %*% beta.y + tau

y <- mu.0 * (1 - z) + mu.1 * z + rnorm(n, 0, 2)

fit <- cibart(y, z, x, n.samples = 100L)
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
