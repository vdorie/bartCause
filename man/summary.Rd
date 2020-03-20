\name{summary.bartcFit}
\alias{summary.bartcFit}
\title{Summary for \code{bartcFit} Objects}
\description{
  Summarizes \code{bartc} fits by calculating target quantities and uncertainty estimates.
}
\usage{
\method{summary}{bartcFit}(object,
                           target = c("pate", "sate", "cate"),
                           ci.style = c("norm", "quant", "hpd"), ci.level = 0.95,
                           pate.style = c("ppd", "var.exp"),
                           \dots)
}
\arguments{
  \item{object}{
    Object of class \code{bartcFit}.
  }
  \item{target}{
    Treatment effect to calculate. One of \code{"pate"} - population average treatment effect,
    \code{"sate"} - sample average treatment effect, and \code{"cate"} - conditional average
    treatment effect.
  }
  \item{ci.style}{
    Means of calculating confidence intervals (posterior credible regions). Options include
    \code{"norm"} - use a normal approximation, \code{"quant"} - the empirical quantites of the
    posterior samples, and \code{"hpd"} - region of highest posterior density.
  }
  \item{ci.level}{
    Level of confidence for intervals.
  }
  \item{pate.style}{
    For target \code{"pate"}, calculate uncertainty by using \code{"ppd"} - the posterior
    predictive distribution or \code{"var.exp"} a variance expansion formula.
  }
  \item{\dots}{ 
    Not used at moment, but present to match \code{summary} generic signature.
  }
}
\details{
  \code{summary} produces a numeric and qualitative summary of a \code{\link{bartc}} fit.
  
  \strong{Target Types}
  
  The SATE and PATE involve calculating predicted response values under different treatment conditions.
  When using \code{\link[=bartc-generics]{extract}} or \code{fitted}, these values are simulated directly
  from the posterior predictive distribution. However, since these quantities all have the same expected
  value, in order to provide consistent results \code{summary} only uses those simulations to derive
  credible intervals. Thus the estimates for CATE, SATE, and PATE all should be reported as the same
  but with increasing degrees of uncertainty.
  
  \strong{Grouped Data}
  
  If a model is fit with a supplied grouping variable and \code{group.effects = TRUE}, the estimates
  will also be reported within groups. When possible, the last line corresponds to the population. Within
  group estimates for resposne methods such as \code{"tmle"} cannot easily be extrapolated to the
  population at large - the means will combine based on the sample sizes but the uncertainty estimates
  will lack correlations.
}
\value{
  An object of class \code{bartcFit.summary} equivalent to a list with named items:
  
  \item{\code{call}}{how \code{bartc} was called}
  \item{\code{method.rsp}}{character string specifying the method used to fit the response surface}
  \item{\code{method.trt}}{character string specifying the method used to fit the treatment assignment
                           mechanism}
  \item{\code{ci.info}}{a named list with items \code{target}, \code{ci.style}, \code{ci.level}, and
                        \code{pate.style} as passed to \code{summary}}
  \item{\code{n.obs}}{total number of observations}
  \item{\code{n.samples}}{number of samples within any one chain}
  \item{\code{n.chains}}{total number of chains}
  \item{\code{commonSup.rule}}{common support rule used when fitting object to produce estimates}
  \item{\code{estimates}}{a \code{data.frame} with columns \code{"estimate"} - the target,
                          \code{"sd"} - standard deviation of posterior of estimate,
                          \code{"ci.lower"} - lower bound of credible region,
                          \code{"ci.upper"} - upper bound of credible region,
                          and optionally \code{"n.cut"} - how many observations were dropped using 
                          common support rule}
}
\author{
Vincent Dorie: \email{vdorie@gmail.com}.
}
