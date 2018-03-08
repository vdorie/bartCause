\name{bartc-plot}
\alias{bartc-plot}
\alias{plot_sigma}
\alias{plot_est}
\alias{plot_indiv}
\alias{plot_support}
\title{Plot methods for \code{bartc}}
\description{
  Visual exploratory data analysis and model fitting diagnostics for causal inference models fit
  using the \code{\link{bartc}} function.
}
\usage{
plot_sigma(x, main = "Traceplot sigma", 
           xlab = "iteration", ylab = "sigma",
           lty = 1:x$n.chains,
           ...)

plot_est(x, main = paste("Traceplot", x$estimand),
         xlab = "iteration", ylab = x$estimand,
         lty = 1:x$n.chains,
         ...)

plot_indiv(x, main = "Histogram Individual Effects",
           xlab = "treatment effect",
           breaks = 20,
           ...)

plot_support(x, main = "Common Support Scatterplot",
             xvar = "pca.1", yvar = "pca.2",
             xlab = NULL, ylab = NULL,
             pch.trt = 21, bg.trt = "black",
             pch.ctl = pch.trt, bg.ctl = NA,
             pch.sup = pch.trt, bg.sup = NA, col.sup = "red", cex.sup = 1.5,
             legend.x = "topleft", legend.y = NULL,
             ...)
}
\arguments{
  \item{x}{
    Object of class \code{bartcFit}.
  }
  \item{main}{
    Character title of plot.
  }
  \item{xlab}{
    Character label of \eqn{x} axis. For \code{plot_support}, if \code{NULL} a default will be used.
  }
  \item{ylab}{
    Character label of \eqn{y} axis. For \code{plot_support}, if \code{NULL} a default will be used.
  }
  \item{lty}{
    For line plots (\code{plot.sigma}, \code{plot.est}), models use the values of \code{lty} to
    visually distinguish each chain.
  }
  \item{breaks}{
    Argument to code{\link{hist}}.
  }
  \item{xvar}{
    Variable for use on \eqn{x} axis. Can be of the form "pca.XX", where "XX" is the number of a
    principal axes ordered by eigenvalues, the value "css" for the common support statistic,
    "p.score" for the propensity score, the number or name of a column used to fit the response
    model, or a given vector.
  }
  \item{yvar}{
    Variable for use on the \eqn{y} axis, of the same form as \code{xvar}.
  }
  \item{pch.trt}{
    \code{pch} point value used when plotting treatment observations.
  }
  \item{bg.trt}{
    \code{bg} background value used when plotting treatment observations.
  }
  \item{pch.ctl}{
    \code{pch} point value used when plotting control observations.
  }
  \item{bg.ctl}{
    \code{bg} background value used when plotting treatment observations.
  }
  \item{pch.sup}{
    \code{pch} point value used when plotting suppressed observations.
  }
  \item{bg.sup}{
    \code{bg} background value used when plotting suppressed observations.
  }
  \item{col.sup}{
    \code{col} color value used when plotting suppressed observations.
  }
  \item{cex.sup}{
    \code{cex} size value used when plotting suppressed observations.
  }
  \item{legend.x}{
    x value passed to \code{\link{legend}}. If \code{NULL}, legend plotting is skipped.
  }
  \item{legend.y}{
    Optional y value passed to \code{legend}
  }
  \item{...}{
    Optional graphical parameters.
  }
}
\details{
  Produces various plots using objects fit by \code{\link{bartc}}. \code{plot_sigma} and \code{plot_est}
  are traditional parameter trace plots that can be used to diagnose the convergence of the posterior
  sampler. If the \code{bartc} model is fit with \code{n.chains} greater than one, by default each chain
  will be plotted with its own line type.
  
  \code{plot_indiv} produces a simple histogram of the distribution of the estimates of the individual
  effects, taken as the average of their posterior samples.
  
  \code{plot_support} is used to visualize the common support diagnostic in the form of a scatterplot.
  Points that the diagnostic excludes are outlined in red. The contents of the \eqn{x} and \eqn{y} axes
  are controlled by the \code{xvar} and \code{yvar} arguments respectively, and those arguments permit
  constructs such as \code{"pca.1"} or \code{"css"}. See the arguments section for possible options.
}
\value{
  None, although plotting occurs as a side -effect.
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
  \code{\link{bartc}}
}
\examples{
\dontrun{
fit <- bartc(...)

plot_sigma(fit)

plot_support(fit, xvar = "x1", yvar = "css", legend.x = NULL)
}
}
\keyword{causal}
\keyword{bart}
\keyword{plot}

