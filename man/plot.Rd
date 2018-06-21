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
             xvar = "tree.1", yvar = "tree.2",
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
    Variable for use on \eqn{x} axis.  Can be one of \code{"tree.XX"}, \code{"pca.XX"}, \code{"css"}, 
    \code{"p.score"}, \code{"y"}, \code{"y0"}, \code{"y1"}, \code{"indiv.diff"}, \code{"p.score"},
    \code{"p.weights"}, the number of name of a column used to fit the response model, or a given
    vector. See below for details.
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
  are controlled by the \code{xvar} and \code{yvar} arguments respectively and can be of the form:

\itemize{
  \item{\code{"tree.XX"}}{ - Variable number \code{"XX"} as selected by variable importance in predicting
    individual treatment effects using a tree fit by \code{\link{rpart}}. \code{"XX"} must be a number
    not exceeding the number of continuous variables used to fit the response model.}
  \item{\code{"pca.XX"}}{ - \code{"XX"}th principal component axis.}
  \item{\code{"css"}}{ - The common support statistic.}
  \item{\code{"y"}}{ - Observed response variable.}
  \item{\code{"y0"}}{ - Predicted values for the response under the control as obtained by \code{\link{fitted}}.}
  \item{\code{"y1"}}{ - Predicted values for the response under the treatment  \code{\link{fitted}}.}
  \item{\code{"indiv.diff"}}{ - Individual treatment effect estimates, or \eqn{\hat{y}(1) - \hat{y}(0)}.}
  \item{\code{"p.score"}}{ - Propensity score used to fit the response model (if applicable).}
  \item{\code{"p.weights"}}{ - Weights used when taking average of individual differences for response method
    \code{"p.weight"}.}
  \item{predictor columns}{ - Given by name or number.}
  \item{supplied vector}{ - Any vector of length equal to the number of observations.}
}
}
\value{
  None, although plotting occurs as a side-effect.
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

plot_support(fit, xvar = "tree.1", yvar = "css", legend.x = NULL)
}
}
\keyword{causal}
\keyword{bart}
\keyword{plot}

