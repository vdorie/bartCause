\name{bcf}
\alias{bcf}
\alias{bartBCF}
\alias{print.bartBCF}
\alias{fitted.bartBCF}
\alias{extract.bartBCF}
\alias{predict.bartBCF}
\alias{residuals.bartBCF}
\title{Bayesian Causal Forests}
\description{
  Fits a Bayesian causal forest: a prognostic surface \eqn{\mu(x)}{mu(x)} and a
  treatment surface \eqn{\tau(x)}{tau(x)}, each its own ensemble of trees,
  combined as \eqn{a\mu(x_i) + b_{z_i}\tau(x_i)}{a mu(x_i) + b_z_i tau(x_i)}
  with amplitudes drawn alongside them. Both counterfactual surfaces, both
  forests' per-draw fitted values, the amplitudes and each forest's per-draw
  variable counts are reported, and no test matrix is involved.
}
\usage{
bcf(formula, data, subset, weights, offset,
    treatment,
    moderators = NULL,
    p.score = NULL,
    n.trees = 200L, base = 0.95, power = 2.0,
    n.trees.treatment = 50L,
    treatment.base = 0.25, treatment.power = 3,
    sd.control = NULL, sd.moderate = 1,
    b.prior.variance = 0.5,
    update.a = TRUE, update.b = TRUE,
    mu.interactions = NULL, tau.interactions = NULL,
    tau.blocks = NULL,
    n.samples = 500L, n.burn = 500L, n.chains = 4L,
    n.threads = dbarts::guessNumCores(), combineChains = TRUE,
    keepSampler = FALSE, verbose = TRUE, seed = NA_integer_, \dots)

\method{print}{bartBCF}(x, digits = max(3L, getOption("digits") - 3L), \dots)

\method{fitted}{bartBCF}(object,
       type = c("mu.obs", "mu.cf", "mu.1", "mu.0", "icate", "mu", "tau",
                "glue", "sigma", "varcount"),
       \dots)

\method{extract}{bartBCF}(object,
        type = c("mu.obs", "mu.cf", "mu.1", "mu.0", "icate", "mu", "tau",
                 "glue", "sigma", "varcount"),
        forest = NULL,
        combineChains = object[["combineChains"]],
        \dots)

\method{predict}{bartBCF}(object, newdata, type = c("mu", "tau"), \dots)

\method{residuals}{bartBCF}(object, \dots)
}
\arguments{
  \item{formula}{
    A formula for the response and the confounders, as in
    \code{\link[dbarts]{dbartsData}}, or a matrix of confounders. The treatment
    is added to the design by \code{bcf} itself and must not appear in the
    formula.
  }
  \item{data}{
    A data frame or list in which to look up the formula's variables, or, when
    \code{formula} is a matrix, the response vector.
  }
  \item{subset}{
    An optional vector of rows to fit. Logical masks, positive and negative
    integer indices, and character row names are all accepted; the treatment
    basis is built from the full-length treatment and restricted to the same
    rows.
  }
  \item{weights}{
    An optional vector of observation weights.
  }
  \item{offset}{
    An optional vector added to the linear predictor. It cancels out of the
    counterfactual difference, so it need not be supplied to obtain one.
  }
  \item{treatment}{
    A binary treatment variable coded 0 and 1, or a reference to one in
    \code{data}. Both arms must be present after subsetting.
  }
  \item{moderators}{
    An optional character vector of design columns the treatment forest may
    split on. Defaults to every column except the treatment and the propensity
    score. Naming either of those is an error rather than a silent drop.
  }
  \item{p.score}{
    An optional propensity score, added to the design as a
    \emph{prognostic-forest-only} column.
  }
  \item{n.trees, base, power}{
    Number of trees and the tree structure prior for the prognostic forest.
  }
  \item{n.trees.treatment, treatment.base, treatment.power}{
    The same, for the treatment forest, whose default prior is shallower.
  }
  \item{sd.control}{
    Scale of the prognostic forest's amplitude prior, or \code{NULL} for the
    family default.
  }
  \item{sd.moderate}{
    Scale of the treatment forest's leaf prior.
  }
  \item{b.prior.variance}{
    Prior variance of the treatment forest's amplitudes \eqn{(b_0, b_1)}{(b0,
    b1)}.
  }
  \item{update.a, update.b}{
    Whether the prognostic and treatment amplitudes are drawn; when
    \code{FALSE} the corresponding amplitude is fixed at its prior center.
  }
  \item{mu.interactions, tau.interactions}{
    Optional \code{\link[dbarts]{interactions}} constraints, one per forest.
  }
  \item{tau.blocks}{
    An optional \code{\link[dbarts]{blocks}} partition for the treatment forest.
    There is no prognostic-forest equivalent; see Details.
  }
  \item{n.samples, n.burn, n.chains, n.threads}{
    Posterior samples kept per chain, burn-in iterations discarded, number of
    chains, and number of threads.
  }
  \item{combineChains}{
    The default the accessors use: when \code{TRUE} they fold the chain margin
    into the sample margin. The fit always stores per-chain arrays.
  }
  \item{keepSampler}{
    Whether to keep the \code{\link[dbarts]{dbartsSampler-class}} object on the
    result, which allows the current-state per-forest variable counts to be
    read.
  }
  \item{verbose}{
    Whether the sampler prints its progress.
  }
  \item{seed}{
    An optional integer seed for the sampler's own random number generators,
    making a fit reproducible independently of the number of threads.
  }
  \item{x}{
    An object of class \code{bartBCF}.
  }
  \item{object}{
    An object of class \code{bartBCF}.
  }
  \item{digits}{
    Number of significant digits to print.
  }
  \item{type}{
    Which quantity to return; see Details.
  }
  \item{forest}{
    For \code{type = "varcount"}, which forest's counts to return -
    \code{"mu"}, \code{"tau"}, or an index. The default returns both, as a
    named list.
  }
  \item{newdata}{
    Not presently used; see Details.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link[dbarts]{dbarts}} or
    \code{\link[dbarts]{dbartsControl}}, such as \code{family} or
    \code{n.thin}.
  }
}
\details{
  A Bayesian causal forest fits two ensembles at once. The prognostic forest
  \eqn{\mu}{mu} is an ordinary sum of trees; the treatment forest \eqn{\tau}{tau}
  reaches the model through a basis of treatment indicators, so each draw
  carries a pair of amplitudes \eqn{(b_0, b_1)}{(b0, b1)} - one per arm - and
  the fitted value is
  \eqn{a\mu(x_i) + b_{z_i}\tau(x_i)}{a mu(x_i) + b_z_i tau(x_i)}. Swapping
  \eqn{b_{z_i}}{b_z_i} for \eqn{b_{1 - z_i}}{b_(1-z_i)} gives the
  counterfactual surface directly, with no test matrix, no out-of-sample
  prediction and no second pass through the trees.

  \strong{The moderator exclusion.} The treatment enters the model through the
  treatment forest's basis, so both forests are masked out of the treatment
  column: a prognostic forest free to split on it would absorb the very effect
  the treatment forest exists to estimate. The propensity score, when supplied,
  is a \emph{prognostic-forest} column only and is masked out of the treatment
  forest. \code{moderators} overrides the treatment forest's split set, but
  naming the treatment or the propensity score is refused rather than silently
  dropped.

  \strong{Variable counts are structurally zero under those masks.} The
  treatment forest's counts are zero on the treatment column and on the
  propensity score column, and the prognostic forest's are zero on the treatment
  column, in every draw. A zero there is the mask, not a finding about the data.

  \strong{Binary responses.} A 0/1 response fits a probit model (or a logistic
  one under \code{family = "logistic"}). The two forests are combined on the
  latent scale and the link is applied last, to both surfaces, so
  \code{"mu.obs"} and \code{"mu.cf"} are probabilities while \code{"mu"},
  \code{"tau"}, \code{"glue"} and \code{"icate"} remain on the latent scale.
  A binary fit has no \code{sigma} or \code{first.sigma} element at all.

  \code{fitted} returns posterior means and \code{extract} the draws. The types
  are:

  \itemize{
    \item \code{"mu.obs"} - the expected value under the observed treatment
          condition
    \item \code{"mu.cf"} - the expected value under the counterfactual condition
    \item \code{"mu.0"}, \code{"mu.1"} - the expected values under control and
          under treatment
    \item \code{"icate"} - the individual conditional average treatment effect,
          \eqn{s (b_1 - b_0) \tau(x_i)}{s (b1 - b0) tau(x_i)} for the response
          scale \eqn{s}{s}. On a continuous response this equals
          \code{"mu.1"} minus \code{"mu.0"} exactly; on a binary one it is the
          latent-scale contrast.
    \item \code{"mu"}, \code{"tau"} - each forest's own fitted values, on the
          \emph{internal} scale the amplitudes multiply. (Unrelated to the
          \code{tau} channel of \code{\link[dbarts]{dbarts}}, which is a
          different quantity.)
    \item \code{"glue"} - the amplitudes \eqn{(a, b_0, b_1)}{(a, b0, b1)} each
          draw combines the forests through, named \code{"a"}, \code{"b.0"} and
          \code{"b.1"}
    \item \code{"sigma"} - the residual standard deviation; an error on a binary
          fit
    \item \code{"varcount"} - per-forest, per-draw predictor split counts, a
          named list unless \code{forest} names one
  }

  \strong{What a Bayesian causal forest cannot do here.} \code{predict} is
  defined and refuses: out-of-sample \eqn{\mu(x)}{mu(x)} and \eqn{\tau(x)}{tau(x)}
  need per-forest saved-tree replay, which \code{dbarts} does not expose, and
  the blended test surface they would be combined through is refused by the
  sampler at creation. For the same reason a missing response cannot be
  recovered, so one is refused rather than dropped. Crossvalidation is refused:
  the crossvalidation engine has no multi-forest form. A varying intercept
  (\code{group.by} with \code{use.ranef = TRUE} through \code{\link{bartc}}) and
  semiparametric \code{parametric} equations are both fit by \code{stan4bart},
  which has no multi-forest form, so they are refused.
  Finally, a prognostic-forest \code{blocks} partition is refused: a block
  partition covers the full design and would have to include the treatment
  column the prognostic forest is masked out of. Use \code{tau.blocks}.
}
\value{
  \code{bcf} returns an object of class \code{bartBCF}, a list whose posterior
  quantities are stored per chain as \code{n.chains x n.samples x n.obs} arrays,
  losing the leading margin when \code{n.chains} is one. Items are:

  \item{\code{forests}}{named list of length two, \code{mu} and \code{tau},
                        each forest's per-draw fitted values on the internal
                        scale}
  \item{\code{glue}}{the per-draw amplitudes, \code{n.chains x n.samples x 3},
                     with \code{"a"}, \code{"b.0"} and \code{"b.1"} named on the
                     last margin}
  \item{\code{mu.hat.obs}, \code{mu.hat.cf}}{the surfaces under the observed and
                     counterfactual treatment conditions}
  \item{\code{varcount}}{named list of length two, keyed as \code{forests} is,
                     each \code{n.chains x n.samples x n.predictors} with the
                     predictor names on the last margin}
  \item{\code{sigma}, \code{first.sigma}}{the residual standard deviation for
                     the kept and the burn-in draws, each an
                     \code{n.chains x n.samples} matrix at every chain count
                     (so \code{extract(type = "sigma")} of a one-chain fit is a
                     \code{1 x n.samples} matrix, where a \code{bart} fit would
                     give a vector). \code{first.sigma} has zero columns when
                     \code{n.burn} is 0, and both are absent on a binary fit}
  \item{\code{y}, \code{trt}}{the training response and treatment}
  \item{\code{name.trt}, \code{name.p.score}}{the resolved design column names;
                     \code{name.p.score} is \code{NULL} when no score was given}
  \item{\code{family}}{the resolved response family}
  \item{\code{response.scale}, \code{response.shift}}{the linear map from the
                     internal scale to the response scale}
  \item{\code{data}}{the \code{\link[dbarts]{dbartsData}} object fit}
  \item{\code{n.trees}}{tree counts, one per forest, named as \code{forests} is}
  \item{\code{n.chains}, \code{n.samples}, \code{n.burn}}{chains run, posterior
                     draws kept per chain, and burn-in iterations discarded}
  \item{\code{n.obs}}{observations fit, after any \code{subset}}
  \item{\code{combineChains}}{the value passed to \code{bcf}, which the
                     accessors take as their default}
  \item{\code{call}}{the matched call, or the response fitter's when the fit
                     came from \code{bartc}}
  \item{\code{fit}}{the sampler, only under \code{keepSampler = TRUE}}
}
\author{
Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
\code{\link{bartc}}, \code{\link[dbarts]{dbarts}}
}
\examples{
## fit a simple linear model
n <- 100L
beta.z <- c(.75, -0.5,  0.25)
beta.y <- c(.5,   1.0, -1.5)

set.seed(725)
x <- matrix(rnorm(3 * n), n, 3)
tau <- rgamma(1L, 0.25 * 16 * rgamma(1L, 1 * 32, 32), 16)

p.score <- as.vector(pnorm(x \%*\% beta.z))
z <- rbinom(n, 1, p.score)

y <- as.vector(x \%*\% beta.y) + tau * z + rnorm(n, 0, 2)

df <- data.frame(y = y, z = z, x1 = x[,1], x2 = x[,2], x3 = x[,3])

# low parameters only for example
fit <- bcf(y ~ x1 + x2 + x3, data = df, treatment = z, p.score = p.score,
           n.samples = 100L, n.burn = 15L, n.chains = 2L, verbose = FALSE)
fit

# the individual conditional average treatment effects
mean(fitted(fit, "icate"))

# the treatment forest never splits on the propensity score or the treatment
fitted(fit, "varcount")$tau
}
\keyword{causal}
\keyword{bart}
