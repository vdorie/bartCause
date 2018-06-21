\name{bartc}
\alias{bartc}
\title{Causal Inference using Bayesian Additive Regression Trees}
\description{
  Fits a collection of treatment and response models using the Bayesian Additive Regresssion Trees
  (BART) algorithm, producing estimates of treatment effects.
}
\usage{
bartc(response, treatment, confounders, data, subset, weights,
      method.rsp = c("bart", "p.weight", "tmle"),
      method.trt = c("none", "glm", "bart", "bart.xval"),
      estimand   = c("ate", "att", "atc"),
      group.by = NULL,
      commonSup.rule = c("none", "sd", "chisq"),
      commonSup.cut  = c(NA_real_, 1, 0.05),
      p.scoreAsCovariate = TRUE, use.rbart = FALSE,
      keepCall = TRUE, verbose = TRUE,
      \dots)
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
    will be searched when supplied.
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
    the average of the individual treatment effect estimates, \code{"p.weight"} - fit the resposne surface
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
  \item{commonSup.rule}{
    Rule for exclusion of observations lacking in common support. Options are \code{"none"} - no
    suppression, \code{"sd"} - exlude units whose predicted counterfactual standard deviation is extreme
    compared to the maximum standard deviation under those units' observed treatment condition, where
    extreme refers to the distribution of all standard deviations of observed treatment conditions,
    \code{"chisq"} - exclude observations according to ratio of the variance of posterior predicted
    counterfactual to the posterior variance of the observed condition, having a Chi Squared distribution
    with one degree of freedom under the null hypothesis of have equal distributions.
  }
  \item{commonSup.cut}{
    Cuttoffs for \code{commonSup.rule}. Ignored for \code{"none"}, when \code{commonSup.rule} is
    \code{"sd"}, refers to how many standard deviations of the distribution of posterior variance for
    counterfactuals an observation can be above the maximum of posterior variances for that treatment
    condition. When \code{commonSup.rule} is \code{"chisq"}, is the \eqn{p} value used for rejection
    of the hypothesis of equal variances.
  }
  \item{p.scoreAsCovariate}{
    A logical such that when \code{TRUE}, the propensity score is added to the response model as a
    covariate.
  }
  \item{use.rbart}{
    Logical specifying the use of \code{\link[=rbart_vi]{rbart}} for when \code{group.by} is supplied.
    When \code{TRUE}, the grouping variable is added as a random intercept to the response model.
  }
  \item{keepCall}{
    A logical such that when \code{FALSE}, the call to \code{bartc} is not kept. This can reduce the 
    amount of information printed by \code{\link[=bartc-generics]{summary}} when passing in data as literals.
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
  \code{bartc} represents a collection of methods that primarily use the Bayesian Additive Regression Trees
  (BART) algorithm to estimate causal treatment effects with binary treatment variables and continuous
  outcomes. This requires models to be fit to the response surface (distribution of the response as a
  function of treatment and confounders, \eqn{p(Y(1), Y(0) | X)} and optionally for treatment assignment
  mechanism (probability of receiving treatment, i.e. propensity score, \eqn{Pr(Z = 1 | X)}). The response surface
  model is used to impute counterfactuals, which may then be adjusted together with the propensity score to
  produce estimates of effects.
  
  Similar to \code{\link{lm}}, models can be specified symbolically. When the \code{data} term is present,
  it will be added to the search path for the \code{response}, \code{treatment}, and \code{confounders}
  variables. The confounders must be specified devoid of any "left hand side", as they appear in both of the
  models.
  
  \strong{Response Surface}
  
  The response surface methods included are:
  
  \itemize{
    \item \code{"bart"} - use BART to fit the response surface and produce individual estimates
          \eqn{\hat{Y}(1)_i}{Y(1)^hat_i} and \eqn{\hat{Y}(0)_i}{Y(0)^hat_i}. Treatment effect estimates are
          obtained by averaging the difference of these across the population of interest.
    \item \code{"p.weight"} - individual effects are esimated as in \code{"bart"}, but treatment effect estimates
          are obtained by using a propensity score weighted average. For the average treatment effect on the
          treated, these weights are \eqn{p(z_i | x_i) / (\sum z / n)}. For ATC, replace \eqn{z} with \eqn{1 - z}.
          For ATE, \code{"p.weight"} is equal to \code{"bart"}.
    \item \code{"tmle"} - individual effects are esimated as in \code{"bart"} and a weighted average is taken
          as in \code{"p.weight"}, however the response surface estimates and propensity scores are corrected
          by using the Targetted Minimum Loss based Estimation method.
  }
  
  \strong{Treatment Assignment}
  
  The treatment assignment models are:
  
  \itemize{
    \item \code{"none"} - no modeling is doing. Only applies when using response method \code{"bart"} and
          \code{p.scoreAsCovariate} is \code{FALSE}.
    \item \code{"glm"} - fit a binomial generalized linear model with logistic link and confounders included
          as linear terms.
    \item \code{"bart"} - fit a binary BART directly to the treatment using all the confounders.
    \item \code{"bart.xval"} - use the \code{\link{xbart}} function to perform cross validation on the node
          prior sensitivity (for that method, \code{k}) before fitting a final model.
    \item Finally, a vector or matrix of propensity scores can be supplied. Propensity score matrices should
          have a number of rows equal to the number of observations in the data and a number of columns
          equal to the number of posterior samples.
  }
  
  \strong{Common Support Rules}
  
  Common support, or that the probability of receiving all treatment conditions is non-zero within every area of
  the covariate space (\eqn{P(Z = 1 | X = x) > 0} for all \eqn{x} in the inferential sample), can be enforced by
  excluding observations with high posterior uncertainty. \code{bartc} supports two common support rules through
  \code{commonSup.rule} argument:
  
  \itemize{
    \item \code{"sd"} - observations are cut from the inferential sample if:
      \eqn{s_i^{f(1-z)} > m_z + a \times sd(s_j^{f(z)}}{s_i^f(1-z) > m_z + a * sd(s_j^f(z))},
      where \eqn{s_i^{f(1-z)}}{s_i^f(1-z)} is the posterior
      standard deviation of the predicted counterfactual for observation \eqn{i}, \eqn{s_j^f(z)} is the posterior
      standard deviation of the prediction for the observed treatment condition of objservation \eqn{j}, 
      \eqn{sd(s_j^{f(z)}}{sd(s_j^f(z))} is the empirical standard deviation of those quantities, and
      \eqn{m_z = max_j \{s_j^{f(z)}\}}{m_z = max_j s_j^f(z)} for all \eqn{j} in the same treatment group,
      i.e. \eqn{Z_j = z}. \eqn{a} is a constant to be passed in using \code{commonSup.cut} and its default is 1.
    \item \code{"chisq"} - observations are cut from the inferential sample if:
      \eqn{(s_i^{f(1-z)} / s_i^{f(z)})^2 > q_\alpha}{s_i^f(1-z) / s_i^f(z))^2 > q_\alpha}, where
      \eqn{s_i} are as above and \eqn{q_\alpha}, is the upper \eqn{\alpha} percentile of a \eqn{\chi^2}
      distribution with one degree of freedom, corresponding to a null hypothesis of equal variance. The default
      for \eqn{\alpha} is 0.05, and it is specified using the \code{commonSup.cut} parameter.
  }
  
  \strong{Missing Data}
  
  Missingness is allowed only in the response. If some response values are \code{NA}, the BART models will be
  trained just for where data are available and those values will be used to make predictions for the missing
  observations. Missing observations are not used when calculating statistics for assessing common support,
  although they may still be excluded on those grounds. Further, missing observations may not be compatible
  with response method \code{"tmle"}.
}
\value{
  \code{bartc} returns an object of class \code{artcFit}. Information about the object can be derived
  by using methods \code{\link[=bartc-generics]{summary}}, \code{\link{plot_sigma}}, \code{\link{plot_est}},
  \code{\link{plot_indiv}}, and \code{\link{plot_support}}. Numerical quantities are recovered with the
  \code{\link[=bartc-generics]{fitted}} and \code{\link{extract}} generics.
  
  Objects of class \code{bartcFit} are lists containing items:
  
  \item{\code{method.rsp}}{character string specifying the method used to fit the response surface}
  \item{\code{method.trt}}{character string specifying the method used to fit the treatment assignment
                           mechanism}
  \item{\code{estimand}}{character string specifying the targetted causal effect}
  \item{\code{fit.rsp}}{object containing the fitted response model}
  \item{\code{data.rsp}}{\code{\link{dbartsData}} object used when fitting the response model}
  \item{\code{fit.trt}}{object containing the fitted treatment model}
  \item{\code{group.by}}{optional factor vector containing the groups in which treatment effects are
                         estimated}
  \item{\code{samples.est}}{matrix or array of posterior samples of the treatment effect estimate}
  \item{\code{samples.indiv.diff}}{matrix or array of posterior samples of the individual treatment effects}
  \item{\code{p.score}}{the vector of propensity scores used as a covariate in the response model,
                        when applicable}
  \item{\code{samples.p.score}}{matrix or array of posterior samples of the propensity score, when applicable}
  \item{\code{name.trt}}{character string giving the name of the treatment variable in the data 
                         of \code{fit.rsp}}
  \item{\code{trt}}{vector of treatment assignments}
  \item{\code{call}}{how \code{bartc} was called}
  \item{\code{n.chains}}{number of independent posterior sampler chains in response model}
  \item{\code{commonSup.rule}}{common support rule used for suppressing observations}
  \item{\code{commonSup.cut}}{common support parameter used to set cutoff when suppression observations}
  \item{\code{sd.obs}}{vector of standard deviations of individual posterior predictors for observed treatment
                       conditions}
  \item{\code{sd.cf}}{vector of standard deviations of individual posterior predictors for counterfactuals}
  \item{\code{commonSup.sub}}{logical vector expressing which observations are used when estimating treatment
                              effects}
}
\references{
Chipman, H., George, E. and McCulloch R. (2006)
  Bayesian Ensemble Learning.
  \emph{Proceedings of the 19th International Conference on Neural Information Processing Systems}, 265--272.
  Cambridge, MA: MIT Press.
  \url{http://papers.nips.cc/paper/3084-bayesian-ensemble-learning.pdf}.

Hill, J. L. (2011)
  Bayesian Nonparametric Modeling for Causal Inference.
  \emph{Journal of Computational and Graphical Statistics} \bold{20(1)}, 217--240.
  Taylor & Francis.
  \url{https://doi.org/10.1198/jcgs.2010.08162}.

Hill, J. L. and Su Y. S. (2013)
  Assessing Lack of Common Support in Causal Inference Using Bayesian Nonparametrics: Implications for Evaluating the Effect of Breastfeeding on Children's Cognitive Outcomes
  \emph{The Annals of Applied Statistics} \bold{7(3)}, 1386--1420.
  \url{https://www.jstor.org/stable/23566478}.
}
\author{
Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
\code{\link{bart2}}
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

fit <- bartc(y, z, x, n.samples = 100L)
summary(fit)

## example to show refitting under the common support rule
fit2 <- refit(fit, commonSup.rule = "sd")
fit3 <- bartc(y, z, x, subset = fit2$commonSup.sub, n.samples = 100L)
}
\keyword{causal}
\keyword{bart}
\keyword{tmle}
