\name{NEWS}
\title{bartCause News}
\encoding{UTF-8}
\section{CHANGES IN VERSION 1.1-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New \code{bcf} function fits a Bayesian causal forest, with
            separate prognostic and treatment ensembles. It reports both
            counterfactual surfaces and each forest's fitted values without a
            test matrix, and returns an object of class \code{bartBCF} with
            \code{print}, \code{fitted}, \code{extract}, \code{residuals} and
            \code{predict} methods.
      \item \code{bartc} accepts \code{method.rsp = "bcf"}, fitting the response
            surface with a Bayesian causal forest and estimating treatment
            effects as for \code{"bart"}. \code{crossvalidate} is not
            available with it.
      \item \code{group.by} can be supplied together with \code{parametric};
            the varying intercept is added to the parametric equation.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Requires \pkg{dbarts} >= 1.0-0. Its sampler differs from earlier
            versions, so fits, and results that depend on them, will not
            match those from 1.0-10 for the same seed.
      \item \code{group.by} with \code{use.ranef = TRUE} is now fit by
            \pkg{stan4bart}, which adds a \code{(1 | group.by)} term to the
            model, for the response surface and for the BART propensity model
            alike; \pkg{dbarts} no longer fits grouped random intercepts.
            \pkg{stan4bart} must be installed for this. The \pkg{dbarts}
            sampling arguments (\code{n.samples}, \code{n.burn},
            \code{n.chains}, \code{n.thin}, \code{n.threads}, \code{n.trees},
            \code{k}, \code{power}, \code{base}) are translated to their
            \pkg{stan4bart} equivalents rather than left at
            \pkg{stan4bart}'s defaults. \code{use.ranef = FALSE} is unchanged.
      \item Semiparametric response models, those with \code{parametric} or
            with \code{group.by} and \code{use.ranef = TRUE}, now refuse a
            response with missing values. In 1.0-10, \code{group.by} with
            \code{use.ranef = TRUE} fit such data and estimated counterfactuals
            for the missing rows; drop the incomplete rows, or use
            \code{use.ranef = FALSE}, which enters the grouping factor as a
            fixed effect in a \pkg{dbarts} fit and still accepts them.
      \item With \code{weights}, \code{method.rsp = "tmle"} now uses the
            \pkg{tmle} package, version 2.0.0 or later, which is required for
            this; a built-in estimator was used in 1.0-10. Weighted TMLE
            estimates and standard errors change. Weights that are all equal
            give the unweighted fit to rounding, and \pkg{tmle}'s own refit of the
            treatment mechanism, which it does for the effect on the treated
            and on the controls when it trims controls, uses a glm and a gam
            only, as \pkg{dbarts} does not fit a probit model under weights
            other than 0 and 1. Without \code{weights} nothing changes.
      \item BART propensity models are always probit and cannot be weighted;
            when \code{weights} are supplied they are ignored in the propensity
            fit, with a message, and enter only the treatment effect
            estimators.
      \item \code{bartc} warns, once for each model it fits, naming any
            argument that did not reach the sampler and why, where it used to
            drop the argument unseen: one that \pkg{dbarts} or
            \pkg{stan4bart} has no argument for, \code{sigdf} and
            \code{sigquant} for a binary response or a model fit by
            \pkg{stan4bart}, which draws the residual standard deviation
            itself, and \code{rngKind} and \code{rngNormalKind}, which
            \pkg{dbarts} 1.0-0 no longer has, and \code{rngSeed} given with
            \code{seed}, which governs as it did in 1.0-10.
            \code{family} and \code{seed} are the \pkg{dbarts} 1.0-0
            arguments that replace \code{sigdf}, \code{sigquant} and
            \code{rngSeed}, which are still taken; see \code{?bartc} for
            giving both. \code{power}, \code{base} and \code{split.probs}
            are \pkg{dbarts}' own arguments and are passed through.
            \code{tree.prior} is not an argument of \code{bartc}, and is
            refused with the \code{power}, \code{base} and
            \code{split.probs} it translates to.
      \item \code{bartc(crossvalidate = TRUE)} given \code{power},
            \code{base} or \code{split.probs} now cross-validates under
            them, where the development builds dropped them from the
            cross-validation; its chosen \code{k} and draws change.
      \item A fit saved by 1.0-10 or earlier cannot be used by \code{predict}
            in this version, whose \pkg{dbarts} cannot replay its saved trees;
            \code{predict} says so and asks for a refit, where it used to stop
            on a missing column. \code{print}, \code{summary}, \code{extract},
            \code{fitted} and \code{refit} still read such a fit.
      \item \code{extract} and \code{predict} with \code{combineChains = FALSE}
            now keep a leading chain dimension of length one for a
            single-chain fit, as \pkg{dbarts} does: scalar quantities are a 1
            by n.samples matrix and per-observation quantities a 1 by
            n.samples by n.obs array. Previously the dimension was dropped.
      \item Without \code{weights}, the \code{"att"} and \code{"atc"}
            \code{"p.weight"} estimate now divides by the sum of the propensity
            scores rather than by the treated (control) count, as the weighted
            estimate always did, so no weights and equal weights agree.
            Estimates move when the propensity score is a BART fit, by a
            small part of a posterior standard deviation for the pooled
            estimate and more for small groups with \code{group.effects}. For
            a \code{"glm"} one they move only through the bounding of extreme
            scores, when a common support rule drops rows, as the kept rows'
            scores then no longer sum to the kept treated count, and, for group
            estimates, when the treatment model has a varying intercept. The
            \code{p.weights} returned now sum to one in each draw, where
            1.0-10 returned values that averaged about one.
    }
  }
  \subsection{BUG-FIXES}{
    \itemize{
      \item The standard errors of weighted \code{method.rsp = "tmle"} fits,
            which came from the built-in estimator, were not converted back to
            the scale of the response, and \code{weights} together with
            \code{posteriorOfTMLE = FALSE} failed with an error.
      \item The standard errors of \code{"p.weight"} fits, and of the
            built-in \code{"tmle"} estimator used when \pkg{tmle} is not
            installed, were on a 0 to 1 scale and are now in the units of the
            response, as the estimates are. The \code{sd} and the intervals
            from \code{summary} for those fits are wider than before, by about
            40 to 50 percent in a simulation with a response of range 38.
      \item The built-in \code{method.rsp = "tmle"} estimator, used when
            \pkg{tmle} is not installed and there are no \code{weights},
            failed with \code{posteriorOfTMLE = FALSE}.
      \item The \code{"p.weight"} standard error of an \code{"att"} or
            \code{"atc"} fit divides each draw by the sum of its scores, as the
            estimate does, with and without \code{weights}. With \code{weights}
            and several draws it summed over every draw and was too small by
            the number of draws.
      \item \code{weights} are now taken for the rows of \code{subset} in
            \code{"p.weight"} and \code{"tmle"} response fits, where those of
            the first rows of the data were used, and with missing responses
            the response is paired with its own rows, weights and scores in
            both.
      \item \code{refit} of an \code{"att"} or \code{"atc"} fit failed with a
            subscript error.
      \item Multi-threaded \code{method.rsp = "tmle"} workers are seeded from
            the session's random number stream, and the cluster is closed if a
            worker fails.
      \item \code{extract} and \code{fitted} of \code{type = "p.weights"}
            no longer fail when the propensity score is a \code{"glm"} fit;
            for \code{"att"} and \code{"atc"} they no longer return \code{Inf}
            when the fit stores no propensity-score bounds; and with
            \code{weights} they weight each observation and sum to one within
            each draw, where the weights were recycled across draws.
            \code{"ate"} with \code{weights} returns the weights for every
            draw, where it errored.
      \item \code{extract(type = "sigma", combineChains = TRUE)} returns a
            plain vector, not a matrix, for a single-chain \pkg{stan4bart}
            fit.
      \item \code{refit} no longer errors on a semiparametric fit.
      \item \code{predict} takes the propensity score's column name from the
            fit. Previously a confounder whose name began with the score's
            (e.g. \code{psps} against a score \code{ps}) could be mistaken for
            it, and predicted scores were written over that confounder.
      \item \code{method.rsp = "tmle"} gave a different estimate with
            \code{combineChains = TRUE} than with chains kept apart, because
            the posterior draws were visited in a different order and
            \code{tmle::tmle} uses R's random number stream. Both now use the
            same order.
      \item \code{predict} with \code{type = "y"} failed, as it sampled the
            posterior predictive draws from an undefined object rather than
            from the predicted means.
      \item \code{refit} on a fit with \code{method.rsp = "p.weight"} or
            \code{"tmle"} errored, as it requested an invalid quantity from
            \code{extract}; with \code{group.effects = TRUE} the
            \code{"tmle"} branch also referenced undefined variables and
            mis-passed its \code{maxIter} and \code{n.threads} settings.
      \item \code{bartc} with a \code{subset} and a propensity score fit on
            that subset (with \code{data} supplied) assigned the score to the
            wrong rows of the response model's data, silently when the lengths
            divided evenly and with an error otherwise; it is now placed at the
            subset's rows.
      \item \code{bartc(..., p.score = <vector>)} matched the
            \code{p.scoreAsCovariate} argument by partial matching and failed
            with an unhelpful length error; the message now says so and points
            to \code{method.trt = <vector>}.
      \item \code{extract} on a fit dropped into R's debugger when the
            result's dimensions could not be read, rather than signalling an
            error.
      \item \code{method.rsp = "tmle"} with \code{estimand = "ate"} no longer
            has \pkg{tmle} also estimate the effects on the treated and on the
            controls for each posterior draw. That refit the treatment
            mechanism with \pkg{SuperLearner} whenever the propensity scores
            were skewed, about 3 seconds a draw, hours at the default number of
            draws; the estimates and standard errors are unchanged.
      \item \code{args.rsp} and \code{args.trt} can be a variable, or any
            expression that evaluates to a list, where only a list written out
            in the call was accepted and a variable failed with an error that
            it was not subsettable.
      \item A prior written in the call as \code{k}, such as
            \code{bartc(..., k = chi(1.25, Inf))}, failed with an error that
            \code{chi} could not be found, as R matched \code{k} to
            \code{keepCall} and \code{bartc} evaluated it, and with
            \code{method.rsp = "tmle"} or \code{combineChains} the fitters
            evaluated it again; it now reaches the \pkg{dbarts} or
            \pkg{stan4bart} fit as written. \code{k} given with another
            argument of an equal value, such as \code{n.trees}, also failed.
      \item \code{subset} written as an expression of the columns of
            \code{data}, such as \code{subset = x1 > 0}, is read within
            \code{data}, as a model frame reads it. It failed with an error
            that the column was not found with \code{group.by}.
      \item \code{subset} with a response that has missing values, fit by
            \pkg{dbarts}, failed with an error; the rows modeled are now those
            of the subset, the missing ones among them being given
            counterfactual estimates as without a subset.
      \item A response with missing values, fit by \pkg{dbarts} without
            the propensity score as a covariate (\code{method.trt = "none"}
            or \code{p.scoreAsCovariate = FALSE}), replaced the response
            column of \code{data} in the calling frame with zeros and
            missing values.
    }
  }
}
\section{CHANGES IN VERSION 1.0-10}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Generics that accept \code{...} arguments but do not use them will
            now issue warnings when called with unknown arguments. Suggestion
            thanks to Joshua Bon.
    }
  }
}
\section{CHANGES IN VERSION 1.0-8}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Compatible with TMLE >= 2.0.1.
    }
  }
}
\section{CHANGES IN VERSION 1.0-7}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Treatment and propensity score (when applicable) are now added
            to parametric model. Thanks to Joshua Bon for noticing their
            absence.
    }
  }
  \subsection{BUG-FIXES}{
    \itemize{
      \item \code{summary} now correctly computes SATE.
      \item \code{summary} now correctly computes sub-group averages.
      \item Propensity score samples correctly stored for stan4bart models.
    }
  }
}
\section{CHANGES IN VERSION 1.0-5}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{bartc} now accepts a \code{seed} argument that yields
             reproducible results, even when running with multiple threads
             and multiple chains.
      \item New \code{parametric} argument to \code{bartc} allows for arbitrary
            linear and multilevel models. Requires the \code{stan4bart} package
            to be installed.
    }
  }
  \subsection{BUG-FIXES}{
    \itemize{
      \item Now works with response models in \code{"tmle"} and
            \code{"p.weight"} when \code{p.scoreAsCovariate} is \code{FALSE}.
      \item Correctly summarizes the number of samples when \code{n.chains = 1}.
      \item \code{sumamry} now correctly computes variance expansion for PATE with
            weighted observations.
      \item \code{summary} correctly drops observations from averages when
            common support rules are being used. Bug report thanks to Marcela
            Veselkova.
    }
  }
}
\section{CHANGES IN VERSION 1.0-4}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{predict} now accepts as options \code{"mu"} and 
            \code{"y"} for when the user-supplied \code{newdata}
            object already has a valid treatment assignment.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Default \code{type} argument for \code{predict} is
            now \code{"mu"}.
      \item \code{bartc} argument \code{crossvalidateBinary}
            renamed to \code{crossvalidate}. Accepts one of \code{TRUE},
            \code{FALSE}, \code{"trt"}, or \code{"rsp"}. This feature
            is experimental and discouraged from use.
      \item Confounder expressions that can evaluate successfully in the calling
            environment but yield a warning are now treated as if they threw
            an error.
      \item \code{summary} for models fit with \code{group.effects} as
            \code{TRUE} will use the level names of the grouping factor
            when reporting the sub-group estimates.
    }
  }
  \subsection{BUG-FIXES}{
    \itemize{
    	\item \code{"ite"}, \code{"sate"}, \code{"y.0"}, and
            \code{"y.1"} for generics now returns correctly.
            Report thanks to Bryan Keller.
    }
  }
}
\section{CHANGES IN VERSION 1.0-3}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Confounders can be passed in as literals when using
            \code{data} argument. This can be of assistance in
            using complicated formulae with many variables.
    }
  }
  \subsection{BUG-FIXES}{
    \itemize{
    	\item Passing check with \code{tmle} version >= 1.5.0.
      \item Fixed effects should now execute correctly.
    }
  }
}
\section{CHANGES IN VERSION 1.0}{
  \itemize{
  	\item Initial public release.
  }
}
