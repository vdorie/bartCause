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
      \item BART propensity models are always probit and cannot be weighted;
            when \code{weights} are supplied they are ignored in the propensity
            fit, with a message, and enter only the treatment effect
            estimators.
      \item \code{extract} and \code{predict} with \code{combineChains = FALSE}
            now keep a leading chain dimension of length one for a
            single-chain fit, as \pkg{dbarts} does: scalar quantities are a 1
            by n.samples matrix and per-observation quantities a 1 by
            n.samples by n.obs array. Previously the dimension was dropped.
      \item Without \code{weights}, the \code{"att"} and \code{"atc"}
            \code{p.weights} are now scaled so that their sum with the
            individual effects is the \code{"p.weight"} estimate, as the
            weighted and \code{"ate"} ones already were. The reported values
            are smaller by a factor of the sample size than in 1.0-10.
    }
  }
  \subsection{BUG-FIXES}{
    \itemize{
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
