\name{NEWS}
\title{bartCause News}
\encoding{UTF-8}
\section{CHANGES IN VERSION 1.0-10}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New \code{bcf} function fits a Bayesian causal forest, with
            separate prognostic and treatment ensembles combined through
            per-draw amplitudes. Both counterfactual surfaces, each forest's
            per-draw fitted values and per-forest variable counts are reported
            without a test matrix. Returns an object of class \code{bartBCF},
            with \code{print}, \code{fitted}, \code{extract}, \code{residuals}
            and \code{predict} methods.
      \item \code{bartc} accepts \code{method.rsp = "bcf"}, fitting the response
            surface with a Bayesian causal forest and estimating treatment
            effects as for \code{"bart"}.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Generics that accept \code{...} arguments but do not use them will
            now issue warnings when called with unknown arguments. Suggestion
            thanks to Joshua Bon.
      \item \code{group.by} with \code{use.ranef = TRUE} is now fit by the
            \code{stan4bart} package, which adds a \code{(1 | group.by)} term
            to the model formula, for the response surface and for the BART
            propensity model alike; \code{dbarts} no longer fits grouped
            random intercepts. \code{stan4bart} is a suggested package, so the
            route errors informatively when it is absent;
            \code{use.ranef = FALSE}, which enters the grouping factor as
            another predictor, needs no additional package. The \code{dbarts}
            sampling arguments are translated onto \code{stan4bart}'s -
            \code{n.samples} and \code{n.burn} onto \code{iter} (which counts
            warmup) and \code{warmup}, \code{n.chains} onto \code{chains},
            \code{n.thin} onto \code{skip}, \code{n.threads} onto \code{cores},
            and \code{n.trees}, \code{k}, \code{power} and \code{base} into
            \code{bart_args} - so these fits sample as they were asked to
            rather than falling back on \code{stan4bart}'s own defaults.
            Arguments written in \code{stan4bart}'s names win where both are
            given. Chains run in a cluster only when there is more than one and
            the fit is long enough to repay building it, and such a fit is
            seeded off R's generator so that \code{set.seed} reproduces it.
      \item \code{group.by} can now be supplied together with
            \code{parametric}: the varying intercept is added to the parametric
            equation rather than refused.
      \item A missing response is now refused for semiparametric response
            surfaces - a \code{parametric} equation, or \code{group.by} with
            \code{use.ranef = TRUE}. The counterfactual comes from flipping
            the treatment in the fitted model frame, which holds only the
            complete cases, so a missing response has no counterfactual to
            recover.
      \item Every BART propensity model is a probit, and a weighted probit has
            no tractable latent form, so the \code{stan4bart}-backed propensity
            models are now fit unweighted with a message, as the
            \code{dbarts}-backed ones already were.
    }
  }
  \subsection{BUG-FIXES}{
    \itemize{
      \item \code{refit} now reads the response and the weights off a
            semiparametric fit instead of a \code{dbartsData} object, and no
            longer errors on one.
      \item \code{predict} now takes the propensity score's column name from
            the fit instead of re-deriving it from the response model's
            predictor names. A confounder whose name began with the score's
            own (e.g. a confounder \code{psps} against a score \code{ps})
            could be resolved as the score, in which case predicted scores
            were written over that confounder's values and the real score
            column was left as supplied.
      \item \code{extract(type = "sigma")} now returns posterior draws in
            chain-major order, matching every other extracted quantity, and
            per-chain sigma is now correctly assigned to its own chain;
            previously multi-chain fits interleaved chains, which also
            mis-paired sigma with the posterior predictive draws.
      \item \code{bcf}'s default \code{n.threads} is now capped at
            \code{n.chains}, matching \code{dbarts}; previously every
            multi-core, default-argument \code{bcf} fit (directly or via
            \code{bartc(method.rsp = "bcf")}) warned that \code{n.threads}
            exceeded \code{n.chains}.
      \item \code{bartc}'s \code{n.chains} is now read from the response
            fitter's own record instead of guessed from the dimensions of
            \code{mu.hat.obs}; \code{combineChains = TRUE} collapses those
            dimensions, which previously made \code{n.chains} read back as 1
            regardless of how many chains were actually run. The
            \code{"p.weight"} and \code{"tmle"} estimate paths now take the
            chain count from the same place rather than from those dimensions.
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
