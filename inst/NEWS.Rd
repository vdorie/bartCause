\name{NEWS}
\title{bartCause News}
\encoding{UTF-8}
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
