\name{freq.tb-class} 
\docType{class}
\alias{freq.tb-class}

\title{Class "freq.tb" }
\description{
  Objects representing frequency counts of all configurations of
  discrete variables existing in a data set.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("freq.tb", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{table}:}{ an \eqn{m} by \eqn{p+1} matrix, where in each
      row, the first \eqn{p} integers give a configuration and the
      \eqn{p+1}-th element records the frequency count. }
    \item{\code{levels}:}{ a \eqn{p}-vector which gives the total
      levels for each variable. }
  }
}
\author{ Zongming Ma and Xiangrui Meng}
\section{Warning }{The user is responsible for giving the correct \code{levels}!}
\examples{
showClass("freq.tb")
}
\keyword{classes}
