\name{sep.pair-class}  
\docType{class}
\alias{sep.pair-class}
\title{Class "sep.pair" }
\description{
  Objects representing separation pairs in a graphical model, with slots
  \code{u} and \code{v}
  for the pair of vertices and slot \code{s} for their separator, which
  indicates the relation \code{u} and \code{v} are separated by \code{s}
  (on some graph \eqn{G}).
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("sep.pair", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{u}, \code{v}:}{Objects of class \code{"character"},
      containing the pair of vertices.}
    \item{\code{s}:}{Object of class \code{"character"}, containing the
      separator of \code{u} and \code{v}. }
  }
}
\section{Methods}{
  \describe{
    \item{all.equal}{\code{signature(target = "sep.pair", current = "sep.pair")}:
      see \code{\link{all.equal}}. }
    \item{show}{\code{signature(object = "sep.pair")}.}
	 }
}
\examples{
showClass("sep.pair")
}
\keyword{classes}
