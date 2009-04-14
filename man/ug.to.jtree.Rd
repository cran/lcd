\name{ug.to.jtree} 
\alias{ug.to.jtree}
\title{ Junction tree construction for undirected graph }
\description{
  Constructs a junction tree for an undirected graph.
}
\usage{
ug.to.jtree(ugamat) 
}
\arguments{
  \item{ugamat}{ the adjacency matrix of the undirected graph. }
}
\details{
  This function and its helpers implement the junction tree construction
  algorithm described in detail in Section 4.3 and 4.4 of Cowell, et al
  (1999). 
}
\value{
  An object of class \code{sep.tree}.
}
\references{
  Cowell, R. G., Dawid, A. P., Lauritzen, S. L. and Spiegelhalter,
  D. J. (1999) \emph{Probabilistic Networks and Expert
    Systems}. Springer-Verlag, New York.
}
\author{ Zongming Ma and Xiangrui Meng}
\seealso{\code{\link{sep.tree-class}}.}
\keyword{ graphs }

