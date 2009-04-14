\name{comp.pat} 
\alias{comp.pat}
\title{ Pattern comparison }
\description{
  Compares a (learned) chain graph pattern to the (supposed) true
  pattern. The two patterns should have the same vertex set in order for
  the function to return a meaningful result.
}
\usage{
comp.pat(truepat, pat)
}
\arguments{
  \item{truepat}{ the adjacency matrix of the true pattern. }
  \item{pat}{ the adjacency matrix of the pattern to be compared with
    the true one. }
}
\value{
  \item{a.total }{total number of complex arrows on the true pattern.}
  \item{a.missing}{number of true complex arrows missing in the pattern to be compared.}
  \item{a.extra }{number of spurious complex arrows present in the
    pattern to be compared.}
  \item{shd }{structural Hamming distance from \code{pat} to \code{truepat}.}
}
\note{
  Structural Hamming distance is defined as the total number of
  operations needed to convert one graph to the other. Each operation
  must be one of the following: (1) add or delete an undirected edge,
  or (2) add, remove or reverse an orientation of an edge.
}
\seealso{ \code{\link{comp.skel}}. }
\references{
  Tsamardinos, I., Brown, L., and Aliferis, C. (2006). The max-min hill-climbing Bayesian network structure
  learning algorithm. \emph{Mach. Learn.}, \bold{65}(1):31-78.
}
\author{Zongming Ma and Xiangrui Meng}
\keyword{ graphs }

