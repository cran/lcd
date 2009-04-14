\name{as.freq.tb}
\alias{as.freq.tb}
\title{ Frequency table transformation }
\description{
  Transforms a discrete data matrix into a raw frequency table.
}
\usage{
as.freq.tb(mat)
}
\arguments{
  \item{mat}{ a data matrix, with each row corresponding to an observation. }
} 
\details{
  Each element in the \code{mat} is expected to an integer or a factor.
}
\value{
  An object of class \code{freq.tb}
}
\author{Zongming Ma and Xiangrui Meng}
\seealso{ \code{\link{compress.freq.tb}}. }
\keyword{ manip }
\keyword{ multivariate }
