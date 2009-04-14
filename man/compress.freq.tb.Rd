\name{compress.freq.tb}
\alias{compress.freq.tb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Frequency table compression }
\description{
  Compresses a frequency table to a subset of the variables.
}
\usage{
compress.freq.tb(tb, subset = colnames(tb@table)[-ncol(tb@table)]) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tb}{ the frequency table to be compressed, an object of class
    \code{freq.tb}. }
  \item{subset}{ a vector of strings indicating the subset of variables
    to be compressed on. The default is the whole variable set. }
}
\details{
  The actual computation routine is implemented in \code{C++}.
}
\value{
  The compressed table, an object of class \code{freq.tb}.
}
\author{ Zongming Ma and Xiangrui Meng }
\keyword{ manip }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
