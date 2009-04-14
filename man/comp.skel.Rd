\name{comp.skel}
\alias{comp.skel}
\title{ Skeleton comparison }
\description{
  Compares a graph skeleton to the (supposed) true skeleton. The two
  skeletons should have the same set of vertices in order for the
  function to return a meaningful result.
}
\usage{
comp.skel(trueskel, skel)
}
\arguments{ 
  \item{trueskel}{ the adjacency matrix of the true skeleton. }
  \item{skel}{ the adjacency matrix of the skeleton to be compared with
    the true one. }
}
\value{
  \item{e.total }{total number of edges on the true skeleton.}
  \item{e.missing }{number of true edges missing in the skeleton to be compared.}
  \item{e.extra }{number of spurious edges present in the skeleton to be compared.}
}
\seealso{ \code{\link{comp.pat}}. }
\keyword{ graphs } 

