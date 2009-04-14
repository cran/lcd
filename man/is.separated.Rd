\name{is.separated}
\alias{is.separated}
\title{ c-separation on the chain graph }
\description{
  Checks whether two vertices \code{u} and \code{v} are
  \eqn{c}-separated by a set \code{sep} on the given chain graph.
}
\usage{
is.separated(u, v, sep, amat)
}
\arguments{
  \item{u}{ one vertex under investigation. }
  \item{v}{ the other vertex under investigation. }
  \item{sep}{ the candidate separator. }
  \item{amat}{ the adjacency matrix of a chain graph. }
}
\details{
  The actual function uses the moralization criterion instead of the
  separation criterion. They are equivalent and the former is easier to
  implement.
}
\value{
  A boolean value indicating whether separation holds or not.
}
\references{
  Lauritzen, S. L. (1996). \emph{Graphical Models.} Clarendon Press,
  Oxford.

  Studeny, M. and Bouckaert, R. R. (1998). On chain graph models for description
  of conditional independence structures. \emph{Annals of Statistics}
  \bold{26} 1434-1495. 
   
}
\author{ Zongming Ma and Xiangrui Meng }
\examples{
is.separated("a", "d", c("b","c"), toy.graph)
} 
\keyword{ graphs }

