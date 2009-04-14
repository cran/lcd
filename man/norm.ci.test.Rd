\name{norm.ci.test} 
\alias{norm.ci.test}
\title{ Conditional independence test for multivariate normal data }
\description{
  Performs a conditional independence test for multivariate normal
  variables \code{u} and \code{v} with conditioning set \code{cond}.
}
\usage{
norm.ci.test(cov, n, u, v, cond = c())
}
\arguments{
  \item{cov}{ the sample covariance matrix. }
  \item{n}{ number of samples used to compute \code{cov}. }
  \item{u}{ name of one variable under investigation. }
  \item{v}{ name of the other variable under investigation. }
  \item{cond}{ a vector of variable names indicating the conditioning
    set in the test. }
}
\details{
  The function performs a likelihood ratio test as described in
  Whittaker (1990), Chapter 6. 
}
\value{
  \item{deviance}{ the Chisq-statistic computed under the null hypothesis. }
  \item{df}{ the degrees of freedom for the Chisq-statistic. }
  \item{p.value}{ the \eqn{p}-value for the test. }
}
\references{
  Whittaker, J. (1990). \emph{Graphical Models in Applied Mathematical
    Multivariate Statistics.} John Wiley and Sons, Chichester, England. 
  }
\author{  Zongming Ma and Xiangrui Meng}

\section{Note}{ Some of the \code{C} routines used were originally written by
  Robert Castelo <robert.castelo@upf.edu>, Alberto Roverato
  <alberto.roverato@unibo.it> in the \code{qp} package. The authors
  added the part for calculate the likelihood ratio test-statistic and
  modified the \code{R} wrapper function. 
}
\keyword{ htest }
\keyword{ multivariate }
