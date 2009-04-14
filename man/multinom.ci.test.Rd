\name{multinom.ci.test}
\alias{multinom.ci.test}
\title{ Conditional independence test for multinomial data }
\description{
  Performs a conditional independence test for variables \code{u} and
  \code{v} with conditioning set \code{cond}.
}
\usage{
multinom.ci.test(tb, u, v, cond = c())
} 
\arguments{
  \item{tb}{ a frequency table summarizing the data, an object of class
    \code{freq.tb}. }
  \item{u}{ name of one variable under investigation. }
  \item{v}{ name of the other variable under investigation. }
  \item{cond}{
    a vector of variable names indicating the conditioning
    set in the test.
  }
}
\details{
  The test is based on the deviance test in log-linear model for count
  data. The null hypothesis is that all interaction terms involving both
  \code{u} and \code{v} vanish.
}
\value{
  \item{deviance}{ the deviance of the saturated model from the model
    under null hypothesis.}
  \item{df}{ the degrees of freedom for the deviance.}
  \item{p.value}{ the approximate \eqn{p}-value for the test based on
    Chi-square approximation.}
}
\references{
  Agresti, A. (2002). \emph{Categorical Data Analysis.} 2nd Ed. John
  Wiley and Sons. Hoboken, NJ. 
}
\author{ Zongming Ma and Xiangrui Meng }
\keyword{ htest }
\keyword{ multivariate }
