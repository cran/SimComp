\name{print.SimTest}
\alias{print.SimTest}
\title{ Print function for SimTest-objects }
\description{
A short print out of the results of \code{SimTestDiff} and \code{SimTestRat},
respectively.
}
\usage{
\method{print}{SimTest}(x, digits = 4, ...)
}
\arguments{
  \item{x}{ an object of class "SimTest" as obtained by calling
            \code{SimTestDiff} or \code{SimTestRat} }
  \item{digits}{ digits for rounding the results }
  \item{\dots}{ arguments to be passed to \code{print} }
}
\value{
A print out containing the margins, estimates, test statistics, raw and
adjusted \emph{p}-values computed by \code{SimTestDiff} or
\code{SimTestRat}, respectively.
}
\author{ Mario Hasler }
\seealso{ \code{\link{print.SimCi}} }
\keyword{ print }
