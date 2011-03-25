\name{summary.SimTest}
\alias{summary.SimTest}
\title{ Summary function for SimTest-objects }
\description{
A detailed print out of the results of \code{SimTestDiff} and \code{SimTestRat},
respectively.
}
\usage{
\method{summary}{SimTest}(object, digits = 4, ...)
}
\arguments{
  \item{object}{ an object of class "SimTest" as obtained by calling
                 \code{SimTestDiff} or \code{SimTestRat} }
  \item{digits}{ digits for rounding the results }
  \item{\dots}{ arguments to be passed to print }
}
\value{
A print out containing the estimates, test statistics, raw and
adjusted \emph{p}-values, estimated covariance correlation matrices of the data
and of the comparisons computed by \code{SimTestDiff} or \code{SimTestRat},
respectively.
}
\author{ Mario Hasler }
\seealso{ \code{\link{summary.SimCi}} }
\keyword{ print }
