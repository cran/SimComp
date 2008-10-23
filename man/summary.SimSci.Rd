\name{summary.SimSci}
\alias{summary.SimSci}
\title{ Summary function for SimSci-objects }
\description{
A detailed print out of the results of \code{SimSciDiff} and \code{SimSciRat},
respectively.
}
\usage{
\method{summary}{SimSci}(object, digits = 4, ...)
}
\arguments{
  \item{object}{ an object of class "SimSci" as obtained by calling
                 \code{SimSciDiff} or \code{SimSciRat} }
  \item{digits}{ digits for rounding the results }
  \item{\dots}{ arguments to be passed to print }
}
\value{
A print out containing the estimates, raw and simultaneous confidence
intervals, estimated correlation matrices of the data and of the
comparisons computed by \code{SimSciDiff} or \code{SimSciRat}, respectively.
}
\author{ Mario Hasler }
\seealso{ \code{\link{summary.SimTest}} }
\keyword{ print }
