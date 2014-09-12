\name{summary.SimCi}
\alias{summary.SimCi}
\title{ Summary function for SimCi-objects }
\description{
A detailed print out of the results of \code{SimCiDiff} and \code{SimCiRat},
respectively.
}
\usage{
\method{summary}{SimCi}(object, digits = 4, ...)
}
\arguments{
  \item{object}{ an object of class "SimCi" as obtained by calling
                 \code{SimCiDiff} or \code{SimCiRat} }
  \item{digits}{ digits for rounding the results }
  \item{\dots}{ arguments to be passed to \code{print} }
}
\value{
A print out containing the estimates, raw and simultaneous confidence
intervals, estimated covariance and correlation matrices of the data and
of the comparisons computed by \code{SimCiDiff} or \code{SimCiRat},
respectively.
}
\author{ Mario Hasler }
\seealso{ \code{\link{summary.SimTest}} }
\keyword{ print }
