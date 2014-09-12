\name{print.SimCi}
\alias{print.SimCi}
\title{ Print function for SimCi-objects }
\description{
A short print out of the results of \code{SimCiDiff} and \code{SimCiRat},
respectively.
}
\usage{
\method{print}{SimCi}(x, digits = 4, ...)
}
\arguments{
  \item{x}{ an object of class "SimCi" as obtained by calling
            \code{SimCiDiff} or \code{SimCiRat} }
  \item{digits}{ digits for rounding the results }
  \item{\dots}{ arguments to be passed to \code{print} }
}
\value{
A print out containing the estimates, raw and simultaneous confidence
intervals computed by \code{SimCiDiff} or \code{SimCiRat}, respectively.
}
\author{ Mario Hasler }
\seealso{ \code{\link{print.SimTest}} }
\keyword{ print }
