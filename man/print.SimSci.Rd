\name{print.SimSci}
\alias{print.SimSci}
\title{ Print function for SimSci-objects }
\description{
A short print out of the results of \code{SimSciDiff} and \code{SimSciRat},
respectively.
}
\usage{
\method{print}{SimSci}(x, digits = 4, ...)
}
\arguments{
  \item{x}{ an object of class "SimSci" as obtained by calling
            \code{SimSciDiff} or \code{SimSciRat} }
  \item{digits}{ digits for rounding the results }
  \item{\dots}{ arguments to be passed to print }
}
\value{
A print out containing the estimates, raw and simultaneous confidence
intervals computed by \code{SimSciDiff} or \code{SimSciRat}, respectively.
}
\author{ Mario Hasler }
\seealso{ \code{\link{print.SimTest}} }
\keyword{ print }
