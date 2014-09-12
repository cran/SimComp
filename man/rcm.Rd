\name{rcm}
\alias{rcm}
\title{ Random Correlation Matrices }
\description{
Correlation matrices with random off-diagonal elements.
}
\usage{
rcm(nrow = NULL, ncol = NULL)
}
\arguments{
  \item{nrow}{ the desired number of rows }
  \item{ncol}{ the desired number of columns }
}
\details{
As a correlation matrix is symmetric, only one of \code{nrow} or \code{ncol}
must be specified.
}
\value{
A symmetric correlation matrix with random elements.
}
\references{
Holmes, R.B. (1991): On random correlation matrices.
\emph{Siam Journal on Matrix Analysis and Applications} 12, 239--272.
}
\author{ Kornelius Rohmeyer and Mario Hasler }
\seealso{ \code{\link{ermvnorm}} }
\examples{
# Example 1:
# A correlation matrix representing three randomly correlated endpoints.

set.seed(1234)
rcm(nrow=3)

# Example 2:
# A correlation matrix representing five randomly correlated endpoints.

set.seed(5678)
rcm(ncol=5)
}
\keyword{ math }
\keyword{ misc }
