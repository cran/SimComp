\name{SimCiDiff}
\alias{SimCiDiff}
\title{ Simultaneous Confidence Intervals for Differences of Means of Multiple Endpoints }
\description{
Simultaneous confidence intervals for general contrasts (linear functions) of
normal means (e.g., "Dunnett", "Tukey", "Williams" ect.) when there is more than
one primary response variable (endpoint). The procedure of Hasler (2008) is applied
for differences of means of normally distributed data. The covariance matrices
(containing the covariances between the endpoints) may be assumed to be equal or
possibly unequal for the different groups.
}
\usage{
SimCiDiff(data, grp, resp = NULL, type = "Dunnett", base = 1, ContrastMat = NULL,
           alternative = "two.sided", covar.equal = FALSE, conf.level = 0.95)
}
\arguments{
  \item{data}{ a data frame containing a grouping variable and the endpoints as
               columns }
  \item{grp}{ a character string with the name of the grouping variable }
  \item{resp}{ a vector of character strings with the names of the endpoints; if
               \code{resp=NULL} (default), all column names of the data frame
               without the grouping variable are chosen automatically }
  \item{type}{ a character string, defining the type of contrast, with the following
               options:
               \itemize{
               \item "Dunnett": many-to-one comparisons
               \item "Tukey": all-pair comparisons
               \item "Sequen": comparisons of consecutive groups
               \item "AVE": comparison of each group with average of all others
               \item "GrandMean": comparison of each group with grand mean of all
                                  groups
               \item "Changepoint": differences of averages of groups of higher
                                    order to averages of groups of lower order
               \item "Marcus": Marcus contrasts
               \item "McDermott": McDermott contrasts
               \item "Williams": Williams trend tests
               \item "UmbrellaWilliams": Umbrella-protected Williams trend tests
               }
               note that \code{type} is ignored if \code{ContrastMat} is specified
               by the user (see below) }
  \item{base}{ a single integer specifying the control group for Dunnett contrasts,
               ignored otherwise }
  \item{ContrastMat}{ a contrast matrix, where columns correspond to groups and rows
                      correspond to contrasts }
  \item{alternative}{ a character string specifying the alternative hypothesis,
                      must be one of \code{"two.sided"} (default), \code{"greater"}
                      or \code{"less"} }
  \item{covar.equal}{ a logical variable indicating whether to treat the covariance
                      matrices (containing the covariances between the endpoints)
                      for the different groups as being equal;
                      if \code{TRUE} then the pooled covariance matrix is used,
                      otherwise the Satterthwaite approximation to the degrees of
                      freedom is used according to Hasler and Hothorn (2008)}
  \item{conf.level}{ a numeric value defining the simultaneous confidence level }
}
\details{
The interest is in simultaneous confidence intervals for several linear combinations
(contrasts) of treatment means in a one-way ANOVA model, and simultaneously for
multiple endpoints. For example, corresponding intervals for the all-pair comparison
of Tukey (1953) and the many-to-one comparison of Dunnett (1955) are implemented,
but allowing for multiple endpoints. Also, the user is free to create other
interesting problem-specific contrasts. An approximate multivariate \emph{t}-
distribution is used to calculate lower and upper limits (see Hasler, 2008).
Simultaneous tests based on these intervals control the familywise error rate in the
strong sense. The covariance matrices of the treatment groups (containing the
covariances between the endpoints) can be assumed to be equal
(\code{covar.equal=TRUE}) or unequal (\code{covar.equal=FALSE}). If being equal, the
pooled covariance matrix is used, otherwise the Satterthwaite approximation to the
degrees of freedom is used according to Hasler and Hothorn (2008). Unequal
covariance matrices occure if either variances or correlations of some endpoints
differ depending on the treatment groups.
}
\value{
An object of class SimCi containing:
\item{estimate}{ a matrix of estimated differences }
\item{lower.raw}{ a matrix of raw (unadjusted) lower limits }
\item{upper.raw}{ a matrix of raw (unadjusted) upper limits }
\item{lower}{ a matrix of lower limits adjusted for multiplicity }
\item{upper}{ a matrix of upper limits adjusted for multiplicity }
\item{CorrMatDat}{ either the estimated common correlation matrix of the data
                   (\code{covar.equal=TRUE}) or the list of the different (one for
                   each treatment) estimated correlation matrices of the data
                   (\code{covar.equal=FALSE}) }
\item{CorrMatComp}{ the estimated correlation matrix to be used for the multivariate
                    \emph{t}-distribution }
\item{degr.fr}{ either a single degree of freedom (\code{covar.equal=TRUE}) or a
                matrix of degrees of freedom (\code{covar.equal=FALSE}) }
}
\note{
All measurement objects of each treatment group must have values for each endpoint.
If there are missing values then the procedure stops. If \code{covar.equal=TRUE},
then the number of endpoints must not be greater than the total sample size minus
the number of treatment groups. If \code{covar.equal=FALSE}, the number of endpoints
must not be greater than the minimal sample size minus 1. Otherwise the procedure
stops.

All the intervals have the same direction for all comparisons and endpoints
(\code{alternative="..."}). In case of doubts, use \code{"two.sided"}.
}
\references{
Hasler, M. (2008): Extensions of Multiple Contrast Tests. \emph{PhD Thesis},
Gottfried-Wilhelm-Leibniz-Universitaet Hannover.

Hasler, M. and Hothorn, L.A. (2008): Multiple contrast tests in the presence
of heteroscedasticity. \emph{Biometrical Journal} 50, 793-800.
}
\author{ Mario Hasler }
\seealso{ \code{\link{SimCiRat}}, \code{\link{SimTestDiff}},
          \code{\link{SimTestRat}},  }
\examples{
data(iris)

interv <- SimCiDiff(data=iris, grp="Species")
summary(interv)
}
\keyword{ htest }