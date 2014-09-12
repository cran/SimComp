\name{plot.SimCi}
\alias{plot.SimCi}
\title{ Plot function for SimCi-objects }
\description{
A plot of the results of \code{SimCiDiff} and \code{SimCiRat},
respectively.
}
\usage{
\method{plot}{SimCi}(x, xlim, xlab, ylim, ...)
}
\arguments{
  \item{x}{ an object of class "SimCi" as obtained by calling
           \code{SimCiDiff} or \code{SimCiRat} }
  \item{xlim}{ a numeric vector of length 2, giving the x coordinate range }
  \item{xlab}{ a title for the x axis }
  \item{ylim}{ a numeric vector of length 2, giving the y coordinate range }
  \item{\dots}{ arguments to be passed to \code{plot} }
}
\value{ A plot of the confidence intervals of a "SimCi" object. }
\author{ Christof Kluss and Mario Hasler }
\seealso{ \code{\link{SimCiDiff}}, \code{\link{SimCiRat}} }
\examples{
# Example 1:
# Simultaneous confidence intervals related to a comparison of the groups
# B and H against the standard S, on endpoint Thromb.count, assuming unequal
# variances for the groups. This is an extension of the well-known Dunnett-
# intervals to the case of heteroscedasticity.

data(coagulation)

interv1 <- SimCiDiff(data=coagulation, grp="Group", resp="Thromb.count",
  type="Dunnett", base=3, alternative="greater", covar.equal=FALSE)
interv1
plot(interv1)

# Example 2:
# Simultaneous confidence intervals related to a comparisons of the groups
# B and H against the standard S, simultaneously on all endpoints, assuming
# unequal covariance matrices for the groups. This is an extension of the well-
# known Dunnett-intervals to the case of heteroscedasticity and multiple
# endpoints.

data(coagulation)

interv2 <- SimCiDiff(data=coagulation, grp="Group", resp=c("Thromb.count","ADP","TRAP"),
  type="Dunnett", base=3, alternative="greater", covar.equal=FALSE)
summary(interv2)
par(mfrow=c(1,3)); plot(interv2)

# Example 3:
# Simultaneous confidence intervals for ratios of means, related to an all-pair
# comparison of the groups B, H and S, simultaneously on all endpoints, assuming unequal
# covariance matrices for the groups. This is an extension of the well-known Tukey-
# intervals to the case of heteroscedasticity and multiple endpoints, and in terms of
# ratios of means instead of differences.

data(coagulation)

interv3 <- SimCiRat(data=coagulation, grp="Group", resp=c("Thromb.count","ADP","TRAP"),
  type="Tukey", alternative="two.sided", covar.equal=FALSE)
summary(interv3)
par(mfrow=c(3,1)); plot(interv3)
}
\keyword{ print }
