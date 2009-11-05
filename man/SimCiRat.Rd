\name{SimCiRat}
\alias{SimCiRat}
\title{ Simultaneous Confidence Intervals for Ratios of Means of Multiple Endpoints }
\description{
Simultaneous confidence intervals for ratios of contrasts (linear functions) of
normal means (e.g., "Dunnett", "Tukey", "Williams" ect.) when there is more than
one primary response variable (endpoint). The procedure of Hasler (2009) is applied
for ratios of means of normally distributed data. The covariance matrices
(containing the covariances between the endpoints) may be assumed to be equal or
possibly unequal for the different groups. For the case of only a single endpoint,
the procedure reduces to the PI procedure of Hasler and Hothorn (2008).
}
\usage{
SimCiRat(data, grp, resp = NULL, type = "Dunnett", base = 1, Num.Contrast = NULL,
          Den.Contrast = NULL, alternative = "two.sided", covar.equal = FALSE,
          conf.level = 0.95)
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
               \item "Dunnett": many-to-one comparisons, with control in the
                                denominator
               \item "Tukey": all-pair comparisons
               \item "Sequen": comparisons of consecutive groups, where the group
                               with lower order is the denomniator
               \item "AVE": comparison of each group with average of all others,
                              where the average is taken as denominator
               \item "GrandMean": comparison of each group with grand mean of all
                                  groups, where the grand mean is taken as
                                    denominator
               \item "Changepoint": ratios of averages of groups of higher order
                                    divided by averages of groups of lower order
               \item "Marcus": Marcus contrasts as ratios
               \item "McDermott": McDermott contrasts as ratios
               \item "Williams": Williams contrasts as ratios
               \item "UmbrellaWilliams": Umbrella-protected Williams contrasts as
                                         ratios
               }
               note that \code{type} is ignored if \code{Num.Contrast} and
               \code{Den.Contrast} are specified by the user (see below) }
  \item{base}{ a single integer specifying the control (i.e. denominator) group for
               Dunnett contrasts, ignored otherwise }
  \item{Num.Contrast}{ a numerator contrast matrix, where columns correspond to
                       groups and rows correspond to contrasts }
  \item{Den.Contrast}{ a denominator contrast matrix, where columns correspond to
                       groups and rows correspond to contrasts }
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
The interest is in simultaneous confidence intervals for several ratios of linear
combinations (contrasts) of treatment means in a one-way ANOVA model, and
simultaneously for multiple endpoints. For example, corresponding intervals for the
all-pair comparison of Tukey (1953) and the many-to-one comparison of Dunnett (1955)
for ratios of means are implemented, but allowing for multiple endpoints. Also, the
user is free to create other interesting problem-specific contrasts. An approximate
multivariate \emph{t}-distribution is used to calculate lower and upper limits (see
Hasler, 2009). Simultaneous tests based on these intervals control the familywise
error rate in the strong sense. The covariance matrices of the treatment groups
(containing the covariances between the endpoints) can be assumed to be equal
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

In contrast to simultaneous confidence intervals for differences, the correlation
matrix for the multivariate \emph{t}-distribution depends on the unknown ratios. The
same problem also arises for the degrees of freedom if the covariance matrices for
the different groups are assumed to be unequal (\code{covar.equal=FALSE}). Both
problems can be handled by a plug-in approach, see the references therefore.
}
\references{
Hasler, M. (2009): Extensions of Multiple Contrast Tests. \emph{PhD Thesis},
Gottfried-Wilhelm-Leibniz-Universitaet Hannover.

Hasler, M. and Hothorn, L.A. (submitted): A Dunnett-type
Procedure for Multiple Endpoints

Hasler, M. and Hothorn, L.A. (2008): Multiple contrast tests in the presence
of heteroscedasticity. \emph{Biometrical Journal} 50, 793-800.

Dilba, G., Bretz, F., and Guiard, V. (2006): Simultaneous confidence sets and
confidence intervals for multiple ratios. \emph{Journal of Statistical Planning and
Inference} 136, 2640-2658.
}
\author{ Mario Hasler }
\seealso{ \code{\link{SimCiDiff}}, \code{\link{SimTestRat}},
          \code{\link{SimTestDiff}},  }
\examples{
# Example 1:
# Simultaneous confidence intervals for ratios of means, related to a
# Dunnett-test for the groups B and H against the standard S, on the (single)
# endpoint Thromb.count, assuming unequal variances for the groups. These are the
# well-known Dunnett-intervals but in the presence of heteroscedasticity and for
# ratios of means.

data(coagulation)

interv1 <- SimCiRat(data=coagulation, grp="Group", resp="Thromb.count", type="Dunnett", base=3, alternative="greater", covar.equal=FALSE)
interv1

# Example 2:
# Simultaneous confidence intervals for ratios of means, related to a
# Dunnett-test for the groups B and H against the standard S, simultaneously on
# all endpoints, assuming unequal covariance matrices for the groups.

data(coagulation)

interv2 <- SimCiRat(data=coagulation, grp="Group", resp=c("Thromb.count","ADP","TRAP"), type="Dunnett", base=3, alternative="greater", covar.equal=FALSE)
summary(interv2)
}
\keyword{ htest }
