\name{SimCiDiff}
\alias{SimCiDiff}
\alias{SimCiDiff.default}
\alias{SimCiDiff.formula}
\title{ Simultaneous Confidence Intervals for General Contrasts (Differences) of 
Means of Multiple Endpoints }
\description{
Simultaneous confidence intervals for general contrasts (linear functions) of 
normal means (e.g., "Dunnett", "Tukey", "Williams" ect.), and for single or 
multiple endpoints (primary response variables) simultaneously. The procedure of 
Hasler and Hothorn (2011) <doi:10.2202/1557-4679.1258> is applied for 
differences of means of normally distributed data. The variances/ covariance 
matrices of the treatment groups (containing the covariances between the 
endpoints) may be assumed to be equal or possibly unequal for the different 
groups (Hasler, 2014 <doi:10.1515/ijb-2012-0015>). For the case of only a single 
endpoint and unequal covariance matrices (variances), the procedure coincides 
with the PI procedure of Hasler and Hothorn (2008) <doi:10.1002/bimj.200710466>.
}
\usage{
\method{SimCiDiff}{default}(data, grp, resp = NULL, na.action = "na.error", type = "Dunnett", 
  base = 1, ContrastMat = NULL, alternative = "two.sided", covar.equal = FALSE, 
  conf.level = 0.95, CorrMatDat = NULL, ...)
\method{SimCiDiff}{formula}(formula, ...)
}
\arguments{
  \item{data}{ a data frame containing a grouping variable and the endpoints as
               columns }
  \item{grp}{ a character string with the name of the grouping variable }
  \item{resp}{ a vector of character strings with the names of the endpoints; if
               \code{resp=NULL} (default), all column names of the data frame
               without the grouping variable are chosen automatically }
  \item{formula}{ a formula specifying a numerical response and a grouping factor 
                  (e.g. \kbd{response ~ treatment})  }
  \item{na.action}{ a character string indicating what should happen when the data 
                    contain \code{NAs};
                    if \code{na.action="na.error"} (default) the procedure stops 
                    with an error message; if \code{na.action="multi.df"} multiple 
                    marginal degrees of freedom are used to adjust for the missing 
                    values problem (not yet available) }
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
  \item{covar.equal}{ a logical variable indicating whether to treat the variances/ 
                      covariance matrices of the treatment groups (containing the 
                      covariances between the endpoints) as being equal;
                      if \code{TRUE} then the pooled variance/ covariance matrix is 
                      used, otherwise the Satterthwaite approximation to the degrees 
                      of freedom is used }
  \item{conf.level}{ a numeric value defining the simultaneous confidence level }
  \item{CorrMatDat}{ a correlation matrix of the endpoints, if \code{NULL} (default) 
                     it is estimated from the data }
  \item{\dots}{ arguments to be passed to SimCiDiff.default }
}
\details{
The interest is in simultaneous confidence intervals for several linear combinations
(contrasts) of treatment means in a one-way ANOVA model, and for single or 
multiple endpoints simultaneously. For example, corresponding intervals for the all-
pair comparison of Tukey (1953) and the many-to-one comparison of Dunnett (1955) are 
implemented, but allowing for heteroscedasticity and multiple endpoints. The user is 
also free to create other interesting problem-specific contrasts. Approximate 
multivariate \emph{t}-distributions are used to calculate lower and upper limits 
(Hasler and Hothorn, 2011 <doi:10.2202/1557-4679.1258>). Simultaneous tests based on 
these intervals control the familywise error rate in admissible ranges and in the 
strong sense. The variances/ covariance matrices of the treatment groups (containing 
the covariances between the endpoints) can be assumed to be equal 
(\code{covar.equal=TRUE}) or unequal (\code{covar.equal=FALSE}). If being equal, the 
pooled variance/ covariance matrix is used, otherwise approximations to the degrees 
of freedom (Satterthwaite, 1946) are used (Hasler, 2014 <doi:10.1515/ijb-2012-0015>; 
Hasler and Hothorn, 2008 <doi:10.1002/bimj.200710466>). Unequal covariance matrices 
occure if variances or correlations of some endpoints differ depending on the 
treatment groups.
}
\value{
An object of class SimCi containing:
\item{estimate}{ a matrix of estimated differences }
\item{lower.raw}{ a matrix of raw (unadjusted) lower limits }
\item{upper.raw}{ a matrix of raw (unadjusted) upper limits }
\item{lower}{ a matrix of lower limits adjusted for multiplicity }
\item{upper}{ a matrix of upper limits adjusted for multiplicity }
\item{CorrMatDat}{ if not prespecified by \code{CorrMatDat}, either the estimated 
                   common correlation matrix of the endpoints 
                   (\code{covar.equal=TRUE}) or a list of different (one for each 
                   treatment) estimated correlation matrices of the endpoints 
                   (\code{covar.equal=FALSE}) }
\item{CorrMatComp}{ the estimated correlation matrix of the comparisons }
\item{degr.fr}{ a matrix of degrees of freedom }
}
\note{
By default (\code{na.action="na.error"}), the procedure stops if there are 
missing values. A new experimental version for missing values is used if 
\code{na.action="multi.df"}. If \code{covar.equal=TRUE}, the number of endpoints 
must not be greater than the total sample size minus the number of treatment 
groups. If \code{covar.equal=FALSE}, the number of endpoints must not be greater 
than the minimal sample size minus 1. Otherwise the procedure stops.

All intervals have the same direction for all comparisons and endpoints
(\code{alternative="..."}). In case of doubt, use \code{"two.sided"}.
}
\references{
Hasler, M. (2014): Multiple contrast tests for multiple endpoints in the presence of 
heteroscedasticity. \emph{The International Journal of Biostatistics} 10, 17--28, 
<doi:10.1515/ijb-2012-0015>.

Hasler, M. and Hothorn, L.A. (2011): A Dunnett-type procedure for multiple endpoints.
\emph{The International Journal of Biostatistics} 7, Article 3, <doi:10.2202/1557-4679.1258>.

Hasler, M. and Hothorn, L.A. (2008): Multiple contrast tests in the presence of
heteroscedasticity. \emph{Biometrical Journal} 50, 793--800, <doi:10.1002/bimj.200710466>.
}
\author{ Mario Hasler }
\seealso{ \code{\link{SimTestDiff}}, \code{\link{SimTestRat}}, 
          \code{\link{SimCiRat}} }
\examples{
# Example 1:
# Simultaneous confidence intervals related to a comparison of the groups
# B and H against the standard S, for endpoint Thromb.count, assuming unequal
# variances for the groups. This is an extension of the well-known Dunnett-
# intervals to the case of heteroscedasticity.

data(coagulation)

interv1 <- SimCiDiff(data=coagulation, grp="Group", resp="Thromb.count",
  type="Dunnett", base=3, alternative="greater", covar.equal=FALSE)
interv1
plot(interv1)

# Example 2:
# Simultaneous confidence intervals related to a comparisons of the groups
# B and H against the standard S, simultaneously for all endpoints, assuming
# unequal covariance matrices for the groups. This is an extension of the well-
# known Dunnett-intervals to the case of heteroscedasticity and multiple
# endpoints.

data(coagulation)

interv2 <- SimCiDiff(data=coagulation, grp="Group", resp=c("Thromb.count","ADP","TRAP"),
  type="Dunnett", base=3, alternative="greater", covar.equal=FALSE)
summary(interv2)
plot(interv2)
}
\keyword{ htest }
