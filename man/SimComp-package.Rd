\name{SimComp-package}
\alias{SimComp-package}
\alias{SimComp}
\docType{package}
\title{
Simultaneous Comparisons for Multiple Endpoints
}
\description{
This package provides simultaneous tests and confidence intervals
for one-way experimental designs with one or many normally
distributed, primary response variables (endpoints). Means of
several groups or dose levels can be compared
\itemize{
  \item by arbitrary contrasts, like the Dunnett or the Tukey test,
  \item for balanced or unbalanced sample sizes,
  \item for a single endpoint or for many endpoints simultaneously,
  \item for homogeneous or heterogeneous variances/ covariance
    matrices of the groups, and
  \item in terms of differences or ratios.
}
Exact or approximate multivariate \emph{t}-distributions,
respectively, are used for quantiles or \emph{p}-values.

For example, the well-known conventional all-pair comparison of
Tukey (1953) can be performed by specifying only a single endpoint
and homogeneous (co-)variances. On the other hand, it's also
possible to do the same, but for many endpoints simultaneously,
with heterogeneous covariance matrices and in terms of ratios.

For multiple comparisons of means of heteroscedastic data, see
Hasler and Hothorn (2008). The test procedure for multiple
endpoints is described by Hasler and Hothorn (2011,2012).
}
\details{
\tabular{ll}{
Package: \tab SimComp\cr
Type: \tab Package\cr
Version: \tab 1.9.0\cr
Date: \tab 2014-02-28\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}

Index:
\itemize{
  \item \code{SimCiDiff}: Simultaneous Confidence Intervals for
    Differences of Means of Multiple Endpoints
  \item \code{SimCiRat}: Simultaneous Confidence Intervals for
    Ratios of Means of Multiple Endpoints
  \item \code{SimTestDiff}: Simultaneous Tests for Differences of
    Means of Multiple Endpoints
  \item \code{SimTestRat}: Simultaneous Tests for Ratios of Means
    of Multiple Endpoints
  \item \code{coagulation}: Data from a clinical study of three
    sets of extracorporeal circulation in heart-lung machines
  \item \code{ermvnorm}: Multivariate Normal Random Numbers with
    Exact Parameters
}
}
\author{
Mario Hasler

Maintainer: Mario Hasler <hasler@email.uni-kiel.de>
}
\references{
Hasler, M. and Hothorn, L.A. (2012): A multivariate Williams-type trend procedure.
\emph{Statistics in Biopharmaceutical Research} 4, 57--65.

Hasler, M. and Hothorn, L.A. (2011): A Dunnett-type procedure for multiple endpoints.
\emph{The International Journal of Biostatistics} 7, Article 3.

Hasler, M. and Hothorn, L.A. (2008): Multiple contrast tests in the presence of
heteroscedasticity. \emph{Biometrical Journal} 50, 793-800.

Dilba, G. et al. (2006): Simultaneous confidence sets and confidence intervals for
multiple ratios. \emph{Journal of Statistical Planning and Inference} 136, 2640-2658.
}
\keyword{ package }
\seealso{ \code{\link[mratios:mratios-package]{mratios}} }
\examples{
# Example 1:
# A Dunnett-test for the groups B and H against the standard S, on
# the (single) endpoint Thromb.count, assuming unequal variances for
# the groups. This is the well-known Dunnett-test but in the
# presence of heteroscedasticity.

data(coagulation)

comp1 <- SimTestDiff(data=coagulation, grp="Group", resp="Thromb.count",
  type="Dunnett", base=3, alternative="greater", covar.equal=FALSE)
comp1

# Example 2:
# A Dunnett-test for the groups B and H against the standard S,
# simultaneously on all endpoints, assuming unequal covariance
# matrices for the groups.

data(coagulation)

comp2 <- SimTestDiff(data=coagulation, grp="Group", resp=c("Thromb.count","ADP","TRAP"),
  type="Dunnett", base=3, alternative="greater", covar.equal=FALSE)
summary(comp2)
}
