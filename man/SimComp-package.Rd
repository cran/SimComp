\name{SimComp-package}
\alias{SimComp-package}
\alias{SimComp}
\docType{package}
\title{
Simultaneous Comparisons for Multiple Endpoints
}
\description{
This package provides simultaneous tests and confidence
intervals for one-way experimental designs with one or
many normally distributed, primary response variables
(endpoints). Means of several groups or dose levels can
be compared
\itemize{
  \item by arbitrary contrasts, like the Dunnett or the
    Tukey test,
  \item for balanced or unbalanced sample sizes,
  \item for homogeneous or heterogeneous group
    variances,
  \item for a single endpoint or for many endpoints
    simultaneously, and
  \item in terms of differences or ratios.
  }
If many endpoints are to be analyzed, the covariance
matrices (containing the covariances between the
endpoints) may be assumed to be equal or unequal for
the different groups. Exact or approximate multivariate
t-distributions, respectively, are used for quantiles
or p-values.

For example, the well-known conventional all-pair
comparison of Tukey (1953) can be performed by
specifying only a single endpoint and homogeneous (co-)
variances. On the other hand, it's also possible to the
same, but for many endpoints simultaneously, with
heterogeneous covariance matrices and in terms of
ratios.

For multiple comparisons of means of heteroscedastic
data, see Hasler and Hothorn (2008). The test procedure
for multiple endpoints is described by Hasler (2009).
}
\details{
\tabular{ll}{
Package: \tab SimComp\cr
Type: \tab Package\cr
Version: \tab 1.4.2\cr
Date: \tab 2009-10-22\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Mario Hasler

Maintainer: Mario Hasler <hasler@email.uni-kiel.de>
}
\references{
Hasler, M. (2009): Extensions of Multiple Contrast Tests.
\emph{PhD Thesis}, Gottfried-Wilhelm-Leibniz-Universitaet
Hannover.

Hasler, M. and Hothorn, L.A. (submitted): A Dunnett-type
Procedure for Multiple Endpoints

Hasler, M. and Hothorn, L.A. (2008): Multiple contrast
tests in the presence of heteroscedasticity.
\emph{Biometrical Journal} 50, 793-800.
}
\keyword{ package }
\seealso{ \code{\link[mratios:mratios-package]{mratios}} }
\examples{
data(iris)

comp <- SimTestRat(data=iris, grp="Species", alternative="greater", Margin=0.8)
summary(comp)
}
