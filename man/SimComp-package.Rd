\name{SimComp-package}
\alias{SimComp-package}
\alias{SimComp}
\docType{package}
\title{
Simultaneous Comparisons for Multiple Endpoints
}
\description{
Simultaneous tests and confidence intervals for general
linear hypotheses when there is more than one primary 
response variable (endpoint). The procedure of Hasler
(2009) is applied for differences or ratios of means of
normally distributed data. The covariance matrices
(containing the covariances between the endpoints) may
be assumed to be equal or possibly unequal for the
different groups.
}
\details{
\tabular{ll}{
Package: \tab SimComp\cr
Type: \tab Package\cr
Version: \tab 1.4.1\cr
Date: \tab 2009-10-14\cr
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

Hasler, M. and Hothorn, L.A. (submitted): Analysis of
randomized clinical dose finding trials with correlated
multiple endpoints using a multivariate Williams-type
trend procedure

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
