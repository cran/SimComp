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
(2008) is applied for differences or ratios of means of
normally distributed data. The covariance matrices
(containing the covariances between the endpoints) may
be assumed to be equal or possibly unequal for the
different groups.
}
\details{
\tabular{ll}{
Package: \tab SimComp\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-10-23\cr
License: \tab GPL\cr
}
}
\author{
Mario Hasler

Maintainer: Mario Hasler <hasler@email.uni-kiel.de>
}
\references{
Hasler, M. (2008): Extensions of Multiple Contrast Tests.
\emph{PhD Thesis}, Gottfried-Wilhelm-Leibniz-Universitaet
Hannover.

Hasler, M. and Hothorn, L.A. (2008): Multiple contrast
tests in the presence of heteroscedasticity.
\emph{Biometrical Journal} 50, 793-800.
}
\keyword{ package }
\seealso{ \code{\link[<multtest>:<multtest>-package]{multtest}},
          \code{\link[<multcomp>:<multcomp>-package]{multcomp}},
          \code{\link[<mratios>:<mratios>-package]{mratios}} }
\examples{
data(iris)

comp <- SimTestRat(data=iris, grp="Species", alternative="greater", Margin=0.8)
summary(comp)
}
