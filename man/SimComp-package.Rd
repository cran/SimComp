\name{SimComp-package}
\alias{SimComp-package}
\alias{SimComp}
\docType{package}
\title{
\packageTitle{SimComp}
}
\description{
\packageDescription{SimComp}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{SimComp}
\packageIndices{SimComp}

}
\author{
\packageAuthor{SimComp}

Maintainer: \packageMaintainer{SimComp}

Thanks to: Frank Schaarschmidt, Gemechis Djira Dilba, Kornelius Rohmeyer
}
\references{
Hasler, M. and Hothorn, L.A. (2018): Multi-arm trials with multiple primary endpoints and
missing values. \emph{Statistics in Medicine} 37, 710--721, <doi:10.1002/sim.7542>.

Hasler, M. (2014): Multiple contrast tests for multiple endpoints in the presence of 
heteroscedasticity. \emph{The International Journal of Biostatistics} 10, 17--28, 
<doi:10.1515/ijb-2012-0015>.

Hasler, M. and Hothorn, L.A. (2012): A multivariate Williams-type trend procedure.
\emph{Statistics in Biopharmaceutical Research} 4, 57--65, <doi:10.1080/19466315.2011.633868>.

Hasler, M. and Hothorn, L.A. (2011): A Dunnett-type procedure for multiple endpoints.
\emph{The International Journal of Biostatistics} 7, Article 3, <doi:10.2202/1557-4679.1258>.

Hasler, M. and Hothorn, L.A. (2008): Multiple contrast tests in the presence of
heteroscedasticity. \emph{Biometrical Journal} 50, 793--800, <doi:10.1002/bimj.200710466>.

Dilba, G. et al. (2006): Simultaneous confidence sets and confidence intervals for
multiple ratios. \emph{Journal of Statistical Planning and Inference} 136, 2640--2658,
<doi:10.1016/j.jspi.2004.11.009>.
}
\keyword{ package }
\examples{
# Example 1:
# A comparison of the groups B and H against the standard S, for endpoint
# Thromb.count, assuming unequal variances for the groups. This is an
# extension of the well-known Dunnett-test to the case of heteroscedasticity.

data(coagulation)

comp1 <- SimTestDiff(data=coagulation, grp="Group", resp="Thromb.count",
  type="Dunnett", base=3, alternative="greater", covar.equal=FALSE)
comp1

# Example 2:
# A comparison of the groups B and H against the standard S, simultaneously
# for all endpoints, assuming unequal covariance matrices for the groups. This is
# an extension of the well-known Dunnett-test to the case of heteroscedasticity
# and multiple endpoints.

data(coagulation)

comp2 <- SimTestDiff(data=coagulation, grp="Group", resp=c("Thromb.count","ADP","TRAP"),
  type="Dunnett", base=3, alternative="greater", covar.equal=FALSE)
summary(comp2)
}
