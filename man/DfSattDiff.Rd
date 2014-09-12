\name{DfSattDiff}
\alias{DfSattDiff}
\title{ Degrees of Freedom Accoding to Satterthwaite (1946) for Differences of Means }
\description{
Degrees of freedom accoding to Satterthwaite (1946) for (multivariate)
\emph{t}-distributions related to multiple contrast tests or corresponding
simultaneous confidence intervals for differences of means. For contrasts representing
a two-sample \emph{t}-test, the degree of freedom coincides with the one of Welch
(1938).
}
\usage{
DfSattDiff(n, sd, type = "Dunnett", base = 1, ContrastMat = NULL)
}
\arguments{
  \item{n}{ a vector of numbers of observations }
  \item{sd}{ a vector of standard deviations }
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
}
\details{
The calculation of critical values or (adjusted) \emph{p}-values related to multiple
contrast tests or corresponding simultaneous confidence intervals is based on a
multivariate \emph{t}-distribution. For homoscedastic data, the respective degree of
freedom only depends on the total sample size and the number of groups. A simple and
well-known special case ist the usual \emph{t}-test. If the data are heteroscedastic,
however, the degree of freedom of a usual \emph{t}-test must be decreased according
to Welch (1938) to come to an approximate solution. Degrees of freedom according to
Satterthwaite (1946) refer to any linear combinations (contrasts) of normal means.
They are applied, for example, when doing multiple contrast tests for heteroscedastic
data according to Hasler and Hothorn (2008) or Hasler (2014). Like Welch (1938), 
Satterthwaite (1946) approximated the degree of freedom by matching first and second
moments. The resulting degree of freedom then depends on the contrast and on the
sample sizes and sample variances per group, respectively.
}
\value{
A vector of degrees of freedom.
}
\note{
The commands \code{SimTestDiff()} and \code{SimCiDiff()} use these degrees of
freedom automatically if \code{covar.equal=FALSE} (default). You don't need to apply
\code{DfSattDiff()} additionally.
}
\references{
Hasler, M. (2014): Multiple contrast tests for multiple endpoints in the presence of 
heteroscedasticity. \emph{The International Journal of Biostatistics} 10, 17--28.

Hasler, M. and Hothorn, L.A. (2008): Multiple contrast tests in the presence of
heteroscedasticity. \emph{Biometrical Journal} 50, 793--800.

Satterthwaite, F.E. (1946): An approximate distribution of estimates of variance
components. \emph{Biometrics} 2, 110--114.

Welch, B.L. (1938): The significance of the difference between two means when the
population variances are unequal. \emph{Biometrika} 29, 350--362.
}
\author{ Mario Hasler }
\seealso{ \code{\link{DfSattRat}} }
\examples{
# Example 1:
# Degrees of freedom for a comparison of group two and three against group one, assuming
# unequal standard deviations for the groups. This is an extension for the well-known
# Dunnett-test to the case of heteroscedasticity.

# Either by specifying the type of contrast:
DfSattDiff(n=c(10,6,6), sd=c(1,3,6), type="Dunnett", base=1)

# Or by specifying the contrast matrix:
DfSattDiff(n=c(10,6,6), sd=c(1,3,6), ContrastMat=rbind(c(-1,1,0),c(-1,0,1)))

# Example 2:
# Degrees of freedom for an all-pair comparison of the groups B, H and S on endpoint ADP,
# assuming unequal standard deviations for the groups. This is an extension for the well-
# known Tukey-test to the case of heteroscedasticity.

data(coagulation)

DfSattDiff(n=tapply(X=coagulation$ADP, INDEX=coagulation$Group, FUN=length),
  sd=tapply(X=coagulation$ADP, INDEX=coagulation$Group, FUN=sd),
  type="Tukey")

# These are the same degrees of freedom as used automatically by command
# \code{SimTestDiff()}:
test <- SimTestDiff(data=coagulation, grp="Group", resp="ADP", type="Tukey",
  covar.equal=FALSE)
test$degr.fr
}
\keyword{ math }
\keyword{ misc }
