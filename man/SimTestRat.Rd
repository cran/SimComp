\name{SimTestRat}
\alias{SimTestRat}
\title{ Simultaneous Tests for Ratios of Means of Multiple Endpoints }
\description{
Simultaneous tests for ratios of contrasts (linear functions) of normal means (e.g.,
"Dunnett", "Tukey", "Williams" ect.) when there is more than one primary response
variable (endpoint). The procedure of Hasler and Hothorn (2011) is applied for
ratios of means of normally distributed data. The covariance matrices (containing
the covariances between the endpoints) may be assumed to be equal or possibly
unequal for the different groups. For the case of only a single endpoint and
unequal covariance matrices (variances), the procedure coincides with the PI
procedure of Hasler and Hothorn (2008).
}
\usage{
SimTestRat(data, grp, resp = NULL, type = "Dunnett", base = 1, Num.Contrast = NULL,
           Den.Contrast = NULL, alternative = "two.sided", Margin = NULL,
           covar.equal = FALSE)
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
  \item{Margin}{ a single numeric value, or a numeric vector corresponding to
                 endpoints, or a matrix where columns correspond to endpoints and
                 rows correspond to contrasts, default is 1 }
  \item{covar.equal}{ a logical variable indicating whether to treat the covariance
                      matrices (containing the covariances between the endpoints)
                      for the different groups as being equal;
                      if \code{TRUE} then the pooled covariance matrix is used,
                      otherwise the Satterthwaite approximation to the degrees of
                      freedom is used according to Hasler and Hothorn (2008)}
}
\details{
The interest is in simultaneous tests for several ratios of linear combinations
(contrasts) of treatment means in a one-way ANOVA model, and simultaneously for
multiple endpoints. For example, the all-pair comparison of Tukey (1953) and the
many-to-one comparison of Dunnett (1955) for ratios of means are implemented, but
allowing for multiple endpoints. Also, the user is free to create other interesting
problem-specific contrasts. An approximate multivariate \emph{t}-distribution is
used to calculate (adjusted) \emph{p}-values (see Hasler and Hothorn, 2011). This
approach controls the familywise error rate in an admissible range and in the
strong sense. The covariance matrices of the treatment groups (containing the
covariances between the endpoints) can be assumed to be equal
(\code{covar.equal=TRUE}) or unequal (\code{covar.equal=FALSE}). If being equal,
the pooled covariance matrix is used, otherwise Satterthwaite approximations to
the degrees of freedom are used according to Hasler and Hothorn (2008). Unequal
covariance matrices occure if variances or correlations of some endpoints differ
depending on the treatment groups.
}
\value{
An object of class SimTest containing:
\item{estimate}{ a matrix of estimated ratios }
\item{statistic}{ a matrix of the calculated test statistics }
\item{p.val.raw}{ a matrix of raw \emph{p}-values }
\item{p.val.adj}{ a matrix of \emph{p}-values adjusted for multiplicity }
\item{CorrMatDat}{ either the estimated common correlation matrix of the data
                   (\code{covar.equal=TRUE}) or the list of the different (one for
                   each treatment) estimated correlation matrices of the data
                   (\code{covar.equal=FALSE}) }
\item{CorrMatComp}{ the estimated correlation matrix to be used for the multivariate
                    \emph{t}-distribution }
\item{degr.fr}{ either a single degree of freedom (\code{covar.equal=TRUE}) or a
                vector of degrees of freedom (\code{covar.equal=FALSE}) related
                to the comparisons }
}
\note{
All measurement objects of each treatment group must have values for each endpoint.
If there are missing values then the procedure stops. If \code{covar.equal=TRUE},
then the number of endpoints must not be greater than the total sample size minus
the number of treatment groups. If \code{covar.equal=FALSE}, the number of endpoints
must not be greater than the minimal sample size minus 1. Otherwise the procedure
stops.

All hypotheses are tested with the same test direction for all comparisons and
endpoints (\code{alternative="..."}). In case of doubt, use \code{"two.sided"}.

If \code{Margin} is a single numeric value or a numeric vector, then the same
value(s) are used for the remaining comparisons or endpoints. If \code{Margin} is
not specified, the default is 1.
}
\references{
Hasler, M. and Hothorn, L.A. (2011): A Dunnett-type procedure for multiple endpoints.
\emph{The International Journal of Biostatistics} 7, Article 3.

Hasler, M. and Hothorn, L.A. (2008): Multiple contrast tests in the presence of
heteroscedasticity. \emph{Biometrical Journal} 50, 793-800.

Dilba, G. et al. (2006): Simultaneous confidence sets and confidence intervals for
multiple ratios. \emph{Journal of Statistical Planning and Inference} 136, 2640-2658.
}
\author{ Mario Hasler }
\seealso{ \code{\link{SimTestDiff}}, \code{\link{SimCiRat}},
          \code{\link{SimCiDiff}},  }
\examples{
# Example 1:
# A Dunnett-test for the groups B and H against the standard S, on the (single)
# endpoint Thromb.count, assuming unequal variances for the groups, and in terms of
# ratios. This is the well-known Dunnett-test but in the presence of
# heteroscedasticity and for ratios of means.

data(coagulation)

comp1 <- SimTestRat(data=coagulation, grp="Group", resp="Thromb.count",
  type="Dunnett", base=3, alternative="greater", covar.equal=FALSE)
comp1

# Example 2:
# A Dunnett-test for the groups B and H against the standard S, simultaneously on
# all endpoints, assuming unequal covariance matrices for the groups, and in terms
# of ratios.

data(coagulation)

comp2 <- SimTestRat(data=coagulation, grp="Group", resp=c("Thromb.count","ADP","TRAP"),
  type="Dunnett", base=3, alternative="greater", covar.equal=FALSE)
summary(comp2)
}
\keyword{ htest }
