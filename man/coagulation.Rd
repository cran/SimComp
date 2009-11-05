\name{coagulation}
\alias{coagulation}
\docType{data}
\title{
Data from a clinical study for the comparison of three sets of extracorporeal circulation
in heart-lung machines
}
\description{
Three sets of extracorporeal circulation in heart-lung machines: treatments H and B, and
standard S. Twelve (S and H each) and eleven (B) male adult patients. The analysis is
based on a set of laboratory parameters restricted to the blood coagulation system,
characterized by three primary endpoints (each as quotient from post- and pre-surgery
values). Higher values indicate a better treatment effect. For more details, see Kropf
et. al, 2000.
}
\usage{data(coagulation)}
\format{
  A data frame with 35 observations on the following 5 variables.
  \describe{
    \item{\code{Patient}}{a numeric vector, the patients' number}
    \item{\code{Thromb.count}}{a numeric vector}
    \item{\code{ADP}}{a numeric vector}
    \item{\code{TRAP}}{a numeric vector}
    \item{\code{Group}}{a factor with levels \code{B} \code{H} \code{S}
      specifying the treatments, where \code{S} is the standard}
  }
}
\source{
Kropf S, Hommel G, Schmidt U, Brickwedel J, and Jepsen MS. Multiple comparisons of
treatments with stable multivariate tests in a two-stage adaptive design, including a test
for non-inferiority. Biometrical Journal 2000; 42(8):951-965
}
\examples{
data(coagulation)
str(coagulation)
}
\keyword{datasets}
