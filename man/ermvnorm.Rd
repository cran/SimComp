\name{ermvnorm}
\alias{ermvnorm}
\title{ Multivariate Normal Random Numbers with Exact Parameters }
\description{
Random numbers of the multivariate normal distribution with EXACT mean vector, EXACT
variance vector and approximate correlation matrix. This function is based on the
function \code{rmvnorm} of the package mvtnorm.
}
\usage{
ermvnorm(n, mean, sd, corr = diag(rep(1, length(mean))),
  mnt = 10000)
}
\arguments{
  \item{n}{ number of observations }
  \item{mean}{ mean vector }
  \item{sd}{ vector of standard deviations }
  \item{corr}{ correlation matrix }
  \item{mnt}{ maximum number of tries for the computation }
}
\details{
Unfortunately, it's very common to present only summary statistics in the literature
when evaluating real data. This makes it hard to retrace or to verify the related
statistical evaluation. Also, the use of such data as an example for other
statistical tests is not possible. For that reason, \code{ermvnorm} allows to
reproduce data by simulation. In contrast to \code{rmvnorm} of the package mvtnorm,
the function \code{ermvnorm} produces random numbers that have EXACTLY the same
paramer values as specified by \code{mean} and \code{sd}. The correlation matrix
\code{corr} is met only approximately.

The simple idea behind \code{ermvnorm} is to apply \code{rmvnorm} of the package
mvtnorm, but only for the first n-2 random numbers. The remaining 2 numbers are
obtained by solving a quadratic equation to achieve the specified values for the
mean vector and for the vector of standard deviations. Depending on the n-2 random
numbers, the underlying quadratic equation can possibly have no solution. In this
case, \code{ermvnorm} creates a new set of n-2 random numbers until a valid data
set is obtained, or until the maximum number of tries \code{mnt} is reached.
}
\value{
A matrix of random numbers with dimension n * length(mean).
}
\references{
Hothorn, T. et al. (2001): On Multivariate t and Gauss Probabilities in R.
\emph{R News} 1, 27-29.
}
\author{
Mario Hasler and Gemechis Djira Dilba
}
\note{
This function is to be used only with caution. Normally, random numbers with exact
mean and standard deviation are not intended to be used. For example, simulations
concerning type I error or power of statistical tests cannot be based on
\code{ermvnorm}.
}
\seealso{ \code{\link{rmvnorm}} }
\examples{
# Example 1:
# A dataset representing two correlated endpoints.

set.seed(1234)
dataset1 <- ermvnorm(n=10, mean=c(10,120), sd=c(1,10),
  corr=rbind(c(1,0.7),c(0.7,1)))
dataset1
mean(dataset1[,1]); mean(dataset1[,2])
sd(dataset1[,1]); sd(dataset1[,2])
round(cor(dataset1),3)
pairs(dataset1)

# Example 2:
# A dataset representing three uncorrelated endpoints.

set.seed(5678)
dataset2 <- ermvnorm(n=20, mean=c(1,12,150), sd=c(0.5,2,20))
dataset2
mean(dataset2[,1]); mean(dataset2[,2]); mean(dataset2[,3])
sd(dataset2[,1]); sd(dataset2[,2]); sd(dataset2[,3])
pairs(dataset2)
}
\keyword{ datagen }
