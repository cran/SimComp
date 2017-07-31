SimCiDiffHom <-
function(trlist, grp, ntr, nep, ssmat, ContrastMat, ncomp, alternative, 
                         conf.level, meanmat, CorrMatDat) {


estimate <- ContrastMat%*%meanmat

defr <- matrix(sum(ssmat[,1])-ntr, nrow=ncomp, ncol=nep)            # matrix of dfs

CovMatDat <- Reduce("+",
                    lapply(trlist, function(x) (nrow(x)-1)*cov(x))
                   )/defr[1,1]                                      # covariance matrix of the data
if (is.null(CorrMatDat)) {
  CorrMatDat <- cov2cor(CovMatDat)                                  # correlation matrix of the data
} else {
  sdmat <- sqrt( diag( diag(CovMatDat),nrow=nep ) )                 # sds on the diagonal
  CovMatDat <- sdmat%*%CorrMatDat%*%sdmat
}

M <- diag(1/ssmat[,1])
CorrMatCont <- diag( 1/sqrt(diag( ContrastMat %*% M %*% t(ContrastMat) )),nrow=ncomp ) %*% 
                                 (ContrastMat %*% M %*% t(ContrastMat)) %*%
               diag( 1/sqrt(diag( ContrastMat %*% M %*% t(ContrastMat) )),nrow=ncomp )
R <- kronecker(CorrMatCont, CorrMatDat, make.dimnames=TRUE)

lower <- upper <- lower.raw <- upper.raw <- matrix(nrow=ncomp,ncol=nep)

if (alternative=="greater") {
  lo1malqu <- qmvt(conf.level,tail="lower.tail",df=as.integer(defr[1,1]),corr=R)$quantile
  univarqu <- qt(p=conf.level, df=as.integer(defr[1,1]))
  for (z in 1:ncomp) {
    for (i in 1:nep) {
      upper[z,i] <- upper.raw[z,i] <- Inf
      lower[z,i]     <- t(ContrastMat[z,])%*%meanmat[,i] -
                        lo1malqu * sqrt( diag(CovMatDat)[i]*( t(ContrastMat[z,])%*%M%*%ContrastMat[z,] ) )
      lower.raw[z,i] <- t(ContrastMat[z,])%*%meanmat[,i] -
                        univarqu * sqrt( diag(CovMatDat)[i]*( t(ContrastMat[z,])%*%M%*%ContrastMat[z,] ) )
    }
  }
}
if (alternative=="less") {
  up1malqu <- qmvt(conf.level,tail="upper.tail",df=as.integer(defr[1,1]),corr=R)$quantile
  univarqu <- qt(p=1-conf.level, df=as.integer(defr[1,1]))
  for (z in 1:ncomp) {
    for (i in 1:nep) {
      upper[z,i]     <- t(ContrastMat[z,])%*%meanmat[,i] -
                        up1malqu * sqrt( diag(CovMatDat)[i]*( t(ContrastMat[z,])%*%M%*%ContrastMat[z,] ) )
      upper.raw[z,i] <- t(ContrastMat[z,])%*%meanmat[,i] -
                        univarqu * sqrt( diag(CovMatDat)[i]*( t(ContrastMat[z,])%*%M%*%ContrastMat[z,] ) )
      lower[z,i] <- lower.raw[z,i] <- -Inf
    }
  }
}
if (alternative=="two.sided") {
  ts1malqu <- qmvt(conf.level,tail="both.tails",df=as.integer(defr[1,1]),corr=R)$quantile
  univarqu <- qt(p=1-(1-conf.level)/2, df=as.integer(defr[1,1]))
  for (z in 1:ncomp) {
    for (i in 1:nep) {
      upper[z,i]     <- t(ContrastMat[z,])%*%meanmat[,i] +
                        ts1malqu * sqrt( diag(CovMatDat)[i]*( t(ContrastMat[z,])%*%M%*%ContrastMat[z,] ) )
      upper.raw[z,i] <- t(ContrastMat[z,])%*%meanmat[,i] +
                        univarqu * sqrt( diag(CovMatDat)[i]*( t(ContrastMat[z,])%*%M%*%ContrastMat[z,] ) )
      lower[z,i]     <- t(ContrastMat[z,])%*%meanmat[,i] -
                        ts1malqu * sqrt( diag(CovMatDat)[i]*( t(ContrastMat[z,])%*%M%*%ContrastMat[z,] ) )
      lower.raw[z,i] <- t(ContrastMat[z,])%*%meanmat[,i] -
                        univarqu * sqrt( diag(CovMatDat)[i]*( t(ContrastMat[z,])%*%M%*%ContrastMat[z,] ) )
    }
  }
}

list(estimate=estimate, lower.raw=lower.raw, upper.raw=upper.raw, lower=lower, upper=upper,
     CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defr, 
     ContrastMat=ContrastMat, alternative=alternative, conf.level=conf.level)


}
