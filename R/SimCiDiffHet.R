SimCiDiffHet <-
function(trlist, grp, ntr, nep, ssmat, ContrastMat, ncomp, alternative, 
                         conf.level, meanmat, CorrMatDat) {


varmat  <- do.call(rbind, lapply(trlist, function(x) apply(x, 2, var)))
estimate <- ContrastMat%*%meanmat

defrmat <- matrix(nrow=ncomp, ncol=nep)
for (j in 1:nep) {
  defrmat[,j] <- DfSattDiff(n=ssmat[,1],sd=sqrt(varmat[,j]),ContrastMat=ContrastMat)
}
defrmat[defrmat<2] <- 2                                             # to be well-defined
defr <- matrix(apply(defrmat,1,min), nrow=ncomp, ncol=nep)          # matrix of dfs, minimum per row/contrast

CovMatDat  <- lapply(trlist, cov)                                   # list of covariance matrices of the data
if (is.null(CorrMatDat)) {
  CorrMatDat <- lapply(CovMatDat,cov2cor)                           # list of correlation matrices of the data
} else {
  sdmat <- lapply(CovMatDat, function(x) sqrt( diag( diag(x),nrow=nep ) )) # sds on the diagonal
  CovMatDat <- lapply(sdmat, function(x) x%*%CorrMatDat%*%x)        # final list of covariance matrices
}

M <- diag(1/ssmat[,1])
R <- NULL
for (z in 1:ncomp) {
  Rrow <- NULL
  for (w in 1:ncomp) {
    Rpart <- matrix(nrow=nep,ncol=nep)
    for (i in 1:nep) {
      for (h in 1:nep) {
        Rpart[i,h] <- ( t(ContrastMat[z,])%*%diag( sapply(CovMatDat, function(x) x[i,h]) )%*%M%*%
                         (ContrastMat[w,]) ) /
                      sqrt( ( t(ContrastMat[z,])%*%diag(varmat[,i])%*%M%*%(ContrastMat[z,]) ) *
                            ( t(ContrastMat[w,])%*%diag(varmat[,h])%*%M%*%(ContrastMat[w,]) ) )
      }
    }
    Rrow <- cbind(Rrow,Rpart)
  }
  R <- rbind(R, Rrow)                                               # correlation matrix for multi-t
}
diag(R) <- 1

lower <- upper <- lower.raw <- upper.raw <- matrix(nrow=ncomp,ncol=nep)

if (alternative=="greater") {
  for (z in 1:ncomp) {
    for (i in 1:nep) {
      lo1malqu <- qmvt(conf.level,tail="lower.tail",df=as.integer(defr[z,i]),corr=R)$quantile
      univarqu <- qt(p=conf.level, df=defrmat[z,i])
      upper[z,i] <- upper.raw[z,i] <- Inf
      lower[z,i]     <- t(ContrastMat[z,])%*%meanmat[,i] -
                        lo1malqu * sqrt( t(ContrastMat[z,])%*%diag(varmat[,i])%*%M%*%ContrastMat[z,] )
      lower.raw[z,i] <- t(ContrastMat[z,])%*%meanmat[,i] -
                        univarqu * sqrt( t(ContrastMat[z,])%*%diag(varmat[,i])%*%M%*%ContrastMat[z,] )
    }
  }
}
if (alternative=="less") {
  for (z in 1:ncomp) {
    for (i in 1:nep) {
      up1malqu <- qmvt(conf.level,tail="upper.tail",df=as.integer(defr[z,i]),corr=R)$quantile
      univarqu <- qt(p=1-conf.level, df=defrmat[z,i])
      upper[z,i]     <- t(ContrastMat[z,])%*%meanmat[,i] -
                        up1malqu * sqrt( t(ContrastMat[z,])%*%diag(varmat[,i])%*%M%*%ContrastMat[z,] )
      upper.raw[z,i] <- t(ContrastMat[z,])%*%meanmat[,i] -
                        univarqu * sqrt( t(ContrastMat[z,])%*%diag(varmat[,i])%*%M%*%ContrastMat[z,] )
      lower[z,i] <- lower.raw[z,i] <- -Inf
    }
  }
}
if (alternative=="two.sided") {
  for (z in 1:ncomp) {
    for (i in 1:nep) {
      ts1malqu <- qmvt(conf.level,tail="both.tails",df=as.integer(defr[z,i]),corr=R)$quantile
      univarqu <- qt(p=1-(1-conf.level)/2, df=defrmat[z,i])
      upper[z,i]     <- t(ContrastMat[z,])%*%meanmat[,i] +
                        ts1malqu * sqrt( t(ContrastMat[z,])%*%diag(varmat[,i])%*%M%*%ContrastMat[z,] )
      upper.raw[z,i] <- t(ContrastMat[z,])%*%meanmat[,i] +
                        univarqu * sqrt( t(ContrastMat[z,])%*%diag(varmat[,i])%*%M%*%ContrastMat[z,] )
      lower[z,i]     <- t(ContrastMat[z,])%*%meanmat[,i] -
                        ts1malqu * sqrt( t(ContrastMat[z,])%*%diag(varmat[,i])%*%M%*%ContrastMat[z,] )
      lower.raw[z,i] <- t(ContrastMat[z,])%*%meanmat[,i] -
                        univarqu * sqrt( t(ContrastMat[z,])%*%diag(varmat[,i])%*%M%*%ContrastMat[z,] )
    }
  }
}

list(estimate=estimate, lower.raw=lower.raw, upper.raw=upper.raw, lower=lower, upper=upper,
     CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defr, 
     ContrastMat=ContrastMat, alternative=alternative, conf.level=conf.level)


}
