SimTestRatHet <-
function(trlist, grp, ntr, nep, ssmat, Num.Contrast, Den.Contrast, ncomp, alternative, 
                          Margin, meanmat, CorrMatDat) {


varmat  <- do.call(rbind, lapply(trlist, function(x) apply(x, 2, var)))
if (any(meanmat<0)) {
  warning("At least one sample mean is negative; check whether the test direction", "\n",
      "is still correct", "\n")
}
estimate <- Num.Contrast%*%meanmat/(Den.Contrast%*%meanmat)

defrmat <- matrix(nrow=ncomp, ncol=nep)
for (j in 1:nep) {
  defrmat[,j] <- DfSattRat(n=ssmat[,1],sd=sqrt(varmat[,j]),Num.Contrast=Num.Contrast,
                           Den.Contrast=Den.Contrast,Margin=Margin[,j])
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
        Rpart[i,h] <- ( t(Num.Contrast[z,]-Margin[z,i]*Den.Contrast[z,])%*%
                        diag( sapply(CovMatDat, function(x) x[i,h]) )%*%M%*%
                         (Num.Contrast[w,]-Margin[w,h]*Den.Contrast[w,]) ) /
                      sqrt( ( t(Num.Contrast[z,]-Margin[z,i]*Den.Contrast[z,])%*%diag(varmat[,i])%*%M%*%
                               (Num.Contrast[z,]-Margin[z,i]*Den.Contrast[z,]) ) *
                            ( t(Num.Contrast[w,]-Margin[w,h]*Den.Contrast[w,])%*%diag(varmat[,h])%*%M%*%
                               (Num.Contrast[w,]-Margin[w,h]*Den.Contrast[w,]) ) )
      }
    }
    Rrow <- cbind(Rrow,Rpart)
  }
  R <- rbind(R, Rrow)                                               # correlation matrix for multi-t
}
diag(R) <- 1

statistic <- matrix(nrow=ncomp, ncol=nep)                           # matrix of test stats
for (z in 1:ncomp) {
  for (i in 1:nep) {
    statistic[z,i] <- ( t(Num.Contrast[z,]-Margin[z,i]*Den.Contrast[z,])%*%meanmat[,i] ) /
                      sqrt( t(Num.Contrast[z,]-Margin[z,i]*Den.Contrast[z,])%*%diag(varmat[,i])%*%M%*%
                             (Num.Contrast[z,]-Margin[z,i]*Den.Contrast[z,]) )
  }
}
p.val <- SimTestP(ncomp=ncomp,nep=nep,alternative=alternative,statistic=statistic,
                  defr.mul=defr,defr.uni=defrmat,R=R)
p.val.adj <- p.val$p.val.adj; p.val.raw <- p.val$p.val.raw

list(estimate=estimate, statistic=statistic, p.val.raw=p.val.raw, p.val.adj=p.val.adj,
     CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defr,
     Num.Contrast=Num.Contrast, Den.Contrast=Den.Contrast, alternative=alternative)


}
