SimTestRatHom <-
function(trlist, grp, ntr, nep, ssmat, Num.Contrast, Den.Contrast, ncomp, alternative, 
                          Margin, meanmat, CorrMatDat) {


if (any(meanmat<0)) {
  warning("At least one sample mean is negative; check whether the test direction", "\n",
      "is still correct", "\n")
}
estimate <- Num.Contrast%*%meanmat/(Den.Contrast%*%meanmat)

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
R <- NULL
for (z in 1:ncomp) { 
  Rrow <- NULL
  for (w in 1:ncomp) {
    Rpart <- matrix(nrow=nep,ncol=nep)
    for (i in 1:nep) {
      for (h in 1:nep) {
        Rpart[i,h] <- CorrMatDat[i,h] * ( t(Num.Contrast[z,]-Margin[z,i]*Den.Contrast[z,])%*%
                                          M%*%(Num.Contrast[w,]-Margin[w,h]*Den.Contrast[w,]) ) /
                      sqrt( (t(Num.Contrast[z,]-Margin[z,i]*Den.Contrast[z,])%*%M%*%
                              (Num.Contrast[z,]-Margin[z,i]*Den.Contrast[z,])) * 
                            (t(Num.Contrast[w,]-Margin[w,h]*Den.Contrast[w,])%*%M%*%
                              (Num.Contrast[w,]-Margin[w,h]*Den.Contrast[w,])) )
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
                      sqrt( diag(CovMatDat)[i] * ( t(Num.Contrast[z,]-Margin[z,i]*Den.Contrast[z,])%*%M%*%
                                                    (Num.Contrast[z,]-Margin[z,i]*Den.Contrast[z,]) ) )
  }
}
p.val <- SimTestP(ncomp=ncomp,nep=nep,alternative=alternative,statistic=statistic,
                  defr.mul=defr,defr.uni=defr,R=R)
p.val.adj <- p.val$p.val.adj; p.val.raw <- p.val$p.val.raw

list(estimate=estimate, statistic=statistic, p.val.raw=p.val.raw, p.val.adj=p.val.adj,
     CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defr,
     Num.Contrast=Num.Contrast, Den.Contrast=Den.Contrast, alternative=alternative)


}
