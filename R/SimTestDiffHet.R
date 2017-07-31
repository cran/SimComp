SimTestDiffHet <-
function(trlist, grp, ntr, nep, ssmat, ContrastMat, ncomp, alternative, 
                           Margin, meanmat, CorrMatDat) {


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

statistic <- matrix(nrow=ncomp, ncol=nep)                           # matrix of test stats
for (z in 1:ncomp) {
  for (i in 1:nep) {
    statistic[z,i] <- ( (t(ContrastMat[z,])%*%meanmat[,i]) - Margin[z,i] ) /
                      sqrt( t(ContrastMat[z,])%*%diag(varmat[,i])%*%M%*%(ContrastMat[z,]) )
  }
}
p.val <- SimTestP(ncomp=ncomp,nep=nep,alternative=alternative,statistic=statistic,
                  defr.mul=defr,defr.uni=defrmat,R=R)
p.val.adj <- p.val$p.val.adj; p.val.raw <- p.val$p.val.raw

list(estimate=estimate, statistic=statistic, p.val.raw=p.val.raw, p.val.adj=p.val.adj,
     CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defr,
     ContrastMat=ContrastMat, alternative=alternative)


}
