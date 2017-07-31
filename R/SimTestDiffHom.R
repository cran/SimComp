SimTestDiffHom <-
function(trlist, grp, ntr, nep, ssmat, ContrastMat, ncomp, alternative, 
                           Margin, meanmat, CorrMatDat) {


estimate <- ContrastMat%*%meanmat

defr <- matrix(sum(ssmat[,1])-ntr, nrow=ncomp, ncol=nep)            # matrix of dfs

CovMatDat <- Reduce("+",                                            # covariance matrix of the data
                    lapply(trlist, function(x) (nrow(x)-1)*cov(x))
                    )/defr[1,1]
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

statistic <- matrix(nrow=ncomp, ncol=nep)                           # matrix of test stats
for (z in 1:ncomp) {
  for (i in 1:nep) {
    statistic[z,i] <- (  (t(ContrastMat[z,])%*%meanmat[,i]) - Margin[z,i]  ) /
                      sqrt( diag(CovMatDat)[i]*( t(ContrastMat[z,])%*%M%*%ContrastMat[z,] ) )
  }
}

p.val <- SimTestP(ncomp=ncomp,nep=nep,alternative=alternative,statistic=statistic,
                  defr.mul=defr,defr.uni=defr,R=R)
p.val.adj <- p.val$p.val.adj; p.val.raw <- p.val$p.val.raw

list(estimate=estimate, statistic=statistic, p.val.raw=p.val.raw, p.val.adj=p.val.adj,
     CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defr,
     ContrastMat=ContrastMat, alternative=alternative)


}
