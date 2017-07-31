SimTestDiffHomMiva <-
function(trlist, grp, ntr, nep, ssmat, ContrastMat, ncomp, alternative, 
                               Margin, meanmat, CorrMatDat) {


estimate <- matrix(nrow=ncomp, ncol=nep)
for (j in 1:nep) {
  estimate[,j] <- ContrastMat[[j]]%*%meanmat[,j] }

varmat <- matrix(nrow=ntr, ncol=nep)                                # variance per treatment and endpoint
for (i in 1:ntr) { for (j in 1:nep) {
  varmat[i,j] <- var( na.omit(trlist[[i]][,j]) ) }}
varvec <- numeric(nep)                                              # variance per endpoint
for (j in 1:nep) {
  varvec[j] <- sum( (ssmat[,j]-1)*varmat[,j] ) / ( sum(ssmat[,j])-ntr )
}

trlist.omit <- lapply(trlist, na.omit)

ssvec.omit <- numeric(ntr)                                          # reduced sample sizes
for (i in 1:ntr) {
  ssvec.omit[i] <- nrow(trlist.omit[[i]])
}

defrvec <- numeric(nep)                                             # df per endpoint
for (j in 1:nep) { defrvec[j] <- sum(ssmat[,j])-ntr }
defr <- matrix(defrvec, nrow=ncomp, ncol=nep, byrow=TRUE)           # matrix of dfs

CovMatDat <- matrix(0,nrow=nep,ncol=nep)                            # covariance matrix of the reduced data
for (i in 1:ntr) {
  CovMatDat <- CovMatDat+(ssvec.omit[i]-1)*cov(trlist.omit[[i]])
}
CovMatDat <- CovMatDat/(sum(ssvec.omit)-ntr)
if (is.null(CorrMatDat)) {
  CorrMatDat <- cov2cor(CovMatDat)                                  # correlation matrix of the data
} else {
  sdmat <- sqrt( diag( diag(CovMatDat),nrow=nep ) )                 # sds on the diagonal
  CovMatDat <- sdmat%*%CorrMatDat%*%sdmat
}

sslist <- list()                                                    # list of sample sizes for pairs of endpoints
for (i in 1:ntr) {
  sslist[[i]] <- matrix(nrow=nep, ncol=nep)
  for (k in 1:nep) { for (m in 1:nep) {
    sslist[[i]][k,m] <- nrow( na.omit(trlist[[i]][,c(k,m)]) )
  }}
}

R <- NULL
for (z in 1:ncomp) { 
  Rrow <- NULL
  for (w in 1:ncomp) {
    Rpart <- matrix(nrow=nep,ncol=nep)
    for (i in 1:nep) {
      for (h in 1:nep) {
        Rpart[i,h] <- CorrMatDat[i,h] * 
                      (t(ContrastMat[[i]][z,])%*%
                       diag( sapply(sslist, function(x) x[i,h])/(ssmat[,i]*ssmat[,h]) )%*%
                       ContrastMat[[h]][w,]) /
                      sqrt( (t(ContrastMat[[i]][z,])%*% diag(1/ssmat[,i]) %*%ContrastMat[[i]][z,]) * 
                             ((ContrastMat[[h]][w,])%*% diag(1/ssmat[,h]) %*%ContrastMat[[h]][w,]) )
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
    statistic[z,i] <- (  (t(ContrastMat[[i]][z,])%*%meanmat[,i]) - Margin[z,i]  ) /
                      sqrt( varvec[i]*( t(ContrastMat[[i]][z,])%*%diag(1/ssmat[,i])%*%
                                        ContrastMat[[i]][z,] ) )
  }
}
p.val <- SimTestP(ncomp=ncomp,nep=nep,alternative=alternative,statistic=statistic,
                  defr.mul=defr,defr.uni=defr,R=R)
p.val.adj <- p.val$p.val.adj; p.val.raw <- p.val$p.val.raw

list(estimate=estimate, statistic=statistic, p.val.raw=p.val.raw, p.val.adj=p.val.adj,
     CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defr,
     ContrastMat=ContrastMat, alternative=alternative)


}
