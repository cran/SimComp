`SimTestRatHet` <-
function(trlist, grp, ntr, nep, ssvec, Num.Contrast, Den.Contrast, alternative, Margin) {


CMat <- Num.Contrast
DMat <- Den.Contrast
ncomp <- nrow(CMat)                                                      # number of comparisons

meanmat <- varmat <- matrix(nrow=ntr, ncol=nep)
for (i in 1:ntr) { for (j in 1:nep) {
  meanmat[i,j]=mean(trlist[[i]][,j]); varmat[i,j]=var(trlist[[i]][,j]) }}
estimate <- CMat%*%meanmat/(DMat%*%meanmat)

defrmat <- matrix(nrow=ncomp, ncol=nep)
for (j in 1:nep) { for (z in 1:ncomp) {
defrmat[z,j]=( (sum((CMat[z,]-Margin[z,j]*DMat[z,])^2*varmat[,j]/ssvec))^2 ) / 
             sum( ( (CMat[z,]-Margin[z,j]*DMat[z,])^4*varmat[,j]^2 ) / ( ssvec^2*(ssvec-1) ) ) }}
defrmat[defrmat<2] <- 2                                                  # to be well-defined

CovMatDat <- CorrMatDat <- list()                                        # list of covariance/ correlation matrices of the data
for (i in 1:ntr) { CovMatDat[[i]]  <- cov(trlist[[i]])
                   CorrMatDat[[i]] <- cov2cor(CovMatDat[[i]]) }

M <- diag(1/ssvec)
R <- NULL
for (z in 1:ncomp) {
  Rrow <- NULL
  for (w in 1:ncomp) {
    Rpart <- matrix(nrow=nep,ncol=nep)
    for (i in 1:nep) { for (h in 1:nep) {
      Rpart[i,h]=( t(CMat[z,]-Margin[z,i]*DMat[z,])%*%
                   diag(unlist( lapply( X=CovMatDat,FUN=function(x){x[i,h]} ) ))%*%M%*%(CMat[w,]-Margin[w,h]*DMat[w,]) ) /
                 sqrt( ( t(CMat[z,]-Margin[z,i]*DMat[z,])%*%diag(varmat[,i])%*%M%*%(CMat[z,]-Margin[z,i]*DMat[z,]) ) *
                       ( t(CMat[w,]-Margin[w,h]*DMat[w,])%*%diag(varmat[,h])%*%M%*%(CMat[w,]-Margin[w,h]*DMat[w,]) ) ) }
    }
    Rrow <- cbind(Rrow,Rpart)
  }
  R <- rbind(R, Rrow)                                                    # correlation matrix for test.stat
}
diag(R)=rep(1,times=nrow(R))

test.stat <- p.val.adj <- p.val.raw <- matrix(nrow=ncomp, ncol=nep)      # matrices of test statistics and p.vals
for (z in 1:ncomp) { for (i in 1:nep) {
  test.stat[z,i]=( t(CMat[z,]-Margin[z,i]*DMat[z,])%*%meanmat[,i] ) /
                 sqrt( t(CMat[z,]-Margin[z,i]*DMat[z,])%*%diag(varmat[,i])%*%M%*%(CMat[z,]-Margin[z,i]*DMat[z,]) )
  if (alternative=="greater") {
    p.val.adj[z,i]=1-pmvt(lower=-Inf,upper=rep(test.stat[z,i],times=ncomp*nep),df=defrmat[z,i],corr=R)[1]
    p.val.raw[z,i]=pt(q=test.stat[z,i],df=defrmat[z,i],lower.tail=FALSE) }
  if (alternative=="less") {
    p.val.adj[z,i]=1-pmvt(lower=rep(test.stat[z,i],times=ncomp*nep),upper=Inf,df=defrmat[z,i],corr=R)[1]
    p.val.raw[z,i]=pt(q=test.stat[z,i],df=defrmat[z,i],lower.tail=TRUE) }
  if (alternative=="two.sided") {
    p.val.adj[z,i]=1-pmvt(lower=rep(-abs(test.stat[z,i]),times=ncomp*nep),upper=rep(abs(test.stat[z,i]),times=ncomp*nep),
                   df=defrmat[z,i],corr=R)[1]
    p.val.raw[z,i]=min(pt(q=abs(test.stat[z,i]),df=defrmat[z,i],lower.tail=FALSE)*2,1) }
}}

list(estimate=estimate, statistic=test.stat, p.val.raw=p.val.raw, p.val.adj=p.val.adj,
     CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defrmat,
     Num.Contrast=CMat, Den.Contrast=DMat, Margin=Margin, alternative=alternative)

}

