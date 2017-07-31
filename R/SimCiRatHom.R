SimCiRatHom <-
function(trlist, grp, ntr, nep, ssmat, Num.Contrast, Den.Contrast, ncomp, alternative, 
                        conf.level, meanmat, CorrMatDat) {


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
        Rpart[i,h] <- CorrMatDat[i,h] * ( t(Num.Contrast[z,]-estimate[z,i]*Den.Contrast[z,])%*%M%*%
                                           (Num.Contrast[w,]-estimate[w,h]*Den.Contrast[w,]) ) /
                      sqrt( (t(Num.Contrast[z,]-estimate[z,i]*Den.Contrast[z,])%*%M%*%
                              (Num.Contrast[z,]-estimate[z,i]*Den.Contrast[z,])) * 
                            (t(Num.Contrast[w,]-estimate[w,h]*Den.Contrast[w,])%*%M%*%
                              (Num.Contrast[w,]-estimate[w,h]*Den.Contrast[w,])) )
      }
    }
    Rrow <- cbind(Rrow,Rpart)
  }
  R <- rbind(R, Rrow)                                               # correlation matrix for multi-t
}
diag(R) <- 1

Azi     <- Bzi     <- Czi     <- Discrimi     <- lower     <- upper     <- matrix(nrow=ncomp,ncol=nep)
Azi.raw <- Bzi.raw <- Czi.raw <- Discrimi.raw <- lower.raw <- upper.raw <- matrix(nrow=ncomp,ncol=nep)
NSD <- 0

if (alternative=="greater") {
  lo1malqu <- qmvt(conf.level,tail="lower.tail",df=as.integer(defr[1,1]),corr=R)$quantile
  univarqu <- qt(p=conf.level, df=as.integer(defr[1,1]))
  for (z in 1:ncomp) {
    for (i in 1:nep) {
      Azi[z,i] <- ( t(Den.Contrast[z,])%*%meanmat[,i] )^2 - 
                         lo1malqu^2 * diag(CovMatDat)[i] * ( t(Den.Contrast[z,])%*%M%*%Den.Contrast[z,] )
      Bzi[z,i] <- -2 * ( (t(Num.Contrast[z,])%*%meanmat[,i]) * (t(Den.Contrast[z,])%*%meanmat[,i]) - 
                         lo1malqu^2 * diag(CovMatDat)[i] * ( t(Num.Contrast[z,])%*%M%*%Den.Contrast[z,] ) )
      Czi[z,i] <- ( t(Num.Contrast[z,])%*%meanmat[,i] )^2 - 
                         lo1malqu^2 * diag(CovMatDat)[i] * ( t(Num.Contrast[z,])%*%M%*%Num.Contrast[z,] )
      Discrimi[z,i] <- Bzi[z,i]^2 - 4*Azi[z,i]*Czi[z,i]
      if ( (Azi[z,i]>0) & (Discrimi[z,i]>=0) ) {
        upper[z,i] <- Inf
        lower[z,i] <- (-Bzi[z,i]-sqrt(Discrimi[z,i])) / (2*Azi[z,i])
      } else {
        upper[z,i] <- Inf; lower[z,i] <- -Inf; NSD <- NSD+1
      }
      Azi.raw[z,i] <- ( t(Den.Contrast[z,])%*%meanmat[,i] )^2 - 
                             univarqu^2 * diag(CovMatDat)[i] * ( t(Den.Contrast[z,])%*%M%*%Den.Contrast[z,] )
      Bzi.raw[z,i] <- -2 * ( (t(Num.Contrast[z,])%*%meanmat[,i]) * (t(Den.Contrast[z,])%*%meanmat[,i]) - 
                             univarqu^2 * diag(CovMatDat)[i] * ( t(Num.Contrast[z,])%*%M%*%Den.Contrast[z,] ) )
      Czi.raw[z,i] <- ( t(Num.Contrast[z,])%*%meanmat[,i] )^2 - 
                             univarqu^2 * diag(CovMatDat)[i] * ( t(Num.Contrast[z,])%*%M%*%Num.Contrast[z,] )
      Discrimi.raw[z,i] <- Bzi.raw[z,i]^2 - 4*Azi.raw[z,i]*Czi.raw[z,i]
      if ( (Azi.raw[z,i]>0) & (Discrimi.raw[z,i]>=0) ) {
        upper.raw[z,i] <- Inf
        lower.raw[z,i] <- (-Bzi.raw[z,i]-sqrt(Discrimi.raw[z,i])) / (2*Azi.raw[z,i])
      } else {
        upper.raw[z,i] <- Inf; lower.raw[z,i] <- -Inf
      }
    }
  }
}
if (alternative=="less") {
  up1malqu <- qmvt(conf.level,tail="upper.tail",df=as.integer(defr[1,1]),corr=R)$quantile
  univarqu <- qt(p=1-conf.level, df=as.integer(defr[1,1]))
  for (z in 1:ncomp) {
    for (i in 1:nep) {
      Azi[z,i] <- ( t(Den.Contrast[z,])%*%meanmat[,i] )^2 - 
                         up1malqu^2 * diag(CovMatDat)[i] * ( t(Den.Contrast[z,])%*%M%*%Den.Contrast[z,] )
      Bzi[z,i] <- -2 * ( (t(Num.Contrast[z,])%*%meanmat[,i]) * (t(Den.Contrast[z,])%*%meanmat[,i]) - 
                         up1malqu^2 * diag(CovMatDat)[i] * ( t(Num.Contrast[z,])%*%M%*%Den.Contrast[z,] ) )
      Czi[z,i] <- ( t(Num.Contrast[z,])%*%meanmat[,i] )^2 - 
                         up1malqu^2 * diag(CovMatDat)[i] * ( t(Num.Contrast[z,])%*%M%*%Num.Contrast[z,] )
      Discrimi[z,i] <- Bzi[z,i]^2 - 4*Azi[z,i]*Czi[z,i]
      if ( (Azi[z,i]>0) & (Discrimi[z,i]>=0) ) {
        upper[z,i] <- (-Bzi[z,i]+sqrt(Discrimi[z,i])) / (2*Azi[z,i])
        lower[z,i] <- -Inf
      } else {
        upper[z,i] <- Inf; lower[z,i] <- -Inf; NSD <- NSD+1
      }
      Azi.raw[z,i] <- ( t(Den.Contrast[z,])%*%meanmat[,i] )^2 - 
                             univarqu^2 * diag(CovMatDat)[i] * ( t(Den.Contrast[z,])%*%M%*%Den.Contrast[z,] )
      Bzi.raw[z,i] <- -2 * ( (t(Num.Contrast[z,])%*%meanmat[,i]) * (t(Den.Contrast[z,])%*%meanmat[,i]) - 
                             univarqu^2 * diag(CovMatDat)[i] * ( t(Num.Contrast[z,])%*%M%*%Den.Contrast[z,] ) )
      Czi.raw[z,i] <- ( t(Num.Contrast[z,])%*%meanmat[,i] )^2 - 
                             univarqu^2 * diag(CovMatDat)[i] * ( t(Num.Contrast[z,])%*%M%*%Num.Contrast[z,] )
      Discrimi.raw[z,i] <- Bzi.raw[z,i]^2 - 4*Azi.raw[z,i]*Czi.raw[z,i]
      if ( (Azi.raw[z,i]>0) & (Discrimi.raw[z,i]>=0) ) {
        upper.raw[z,i] <- (-Bzi.raw[z,i]+sqrt(Discrimi.raw[z,i])) / (2*Azi.raw[z,i])
        lower.raw[z,i] <- -Inf
      } else {
        upper.raw[z,i] <- Inf; lower.raw[z,i] <- -Inf
      }
    }
  }
}
if (alternative=="two.sided") {
  ts1malqu <- qmvt(conf.level,tail="both.tails",df=as.integer(defr[1,1]),corr=R)$quantile
  univarqu <- qt(p=1-(1-conf.level)/2, df=as.integer(defr[1,1]))
  for (z in 1:ncomp) {
    for (i in 1:nep) {
      Azi[z,i] <- ( t(Den.Contrast[z,])%*%meanmat[,i] )^2 - 
                         ts1malqu^2 * diag(CovMatDat)[i] * ( t(Den.Contrast[z,])%*%M%*%Den.Contrast[z,] )
      Bzi[z,i] <- -2 * ( (t(Num.Contrast[z,])%*%meanmat[,i]) * (t(Den.Contrast[z,])%*%meanmat[,i]) - 
                         ts1malqu^2 * diag(CovMatDat)[i] * ( t(Num.Contrast[z,])%*%M%*%Den.Contrast[z,] ) )
      Czi[z,i] <- ( t(Num.Contrast[z,])%*%meanmat[,i] )^2 - 
                         ts1malqu^2 * diag(CovMatDat)[i] * ( t(Num.Contrast[z,])%*%M%*%Num.Contrast[z,] )
      Discrimi[z,i] <- Bzi[z,i]^2 - 4*Azi[z,i]*Czi[z,i]
      if ( (Azi[z,i]>0) & (Discrimi[z,i]>=0) ) {
        upper[z,i] <- (-Bzi[z,i]+sqrt(Discrimi[z,i])) / (2*Azi[z,i])
        lower[z,i] <- (-Bzi[z,i]-sqrt(Discrimi[z,i])) / (2*Azi[z,i])
      } else {
        upper[z,i] <- Inf; lower[z,i] <- -Inf; NSD <- NSD+1
      }
      Azi.raw[z,i] <- ( t(Den.Contrast[z,])%*%meanmat[,i] )^2 - 
                             univarqu^2 * diag(CovMatDat)[i] * ( t(Den.Contrast[z,])%*%M%*%Den.Contrast[z,] )
      Bzi.raw[z,i] <- -2 * ( (t(Num.Contrast[z,])%*%meanmat[,i]) * (t(Den.Contrast[z,])%*%meanmat[,i]) - 
                             univarqu^2 * diag(CovMatDat)[i] * ( t(Num.Contrast[z,])%*%M%*%Den.Contrast[z,] ) )
      Czi.raw[z,i] <- ( t(Num.Contrast[z,])%*%meanmat[,i] )^2 - 
                             univarqu^2 * diag(CovMatDat)[i] * ( t(Num.Contrast[z,])%*%M%*%Num.Contrast[z,] )
      Discrimi.raw[z,i] <- Bzi.raw[z,i]^2 - 4*Azi.raw[z,i]*Czi.raw[z,i]
      if ( (Azi.raw[z,i]>0) & (Discrimi.raw[z,i]>=0) ) {
        upper.raw[z,i] <- (-Bzi.raw[z,i]+sqrt(Discrimi.raw[z,i])) / (2*Azi.raw[z,i])
        lower.raw[z,i] <- (-Bzi.raw[z,i]-sqrt(Discrimi.raw[z,i])) / (2*Azi.raw[z,i])
      } else {
        upper.raw[z,i] <- Inf; lower.raw[z,i] <- -Inf
      }
    }
  }
}

list(estimate=estimate, NSD=NSD, lower.raw=lower.raw, upper.raw=upper.raw, 
     lower=lower, upper=upper, CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, 
     degr.fr=defr, Num.Contrast=Num.Contrast, Den.Contrast=Den.Contrast, 
     alternative=alternative, conf.level=conf.level)


}
