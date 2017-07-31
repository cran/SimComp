SimTestP <-
function(ncomp, nep, alternative, statistic, defr.mul, defr.uni, R) {


p.val.adj <- p.val.raw <- matrix(nrow=ncomp, ncol=nep)              # matrices of p.vals
for (z in 1:ncomp) {
  for (i in 1:nep) {
    if (alternative=="greater") {
      p.val.adj[z,i] <- 1-pmvt(lower=-Inf,
                               upper=rep(statistic[z,i],times=ncomp*nep),
                               df=as.integer(defr.mul[z,i]),corr=R)[1]
      p.val.raw[z,i] <- pt(q=statistic[z,i],df=as.integer(defr.uni[z,i]),lower.tail=FALSE)
    }
    if (alternative=="less") {
      p.val.adj[z,i] <- 1-pmvt(lower=rep(statistic[z,i],times=ncomp*nep),
                               upper=Inf,
                               df=as.integer(defr.mul[z,i]),corr=R)[1]
      p.val.raw[z,i] <- pt(q=statistic[z,i],df=as.integer(defr.uni[z,i]),lower.tail=TRUE)
    }
    if (alternative=="two.sided") {
      p.val.adj[z,i] <- 1-pmvt(lower=rep(-abs(statistic[z,i]),times=ncomp*nep),
                               upper=rep( abs(statistic[z,i]),times=ncomp*nep),
                               df=as.integer(defr.mul[z,i]),corr=R)[1]
      p.val.raw[z,i] <- min(pt(q=abs(statistic[z,i]),df=as.integer(defr.uni[z,i]),
                               lower.tail=FALSE)*2,1)
    }
  }
}

list(p.val.adj=p.val.adj, p.val.raw=p.val.raw)


}
