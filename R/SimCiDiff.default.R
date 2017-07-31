SimCiDiff.default <-
function(data, grp, resp=NULL, na.action="na.error", type="Dunnett", base=1,
                              ContrastMat=NULL, alternative="two.sided", covar.equal=FALSE,
                              conf.level=0.95, CorrMatDat=NULL, ...) {


check.out <- SimCheck(data=data, grp=grp, resp=resp, na.action=na.action, type=type, base=base, 
                      ContrastMat=ContrastMat, Num.Contrast=NULL, Den.Contrast=NULL, 
                      alternative=alternative, Margin=NULL, covar.equal=covar.equal,
                      conf.level=conf.level, CorrMatDat=CorrMatDat)

test <- if (covar.equal==TRUE) {
  if (na.action=="multi.df") stop("Procedure not yet available") else "SimCiDiffHom"
} else {
  if (na.action=="multi.df") stop("Procedure not yet available") else "SimCiDiffHet"
}
out <- do.call(test,
               list(trlist=check.out$trlist, grp=grp, ntr=check.out$ntr, nep=check.out$nep, 
                    ssmat=check.out$ssmat, ContrastMat=check.out$ContrastMat, 
                    ncomp=check.out$ncomp, alternative=alternative, conf.level=conf.level,
                    meanmat=check.out$meanmat, CorrMatDat=check.out$CorrMatDat))
out$type <- check.out$type
out$test.class <- "differences"
out$covar.equal <- covar.equal
out$comp.names <- check.out$comp.namesDiff
out$resp <- check.out$resp
out$na.action <- na.action
rownames(out$degr.fr) <- rownames(out$lower.raw) <- rownames(out$upper.raw) <- rownames(out$lower) <- 
  rownames(out$upper) <- check.out$comp.namesDiff
colnames(out$degr.fr) <- colnames(out$lower.raw) <- colnames(out$upper.raw) <- colnames(out$lower) <- 
  colnames(out$upper) <- check.out$resp
rownames(out$CorrMatComp) <- colnames(out$CorrMatComp) <- rep(check.out$resp,times=check.out$ncomp)
if (na.action=="multi.df") {
  rownames(out$estimate) <- check.out$comp.namesDiff
  colnames(out$estimate) <- check.out$resp
}
if (covar.equal==FALSE) {
  names(out$CovMatDat) <- check.out$tr.names
  if (class(check.out$CorrMatDat)!="UserMatrix") {
    names(out$CorrMatDat) <- check.out$tr.names
  }
}
if ( (class(check.out$CorrMatDat)=="UserMatrix") & (covar.equal==TRUE) ) {
  rownames(out$CovMatDat) <- colnames(out$CovMatDat) <- check.out$resp
}
class(out) <- "SimCi"
return(out)


}
