SimTestDiff.default <-
function(data, grp, resp=NULL, na.action="na.error", type="Dunnett", base=1,
                                ContrastMat=NULL, alternative="two.sided", Margin=0,
                                covar.equal=FALSE, CorrMatDat=NULL, ...) {


check.out <- SimCheck(data=data, grp=grp, resp=resp, na.action=na.action, type=type, base=base, 
                      ContrastMat=ContrastMat, Num.Contrast=NULL, Den.Contrast=NULL, 
                      alternative=alternative, Margin=Margin, covar.equal=covar.equal,
                      conf.level=NULL, CorrMatDat=CorrMatDat)

test <- if (covar.equal==TRUE) {
  if (na.action=="multi.df") "SimTestDiffHomMiva" else "SimTestDiffHom"
} else {
  if (na.action=="multi.df") stop("Procedure not yet available") else "SimTestDiffHet"
}
out <- do.call(test,
               list(trlist=check.out$trlist, grp=grp, ntr=check.out$ntr, 
                    nep=check.out$nep, ssmat=check.out$ssmat, 
                    ContrastMat=check.out$ContrastMat, ncomp=check.out$ncomp, 
                    alternative=alternative, Margin=check.out$Margin,
                    meanmat=check.out$meanmat, CorrMatDat=check.out$CorrMatDat))
out$type <- check.out$type
out$Margin <- check.out$Margin
out$test.class <- "differences"
out$covar.equal <- covar.equal
out$comp.names <- check.out$comp.namesDiff
out$resp <- check.out$resp
out$na.action <- na.action
rownames(out$degr.fr) <- rownames(out$statistic) <- rownames(out$p.val.raw) <- rownames(out$p.val.adj) <- 
  rownames(out$Margin) <- check.out$comp.namesDiff
colnames(out$degr.fr) <- colnames(out$statistic) <- colnames(out$p.val.raw) <- colnames(out$p.val.adj) <- 
  colnames(check.out$Margin) <- check.out$resp
rownames(out$CorrMatComp) <- colnames(out$CorrMatComp) <- rep(check.out$resp,times=check.out$ncomp)
if (na.action=="multi.df") {
  rownames(out$estimate) <- check.out$comp.namesDiff
  colnames(out$estimate) <- check.out$resp
}
if (covar.equal==FALSE) {
  names(out$CovMatDat) <- check.out$tr.names
  if (!inherits(check.out$CorrMatDat, "UserMatrix")) {
    names(out$CorrMatDat) <- check.out$tr.names
  }
}
if ( (inherits(check.out$CorrMatDat, "UserMatrix")) & (covar.equal==TRUE) ) {
  rownames(out$CovMatDat) <- colnames(out$CovMatDat) <- check.out$resp
}
class(out) <- "SimTest"
return(out)


}
