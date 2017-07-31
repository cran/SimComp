SimTestRat.default <-
function(data, grp, resp=NULL, na.action="na.error", type="Dunnett", base=1,
                       Num.Contrast=NULL,Den.Contrast=NULL, alternative="two.sided", Margin=1,
                       covar.equal=FALSE, CorrMatDat=NULL, ...) {


check.out <- SimCheck(data=data, grp=grp, resp=resp, na.action=na.action, type=type, base=base, 
                      ContrastMat=NULL, Num.Contrast=Num.Contrast, Den.Contrast=Den.Contrast, 
                      alternative=alternative, Margin=Margin, covar.equal=covar.equal,
                      conf.level=NULL, CorrMatDat=CorrMatDat)



test <- if (covar.equal==TRUE) {
  if (na.action=="multi.df") stop("Procedure not yet available") else "SimTestRatHom"
} else {
  if (na.action=="multi.df") stop("Procedure not yet available") else "SimTestRatHet"
}
out <- do.call(test,
               list(trlist=check.out$trlist, grp=grp, ntr=check.out$ntr, 
                    nep=check.out$nep, ssmat=check.out$ssmat, 
                    Num.Contrast=check.out$Num.Contrast, Den.Contrast=check.out$Den.Contrast, 
                    ncomp=check.out$ncomp, alternative=alternative, Margin=check.out$Margin,
                    meanmat=check.out$meanmat, CorrMatDat=check.out$CorrMatDat))
out$type <- check.out$type
out$Margin <- check.out$Margin
out$test.class <- "ratios"
out$covar.equal <- covar.equal
out$comp.names <- check.out$comp.namesRat
out$resp <- check.out$resp
out$na.action <- na.action
rownames(out$degr.fr) <- rownames(out$statistic) <- rownames(out$p.val.raw) <- rownames(out$p.val.adj) <- 
  rownames(out$Margin) <- check.out$comp.namesRat
colnames(out$degr.fr) <- colnames(out$statistic) <- colnames(out$p.val.raw) <- colnames(out$p.val.adj) <- 
  colnames(check.out$Margin) <- check.out$resp
rownames(out$CorrMatComp) <- colnames(out$CorrMatComp) <- rep(check.out$resp,times=check.out$ncomp)
if (na.action=="multi.df") {
  rownames(out$estimate) <- check.out$comp.namesRat
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
class(out) <- "SimTest"
return(out)


}
