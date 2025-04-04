SimCiRat.default <-
function(data, grp, resp=NULL, na.action="na.error", type="Dunnett", base=1,
                             Num.Contrast=NULL, Den.Contrast=NULL, alternative="two.sided", covar.equal=FALSE,
                             conf.level=0.95, CorrMatDat=NULL, ...) {


check.out <- SimCheck(data=data, grp=grp, resp=resp, na.action=na.action, type=type, base=base, 
                      ContrastMat=NULL, Num.Contrast=Num.Contrast, Den.Contrast=Den.Contrast, 
                      alternative=alternative, Margin=NULL, covar.equal=covar.equal,
                      conf.level=conf.level, CorrMatDat=CorrMatDat)

test <- if (covar.equal==TRUE) {
  if (na.action=="multi.df") stop("Procedure not yet available") else "SimCiRatHom"
} else {
  if (na.action=="multi.df") stop("Procedure not yet available") else "SimCiRatHet"
}
out <- do.call(test,
               list(trlist=check.out$trlist, grp=grp, ntr=check.out$ntr, 
                    nep=check.out$nep, ssmat=check.out$ssmat, 
                    Num.Contrast=check.out$Num.Contrast, Den.Contrast=check.out$Den.Contrast, 
                    ncomp=check.out$ncomp, alternative=alternative, conf.level=conf.level,
                    meanmat=check.out$meanmat, CorrMatDat=check.out$CorrMatDat))
out$type <- check.out$type
out$test.class <- "ratios"
out$covar.equal <- covar.equal
out$comp.names <- check.out$comp.namesRat
out$resp <- check.out$resp
out$na.action <- na.action
rownames(out$degr.fr) <- rownames(out$lower.raw) <- rownames(out$upper.raw) <- rownames(out$lower) <- 
  rownames(out$upper) <- check.out$comp.namesRat
colnames(out$degr.fr) <- colnames(out$lower.raw) <- colnames(out$upper.raw) <- colnames(out$lower) <- 
  colnames(out$upper) <- check.out$resp
rownames(out$CorrMatComp) <- colnames(out$CorrMatComp) <- rep(check.out$resp,times=check.out$ncomp)
if (na.action=="multi.df") {
  rownames(out$estimate) <- check.out$comp.namesRat
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
class(out) <- "SimCi"
return(out)


}
