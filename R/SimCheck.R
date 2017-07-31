SimCheck <-
function(data, grp, resp, na.action, type, base, 
                     ContrastMat=NULL, Num.Contrast=NULL, Den.Contrast=NULL, 
                     alternative, Margin=NULL, covar.equal, conf.level=NULL,
                     CorrMatDat) {


if (length(grp) > 1) {
  stop("Specify only one grouping variable")
}
tr.names <- levels(data[,grp])
ntr <- length(tr.names)                                             # number of treatments
if (ntr<2) {
  stop("Number of levels of grp must be greater than 1")
}

if (is.null(resp)) {
  resp <- setdiff(names(data), grp)
} 
nep <- length(resp)                                                 # number of endpoints

trlist <- lapply(split(x=data,f=data[,grp]),                        # splitted data
                 function(x) as.matrix(x[,resp]))                   # ... without factor
if (!is.numeric(unlist(trlist))) {
  stop("Response variables must be numeric")
}

ssmat <- matrix(nrow=ntr, ncol=nep)                                 # sample size per treatment and endpoint
for (i in 1:ntr) { 
  ssmat[i,] <- apply(trlist[[i]], 2, function(x) length(na.omit(x))) }
rownames(ssmat) <- tr.names; colnames(ssmat) <- resp

if (is.null(CorrMatDat)) {                                          # CorrMatDat not given
  if (covar.equal==TRUE) {
    if ( any(apply(ssmat, 2, function(x) sum(x-1)<nep)) ) {
      stop("Not enough observations or too many endpoints")
    }
  } else {
    if ( any((ssmat-1)<nep) ) {
      stop("Not enough observations or too many endpoints")
    }
  }
} else {                                                            # CorrMatDat given
  if (any(dim(CorrMatDat)<2)) {
    stop("Number of rows and columns of CorrMatDat must be greater than 1")
  }
  if (nrow(CorrMatDat)!=nep || ncol(CorrMatDat)!=nep || nrow(CorrMatDat)!=ncol(CorrMatDat)) {
    stop("Number of rows and columns of CorrMatDat and number of endpoints must be equal")
  }
  if (any(CorrMatDat>1) || any(CorrMatDat<(-1)) || any(diag(CorrMatDat)!=1)) {
    stop("CorrMatDat is not a correlation matrix")
  }
  rownames(CorrMatDat) <- colnames(CorrMatDat) <- resp
  class(CorrMatDat) <- "UserMatrix"
}

na.action <- match.arg(na.action, choices=c("na.error","multi.df"))
if (na.action=="na.error") {
  if (NA %in% unlist(trlist)) {
    stop("Missing values within the data; na.action=\"multi.df\" could be appropriate")
  }
}

if (is.null(ContrastMat) & is.null(Num.Contrast) & is.null(Den.Contrast)) { # all contrast matrices missing
  type <- match.arg(type, choices=c("Dunnett","Tukey","Sequen","AVE","GrandMean","Changepoint",
                                    "Marcus","McDermott","Williams","UmbrellaWilliams"))
  ContrastMat <- contrMat(n=ssmat[,1], type=type, base=base)
  ContrastMatRat <- contrMatRatio(n=ssmat[,1], type=type, base=base)
  Num.Contrast <- ContrastMatRat$numC
  Den.Contrast <- ContrastMatRat$denC
  ncomp <- nrow(ContrastMat)                                        # number of comparisons
  comp.namesDiff <- rownames(ContrastMat)
  comp.namesRat <- ContrastMatRat$rnames
} else {                                                            # at least one contrast matrix given
  type <- "User defined"
  if (!is.null(ContrastMat)) {                                      # ContrastMat given
    if (ncol(ContrastMat)!=ntr) {
      stop("Number of columns of ContrastMat and number of groups must be equal")
    }
    if (any(apply(ContrastMat==0, 1, all))) {
      stop("At least one row of ContrastMat is a vector with all components", 
           "\n", "equal to zero")
    }
    ncomp <- nrow(ContrastMat)                                      # number of comparisons
    if (is.null(rownames(ContrastMat))) {
      rownames(ContrastMat) <- paste("C", 1:nrow(ContrastMat), sep="")
    }
    colnames(ContrastMat) <- tr.names
    comp.namesDiff <- rownames(ContrastMat)
    comp.namesRat <- NULL
  } else {                                                          # ContrastMat not given
    if (!is.null(Num.Contrast) & is.null(Den.Contrast)) {
      stop("Num.Contrast is specified, but Den.Contrast is missing")
    }
    if (is.null(Num.Contrast) & !is.null(Den.Contrast)) {
      stop("Den.Contrast is specified, but Num.Contrast is missing")
    }
    if (!is.null(Num.Contrast) & !is.null(Den.Contrast)) {
      if (nrow(Num.Contrast)!=nrow(Den.Contrast)) {
        stop("Number of rows of Num.Contrast and Den.Contrast must be equal")
      }
      if (ncol(Num.Contrast)!=ntr | ncol(Den.Contrast)!=ntr) {
        stop("Number of columns of Num.Contrast and Den.Contrast and number of groups must be equal")
      }
      if (any(c(apply(Num.Contrast==0, 1, all),apply(Den.Contrast==0, 1, all)))) {
        stop("At least one row of Num.Contrast or Den.Contrast is a vector with all components", "\n",
             "equal to zero")
      }
      ncomp <- nrow(Num.Contrast)                                   # number of comparisons
      if (is.null(rownames(Num.Contrast)) && is.null(rownames(Den.Contrast))) {
        comp.namesRat <- paste("C", 1:nrow(Num.Contrast), sep="")
      } else {
        if (any(rownames(Num.Contrast)!=rownames(Den.Contrast))) {
          comp.namesRat <- paste(rownames(Num.Contrast), rownames(Den.Contrast), sep="/")
        } else {
          comp.namesRat <- rownames(Num.Contrast)
        }
      }
    }
    colnames(Num.Contrast) <- colnames(Den.Contrast) <- tr.names
    comp.namesDiff <- NULL
  }
}

if (na.action=="multi.df") {
  ContrastMatList <- Num.ContrastList <- Den.ContrastList <- list()
  if (type=="User defined") {
    for (j in 1:nep) {
      ContrastMatList[[j]] <- ContrastMat
      Num.ContrastList[[j]] <- Num.Contrast; Den.ContrastList[[j]] <- Den.Contrast
    }
  } else {
    for (j in 1:nep) {
      ContrastMatList[[j]] <- contrMat(n=ssmat[,j], type=type, base=base)
      ContrastMatRat <- contrMatRatio(n=ssmat[,j], type=type, base=base)
      Num.ContrastList[[j]] <- ContrastMatRat$numC
      Den.ContrastList[[j]] <- ContrastMatRat$denC
    }
  }
  names(ContrastMatList) <- names(Num.ContrastList) <- names(Den.ContrastList) <- resp
  ContrastMat <- ContrastMatList
  Num.Contrast <- Num.ContrastList; Den.Contrast <- Den.ContrastList
}

alternative <- match.arg(alternative, choices=c("two.sided","less","greater"))

if (is.numeric(Margin)) {
  if (length(Margin)==1) {
    Margin <- matrix(Margin,nrow=ncomp,ncol=nep)
  }
  if (is.vector(Margin)) {
    if (length(Margin)!=nep) {
      stop("Margin must be a single numeric value, or a numeric vector with length equal to the 
           number of endpoints,", "\n", "or a matrix with number of rows equal to the number of 
           contrasts and number of columns equal to the number of endpoints")
    } else {
      Margin <- matrix(Margin,nrow=ncomp,ncol=nep,byrow=TRUE)
    }
  }
  if (is.matrix(Margin)) {
    if (nrow(Margin)!=ncomp) {
      stop("Number of rows of Margin and number of contrasts must be equal")
    }
    if (ncol(Margin)!=nep) {
      stop("Number of columns of Margin and number of endpoints must be equal")
    }
  }
}

if (!is.null(conf.level)) {
  if (length(conf.level)>1 || !is.numeric(conf.level) || conf.level>=1 || conf.level<=0 ) {
    stop("conf.level must be a single numeric value between 0 and 1")
  }
}

meanmat <- do.call(rbind, lapply(trlist, function(x) colMeans(x,na.rm=TRUE)))

list(trlist=trlist,tr.names=tr.names,ntr=ntr,resp=resp,nep=nep,ssmat=ssmat,type=type,
     ContrastMat=ContrastMat,Num.Contrast=Num.Contrast,Den.Contrast=Den.Contrast,ncomp=ncomp,
     comp.namesDiff=comp.namesDiff,comp.namesRat=comp.namesRat,Margin=Margin,meanmat=meanmat,
     CorrMatDat=CorrMatDat)


}
