DfSattRat <-
function(n, sd, type="Dunnett", base=1, Num.Contrast=NULL, Den.Contrast=NULL,
                      Margin=1) {


if (length(n)!=length(sd)) {
  stop("Lengths of n and sd must be equal")
}

if ( is.null(Num.Contrast)+is.null(Den.Contrast)<2 ) {              # if at most one matrix is missing
  type <- "User defined"
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
    if (ncol(Num.Contrast)!=length(n) | ncol(Den.Contrast)!=length(n)) {
      stop("Number of columns of Num.Contrast and Den.Contrast and number of groups must be equal")
    }
    if (any(c(apply(Num.Contrast==0, 1, all),apply(Den.Contrast==0, 1, all)))) {
      stop("At least one row of Num.Contrast or Den.Contrast is a vector with all components", "\n",
           "equal to zero")
    }
    if (is.null(rownames(Num.Contrast)) && is.null(rownames(Den.Contrast))) {
      comp.names <- paste("C", 1:nrow(Num.Contrast), sep="")
    } else {
      if (any(rownames(Num.Contrast)!=rownames(Den.Contrast))) {
        comp.names <- paste(rownames(Num.Contrast), rownames(Den.Contrast), 
                       sep="/")
      } else {
        comp.names <- rownames(Num.Contrast)
      }
    }
  }
} else {                                                            # if both matrices are missing
  type <- match.arg(type, choices=c("Dunnett","Tukey","Sequen","AVE","GrandMean","Changepoint",
                                    "Marcus","McDermott","Williams","UmbrellaWilliams"))
  ContrastMatRat <- contrMatRatio(n=n, type=type, base=base)
  Num.Contrast <- ContrastMatRat$numC
  Den.Contrast <- ContrastMatRat$denC
  comp.names <- ContrastMatRat$rnames
}

if (is.numeric(Margin)) {
  if (length(Margin)==1) {
    Margin <- rep(Margin,nrow(Num.Contrast))
  }
  if (length(Margin)!=nrow(Num.Contrast)) {
    stop("Margin must be a single numeric value, or a numeric vector with length equal to the 
         number of contrasts")
  }
}

defrvec <- numeric(nrow(Num.Contrast))
for (z in 1:nrow(Num.Contrast)) {
  defrvec[z] <- ( (sum((Num.Contrast[z,]-Margin[z]*Den.Contrast[z,])^2*sd^2/n))^2 ) / 
                sum( ( (Num.Contrast[z,]-Margin[z]*Den.Contrast[z,])^4*sd^4 ) / ( n^2*(n-1) ) )
}
names(defrvec) <- comp.names

return(defrvec)


}
