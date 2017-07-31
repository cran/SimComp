DfSattDiff <-
function(n, sd, type="Dunnett", base=1, ContrastMat=NULL) {


if (length(n)!=length(sd)) {
  stop("Lengths of n and sd must be equal")
}

if (!is.null(ContrastMat)) {
  if (ncol(ContrastMat)!=length(n)) {
    stop("Number of columns of ContrastMat and length of n must be equal")
  }
  if (any(apply(ContrastMat==0, 1, all))) {
    stop("At least one row of ContrastMat is a vector with all components", 
         "\n", "equal to zero")
  }
  type <- "User defined"
  if (is.null(rownames(ContrastMat))) {
    rownames(ContrastMat) <- paste("C", 1:nrow(ContrastMat), sep="")
  }
} else {
  type <- match.arg(type, choices=c("Dunnett","Tukey","Sequen","AVE","GrandMean","Changepoint",
                                    "Marcus","McDermott","Williams","UmbrellaWilliams"))
  ContrastMat <- contrMat(n=n, type=type, base=base)
}
comp.names <- rownames(ContrastMat)

defrvec <- numeric(nrow(ContrastMat))
for (z in 1:nrow(ContrastMat)) {
  defrvec[z] <- ( (sum((ContrastMat[z,])^2*sd^2/n))^2 ) / 
                sum( ( (ContrastMat[z,])^4*sd^4 ) / ( n^2*(n-1) ) )
}
names(defrvec) <- comp.names

return(defrvec)


}
