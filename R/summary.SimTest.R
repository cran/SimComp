summary.SimTest <-
function(object, digits=4, ...) {


cat("", "\n")
if (object$test.class=="differences") {
  if (!is.list(object$ContrastMat)) {
    cat("Contrast matrix:", "\n")
  } else {
    cat("Contrast matrices:", "\n")
  }
  print(object$ContrastMat)
} else {
  if (!is.list(object$Num.Contrast)) {
    cat("Numerator contrast matrix:", "\n")
    print(object$Num.Contrast)
    cat("", "\n")
    cat("Denominator contrast matrix:", "\n")
    print(object$Den.Contrast)
  } else {
    cat("Numerator contrast matrices:", "\n")
    print(object$Num.Contrast)
    cat("", "\n")
    cat("Denominator contrast matrices:", "\n")
    print(object$Den.Contrast)
  }
}

cat("", "\n")
if (object$covar.equal==TRUE) {
  cat("Estimated common covariance matrix of the data:", "\n")
  print(round(object$CovMatDat, digits), digits=digits)
} else {
  cat("Estimated covariance matrices of the data:", "\n")
  print(lapply(object$CovMatDat, round, digits=digits), digits=digits)
}

cat("", "\n")

if (class(object$CorrMatDat)=="UserMatrix") {
  cat("Specified common correlation matrix of the data:", "\n")
  print(round(object$CorrMatDat, digits), digits=digits)
} else {
  if (object$covar.equal==TRUE) {
    cat("Estimated common correlation matrix of the data:", "\n")
    print(round(object$CorrMatDat, digits), digits=digits)
  } else {
    cat("Estimated correlation matrices of the data:", "\n")
    print(lapply(object$CorrMatDat, round, digits=digits), digits=digits)
  }
}

cat("", "\n")
cat("Estimated correlation matrix of the comparisons:", "\n")
print(round(object$CorrMatComp,digits), digits=digits)

cat("", "\n")
cat("Alternative hypotheses: True", object$test.class)
  if (object$alternative=="greater") cat(" greater than ")
  if (object$alternative=="less") cat(" less than ")
  if (object$alternative=="two.sided") cat(" not equal to ")
  cat("the margins", "\n")

comparison <- rep(object$comp.names, each=length(object$resp))
endpoint <- rep(object$resp, times=length(object$comp.names))
margin <- estimate <- statistic <- degr.fr <- p.value.raw <- p.value.adj <- NULL
for (i in 1:length(object$comp.names)) {
  margin      <- c(margin,      round(object$Margin[i,],   digits))
  estimate    <- c(estimate,    round(object$estimate[i,], digits))
  statistic   <- c(statistic,   round(object$statistic[i,],digits))
  degr.fr     <- c(degr.fr,     round(object$degr.fr[i,],  digits))
  p.value.raw <- c(p.value.raw, round(object$p.val.raw[i,],digits))
  p.value.adj <- c(p.value.adj, round(object$p.val.adj[i,],digits))
}
out <- data.frame(comparison, endpoint, margin, estimate, statistic, degr.fr, p.value.raw, p.value.adj)
cat("", "\n")
print(out, digits=digits)
cat("", "\n")


}
