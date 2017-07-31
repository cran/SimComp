summary.SimCi <-
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
if (object$covar.equal==TRUE) {
  cat("Estimated common correlation matrix of the data:", "\n")
  print(round(object$CorrMatDat, digits), digits=digits)
} else {
  cat("Estimated correlation matrices of the data:", "\n")
  print(lapply(object$CorrMatDat, round, digits=digits), digits=digits)
}

cat("", "\n")
cat("Estimated correlation matrix of the comparisons:", "\n")
print(round(object$CorrMatComp,digits), digits=digits)

comparison <- rep(object$comp.names, each=length(object$resp))
endpoint <- rep(object$resp, times=length(object$comp.names))
estimate <- degr.fr <- lower.raw <- upper.raw <- lower <- upper <- NULL
for (i in 1:length(object$comp.names)) {
  estimate  <- c(estimate,  round(object$estimate[i,], digits))
  degr.fr   <- c(degr.fr,   round(object$degr.fr[i,],  digits))
  lower.raw <- c(lower.raw, round(object$lower.raw[i,],digits))
  upper.raw <- c(upper.raw, round(object$upper.raw[i,],digits))
  lower     <- c(lower,     round(object$lower[i,],    digits))
  upper     <- c(upper,     round(object$upper[i,],    digits))
}
out <- data.frame(comparison, endpoint, estimate, degr.fr, lower.raw, upper.raw, lower, upper)
cat("", "\n")
print(out, digits=digits)
cat("", "\n")

if (object$test.class=="ratios" && object$NSD>0) {
  cat("The mean in", object$NSD, "denominators is not significantly different from zero.", "\n")
  cat("", "\n")
}


}
