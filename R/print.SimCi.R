print.SimCi <-
function(x, digits=4, ...) {


cat("", "\n")
cat("Simultaneous", x$conf.level*100)
cat("% confidence intervals for", x$test.class, "of means of multiple endpoints", "\n")
cat("Assumption: ")
  if (x$covar.equal==TRUE) cat("Homogeneous ") else cat("Heterogeneous ")
  cat("covariance matrices for the groups", "\n")
if (x$na.action=="multi.df") cat("Approach for missing values: multi.df", "\n")

comparison <- rep(x$comp.names, each=length(x$resp))
endpoint <- rep(x$resp, times=length(x$comp.names))
estimate <- degr.fr <- lower.raw <- upper.raw <- lower <- upper <- NULL
for (i in 1:length(x$comp.names)) {
  estimate  <- c(estimate,  round(x$estimate[i,], digits))
  degr.fr   <- c(degr.fr,   round(x$degr.fr[i,],  digits))
  lower.raw <- c(lower.raw, round(x$lower.raw[i,],digits))
  upper.raw <- c(upper.raw, round(x$upper.raw[i,],digits))
  lower     <- c(lower,     round(x$lower[i,],    digits))
  upper     <- c(upper,     round(x$upper[i,],    digits))
}
out <- data.frame(comparison, endpoint, estimate, degr.fr, lower.raw, upper.raw, lower, upper)
cat("", "\n")
print(out, digits=digits)
cat("", "\n")

if (x$test.class=="ratios" && x$NSD>0) {
  cat("The mean in", x$NSD, "denominators is not significantly different from zero.", "\n")
  cat("", "\n")
}

invisible(x)


}
