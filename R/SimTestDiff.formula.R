SimTestDiff.formula <-
function(formula, ...) {
  resp <- all.vars(formula[[2]]) 
  grp  <- all.vars(formula[[3]]) 
  SimTestDiff(grp=grp, resp=if (length(resp)==0) NULL else resp, ...)
}
