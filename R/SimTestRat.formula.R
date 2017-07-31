SimTestRat.formula <-
function(formula, ...) {
  resp <- all.vars(formula[[2]]) 
  grp  <- all.vars(formula[[3]]) 
  SimTestRat(grp=grp, resp=if (length(resp)==0) NULL else resp, ...)
}
