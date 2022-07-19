#expression = conc$time

# nonlinear term definition for gnm() formula
log_sum_exp <- function(expression, inst = NULL){
  list(predictors = list(a1 = 1, c1 = 1, a2 = 1, c2 = 1, a3 = 1, c3 = 1),
       variables = list(substitute(expression)),
       term = function(predLabels, varLabels) {
         paste("log(exp(", predLabels[1],"+",predLabels[2],"*",varLabels[1], ")+
               exp(", predLabels[3],"+",predLabels[4],"*",varLabels[1], ")+
               exp(", predLabels[5],"+",predLabels[6],"*",varLabels[1], "))", sep = "")
       })
}
class(log_sum_exp) <- "nonlin"

