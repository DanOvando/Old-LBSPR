#' This function Predict the value of kslope parameter from the LM and SimPars.  
#' @name PredictKSlope  
#' @title Predict the value of kslope parameter from the LM and SimPars.
#' @param SimPars An object of class \code{list} that contains all parameters
#'   required to run GTG model. 
#' @author Adrian Hordyk 
#' @export

PredictKSlope <- function(SimPars) {
  data(KSlopeMod)
  KSlope <- exp(predict(KSlopeMod, newdata=list(LinfMK=log(SimPars$Linf)-log(SimPars$MK), CVLinf=SimPars$CVLinf, RelL50=(SimPars$L50/SimPars$Linf), Mpow=SimPars$Mpow)))  
  return(KSlope)
}
