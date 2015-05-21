#' This function Predict the value of Mslope parameter from the LM and Pars.  
#' @name PredictMSlope  
#' @title Predict the value of Mslope parameter from the LM and Pars.
#' @param Pars An object of class \code{list} that contains all parameters
#'   required to run GTG model. 
#' @author Adrian Hordyk 
#' @export

PredictMSlope <- function(Pars) {
  data(MSlopeMod)
  Mslope <- exp(predict(MslopeMod, newdata=list(logLinf=log(Pars$Linf), RelL50=Pars$relL50, RelL95=Pars$relL95, Mpow=Pars$Mpow, logMK=log(Pars$MK)))) 
  return(as.numeric(Mslope))
}
