#' Beverton-Holt Stock Recruit Relationship.
#' @name BHRecruitFun
#' @title Calculates the recruitment from the Beverton-Holt SRR model.
#' @param currEggProd Current egg production - total number of eggs or spawning biomass at the 
#' beginning of the year.
#' @param steepness Steepness (h) of the stock-recruit relationship.
#' @param R0 Number of recruits at unfished equilibrium.
#' @param UnfishedEggProd Egg production (or spawning biomass) in the equilibrium unfished state.
#' @param RecDev Log-normally distributed recruitment deviation.
#' @return  To add details later.
#' @author Adrian Hordyk
#' @seealso \code{\link{}}
#' @export

BHRecruitFun  <- function(currEggProd, steepness, R0, UnfishedEggProd, RecDev=1) {
  Init.EggRec	<- UnfishedEggProd/R0
  Alpha <- Init.EggRec * ((1-steepness)/(4*steepness))
  Beta  <- (5*steepness - 1)/(4*steepness*R0)
  Recruits <- currEggProd/(Alpha + Beta * currEggProd)	
  Recruits <- round(Recruits * RecDev,0)
  return(Recruits)
}
  