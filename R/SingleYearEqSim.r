#' Equilibrium model to project a single year of the age-structured GTG model (usually first unfished year)
#' @name SingleYearEq
#' @title Single Year Equilibrium Projection Model
#' @param FparYr Fishing mortality (in monthly units of time). Usually set to 0 for unfished equlibrium.
#' @param SimPars An object of class \code{list} that contains all parameters
#'   required to run GTG model.  Full description of model to be added at later
#'   date.
#' @return  To add details later.
#' @author Adrian Hordyk
#' @seealso \code{\link{}}
#' @export

# Initial Year - Equilibrium model accounting for recruitment by month #
SingleYearEq <- function(FparYr, SimPars) {
  with(SimPars, {
    # Setup Matrices 
    NTimeSteps <- 12
    SpUnFPR <- SpFPR <- CurrEggs <- SpawnBioTS <- rep(NA, NTimeSteps)
    CatchTS <- rep(NA, NTimeSteps)
    NatAgeMonth <- NatAgeMonthUF <- matrix(NA, nrow=length(AgeVec), ncol=NTimeSteps)
    LenDat <- LenDatPopUF <- LenDatPopF <- matrix(0, nrow=length(LenMids), ncol=NTimeSteps)
    TrackRecruits <- matrix(NA, nrow=NTimeSteps, ncol=NGTG)
    
    RecGTGMonth <- TrackRecruits[1,] <- R0 * Probs * RecGTGMonthProb[1]
    MKGTG <- (MK)+ Mslope*(DiffLinfs-Linf)
    MparGTG <- MKGTG * TSkpar
    # Initial Equilibrium Unfished Age Structure - First Month 
    RunGTGs1stMonth <- sapply(1:NGTG, function(GTG) {
	  RecGTG <- Probs[GTG] * RecGTGMonth
      SingleGTGDynamic(AgeVec, LinfGTG=DiffLinfs[GTG], MparGTG=MparGTG[GTG], RecGTGMonth=RecGTG, TSkpar=TSkpar, L50, L95, SL50, SL95, Walpha, Wbeta, FecB, LenMids, LenBins, Mpow, FparYr=FparYr, MeanLinf=Linf, LastNAgeVec=NULL, LastNAgeUFVec=NULL, TrackRecruitsGTG=TrackRecruits[,GTG], TS=1, RecGTGMonthVec=MonthRecProb, RecProb=Probs[GTG], R0=R0, RecGTGMonthProb=RecGTGMonthProb)
    })
    
    FirstMonthEqAgeFs <- sapply(1:NGTG, function (X) RunGTGs1stMonth[,X]$Fished)
    FirstMonthEqAgeUF <- sapply(1:NGTG, function (X) RunGTGs1stMonth[,X]$Unfished)
    NatAgeMonth[,1] <-  apply(FirstMonthEqAgeFs, 1, sum) # Age structure at end of month 1 
	
    NatAgeMonthUF[,1] <- apply(FirstMonthEqAgeUF, 1, sum) # Age structure at end of month 1 
    SpawnBioTS[1] <- sum(apply(sapply(1:NGTG, function (X) RunGTGs1stMonth[,X]$SpawnBioUF),1,sum))
	SpUnFPR[1] <- sum(sapply(1:NGTG, function (X) RunGTGs1stMonth[,X]$SpUnFPR2))
	SpFPR[1] <- sum(sapply(1:NGTG, function (X) RunGTGs1stMonth[,X]$SpFPR2))
    CurrEggs[1] <- sum(sapply(1:NGTG, function (X) RunGTGs1stMonth[,X]$CurrEggs))
    CatchTS[1] <- sum(sapply(1:NGTG, function (X) RunGTGs1stMonth[,X]$Catch))
    
    LastNAgeVec <- FirstMonthEqAgeFs
    LastNAgeUFVec <- FirstMonthEqAgeUF
    # Run First Year (months 2 to 12) 
    for (TS in 2:12) {
      RecGTGMonth <- R0 * RecGTGMonthProb[TS] 
      TrackRecruits[TS,] <- Probs * RecGTGMonth
      
      # Run Each GTG one TS #
      RunGTGs <- sapply(1:NGTG, function(GTG) {
          RecGTG <- Probs[GTG] * RecGTGMonth
          SingleGTGDynamic(AgeVec, LinfGTG=DiffLinfs[GTG], MparGTG=MparGTG[GTG], RecGTGMonth=RecGTG, TSkpar=TSkpar, L50, L95, SL50, SL95, Walpha, Wbeta, FecB, LenMids, LenBins, Mpow, FparYr=FparYr, MeanLinf=Linf, LastNAgeVec=LastNAgeVec[,GTG], LastNAgeUFVec=LastNAgeUFVec[,GTG], TrackRecruitsGTG=TrackRecruits[,GTG], TS, RecGTGMonthVec=MonthRecProb, RecProb=Probs[GTG], R0=R0, RecGTGMonthProb=RecGTGMonthProb)
      })
      LastNAgeVec <- sapply(1:NGTG, function (X) RunGTGs[,X]$Fished)  
      LastNAgeUFVec <- sapply(1:NGTG, function (X) RunGTGs[,X]$Unfished)
      
      spBio <- sapply(1:NGTG, function (X) RunGTGs[,X]$SpawnBioF)
      SpawnBioTS[TS] <- sum(apply(spBio, 1, sum))
	  SpUnFPR[TS] <- sum(sapply(1:NGTG, function (X) RunGTGs[,X]$SpUnFPR2))
	  SpFPR[TS] <- sum(sapply(1:NGTG, function (X) RunGTGs[,X]$SpFPR2))
      CurrEggs[TS] <- sum(sapply(1:NGTG, function(X) RunGTGs[,X]$CurrEggs)) # Age based calc 
      CatchTS[TS] <- sum(sapply(1:NGTG, function (X) RunGTGs[,X]$Catch))
      
      NatAgeMonth[,TS] <- apply(LastNAgeVec, 1, sum)
      NatAgeMonthUF[,TS] <- apply(LastNAgeUFVec, 1, sum)
      print(paste("Month", TS, "of", 12))
    }
	
	# Generate Size Composition 
	
	
	
    output <- NULL
    output$EndYrAgeVecGTG <- LastNAgeVec
    output$EndYrAgeVecUFGTG <- LastNAgeUFVec
    output$NatAgeMonth <- NatAgeMonth
    output$NatAgeMonthUF <- NatAgeMonthUF
    output$TrackRecruits <- TrackRecruits
    output$SBBiomassYr <- sum(SpawnBioTS[1:12])
    output$CurrEggs <- sum(CurrEggs[1:12])
    output$CatchAnnual <- sum(CatchTS[1:12])
	output$SpFPR <- SpFPR
	output$SpUnFPR <- SpUnFPR
	output$SPR <- sum(SpFPR)/sum(SpUnFPR)
    return(output)
  })	
}

