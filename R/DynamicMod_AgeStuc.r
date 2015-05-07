#' Age-structured GTG dynamic projection model.
#' of stock and relative yield. This model is parameterised with the life-history ratios.
#' @name DynamicMod_AgeStuc
#' @title Projection Model
#' @param SimPars An object of class \code{list} that contains all parameters
#'   required to run GTG model.  Full description of model to be added at later
#'   date.
#' @param RecDevVec A vector of log-normally distributed recruitment deviations of length equal to
#' number of projection years (to be added to SimPars later).
#' @return  To add details later.
#' @author Adrian Hordyk
#' @seealso \code{\link{}}
#' @export

DynamicMod_AgeStuc <- function(SimPars, RecDevVec) {
  with(SimPars, { 
    TSFpar <- FM * TSMpar

    # Setup Empty Matrices & Vectors #
    # Time-step 
    SpUnFPR <- SpFPR <- CurrEggs <- SpawnBioTS <- rep(NA, NTimeSteps)
    CatchTS <- rep(NA, NTimeSteps)
    NatAgeMonth <- NatAgeMonthUF <- matrix(NA, nrow=length(AgeVec), ncol=NTimeSteps)
    LenDat <- LenDatPopUF <- LenDatPopF <- matrix(0, nrow=length(LenMids), ncol=NTimeSteps)
    TrackRecruits <- matrix(NA, nrow=NTimeSteps, ncol=NGTG)
    
    # Annual 
    SPR <-  AnnualEgProd <- rep(NA, NYears)
    SBBiomassYr <- SBpr <- AnnualRecVec <- rep(NA, NYears)
    CatchAnnual <- rep(NA, NYears)
    CatchLenComp <- matrix(NA, nrow=length(LenMids), ncol=NYears)
    
    # Run First Year 
	AnnualRec <- R0
    FirstYrFpar <- 0 # Unfished Eq
    RunSingleYr <- SingleYearEq(FparYr=FirstYrFpar, SimPars)
    
    # First 12 TS 
    TrackRecruits[1:12,] <- RunSingleYr$TrackRecruits
    NatAgeMonth[,1:12] <- RunSingleYr$NatAgeMonth
	NatAgeMonthUF[,1:12] <- RunSingleYr$NatAgeMonthUF
    # First Yr #
    SBBiomassYr[1] <- RunSingleYr$SBBiomassYr
	B0 <- SBBiomassYr[1]
    AnnualEgProd[1] <- RunSingleYr$CurrEggs
    AnnualRecVec[1] <- R0
    SPR[1] <- 1
    CatchAnnual[1] <- RunSingleYr$CatchAnnual
    UnfishedEggProd <- RunSingleYr$CurrEggs
    LastNAgeVec <- RunSingleYr$EndYrAgeVecGTG
    LastNAgeUFVec <- RunSingleYr$EndYrAgeVecUFGTG
      
    # Year 2 onwards 
    for (TS in 13:NTimeSteps) {
      if (TS == 13) { TrackMonth <- 0; Yr <- 1 }
      TrackMonth <- TrackMonth + 1
      
      RecYr <- AnnualRecVec[Yr]
      	
      RecGTGMonth <- RecYr * RecGTGMonthProb[TS] # Add variable recruitment here
      TrackRecruits[TS,] <- Probs * RecGTGMonth
      # Run Each GTG one TS #
      RunGTGs <- sapply(1:NGTG, function(GTG) {
          RecGTG <- Probs[GTG] * RecGTGMonth
          SingleGTGDynamic(AgeVec, LinfGTG=DiffLinfs[GTG], MparGTG=MparGTG[GTG], RecGTGMonth=RecGTG, TSkpar=TSkpar, L50, L95, SL50, SL95, Walpha, Wbeta, FecB, LenMids, LenBins, Mpow, FparYr=TSFpar, MeanLinf=Linf, LastNAgeVec=LastNAgeVec[,GTG], LastNAgeUFVec=LastNAgeUFVec[,GTG], TrackRecruitsGTG=TrackRecruits[,GTG], TS, R0=NULL)
      })
      
      LastNAgeVec <- sapply(1:NGTG, function (X) RunGTGs[,X]$Fished)  
      LastNAgeUFVec <- sapply(1:NGTG, function (X) RunGTGs[,X]$Unfished)
      
      # SpFPR[TS] <- sum(sapply(1:NGTG, function(X) RunGTGs[,X]$SpFPR))
      # SpUnFPR[TS] <- sum(sapply(1:NGTG, function(X) RunGTGs[,X]$SpUnFPR))
      SpFPR[TS] <- sum(sapply(1:NGTG, function(X) RunGTGs[,X]$SpFPR2)) # Age based calc
      SpUnFPR[TS] <- sum(sapply(1:NGTG, function(X) RunGTGs[,X]$SpUnFPR2)) # Age based calc 
      CurrEggs[TS] <- sum(sapply(1:NGTG, function(X) RunGTGs[,X]$CurrEggs)) # Age based calc   
      
      NatAgeMonth[,TS] <- apply(LastNAgeVec, 1, sum)
      NatAgeMonthUF[,TS] <- apply(LastNAgeUFVec, 1, sum)
      
      # Spawn Biomass
      spBio <- sapply(1:NGTG, function (X) RunGTGs[,X]$SpawnBioF)
      SpawnBioTS[TS] <- sum(apply(spBio, 1, sum))
      
      CatchTS[TS] <- sum(sapply(1:NGTG, function (X) RunGTGs[,X]$Catch))
      
      # Make Size Comp 
      LenDat[,TS] <- apply(sapply(1:NGTG, function(X) RunGTGs[,X]$LengthCompCatch), 1, sum)
      LenDatPopF[,TS] <- apply(sapply(1:NGTG, function(X) RunGTGs[,X]$LengthCompFished), 1, sum)
      LenDatPopUF[,TS] <- apply(sapply(1:NGTG, function(X) RunGTGs[,X]$LengthCompUF), 1, sum)
      
      # Annual Calculations
      if (TrackMonth == 12) {
        Yr <- Yr + 1 
    	# SPR[Yr] <- sum(SpFPR[(TS-11):TS])/sum(SpUnFPR[(TS-11):TS]) 
    	EggPRFished <- sum(SpFPR[(TS-11):TS])
    	EggPRUNFished <- sum(SpUnFPR[(TS-11):TS]) 
    	SPR[Yr] <- EggPRFished/EggPRUNFished
    	
    	# Length Composition Sample 
    	if (TSFpar > 0) {  
          SampleMonthProb <- SampleProb * SampleSize
		  MonthlySamples <- sapply(1:12, function (X) rmultinom(n=1, size=SampleMonthProb[X], prob=LenDat[,TS-12+X]))
    	  CatchLenComp[,Yr] <- apply(MonthlySamples, 1, sum)
    	}
	
    	# Annual Spawn Biomass 
        AnnualEgProd[Yr] <- sum(CurrEggs[(TS-11):TS])
    	SBBiomassYr[Yr] <- sum(SpawnBioTS[(TS-11):TS])
    	AnnualRecVec[Yr] <- BHRecruitFun(currEggProd= AnnualEgProd[Yr], steepness, R0, UnfishedEggProd, RecDev=RecDevVec[Yr])	
    	TrackMonth <- 0 
    	
    	# Catch 
    	CatchAnnual[Yr] <- sum(CatchTS[(TS-11):TS])
    	
    	# Adjust F for next year 
    	TSFpar <- TSFpar 
      }
      
      print(paste("Month", TS, "of", NTimeSteps))
    }
    output <- NULL
    output$AnnSpBio <- SBBiomassYr
    output$AnnEgProd  <- AnnualEgProd
    output$AnnCatch  <- CatchAnnual
    output$AnnCatchComp  <- CatchLenComp
    output$AnnRec  <- AnnualRecVec
    output$SPR <- SPR
    output$AnnF <- TSFpar
	output$TSMpar <- TSMpar
	output$LenMids <- LenMids
	output$AgeVec <- AgeVec
	output$B0 <- B0
	output$NatAgeMonth <- NatAgeMonth
	output$NatAgeMonthUF <- NatAgeMonthUF
	output$TrackRecruits <- TrackRecruits
	
    return(output)
  })	
}

