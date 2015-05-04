#' Function to simulate age-structured growth-type-group (GTG) model to 
#' generate size equilibrium composition of population and catch, as well as SPR
#' of stock and relative yield.  
#' @name SimMod_AgeEq
#' @title Generate size structure using GTG model and calculate relative YPR and SPR
#' @param SimPars An object of class \code{list} that contains all parameters
#'   required to run GTG model.  Full description of model to be added at later
#'   date.
#' @return  To add details later.
#' @author Adrian Hordyk
#' @export


SimMod_AgeEq <- function(SimPars) {
  with(SimPars, {
    # Mortality for given time-step
    TSMpar <- Mpar / TStep
    TSkpar <- kpar/ TStep
	TSFpar <- FM * TSMpar
    MaxAge <- round(-log(0.001)/TSMpar,0) # Maximum age 
    AgeVec <- 0:MaxAge # in units of time-step
    
	# GTG set-up 
	SDLinf <- CVLinf * Linf
    DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, by=1) # Linfs of GTGs by units of 1 
    NGTG <- length(DiffLinfs)
    Probs <- dnorm(DiffLinfs, Linf, sd=SDLinf)/sum(dnorm(DiffLinfs, Linf, sd=SDLinf)) 
    ToSize <- max(DiffLinfs) 
    ToSize <- ceiling(max(ToSize)/Linc)*Linc # round up 
    LenBins <- seq(from=0, to=ToSize, by=Linc)
    LenMids <- seq(from=LenBins[2]-0.5*Linc, by=Linc, length=length(LenBins)-1)
	
	# Selectivity by size - add MLL later 
    SelectGearLen <- 1.0/(1+exp(-log(19)*(LenMids-SL50)/((SL95)-SL50)))
	# SelectMLLLen <- 1.0/(1+exp(-log(19)*(LenMids-MLL50)/((MLL95)-MLL50)))
    # SelectMLLLen <- SelectMLLLen * SelectGearLen
	
	# Single GTG Model 
    SingleGTGAgeModel <- function(AgeVec, LinfGTG, MparGTG, RecGTG, kpar, L50, L95, SL50, SL95, Walpha, Wbeta, FecB, LenMids, LenBins, Mpow, Fpar, MaxAge, MeanLinf) {  
      LenVec <-  LinfGTG * (1-exp(-kpar*(AgeVec+0.5)))
      SelAge <- 1.0/(1+exp(-log(19)*(LenVec-SL50)/((SL95)-SL50)))
      
      MatAge <- rep(MparGTG, length(AgeVec)) # add Variable M here - need to condition on size
      if (Mpow > 0) MatAge <- MparGTG * (MeanLinf/(LenVec))^Mpow
      FatAge <- SelAge * Fpar # F for each age 
      ZatAge <- FatAge + MatAge # Z for each GTG 
      
      # Age Model #
      NatAgeUnfished <- rep(0, length(AgeVec))
      NatAgeFished <- rep(0, length(AgeVec))
      NatAgeCatch <- rep(0, length(AgeVec))
      for (age in 1: length(AgeVec)) {
        if (age == 1) {
    	  NatAgeUnfished[age] <- RecGTG  # Assume recruitment is continuous process - change this later to make seasonal
    	  NatAgeFished[age] <- RecGTG
    	} 
        if (age > 1) {
    	  NatAgeUnfished[age] <- NatAgeUnfished[age-1] * exp(-MatAge[age-1])  
          NatAgeFished[age] <- NatAgeFished[age-1] * exp(-ZatAge[age-1])
    	}  
    	if (AgeVec[age] == MaxAge) {
    	  NatAgeUnfished[age] <- NatAgeUnfished[age-1] * (exp(-MatAge[age-1])/(1-exp(-MatAge[age])))
          NatAgeFished[age] <- NatAgeFished[age-1] * (exp(-ZatAge[age-1])/(1-exp(-ZatAge[age])))
    	}  
        NatAgeCatch[age] <- FatAge[age]/ZatAge[age] * NatAgeFished[age] * (1-exp(-ZatAge[age]))
      }
     
      # Length Structure # 
      N <- NatAgeCatch
      LengthCompCatch <- LengthCompFished <- LengthCompUF <- rep(0, length(LenMids))
      for (A in 1:length(N)) {
        Ind <- min(which(LenBins > LenVec[A]))
    	LengthCompUF[Ind-1] <- LengthCompUF[Ind-1] + NatAgeUnfished[A]
    	LengthCompFished[Ind-1] <- LengthCompFished[Ind-1] + NatAgeFished[A]
        LengthCompCatch[Ind-1] <- LengthCompCatch[Ind-1] + NatAgeCatch[A]
      } 
    
      # Calc SPR 
      MatVec <- 1.0/(1+exp(-log(19)*(LenMids-L50)/((L95)-L50)))
      WatAge <- Walpha * LenMids ^ Wbeta
      FecVec <- (LenMids ^ FecB) * MatVec
      FecVec <- FecVec/max(FecVec)
    
      SpUnF <- LengthCompUF * FecVec
      SpF <- LengthCompFished * FecVec
      SpUnFPR <- sum(SpUnF)/RecGTG
      SpFPR <- sum(SpF)/RecGTG 
    
      # MatVec <- 1.0/(1+exp(-log(19)*(LenVec-L50)/((L95)-L50)))
      # WatAge <- Walpha * LenVec ^ Wbeta
      # FecVec <- (LenVec ^ FecB) * MatVec
      # FecVec <- FecVec/max(FecVec)  
     
      # SpUnF <- NatAgeUnfished * FecVec
      # SpF <- NatAgeFished * FecVec
      # SpUnFPR <- sum(SpUnF)/RecGTG
      # SpFPR <- sum(SpF)/RecGTG
      
      FitPR <- SpUnFPR
     
      return(list(Age=AgeVec, Unfished=NatAgeUnfished, Fished=NatAgeFished, Catch=NatAgeCatch, LatA=LenVec, SpUnFPR=SpUnFPR, SpFPR=SpFPR, LengthCompUF=LengthCompUF, LengthCompFished=LengthCompFished, LengthCompCatch=LengthCompCatch, FitPR=FitPR))
    }

	MKGTG <- (TSMpar/TSkpar)+ Mslope*(DiffLinfs-Linf)
	# Run population model for each GTG
	RunGTGs <- sapply(1:NGTG, function(X) {
      LinfGTG <- DiffLinfs[X]
      MparGTG <- MKGTG[X] * TSkpar
      RecGTG <- Probs[X] * R0
      SingleGTGAgeModel(AgeVec, LinfGTG, MparGTG, RecGTG, TSkpar, L50, L95, SL50, SL95, Walpha, Wbeta, FecB, LenMids, LenBins, Mpow, TSFpar, MaxAge, Linf)
    })
	
    SPR <- sum(sapply(1:NGTG, function(X) RunGTGs[,X]$SpFPR))/sum(sapply(1:NGTG, function(X) RunGTGs[,X]$SpUnFPR)) # Calc SPR 
    FitPR <- sapply(1:NGTG, function(X) RunGTGs[,X]$FitPR) # Calc fit for each GTG 
	ObjFun <- sum((FitPR - median(FitPR, na.rm=TRUE))^2, na.rm=TRUE) # This needs to be minimised to make fitness approximately equal across GTG - by adjusting Mslope 
	Pen <- 0; if (min(MKGTG) <= 0 ) Pen <- (1/abs(min(MKGTG)))^2 * 1E5 # Penalty for optimising Mslope   
    ObjFun <- ObjFun + Pen
	
	# Make Size Comp 
    LenDat <- apply(sapply(1:NGTG, function(X) RunGTGs[,X]$LengthCompCatch), 1, sum)
    LenDatPopF <- apply(sapply(1:NGTG, function(X) RunGTGs[,X]$LengthCompFished), 1, sum)
    LenDatPopUF <- apply(sapply(1:NGTG, function(X) RunGTGs[,X]$LengthCompUF), 1, sum)
	
   output <- NULL
   output$SPR <- SPR
   output$LenDat <- LenDat
   output$LenDatPopF <- LenDatPopF
   output$LenDatPopUF <- LenDatPopUF
   output$LenMids <- LenMids 
   output$LenBins <- LenBins
   output$ObjFun <- ObjFun
   output$Pen <- Pen   
   
   # add catch etc later
   return(output)
  })
}
