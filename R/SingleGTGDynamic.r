#' Function for running age-structured GTG model. Not meant to be run.
#' @name SingleGTGDynamic
#' @title Single GTG Dynamic Projection Model Function
#' @param This function is called by other functions and parameters do not need to be set.
#' @author Adrian Hordyk
#' @seealso \code{\link{}}
#' @export

SingleGTGDynamic <- function(AgeVec, LinfGTG, MparGTG, RecGTGMonth, TSkpar, L50, L95, SL50, SL95, Walpha, Wbeta, FecB, LenMids, LenBins, Mpow, FparYr, MeanLinf, LastNAgeVec, LastNAgeUFVec, TrackRecruitsGTG, TS, RecGTGMonthVec, RecProb, R0) {  
  LenVec <-  LinfGTG * (1-exp(-TSkpar*(AgeVec+0.5)))
  SelAge <- 1.0/(1+exp(-log(19)*(LenVec-SL50)/((SL95)-SL50)))
  
  MatAge <- rep(MparGTG, length(AgeVec))
  if (Mpow > 0) MatAge <- MparGTG * (MeanLinf/(LenVec))^Mpow
  FatAge <- SelAge * FparYr # F for each age 
  ZatAge <- FatAge + MatAge # Z for each GTG 
  # Age Model #
  MaxAge <- max(AgeVec)
  NatAgeUnfished <- rep(0, length(AgeVec))
  NatAgeFished <- rep(0, length(AgeVec))
  NatAgeCatch <- rep(0, length(AgeVec))
  
  # First Month - Unfished Equilibrium #
  if (TS == 1) {
    rec <- R0 * RecProb
    NatAgeFished[1] <- NatAgeUnfished[1] <- RecGTGMonthVec[1] * rec 
    KK <- 12:2
    for (age in 2:12) {
   	  NatAgeUnfished[age] <- rec * RecGTGMonthVec[KK[age-1]] * exp(-MatAge[age-1]*(age-1))
	  NatAgeFished[age] <- rec * RecGTGMonthVec[KK[age-1]] * exp(-ZatAge[age-1]*(age-1))
	}  
    for (age in 13:length(AgeVec)) { 
	  NatAgeUnfished[age] <- NatAgeUnfished[age-12] * exp(-sum(MatAge[(age-12):(age-1)]))
	  NatAgeFished[age] <- NatAgeFished[age-12] * exp(-sum(ZatAge[(age-12):(age-1)]))
	}  
	NatAgeCatch[age] <- FatAge[age]/ZatAge[age] * NatAgeFished[age] * (1-exp(-ZatAge[age]))
  }
  
  # Other Months
  if (TS > 1) {
    for (age in 1: length(AgeVec)) {
      if (age == 1) {
	    NatAgeUnfished[age] <- RecGTGMonth  
	    NatAgeFished[age] <- RecGTGMonth
	  } 
      if (age > 1) {
	    NatAgeUnfished[age] <- LastNAgeUFVec[age-1] * exp(-MatAge[age-1])  
        NatAgeFished[age] <- LastNAgeVec[age-1] * exp(-ZatAge[age-1])
		# if (age == 25) print(LastNAgeVec[age-1])
	  }  
	  if (AgeVec[age] == MaxAge) {
	    NatAgeUnfished[age] <- LastNAgeUFVec[age-1] * (exp(-MatAge[age-1])) #/(1-exp(-MatAge[age])))
        NatAgeFished[age] <- LastNAgeVec[age-1] * (exp(-ZatAge[age-1])) # /(1-exp(-ZatAge[age])))
	  }  
      NatAgeCatch[age] <- FatAge[age]/ZatAge[age] * NatAgeFished[age] * (1-exp(-ZatAge[age]))
    }
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
  # # SPR by Length - can't calculate transitional SPR with recruitment variability 
  # MatVec <- 1.0/(1+exp(-log(19)*(LenMids-L50)/((L95)-L50)))
  # FecVec <- (LenMids ^ FecB) * MatVec
  # SpUnF <- LengthCompUF * FecVec
  # SpF <- LengthCompFished * FecVec
  # SpUnFPR <- SpFPR <- FitPR <- 0
  # if (TS > MaxAge) {	
    # SpUnFPR <- sum(SpUnF) #RecGTGMonth
    # SpFPR <- sum(SpF) #RecGTGMonth 
  # }
 
  # SPR by Age 
  MatVec <- 1.0/(1+exp(-log(19)*(LenVec-L50)/((L95)-L50)))
  WghtVec <- Walpha * LenVec ^ Wbeta
  FecVec <- (LenVec ^ FecB) * MatVec
  
  SpUnF <- NatAgeUnfished * FecVec
  SpF <- NatAgeFished * FecVec
  if (TS > MaxAge) {	
    RecVec <- rev(TrackRecruitsGTG[(TS-MaxAge):TS])
	ind <- which(RecVec > 0)
    SpUnFPR2 <- sum(SpUnF[ind]/RecVec[ind])
    SpFPR2 <- sum(SpF[ind]/RecVec[ind])
  }	
  
  if (TS <= MaxAge) { 
    RecVec <- c(rep(R0, MaxAge - TS),rev(TrackRecruitsGTG[1:TS]))
	ind <- which(RecVec > 0)
    SpUnFPR2 <- sum(SpUnF[ind]/RecVec[ind])
    SpFPR2 <- sum(SpF[ind]/RecVec[ind])
  }
  
  SpawnBioUF <- NatAgeUnfished * WghtVec * MatVec
  SpawnBioF <- NatAgeFished * WghtVec * MatVec
  CurrEggs <- sum(SpF)
  
  
  # Catch 
  Catch <- sum(NatAgeCatch * WghtVec)
 
  return(list(Age=AgeVec, Unfished=NatAgeUnfished, Fished=NatAgeFished, NatAgeCatch=NatAgeCatch, Catch=Catch, LatA=LenVec, SpUnFPR2=SpUnFPR2, SpFPR2=SpFPR2, LengthCompUF=LengthCompUF, LengthCompFished=LengthCompFished, LengthCompCatch=LengthCompCatch, SpawnBioUF=SpawnBioUF, SpawnBioF=SpawnBioF, CurrEggs=CurrEggs))
}


