#' This function to be run after \code{LoadSimPars} and calculates the derived parameters. Re-run this 
#' function each time of the \code{SimPars} parameters is manually changed.
#' @name SimParsCalc  
#' @title Calculate values for derived simulations parameters. 
#' @param SimPars An object of class \code{list} that contains all parameters
#'   required to run GTG model. This function calculates additional parameters and adds them to \code{SimPars}.
#' @param ModType A character string of either \code{Len} or \code{Age} to indicate if the equilibrium simulation model is length or age structured respectively.
#' @return Returns a list of simulation parameters
#' @author Adrian Hordyk 
#' @export

SimParsCalc <- function(SimPars, ModType="Len", OptMslope=FALSE, PredictMslope=TRUE, OverMslope=NULL) {

  OptimiseFitness <- function(logMslope, SimPars, Function) {
    Mslope <- exp(logMslope) #+ 0.000000001
    SimPars$Mslope <- Mslope
	SimPars$FM <- 0
    return(Function(SimPars=SimPars)$ObjFun)
  }
  
  getFMFun <- function(testFM, setSPR, simpars) {
    simpars$FM <- testFM
    return(sum(setSPR - SimMod_LHR(simpars)$SPR)^2)
  }
  
  with(SimPars, {
    SimPars$SDLinf <- SDLinf<- CVLinf * Linf
    SimPars$kpar <- kpar <- Mpar/MK
    # Adjust rate parameters
	SimPars$TSMpar <- TSMpar <- Mpar / TStep
	SimPars$TSFpar <- TSFpar <- FM * TSMpar
    SimPars$TSkpar <- TSkpar <- kpar/ TStep
	SimPars$MaxAge <- MaxAge <- round(-log(0.001)/TSMpar,0) # Maximum age 
    SimPars$AgeVec <- NULL
	if (is.na(MaxAge) != 1) SimPars$AgeVec <- 0:MaxAge  # in units of time-step
	
    # Projection  
	SimPars$NTimeSteps <- NTimeSteps <- NYears * TStep
	
	# Recruitment 
	SimPars$recK   	<- (4*steepness)/(1-steepness) # Goodyear compensation ratio 
	#Recruitment probability by Month
    MonthRecProb <- dnorm(1:12, mean=MeanMonth, sd=MonthSD)
    SimPars$MonthRecProb <- MonthRecProb <- MonthRecProb/sum(MonthRecProb)
	SimPars$RecGTGMonthProb <- NULL
	if (is.na(sum(MonthRecProb)) != 1)   SimPars$RecGTGMonthProb <- rep(MonthRecProb, NTimeSteps/TStep)
	
	
	# GTG set-up
    if(is.na(GTGLinfBy)) {
       SimPars$DiffLinfs <- DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, length=SimPars$NGTG) # Linfs of GTGs by units of GTGLinfBy 
	   SimPars$GTGLinfBy <- GTGLinfBy <- DiffLinfs[2] - DiffLinfs[1]
	}
    if(is.na(NGTG)) {
	  SimPars$DiffLinfs <- DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, by=GTGLinfBy) # Linfs of GTGs by units of GTGLinfBy 
      SimPars$NGTG <- length(DiffLinfs)
	} 
    	
    SimPars$Probs <- dnorm(DiffLinfs, Linf, sd=SDLinf)/sum(dnorm(DiffLinfs, Linf, sd=SDLinf)) 
    SimPars$ToSize <- ToSize <- max(DiffLinfs) 
    ToSize <- ceiling(max(ToSize)/Linc)*Linc # round up 
    SimPars$LenBins <- LenBins <- seq(from=0, to=ToSize, by=Linc)
    SimPars$LenMids <- LenMids <- seq(from=LenBins[2]-0.5*Linc, by=Linc, length=length(LenBins)-1)
	
	# Selectivity by size - add MLL later 
    SimPars$SelectGearLen <- 1.0/(1+exp(-log(19)*(LenMids-SL50)/((SL95)-SL50)))
    # SelectMLLLen <- 1.0/(1+exp(-log(19)*(LenMids-MLL50)/((MLL95)-MLL50)))
    # SelectMLLLen <- SelectMLLLen * SelectGearLen
	
	# Sampling 
	SampleProb <- dnorm(1:12, mean=SampleMonthMean, sd=SampleMonthSD)
    SimPars$SampleProb <- SampleProb/sum(SampleProb)
	
	
	# Optimise Mslope 
	if (ModType == "Len") Function <- SimMod_LHR				
	if (ModType == "Age") Function <- SimMod_AgeEq
	SimPars$AssessOpt <- FALSE
	
	if (length(OverMslope) < 1 & PredictMslope == FALSE & OptMslope==FALSE) stop("Mslope parameters not set")
	if (length(OverMslope) > 0) {
	  print(paste("Mslope set at ", OverMslope))
	  SimPars$Mslope <- OverMslope
	  OptMslope <- PredictMslope <- FALSE
	} 
    if (PredictMslope) {
	  print("Predicting Mslope")
	  SimPars$Mslope <- Mslope <- PredictMSlope(SimPars)
	  OptMslope <- FALSE
	}
	if (OptMslope) {
	  print("Optimising for Mslope - this may take a short while...")
	  SimPars$Mslope <- Mslope <- exp(optimise(OptimiseFitness, interval=log(c(0.0001, 0.1)), SimPars=SimPars, Function=Function)$minimum)
	  # M per GTG 
	  SimPars$MKGTG <- MKGTG <- MK + Mslope*(DiffLinfs-Linf)
	  if (min(SimPars$MKGTG) <= 0) {
	    print("MKgtg is negative for some GTG. Try change NGTG")
	    # readline("********** WARNING - Press Enter to continue **********")
	  }
	}

	SimPars$MKGTG <- MKGTG <- MK + Mslope*(DiffLinfs-Linf)
    SimPars$MparGTG <- MparGTG <- MKGTG * TSkpar
	
	# Start SPR and F/M
	if (startSPR > 0 & is.na(startSPR) !=1) {
	  print("Optimising F/M for startSPR - this may take a short while...")
	  SimPars$FM <- FM <- optimise(getFMFun, interval=c(0,30), setSPR=startSPR, simpars=SimPars)$minimum
	}
 
  print("Derived simulation parameters successfully calculated") 
  return(SimPars)
  })
}


