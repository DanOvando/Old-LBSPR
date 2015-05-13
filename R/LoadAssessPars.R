#' This function reads in a csv file that contains the parameters required for 
#' running the LB-SPR Assessment model.  See LoadAssessPars ReadMe for details
#' of contents of SimParFile.
#' @name LoadAssessPars
#' @title Load LBSPR Assessment Parameters
#' @param PathtoAssessFile an object of class \code{character} containing the
#'   full path to the Assessment Parameter file
#' @param AssessParFileName an object of class \code{character} containing the
#'   name of Assessment Parameter file
#' @param AssessParExt an object of class \code{character} containing the file 
#'   extension of Assessment Parameter file (default is ".csv" and currently can
#'   not handle anything else)
#' @param ind an optional \code{numeric} value indicating which column contains 
#'   the parameters (default is 1).
#' @param LenMids A vector of the midpoints of the length classes for the observed length data.
#' @return Returns a list of assessment parameters (\code{AssessPars}) that is
#'   used in other functions.
#' @author Adrian Hordyk
#' @export


LoadAssessPars <- function(PathtoAssessFile="~/PathToAssessFile", AssessParFileName="AssessPars", AssessParExt=".csv", ind=1, LenMids, LHR_OverRide=list()) {
  
  OptimiseFitness <- function(logMslope, SimPars, Function) {
    Mslope <- exp(logMslope) #+ 0.000000001
    SimPars$Mslope <- Mslope
    return(Function(SimPars=SimPars)$ObjFun)
  }
  if(AssessParExt == ".csv") {
    Dat <- read.csv(file.path(PathtoAssessFile, paste0(AssessParFileName, AssessParExt)), as.is=TRUE)
    row.names(Dat) <- Dat[,1]
    
    # Load new parameters 
	MK     <- LHR_OverRide$MK
	Linf   <- LHR_OverRide$Linf 
	CVLinf <- LHR_OverRide$CVLinf
	L50    <- LHR_OverRide$L50
	L95    <- LHR_OverRide$L95
	Walpha <- LHR_OverRide$Walpha
	Wbeta  <- LHR_OverRide$Wbeta
	FecB   <- LHR_OverRide$FecB
	Mpow   <- LHR_OverRide$Mpow
	NGTG   <- LHR_OverRide$NGTG

    if (is.null(MK)) 	MK   <- as.numeric(Dat["MK",ind+1])
    if (is.null(Linf)) 	Linf   <- as.numeric(Dat["Linf",ind+1])
    if (is.null(CVLinf))CVLinf <- as.numeric(Dat["CVLinf",ind+1])
    if (is.null(L50)) 	L50    <- as.numeric(Dat["L50",ind+1])
    if (is.null(L95)) 	L95    <- as.numeric(Dat["L95",ind+1])
    if (is.null(Walpha))Walpha <- as.numeric(Dat["Walpha",ind+1])
    if (is.null(Wbeta)) Wbeta  <- as.numeric(Dat["Wbeta",ind+1])
    if (is.null(FecB)) 	FecB   <- as.numeric(Dat["FecB",ind+1])
    if (is.null(Mpow)) 	Mpow   <- as.numeric(Dat["Mpow",ind+1])
    if (is.null(NGTG)) 	NGTG   <- as.numeric(Dat["NGTG",ind+1])

    SDLinf <- as.numeric(CVLinf * Linf)
    MaxSD  <- as.numeric(Dat["MaxSD",ind+1])
    GTGLinfdL <- ((Linf + MaxSD * SDLinf) - (Linf - MaxSD * SDLinf))/(NGTG-1)
  
	DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, length=NGTG)

    SL50Min <- as.numeric(Dat["SL50Min", ind+1])
    if (length(SL50Min) < 1 | is.na(SL50Min)) SL50Min <- 1
    SL50Max <- as.numeric(Dat["SL50Max", ind+1])
    if (length(SL50Max) < 1| is.na(SL50Max)) SL50Max <- 0.95 * Linf
    DeltaMin <- as.numeric(Dat["DeltaMin", ind+1])
    if (length(DeltaMin) < 1| is.na(DeltaMin)) DeltaMin <- 0.01
    DeltaMax <- as.numeric(Dat["DeltaMax", ind+1])
    if (length(DeltaMax) < 1| is.na(DeltaMax)) DeltaMax <- 0.5 * Linf
    
    AssessPars <- list(MK=MK, Linf=Linf, CVLinf=CVLinf, SDLinf=SDLinf, L50=L50, L95=L95, 
                    Walpha=Walpha, Wbeta=Wbeta, FecB=FecB, Mpow=Mpow, 
                    NGTG=NGTG, GTGLinfdL=GTGLinfdL, MaxSD=MaxSD, SL50Min=SL50Min, 
                    SL50Max=SL50Max, DeltaMin=DeltaMin, DeltaMax=DeltaMax, DiffLinfs=DiffLinfs)
	
	AssessPars$LenMids <- LenMids
	AssessPars$Linc <- BY <- LenMids[2] - LenMids[1]
	AssessPars$LenBins <- seq(from=LenMids[1] - 0.5*BY, by=BY, length=length(LenMids)+1)
	AssessPars$AssessOpt <- TRUE

	AssessPars$Mslope <- exp(optimise(OptimiseFitness, interval=log(c(0.000001, 0.1)), SimPars=AssessPars, Function=SimMod_LHR)$minimum) # Run Sim Model and fit optimal MSlope
  }
  
  if(AssessParExt != ".csv") stop("Unrecognized file extension")
  print("Assessment parameters successfully loaded")
  return(AssessPars)
}
