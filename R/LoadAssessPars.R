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


LoadAssessPars <- function(PathtoAssessFile="~/PathToAssessFile", AssessParFileName="AssessPars", AssessParExt=".csv", ind=1, LenMids, IncludeDatFile=FALSE, LenDatType=list("raw", "comp"), LHR_OverRide=list(), OptMslope=FALSE, PredictMslope=TRUE, OverMslope=NULL, datCol=1) {
  
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
	By     <- LHR_OverRide$By

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
	if (is.null(By))  By     <- as.numeric(Dat["DatLinc",ind+1])
	
	relL50 <- L50/Linf
	relL95 <- L95/Linf
	
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
    
    AssessPars <- list(MK=MK, Linf=Linf, CVLinf=CVLinf, SDLinf=SDLinf, L50=L50, L95=L95, relL50=relL50, relL95=relL95, Walpha=Walpha, Wbeta=Wbeta, FecB=FecB, Mpow=Mpow,               NGTG=NGTG, GTGLinfdL=GTGLinfdL, MaxSD=MaxSD, SL50Min=SL50Min, SL50Max=SL50Max, DeltaMin=DeltaMin, DeltaMax=DeltaMax, DiffLinfs=DiffLinfs)
					
	# Read in data file if it exists
	if (!IncludeDatFile & length(LenMids) < 1) stop("Require either length data file or length bin mid-points")
	if (IncludeDatFile) {
	  DataFile <- Dat["DataFile",ind+1]
	  if (is.na(DataFile) | is.null(DataFile) | length(DataFile) <1) stop("Data file not specified")
      WD <- getwd()
	  setwd(PathtoAssessFile)
	  readLenDat <- try(read.csv(DataFile, header=FALSE), silent=TRUE)
	  if (class(readLenDat) == "try-error") {
	    print(paste("File", DataFile, "not found"))
		print("Files in current directory are:")
		print(list.files())
		stop()
	  }

	  if (ncol(readLenDat) == 1 & tolower(LenDatType) == "raw") {
	    rawLenDat <- as.vector(unlist(readLenDat))
		Max <- max(max(rawLenDat) * 1.1, Linf*1.25)
	    LenBins <- seq(from=0, to=Max, by=By) 
	    LenMids <- seq(from=LenBins[1] + 0.5*By, by=By, length=length(LenBins)-1)
		
		LenDat <- as.vector(table(cut(rawLenDat, LenBins)))
		AssessPars$LenDat <- LenDat
	  }
	 if (ncol(readLenDat) == 1 & tolower(LenDatType) != "raw") stop("Length data file only has one column but data type is not specified as 'raw'")
	 if (ncol(readLenDat) > 1 & tolower(LenDatType) != "comp") stop("Length data file only has two columns but data type is not specified as 'comp'") 
	 Years <- NULL
	 Ind <- 1:nrow(readLenDat)
	 if (ncol(readLenDat) >= 2 & tolower(LenDatType) == "comp") {
	   if(is.na(readLenDat[1,1])) { # time-series info?
	     Years <- readLenDat[1,2:ncol(readLenDat)]
		 Ind <- 2:nrow(readLenDat)
	   }
	   LenMids <- readLenDat[Ind,1]
	   By <- LenMids[2] - LenMids[1]
	   X <- 2:length(LenMids)
	   if (all(LenMids[X] - LenMids[X-1] == By) == FALSE) stop("Length classes not equidistant")
	    # add zeros if minimum length mids is not 0
	   LenDat <- readLenDat[Ind,datCol+1]
	   
	   if (LenMids[1] - 0.5*By != 0) {
	     lenSeq <- seq(from=By, to=LenMids[1] -By, by=By) 
		 Zeros <- rep(0, length(lenSeq))
		 LenMids <- c(lenSeq,LenMids)
		 LenDat <- c(Zeros, LenDat)
	   }
	   
	   AssessPars$LenDat <- LenDat
	 }
	}
	
	AssessPars$LenMids <- LenMids
	AssessPars$Linc <- By <- LenMids[2] - LenMids[1]
	AssessPars$LenBins <- seq(from=LenMids[1] - 0.5*By, by=By, length=length(LenMids)+1)
	AssessPars$AssessOpt <- TRUE
	AssessPars$Years <- Years
    
    # Optimise Mslope 
	if (length(OverMslope) < 1 & PredictMslope == FALSE & OptMslope==FALSE) stop("Mslope parameters not set")
	if (length(OverMslope) > 0) {
	  print(paste("Mslope set at ", OverMslope))
	  AssessPars$Mslope <- OverMslope
	  OptMslope <- PredictMslope <- FALSE
	} 
    if (PredictMslope) {
	  print("Predicting Mslope")
	  AssessPars$Mslope <- Mslope <- PredictMSlope(AssessPars)
	  OptMslope <- FALSE
	}
	if (OptMslope) {
	  print("Optimising for Mslope - this may take a short while...")
	  AssessPars$Mslope <- Mslope <- exp(optimise(OptimiseFitness, interval=log(c(0.0001, 0.1)), SimPars=AssessPars, Function=SimMod_LHR)$minimum)
	  # M per GTG 
	  AssessPars$MKGTG <- MKGTG <- MK + Mslope*(DiffLinfs-Linf)
	  if (min(AssessPars$MKGTG) <= 0) {
	    print("MKgtg is negative for some GTG. Try change NGTG")
	    # readline("********** WARNING - Press Enter to continue **********")
	  }
	}
	AssessPars$MKMin <- as.numeric(Dat["MKMin",ind+1])
	AssessPars$MKMax <- as.numeric(Dat["MKMax",ind+1])
	AssessPars$LinfMin <- as.numeric(Dat["LinfMin",ind+1])
	AssessPars$LinfMax <- as.numeric(Dat["LinfMax",ind+1])
	
  }
  
  if(AssessParExt != ".csv") stop("Unrecognized file extension")
  print("Assessment parameters successfully loaded")
  setwd(WD)
  return(AssessPars)
}
