#' This function reads in a csv file that contains the parameters for simulating the LB-SPR model.  See SimData ReadMe for details of
#' contents of SimParFile.  
#' @name LoadSimPars  
#' @title Load LBSPR Simulation Parameters
#' @param PathtoSimFile an object of class \code{character} containing the full path to the Simulation Parameter file
#' @param SimParFileName an object of class \code{character} containing the name of Simulation Parameter file
#' @param SimParExt an object of class \code{character} containing the file extension of Simulation Parameter file (default is ".csv")
#' @param ind an optional \code{numeric} value indicating which column contains the parameters (default is 1)
#' @return Returns a list of simulation parameters
#' @author Adrian Hordyk 
#' @export

LoadSimPars <- function(PathtoSimFile="~/PathToSimFile", SimParFileName="SimData", SimParExt=".csv", ind=1) {
  
  if(SimParExt == ".csv") {
    Dat <- read.csv(file.path(PathtoSimFile, paste0(SimParFileName, SimParExt)))
	row.names(Dat) <- Dat[,1]
		
	# Load new parameters 
	MK     <- Dat["MK",ind+1]
	Linf   <- Dat["Linf",ind+1]
	CVLinf <- Dat["CVLinf",ind+1]
	L50    <- Dat["L50",ind+1]
	L95    <- Dat["L95",ind+1]
	Walpha <- Dat["Walpha",ind+1]
	Wbeta  <- Dat["Wbeta",ind+1]
	FecB   <- Dat["FecB",ind+1]
	Mpow   <- Dat["Mpow",ind+1]
	NGTG   <- Dat["NGTG",ind+1]
	SDLinf <- CVLinf * Linf
	MaxSD  <- Dat["MaxSD",ind+1]
	GTGLinfdL <- ((Linf + MaxSD * SDLinf) - (Linf - MaxSD * SDLinf))/(NGTG-1);	

	Linc   <- Dat["Linc",ind+1]
	R0     <- Dat["R0",ind+1]
	recK   <- Dat["recK",ind+1]
	SL50   <- Dat["SL50",ind+1]
	SL95   <- Dat["SL95",ind+1]
	FM     <- Dat["FM",ind+1]
	SPR    <- Dat["SPR",ind+1]
	MLL	   <- Dat["MLL",ind+1]
	DisMortFrac <- Dat["DisMortFrac",ind+1] 
	SimPars <- list(MK=MK, Linf=Linf, CVLinf=CVLinf, L50=L50, L95=L95, 
					Walpha=Walpha, Wbeta=Wbeta, FecB=FecB, Mpow=Mpow, 
					NGTG=NGTG, GTGLinfdL=GTGLinfdL, MaxSD=MaxSD, Linc=Linc, R0=R0, recK=recK, 
					SL50=SL50, SL95=SL95, FM=FM, SPR=SPR, MLL=MLL, DisMortFrac=DisMortFrac)
	SimPars$kslope <- as.numeric(PredictKSlope(SimPars))
  }
  
  if(SimParExt != ".csv") stop("Unrecognized file extension")
  print("Simulation parameters successfully loaded")
  print(as.data.frame(SimPars))
  return(SimPars)
}





