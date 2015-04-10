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
#' @param sep an object of class \code{character} that specifies if decimals points in csv are indicated by \code{comma} (default) or \code{period} (using read.csv2).
#' @return Returns a list of assessment parameters (\code{AssessPars}) that is
#'   used in other functions.
#' @author Adrian Hordyk
#' @seealso \code{\link{}}
#' @export
#' @examples      
#' \dontrun{
#' 
#' }  

LoadAssessPars <- function(PathtoAssessFile="~/PathToAssessFile", AssessParFileName="AssessPars", AssessParExt=".csv", ind=1, sep="comma") {
  
  if(AssessParExt == ".csv") {
    Dat <- read.csv(file.path(PathtoAssessFile, paste0(AssessParFileName, AssessParExt)))
    row.names(Dat) <- Dat[,1]
    
    # Load new parameters 
    Mk     <- Dat["MK",ind+1]
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
    
    SL50Min <- Dat["SL50Min", ind+1]
    if (length(SL50Min) < 1) SL50Min <- 10
    SL50Max <- Dat["SL50Max", ind+1]
    if (length(SL50Max) < 1) SL50Max <- 0.95 * Linf
    DeltaMin <- Dat["DeltaMin", ind+1]
    if (length(DeltaMin) < 1) DeltaMin <- 0.01
    DeltaMax <- Dat["DeltaMax", ind+1]
    if (length(DeltaMax) < 1) DeltaMax <- 0.1 * Linf
    
    AssessPars <- list(Mk=Mk, Linf=Linf, CVLinf=CVLinf, L50=L50, L95=L95, 
                    Walpha=Walpha, Wbeta=Wbeta, FecB=FecB, Mpow=Mpow, 
                    NGTG=NGTG, GTGLinfdL=GTGLinfdL, MaxSD=MaxSD, SL50Min=SL50Min, 
                    SL50Max=SL50Max, DeltaMin=DeltaMin, DeltaMax=DeltaMax)
  }
  
  if(AssessParExt != ".csv") stop("Unrecognized file extension")
  print("Assessment parameters successfully loaded")
  print(as.data.frame(AssessPars))
  return(AssessPars)
}