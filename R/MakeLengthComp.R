#' \code{MakeLengthCompFun} is used to create a length frequency composition 
#' from raw length data.
#' @name MakeLengthComp
#' @aliases MakeLengthCompFun
#' @title Create Length Frequency Composition
#' @param PathtoLenDat An object of class \code{character} that contains the 
#'   full path to the raw length data file.
#' @param LenDatFileName An object of class \code{character} that contains the 
#'   name of the length data file.
#' @param LenDatExt An object of class \code{character} that contains the 
#'   extension of the length data file (default is ".csv" and currently does not
#'   handle anything else).
#' @param Header An \code{logical} indicating if the first row of the CSV file
#'   is a header (Default is FALSE).
#' @param Linc An \code{numeric} value indicating the width of the length class
#'   (in same units as \code{AssessPars}).
#' @param AssessPars A \code{list} containing the parameters required for the
#'   LB-SPR assessment. Created by \code{LoadAssessPars} function.
#' @param Multi An \code{numeric} object that specifies the maximum length class as a mulitiple of \code{AssessPars$Linf}.  
#' @param sep Decimals are indicated by \code{"period"} (default) or \code{"comma"} 
#' @return A \code{list} containing the Length Bins (\code{LenBins}), the Length
#'   Frequency (\code{LenFreq}) and the midpoints of the length bins 
#'   (\code{LenMids}).
#' @author Adrian Hordyk 
#' @export


MakeLengthComp <- function(PathtoLenDat, LenDatFileName, LenDatExt=".csv", Header=FALSE, Linc=5, AssessPars, Multi=1.25, sep="period") {
  
  if (sep=="period") LenDat <- read.csv(file.path(PathtoLenDat, paste0(LenDatFileName, LenDatExt)), header=Header)
  if (sep=="comma") LenDat <- read.csv2(file.path(PathtoLenDat, paste0(LenDatFileName, LenDatExt)), header=Header)
  # add error checkers etc
  NCol <- ncol(LenDat)
  if (NCol > 2) stop("Too many columns in data file")
  if (NCol == 2) {
    print("Length frequency data - NOT CURRENTLY SUPPORTED")
    return(NULL)
    # LenFreq <- LenDat[,2]
    # LenMids <- LenDat[,1] # add some checks for smallest length bin 
  }
  if (NCol == 1) {
    print("Raw length frequency data")
    # Bin Raw Data 
    LenBins <- seq(from=0, to=AssessPars$Linf*Multi, by=Linc)
    LenFreq <- as.vector(table(cut(unlist(LenDat), LenBins)))
    LenMids <- seq(from=LenBins[2]-0.5*Linc, by=Linc, length=length(LenBins)-1)
  } 
  Output <- NULL 
  Output$LenBins <- LenBins
  Output$LenFreq <- LenFreq
  Output$LenMids <- LenMids
  return(Output)
} 
