#' Function to download the LBSPR binary from online repository and copy to local machine.
#' @name CopyADMBFile  
#' @title Download LBSPR executable. 
#' @param System A \code{character} string of either "Win" for Windows OS compiled executable, or anything else (not blank) for UNIX OS
#' @param DestPath A code{character} string detailing the full path where the binary file is to be downloaded.
#' @author Adrian Hordyk 
#' @export

CopyADMBFile <- function(System="Win", DestPath=NULL) {
  if (System == "Win") {
    destPath <- paste0(DestPath, "/lbspr.exe")
	
	URL <- paste0("https://github.com/AdrianHordyk/LBSPR_ADMB/releases/download/v1.00/lbspr.exe")
    download(URL, destPath, mode="wb")
  }
  if (System != "Win") {
    destPath <- paste0(DestPath, "/lbspr")
	URL <- paste0("https://github.com/AdrianHordyk/LBSPR_ADMB/releases/download/v1.00/lbspr")
    download(URL, destPath, mode="wb")
  }  
}
