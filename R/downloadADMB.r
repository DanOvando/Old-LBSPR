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
	
	URL <- paste0("https://dl.dropboxusercontent.com/u/24856730/LBSPR/release/current/lbspr.exe")
    download(URL, destPath, mode="wb")
  }
  if (System != "Win") {
    destPath <- paste0(DestPath, "/lbspr")
	URL <- paste0("https://dl.dropboxusercontent.com/u/24856730/LBSPR/release/current/lbspr")
    download(URL, destPath, mode="wb")
  }  
}
