#' Function to download the LBSPR binary from GitHub and copy to local machine.
#' @name CopyADMBFile  
#' @title Download LBSPR executable. 
#' @param System A \code{character} string of either "Win" for Windows OS compiled executable, or anything else (not blank) for UNIX OS
#' @param DestPath A code{character} string detailing the full path where the binary file is to be downloaded.
#' @param vers A \code{character} string representing the version. Default value should be the latest release. 
#' @author Adrian Hordyk 
#' @export

CopyADMBFile <- function(System="Win", DestPath=NULL, vers="1.00") {
  if (System == "Win") {
    destPath <- paste0(DestPath, "/lbspr.exe")
	URL <- paste0("https://github.com/AdrianHordyk/LBSPR_ADMB/releases/download/v", vers, "/lbspr.exe")
    download(URL, destPath)
  }
  if (System != "Win") {
    destPath <- paste0(DestPath, "/lbspr")
	URL <- paste0("https://github.com/AdrianHordyk/LBSPR_ADMB/releases/download/v", vers, "/lbspr")
    download(URL, destPath)
  }  
}