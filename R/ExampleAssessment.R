
# # Example Assessment
# tempfunction <- function() {
# ADMBDir <- "E:/Dropbox/Projects/LBSPRModel/"
# setwd("E:/Dropbox/Projects/LBSPRModel/LBSPR")
# 
# # Load Assessment Parameters
# AssessPars <- LoadAssessPars(PathtoAssessFile="data/", AssessParFileName="ExampleAssessPars", AssessParExt=".csv", ind=1) 
# print(t(AssessPars))
# 
# # Load Raw Length Data and create Length Compostion
# LenData <- MakeLengthComp(PathtoLenDat="data/", LenDatFileName="ExampleRawData", LenDatExt=".csv", Header=FALSE, Linc=5, AssessPars, Multi=1.25)
# print(LenData)
# 
# # Run Example Assessment
# RunAssess <- RunLBSPRAssess(AssessPars, LenFreq=LenData$LenFreq, LenMids=LenData$LenMids, ADMBDir, ExName="lbspr", showOutput=FALSE, MaxCount=5) 
# 
# }