
Mpar <- 0.1
Mpow <- 0.1 
MaxAge <- as.integer(-log(0.01)/Mpar)
Ages <- 0:MaxAge
Fpar <- 0.1  

Linf <- 100 
CVLinf <- 0.1 
kpar <- 0.075
SL50 <- 25
SL95 <- 30
L50 <- 65
L95 <- 75 
FecB <- 3 

SDLinf <- CVLinf * Linf
NGTG <- 200 
Linc <- 5
MaxSD <- 3 
Walpha <- 0.0001
Wbeta <- 3 

LenBins <- seq(from=0, by=Linc, to=Linf + MaxSD * SDLinf)
LenMids <- seq(from=LenBins[1] + 0.5*Linc, by=Linc, length=length(LenBins)-1)

kslope <- 0.002019024

tempFun <- function(kslope) {
# Growth-Type-Group Model Setup 
DiffLinfs <- seq(from=Linf - MaxSD * SDLinf, to=Linf + MaxSD * SDLinf, length=NGTG)
GTGLinfdL <- DiffLinfs[2] - DiffLinfs[1]
RecProbs <- dnorm(DiffLinfs, Linf, SDLinf)/sum(dnorm(DiffLinfs, Linf, SDLinf)) # Recruits normally distributed across GTGs

AgeGTG <- matrix(NA, nrow=length(Ages), ncol=NGTG)
AgeGTG[1,] <- RecProbs

AgeLenMatrix <- sapply(seq_along(DiffLinfs), function (X) DiffLinfs[X] * (1-exp(-kpar*Ages)))
AgeLenMatrix[AgeLenMatrix == 0] <- 1
MatureMat <- sapply(seq_along(DiffLinfs), function (X) 1.0/(1+exp(-log(19)*(AgeLenMatrix[,X]-L50)/(L95-L50)))) # Maturity Schedule
FecMat <- (AgeLenMatrix^FecB) * MatureMat
FecMat <- FecMat/max(FecMat)

MVec <- Mpar * (Linf/AgeLenMatrix[,round(NGTG/2,0)])^Mpow # M for mean growth type

MMat <- matrix(NA, nrow=length(MVec), ncol=NGTG)
for (X in 1:length(MVec)){
  for (Y in 1:NGTG) {
    MMat[X,Y] <- MVec[X] + kslope*(DiffLinfs[Y] - Linf)
  }
}
 
for (GTG in 1:NGTG) {
  for (Age in 2:length(Ages)) {
    AgeGTG[Age,GTG] <- AgeGTG[Age-1,GTG] * exp(-MMat[Age-1,GTG])
  }
}

TotalFecund <- FecMat * AgeGTG
Fitness <- apply(TotalFecund, 2, sum)
FitnessPR <- Fitness/RecProbs

# Make Size Comp 
SizeComp <- matrix(0, nrow=length(LenMids), ncol=NGTG)
for (GTG in 1:NGTG) {
  for (Age in 1:length(Ages)) {
    currLen <- AgeLenMatrix[Age,GTG]
	Ind <- min(which(LenBins > currLen))
	SizeComp[Ind-1,GTG] <- SizeComp[Ind-1,GTG] + AgeGTG[Age,GTG]
  }
}
SizeMat <- apply(SizeComp, 1, sum)
# plot(LenMids, SizeMat)

Pen <- 0
if (min(MMat) < 0 ) Pen <- 1E9 
SSQ <- sum((FitnessPR - median(FitnessPR))^2)
print(SSQ)
return(SSQ + Pen)

}

optimise(tempFun, interval=c(0, 0.1))


lines(FitnessPR)

SaveMMat[,1] 
MMat[,1]