#' Function to simulate length-structured growth-type-group (GTG) model to 
#' generate size equilibrium composition of population and catch, as well as SPR
#' of stock and relative yield.  More details to be added.
#' @name GenSPRYPR
#' @title Generate size structure using GTG model
#' @param SimPars An object of class \code{list} that contains all parameters
#'   required to run GTG model.  Full description of model to be added at later
#'   date.
#' @return  To add details later.
#' @author Adrian Hordyk
#' @seealso \code{\link{}}
#' @export
#' @examples      
#' \dontrun{
#' 
#' }  

EqGTGModelStn <- function(kslope, SimPars, Mpar) {
  with(SimPars, {
  Fpar <- Mpar * FM
  kpar <- Mpar/ MK
  
  SDLinf <- CVLinf * Linf # Standard Deviation of Pop. Linf 
  
  # Growth-Type-Group Model Setup 
  DiffLinfs <- seq(from=Linf - MaxSD * SDLinf, to=Linf + MaxSD * SDLinf, length=NGTG)
  GTGLinfdL <- DiffLinfs[2] - DiffLinfs[1]
  RecProbs <- dnorm(DiffLinfs, Linf, SDLinf)/sum(dnorm(DiffLinfs, Linf, SDLinf)) # Recruits normally distributed across GTGs
  
  # Set up Length Bins of Population Model 
  LenBins <- seq(from=0, by=Linc, to=Linf + MaxSD * SDLinf)
  LenMids <- seq(from=LenBins[1] + 0.5*Linc, by=Linc, length=length(LenBins)-1)
  
  MatLen <- 1.0/(1+exp(-log(19)*(LenMids-L50)/(L95-L50))) # Maturity Schedule 
  Weight <- Walpha * LenMids^Wbeta # Weight-at-length 
  FecLen <- MatLen * LenMids^FecB # Relative Fecundity at Length 
  FecLen <- FecLen/max(FecLen) # Divide by maximum to make it more useful for plotting. Can add extra parameter(s) so that fecundity is absolute rather than relative (makes no difference to SPR)
  
  VulLen <- 1.0/(1+exp(-log(19)*(LenBins-SL50)/(SL95-SL50))) # Selectivity-at-Length
  SelLen <- VulLen
  # # Minimum Legal Length 
  # if (length(MLL) > 0 & MLL > 0 ) {
  # Legal <- rep(0, length(LenBins))
  # Legal[LenBins >= MLL] <- 1 
  # SelLen <- VulLen * Legal
  # } 
  
  MVec <- Mpar * (Linf/(LenBins+0.5*Linc))^Mpow # Vector of natural mortality for mean GTG
  MMat <- matrix(MVec, nrow=length(MVec), ncol=NGTG) # Matrix of M for each GTG
  tempFun <- function(X) MVec + kslope*(DiffLinfs[X] - Linf)
  MMat <- sapply(seq_along(DiffLinfs), function (X) tempFun(X))
  
  FVec <- Fpar * SelLen # Vector of F for each size class
  ZMat <- MMat + FVec # Matrix of total mortality for each GTG
  
  
  # Set Up Empty Matrices 
  NPRFished <- NPRUnfished <- matrix(0, nrow=length(LenBins), ncol=NGTG) # number-per-recruit at length
  NatLUnFishedPop <- NatLFishedPop <- NatLUnFishedCatch <- NatLFishedCatch <- FecGTGUnfished <- matrix(0, nrow=length(LenMids), ncol=NGTG) # number per GTG in each length class 
    
  EPR_GTG <- matrix(NA, nrow=NGTG, ncol=2)
  YPR <- rep(NA, NGTG)
 
 # Loop over Growth-Type-Groups 	
  for (GTG in 1:NGTG) {
    NPRUnfished[1, GTG] <- NPRFished[1, GTG] <- RecProbs[GTG]
	options(warn=-1) # turn off warnings
    ALen <- -(log(1-LenBins/DiffLinfs[GTG]))/kpar # Calculate Age for length
	options(warn=0) # turn back on
    ALen[is.na(ALen >= 0) ] <- NA # Get rid of NaNs 
	HighAge <- which.max(ALen)
	LastInc <- ALen[HighAge-1] - ALen[HighAge-2]
	if (ALen[HighAge] != Inf) {
	  LastInc <- ALen[HighAge] - ALen[HighAge-1]
      ALen[HighAge+1] <- ALen[HighAge] + 2 * LastInc # Make a maximum age for last length class 
	}    
    if (ALen[HighAge] == Inf) {
	  LastInc <- ALen[HighAge-1] - ALen[HighAge-2]
	  ALen[HighAge] <- ALen[HighAge-1] + 2 * LastInc # Make a maximum age for last length class 
	}
	for (L in 2:(HighAge)) {
      if (LenBins[L] < DiffLinfs[GTG]) {
        NPRUnfished[L, GTG] <- NPRUnfished[L-1, GTG] * exp(-MMat[L-1, GTG] * (ALen[L+1] - ALen[L]))
        NPRFished[L, GTG] <- NPRFished[L-1, GTG] * exp(-ZMat[L-1, GTG] * (ALen[L+1] - ALen[L]))
      }
    }
    for (L in 1:length(LenMids)) {
      NatLUnFishedPop[L, GTG] <- (NPRUnfished[L,GTG] - NPRUnfished[L+1,GTG])/MMat[L, GTG]
      NatLFishedPop[L, GTG] <- (NPRFished[L,GTG] - NPRFished[L+1,GTG])/ZMat[L, GTG]  
	  FecGTGUnfished[L, GTG] <- NatLUnFishedPop[L, GTG] * FecLen[L]
    }	
	
    NatLUnFishedCatch[,GTG] <- NatLUnFishedPop[, GTG] * 1.0/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50))) 
    NatLFishedCatch[,GTG]  <- NatLFishedPop[, GTG] * 1.0/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50)))
	   
    # Eggs-per-recruit for each GTG 
    EPR_GTG[GTG,1] <- sum(NatLUnFishedPop[, GTG] * FecLen)
    EPR_GTG[GTG,2] <- sum(NatLFishedPop[, GTG] * FecLen)
      
    # YPR 
    YPR[GTG] <- sum(NatLFishedPop[, GTG]  * Weight * 1.0/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50)))) * FM
  }	

  # Calc Unfished Fitness 
  Fit <- apply(FecGTGUnfished, 2, sum, na.rm=TRUE) # Total Fecundity per Group
  FitPR <- Fit/RecProbs # Fitness per-recruit
  ObjFun <- sum((FitPR - median(FitPR, na.rm=TRUE))^2, na.rm=TRUE) # This needs to be minimised to make fitness approximately equal across GTG - by adjusting kslope
  Pen <- 0; if (min(MMat) < 0 ) Pen <- (1/abs(min(MMat)))^2 * 1E5 # Penalty for optimising kslope   
  ObjFun <- ObjFun + Pen
  # print(cbind(kslope, ObjFun, Pen))
  
  # Calc SPR
  EPR0 <- apply(EPR_GTG,2,sum)[1]
  EPRf <- apply(EPR_GTG,2,sum)[2]
  SPR <- EPRf/EPR0  
  
  # Equilibrium Relative Recruitment 
  reca <- recK/EPR0
  recb <- (reca * EPR0 - 1)/(R0*EPR0)
  RelRec <- max(0, (reca * EPRf-1)/(recb*EPRf))
    
  Yield <- sum(YPR) * RelRec

  # Expected Length Structure of Catch 
  ExpectedLenCatchFished <- apply(NatLFishedCatch, 1, sum)/sum(apply(NatLFishedCatch, 1, sum))
  ExpectedLenPopFished <- apply(NatLFishedPop, 1, sum)/sum(apply(NatLFishedPop, 1, sum))
  ExpectedLenCatchUnfished <- apply(NatLUnFishedCatch, 1, sum)/sum(apply(NatLUnFishedCatch, 1, sum))
  ExpectedLenPopUnfished <- apply(NatLUnFishedPop, 1, sum)/sum(apply(NatLUnFishedPop, 1, sum))
  
  Output <- NULL 
  Output$SPR <- SPR
  Output$Yield <- Yield 
  Output$ExpLenCatchFished <- ExpectedLenCatchFished
  Output$ExpLenPopFished <- ExpectedLenPopFished
  Output$ExpLenCatchUnfished <- ExpectedLenCatchUnfished
  Output$ExpLenPopUnfished <- ExpectedLenPopUnfished
  Output$LenBins <- LenBins
  Output$LenMids <- LenMids
  Output$NGTG <- NGTG
  Output$GTGdL <- GTGLinfdL
  Output$DiffLinfs <- DiffLinfs
  Output$RecProbs <- RecProbs
  Output$Weight <- Weight
  Output$FecLen <- FecLen 
  Output$MatLen <- MatLen 
  Output$SelLen <- SelLen
  Output$MMat <- MMat 
  Output$FVec <- FVec 
  Output$ZMat <- ZMat
  Output$ObjFun <- ObjFun 
  Output$Pen <- Pen
  Output$Mpar <- Mpar 
  Output$kpar <- kpar 
  
  return(Output)
  })
}

OptimiseFitness <- function(logkslope, SimPars, Function, Mpar=NULL) {
  Kslope <- exp(logkslope) + 0.0000001
  return(Function(kslope=Kslope, SimPars=SimPars, Mpar=Mpar)$ObjFun)
}  


TestKslopeFun <- function(Function) {
  SimPars$NGTG <- as.integer(runif(1, 20, 60))
  SimPars$Linf <- runif(1, 40, 300) 
  relL50 <- runif(1, 0.4, 0.9)
  SimPars$L50 <- relL50 * SimPars$Linf 
  SimPars$L95 <- SimPars$L50 + 0.1 * SimPars$Linf 
  SimPars$Mpow <- runif(1, 0, 0.2)
  Mpar <- runif(1, 0.05, 0.5)
  SimPars$MK <- runif(1, 0.5, 3)
  
  SimPars$Linc <- 5
  if (SimPars$Linf/5 < 30) SimPars$Linc <- 1
  
  SimPars$FM <- 0 
  Runopt <- optimise(OptimiseFitness, interval=log(c(0.000001, 0.2)), SimPars=SimPars, Function= Function, Mpar=Mpar)
  kslope <- exp(Runopt$minimum)
  RunFunc <- Function(kslope=kslope, SimPars=SimPars, Mpar=Mpar)
  OptSuccess <- TRUE
  if (RunFunc$Pen) OptSuccess <- FALSE
  
  Output <- c(kslope, SimPars$NGTG, SimPars$Linf, SimPars$L50/SimPars$Linf, RunFunc$Mpar, SimPars$Mpow, SimPars$MK, RunFunc$kpar, OptSuccess)
  return(Output)
}

# Life History Ratio Model Test 
N <- 100
Run <- sapply(1:N, function(X) {
  print(X)
  TestKslopeFun(EqGTGModelStn)
 })  
 
Fails <- which(Run[9,] != 1)
if (length(Fails) > 0) Run <- Run[,-Fails]

Run <- t(Run)
colnames(Run) <- c("kslope", "NGTG", "Linf", "L50/Linf", "Mpar", "Mpow", "MK", "kpar", "OptSuccess")

# dput(Run, file="Run2.dat")
Run <- dget(file="Run2.dat")

par(mfrow=c(3,3), cex.lab=1.5, bty="l", mar=c(5,5,2,2), oma=c(1,1,1,1))
plot(Run[,2], Run[,1], xlab="NGTG", ylab="kslope")
plot(Run[,3], Run[,1], xlab="Linf", ylab="kslope")
plot(Run[,4], Run[,1], xlab="L50/Linf", ylab="kslope")
plot(Run[,5], Run[,1], xlab="Mpar", ylab="kslope")
plot(Run[,6], Run[,1], xlab="Mpow", ylab="kslope")
plot(Run[,8], Run[,1], xlab="Kpar", ylab="kslope")
plot(Run[,7], Run[,1], xlab="MK", ylab="kslope")

require(akima); require(fields)
par(mfrow=c(3,5), mar=c(4,5,2,5), oma=c(1,1,1,1))
Breaks <- seq(from=0, to=0.02, length=21)
Mat <- matrix(c(2,3,2,4,2,5,2,6,2,7,3,4,3,5,3,6,3,7,4,5,4,6,4,7,5,6,5,7,6,7), ncol=2, byrow=TRUE)
for (X in 1:nrow(Mat)) {
  x <- Run[,Mat[X,1]]
  y <- Run[,Mat[X,2]]
  z <- Run[,1]
  s <- interp(x,y,z)
  image.plot(s, xlab=colnames(Run)[Mat[X,1]], ylab=colnames(Run)[Mat[X,2]],  nlevel=20, breaks=Breaks, cex.lab=1.5)
}

# install.packages("akima")
# install.packages("fields")

# Life History Ratio Model Test 
N <- 1000
Run <- sapply(1:N, function(X) {
  print(X)
  TestKslopeFun(GenSPRYPR)
 })  
 
Fails <- which(Run[7,] != 1)
if (length(Fails) > 0) Run <- Run[,-Fails]

Run <- t(Run)
colnames(Run) <- c("kslope", "NGTG", "Linf", "L50/Linf",  "Mpow", "MK", "OptSuccess")

dput(Run, file="Run_LHMod.dat")

par(mfrow=c(2,3), cex.lab=1.5, bty="l", mar=c(5,5,2,2), oma=c(1,1,1,1))
plot(Run[,2], Run[,1], xlab="NGTG", ylab="kslope")
plot(Run[,3], Run[,1], xlab="Linf", ylab="kslope")
plot(Run[,4], Run[,1], xlab="L50/Linf", ylab="kslope")
plot(Run[,5], Run[,1], xlab="Mpow", ylab="kslope")
plot(Run[,6], Run[,1], xlab="MK", ylab="kslope")

require(akima); require(fields)
par(mfrow=c(3,5), mar=c(4,5,2,5), oma=c(1,1,1,1))
Breaks <- seq(from=0, to=0.02, length=21)
Mat <- matrix(c(2,3,2,4,2,5,2,6,2,7,3,4,3,5,3,6,3,7,4,5,4,6,4,7,5,6,5,7,6,7), ncol=2, byrow=TRUE)
for (X in 1:nrow(Mat)) {
  x <- Run[,Mat[X,1]]
  y <- Run[,Mat[X,2]]
  z <- Run[,1]
  s <- interp(x,y,z)
  image.plot(s, xlab=colnames(Run)[Mat[X,1]], ylab=colnames(Run)[Mat[X,2]],  nlevel=20, breaks=Breaks, cex.lab=1.5)
}


# Fix up code and re-compile package 

# write code in other directory to run models
 
# confirm two models give same answer 

# test for relationship between kslope and Linf & MK in LHR model 

# formalise results and continue to write manuscript




RndFun <- function(Par, ParVal, SimPars, Function, Mpar) {
  SimPars[Par] <- ParVal
  SimPars$L50 <- 0.6 * SimPars$Linf 
  SimPars$L95 <- SimPars$L50 + 0.1 * SimPars$Linf
  if (Par == "kpar") SimPars$MK <- Mpar/SimPars$kpar 
  Runopt <- optimise(OptimiseFitness, interval=log(c(0.000001, 0.1)), SimPars=SimPars, Function= Function, Mpar=Mpar)
  kslope <- exp(Runopt$minimum)
  RunFunc <- Function(kslope=kslope, SimPars=SimPars, Mpar=Mpar)
  OptSuccess <- TRUE
  if (RunFunc$Pen) OptSuccess <- FALSE
  Output <- data.frame(kslope, Par, ParVal, OptSuccess)
}

SimPars$NGTG <- 50
SimPars$Linf <- 100
SimPars$L50 <- 60
SimPars$L95 <- 70 
SimPars$Mpow <- 0
Mpar <- 0.1
SimPars$MK <- 1.5

Function <- EqGTGModelStn
Par <- "Linf" 
ParValVec <- seq(from=50, to=300, by=25)
Cnt <- 0
parVec <- c(0.1, 0.2, 0.3, 0.4, 0.5)
SaveMat <- matrix(NA, nrow=length(ParValVec), ncol=length(parVec))
for (Mpar in parVec) {
  Cnt <- Cnt + 1 
  VarLinf <- sapply(seq_along(ParValVec), function(X) RndFun(Par, ParValVec[X], SimPars=SimPars, Function=Function, Mpar=Mpar)) 
  VarLinf <- t(VarLinf)
  SaveMat[,Cnt] <- unlist(VarLinf[,1])
  print(Cnt)
}

matplot(x=ParValVec, y=SaveMat)

# NGTGVec <- ParValVec
# SlopeNGTG <- SaveMat

matplot(x=LinfVec, y=SlopeLinf)
matplot(x=kparVec, y=Slopekpar)
matplot(x=MKVec, y=SlopeMK)
matplot(x=MpowVec, y=SlopeMpow)
matplot(x=NGTGVec, y=SlopeNGTG)