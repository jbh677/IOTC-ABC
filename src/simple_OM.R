######################################################
# simple OM for use in ABC estimator #################
######################################################
# R. Hillary, CSIRO O & A 2022 #######################
######################################################

library(FLCore)
library(Rcpp)
library(RcppArmadillo)
source("abc_utils.R")

# compile

sourceCpp("projhy.cpp")
 
######################
# some biology stuff #
######################

nages <- 20
ages <- 1:nages
lbins <- c(seq(20,100,by=10),c(120,140,160,180))
nbins <- length(lbins)-1
mulbins <- mulbins <- 0.5*(lbins[-1]+lbins[-(nbins+1)])
B0 <- 1e+6 # tonnes
hh <- 0.75 # steepness
M <- rep(0.4,nages)
sigmar <- 0.5
initdel <- 1 # starting depletion level
rho <- 0 # recruitment autocorrelation
a1 <- 2 # ref age 1
l1 <- 50 # length at age a1
a2 <- 10 # ref age 2
l2 <- 120 # length at age a2
k <- 0.2
sdl <- 0.15
awl <- 5e-6 # kg units
bwl <- 3
mula <- l1+(l2-l1)*(1-exp(-k*(ages-a1)))/(1-exp(-k*(a2-a1)))
lm50 <- 90
numl <- 9
ml <- mula^numl/(lm50^numl+mula^numl)
#ml <- 1/(1+19^{-(mula-lm50)/(lm95-lm50)})
resb <- get.wt.mat(mula,mulbins,ml,ages,awl,bwl,lm50,numl,sdl)

Biol <- list(B0=B0,
             hh=hh,
             sigmar=sigmar,
             rho=rho,
             initdel=initdel,
             nages=nages,
             M=M,
             mula=mula,
             sdla=sdl,
             wt=resb$wt,
             mat=resb$mat,
             pla=resb$pla)

######################
# some fishery stuff #
######################

nf <- 1 # number of fisheries
self <- matrix(nrow=nf,ncol=nages)
fshare <- rep(1/nf,nf) # relative F share by fishery

# simple double-normal relationship for 1 fleet as example

smax <- 4
sl <- 1.25 
sr <- 25
for(f in 1:nf) {
   
  for(a in 1:nages) self[f,a] <- ifelse(a <= smax,2^{-(a-smax)^2/sl^2},2^{-(a-smax)^2/sr^2})
 
}

Fleet <- list(nf=nf,self=self,fshare=fshare)

############################
# get MSY reference points #
############################

# stochastic or deterministic

msytype <- "deterministic"
hmax <- 0.4
MSYvars <- get.MSY(Biol,Fleet,msytype,hmax)

###################################
# create harvest rate time-series #
###################################

# set up

ny <- 50
yrs <- (2022-ny+1):2022
nits <- 500

# linear increasing from low (0.05) to Hmsy

hy <- seq(0.05,MSYvars$Hmsy,length=ny)
rngseed <- 13
res1 <- projmultiyearh(hy,Biol,Fleet,ny,nits,rngseed)
colnames(res1$dSSB) <- as.character(yrs)
dimnames(res1$C)[[3]] <- as.character(yrs)
dimnames(res1$N)[[2]] <- as.character(ages)
dimnames(res1$N)[[3]] <- as.character(yrs)

boxplot(res1$dSSB,col='magenta',outline=F,ylim=c(0,1.1))

##########################
# generate observed data #
##########################

# CPUE

obsCPUE <- list(f=1,q=1e-6,cpuecv=0.15,yrs.cpue=yrs,unit='numbers')
res.CPUE <- get.CPUE(Biol,Fleet,res1$N,obsCPUE)
boxplot(res.CPUE,col='magenta',outline=F)
CPUE <- FLQuant(dim=c(1,length(obsCPUE$yrs.cpue),1,1,1,nits))
CPUE[1,,1,1,1,] <- t(res.CPUE)
dimnames(CPUE)[[2]] <- as.character(obsCPUE$yrs.cpue)
bwplot(data~as.factor(year),CPUE)

# catch size composition

nycal <- 20 # 20 years of historical CAL
obsCAL <- list(f=1,neff=50,lbins=lbins,yrs.cal=yrs[(length(yrs)-nycal+1):length(yrs)])
res.CAL <- get.CAL(Biol,Fleet,res1$N,obsCAL)
CAL <- FLQuant(quant='length',dim=c(nbins,length(obsCAL$yrs.cal),1,1,1,nits))
dimnames(CAL)[[1]] <- as.character(mulbins)
dimnames(CAL)[[2]] <- as.character(obsCAL$yrs.cal)
CAL[,,1,1,1,] <- aperm(res.CAL,c(2,3,1))
ply <- apply(CAL,-1,function(x){x <- x/sum(x)})
bwplot(data~as.factor(length)|as.factor(year),CAL)
bwplot(data~as.factor(length)|as.factor(year),ply)

############################################################################
# expand the key variables of use from the biological OM (MSY ratios etc.) #
############################################################################

dBmsy <- res1$dSSB/MSYvars$Bmsy
dHmsy <- hy/MSYvars$Hmsy
BtoBmsy <- dSSB <- HtoHmsy <- FLQuant(dim=c(1,length(yrs),1,1,1,nits))
C <- FLQuant(dim=c(1,length(yrs),nf,1,1,nits))
dimnames(BtoBmsy)[[2]] <- dimnames(dSSB)[[2]] <- dimnames(HtoHmsy)[[2]] <- dimnames(C)[[2]] <- as.character(yrs)
BtoBmsy[1,,1,1,1,] <- t(dBmsy)
HtoHmsy[1,,1,1,1,] <- dHmsy
dSSB[1,,1,1,1,] <- t(res1$dSSB)
C[1,,,1,1,] <- aperm(res1$C,c(3,2,1))

# convert to FLStock/FLBiol as well for bettr visualisation

# save it all

save.image("filestore/ALBlike_OM.rda")

