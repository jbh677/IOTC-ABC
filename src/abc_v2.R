################################################################
# ABC example sampler for IOTC v2 ##############################
################################################################
# R. Hillary, I. Mosqueira 2022 ################################
################################################################

library(FLCore)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
source("abc_utils.R")
logit <- function(x){return(log(x/(1-x)))}
ilogit <- function(x){return(1/(1+exp(-x)))}

# compile required C++ code

sourceCpp("abc_pdyn.cpp") # calculate predicted popn and data given parameters

# load up the OM objects

load("filestore/ALBlike_OM.rda")

######################################################
# define the variables required for ABC MCMC sampler #
######################################################

# priors: estimated parameters
# 1. mean BtoBmsy ratio (logN mu and sd) and year range
# 2. mean HtoHmsy ratio (logN mu and sd) and year range
# 3. mean BtoB0 ratio (logN mu and sd) and year range
# 4. B0 (logN mu and sd)
# 5. initial harvest rate (logistic mu and sd)
# 6. selectivity parameters (logN mu and sd)

mubmsy <- log(1)
sdbmsy <- 0.1
muhmsy <- log(1)
sdhmsy <- 0.1
mubdep <- log(0.3)
sdbdep <- 0.5
mulnB0 <- log(1e+6)
sdlnB0 <- 0.5
muhinit <- log(0.2/0.8)
sdhinit <- 0.5
muselpars <- log(c(10,7.5,25))
sdselpars <- c(0.5,0.5,0.5)

# priors: enforced priors (i.e. not freely estimated)
# 1. M (logN mu and sd)
# 2. steepness (modified logistic)
# 3. recruitment variance (inverse gamma)

muM <- log(0.4)
sdM <- 0.05
muh <- Biol$hh
cvh <- 0.05 # gives approximate 95% prior CI of 0.65-0.85
sdh <- cvh*muh
mueta <- log((Biol$hh-0.2)/(1-Biol$hh)) # inverse is (exp(mulh)+0.2)/(1+exp(mulh))
sdeta <- sdh*(4/((5*muh-1)*(1-muh))) # delta method 
musigmar <- 0.5^2
mux <- 1/musigmar
sigmax <-  0.05 
vx <- sigmax^2*mux^4
gamr <- mux^2/vx
psir <- mux/vx

# time ranges required

yrs.abc <- yrs[(length(yrs)-20+1):(length(yrs))] # last 20 years of data
yrs.bmsy <- yrs.abc[(length(yrs.abc)-3+1):length(yrs.abc)] # last 3 years
yrs.hmsy <- yrs.abc[(length(yrs.abc)-3+1):length(yrs.abc)] # last 3 years
yrs.bdep <- yrs.abc[(length(yrs.abc)-3+1):length(yrs.abc)] # last 3 years

##################
# MCMC variables #
##################

# Pick the iteration

nyabc <- length(yrs.abc)
ystart <- yrs.abc[1]
yend <- yrs.abc[length(yrs.abc)]
itest <- 52
Cx <- t(as.matrix(as.vector(iter(window(C,start=ystart,end=yend),itest)),nrow=1))
Ix <- as.vector(iter(window(CPUE,start=ystart,end=yend),itest))
LFx <- as.vector(iter(window(ply,start=ystart,end=yend),itest))
LFx <- matrix(LFx,ncol=length(yrs.abc),nrow=nbins,byrow=FALSE)

# get approximate obs erro for CPUE using loess

idf <- data.frame(t=yrs.abc,y=log(Ix))
ires <- loess(y~t,idf)
plot(yrs.abc,ires$fitted,type='l')
points(yrs.abc,log(Ix))
sd.cpue <- sd(residuals(ires))

# define weighting for LF data discrepancy

lam.lf <- 2

# total number of parameters (rec. vars. + B0 + hinit + sel. pars. + hh + M + sigmar)

npar <- length(yrs.abc)+8

# initial estimated parameter vector

theta <- c(log(1e+6),logit(hy[30]),rep(0,nyabc),log(c(smax,sl,sr)),mueta,muM,sqrt(psir/gamr))

# initial discrepancy value

rex <- abc.mcmc3(theta,'B0hinit',rep(0,npar),Cx,LFx,lam.lf,Ix,sd.cpue)
piold <- rex$objf

# iterations etc.

nits <- 1000
burn <- 100
thin <- 500
theta.mcmc <- matrix(nrow=nits,ncol=npar)

# RW variance

#ewsd <- rep(0.075,npar)
rwsd <- rep(NA,npar)
rwsd[1:2] <- 0.04
rwsd[3:(npar-6)] <- 0.06
rwsd[(npar-5):(npar-3)] <- 0.05
rwsd[(npar-2):(npar-1)] <- 0.03

acp.B0hinit <- acp.recdevs <- acp.selectivity <- acp.steepM <- 0 # acceptance rates

system.time(for(n in 1:(burn+thin*nits)) {
  
  ###############################
  # approximate MCMC algorithm: #
  # 1. B0 and hinit #############
  # 2. recdevs ##################
  # 3. selectivity ##############
  # 4. steepness & M ############
  ###############################
  
  # B0 and hinit
  
  rexold <- rex
  parx <- 'B0hinit'
  rex <- abc.mcmc3(theta,parx,rwsd,Cx,LFx,lam.lf,Ix,sd.cpue)
  uvar <- log(runif(1,0,1))
  accpt <- min(rex$objf-piold,0)
  if(accpt > uvar) {
    
    theta <- rex$par
    piold <- rex$objf
    rexold <- rex
    if(n > burn) acp.B0hinit <- acp.B0hinit+1
    
  }
  
  # recdevs

  parx <- 'recdevs'
  rex <- abc.mcmc3(theta,parx,rwsd,Cx,LFx,lam.lf,Ix,sd.cpue)
  uvar <- log(runif(1,0,1))
  accpt <- min(rex$objf-piold,0)
  if(accpt > uvar) {
    
    theta <- rex$par
    piold <- rex$objf
    rexold <- rex
    if(n > burn) acp.recdevs <- acp.recdevs+1
    
  }
  
  # selectivity
  
  parx <- 'selectivity'
  rex <- abc.mcmc3(theta,parx,rwsd,Cx,LFx,lam.lf,Ix,sd.cpue)
  uvar <- log(runif(1,0,1))
  accpt <- min(rex$objf-piold,0)
  if(accpt > uvar) {
    
    theta <- rex$par
    piold <- rex$objf
    if(n > burn) acp.selectivity <- acp.selectivity+1
    
  }
  
  # steepness & M
  
  #parx <- 'steepM'
  #rex <- abc.mcmc3(theta,parx,rwsd,Cx,LFx,lam.lf,Ix,sd.cpue)
  #uvar <- log(runif(1,0,1))
  #accpt <- min(rex$objf-piold,0)
  #if(accpt > uvar) {
    
    #theta <- rex$par
    #piold <- rex$objf
    #if(n > burn) acp.steepM <- acp.steepM+1
    
  #}
  
  # outputs
  
  if(n > burn & (n-burn) %% thin == 0) theta.mcmc[(n-burn)/thin,] <- theta
  
  # progress
  
  if(n %% 1000 == 0) {
    
    cat("\r",n,"of",burn+thin*nits)  
    flush.console()
    
  }
  
})
cat("\n")
cat("B0hinit acceptance rate: ",round(acp.B0hinit/(thin*nits),4),"\n")
cat("recdevs acceptance rate: ",round(acp.recdevs/(thin*nits),4),"\n")
cat("selectivity acceptance rate: ",round(acp.selectivity/(thin*nits),4),"\n")
cat("steepM acceptance rate: ",round(acp.steepM/(thin*nits),4),"\n")

############################
# reconstruct the dynamics #
############################

recon <- reconstructv2(Biol,Fleet,theta.mcmc,Cx,yrs.abc)

# population dynamic summary relative to true values

# SSB

ssbq <- apply(t(matrix(recon$SSB[1,,1,1,1,],nrow=20,ncol=nits)),2,quantile,c(0.025,0.5,0.975))
strue <- res1$dSSB[itest,(length(yrs)-19):(length(yrs))]*Biol$B0
smax <- max(c(max(strue),max(ssbq[3,])))
plot(yrs.abc,strue,xlab='year',ylab='SSB',type='l',lwd=1.5,col='magenta',ylim=c(0,smax))
lines(yrs.abc,ssbq[2,],lty=1,col='blue',lwd=1.5)
lines(yrs.abc,ssbq[1,],lty=2,col='blue',lwd=1.5)
lines(yrs.abc,ssbq[3,],lty=2,col='blue',lwd=1.5)
legend('bottomleft',lty=c(1,1),col=c("magenta","blue"),legend=c("true","estimated"),bty='n')

# depletion

ssb <- t(matrix(recon$SSB[1,,1,1,1,],nrow=20,ncol=nits))
bb0 <- exp(theta.mcmc[,1])
dssb <- apply(ssb,2,function(x,bb0){x <- x/bb0},bb0)
dssbq <- apply(dssb,2,quantile,c(0.025,0.5,0.975))
dtrue <- res1$dSSB[itest,(length(yrs)-19):(length(yrs))]
ymax <- max(c(max(dtrue),max(dssbq[3,])))
plot(yrs.abc,dtrue,xlab='year',ylab='SSB depletion',type='l',lwd=1.5,col='magenta',ylim=c(0,ymax))
lines(yrs.abc,dssbq[2,],lty=1,col='blue',lwd=1.5)
lines(yrs.abc,dssbq[1,],lty=2,col='blue',lwd=1.5)
lines(yrs.abc,dssbq[3,],lty=2,col='blue',lwd=1.5)
legend('bottomleft',lty=c(1,1),col=c("magenta","blue"),legend=c("true","estimated"),bty='n')

# recruitment 

rec <- t(matrix(recon$N[1,,1,1,1,],nrow=20,ncol=nits))
drecq <- apply(rec,2,quantile,c(0.025,0.5,0.975))
rtrue <- res1$N[itest,1,(length(yrs)-19):(length(yrs))]
ymax <- max(c(max(rtrue),max(drecq[3,])))
plot(yrs.abc,rtrue,xlab='year',ylab='Recruitment',type='l',lwd=1.5,col='magenta',ylim=c(0,ymax))
lines(yrs.abc,drecq[2,],lty=1,col='blue',lwd=1.5)
lines(yrs.abc,drecq[1,],lty=2,col='blue',lwd=1.5)
lines(yrs.abc,drecq[3,],lty=2,col='blue',lwd=1.5)
legend('bottomleft',lty=c(1,1),col=c("magenta","blue"),legend=c("true","estimated"),bty='n')

# harvest rate MSY ratio

hmsytrue <- get.MSY(Biol,Fleet,'deterministic',0.8)$Hmsy
hhy <- t(matrix(recon$H[1,,1,1,1,],nrow=20,ncol=nits))
Hmsy <- recon$Hmsy
hhmsy <- apply(hhy,2,function(x,Hmsy){x <- x/Hmsy},Hmsy)
dhyq <- apply(hhmsy,2,quantile,c(0.025,0.5,0.975))
htrue <- hy[(length(yrs)-19):(length(yrs))]/hmsytrue
ymax <- max(c(max(htrue),max(dhyq[3,])))
plot(yrs.abc,htrue,xlab='year',ylab='H-to-Hmsy',type='l',lwd=1.5,col='magenta',ylim=c(0,ymax))
lines(yrs.abc,dhyq[2,],lty=1,col='blue',lwd=1.5)
lines(yrs.abc,dhyq[1,],lty=2,col='blue',lwd=1.5)
lines(yrs.abc,dhyq[3,],lty=2,col='blue',lwd=1.5)
legend('bottomleft',lty=c(1,1),col=c("magenta","blue"),legend=c("true","estimated"),bty='n')

####################
# fits to the data #
####################

# LF data

lfhat <- array(recon$LF[,,1,1,1,],dim=c(12,20,nits))
yy <- as.character(yrs.abc)
df <- expand.grid(obs=NA,med=NA,lq=NA,uq=NA,y=yy,l=mulbins)
for(y in 1:length(yy)) {
  
  lx <- lfhat[,y,]
  qlx <- apply(lx,1,quantile,c(0.025,0.5,0.975))
  lobs <- LFx[,y]
  df[df$y == yy[y],'obs'] <- lobs
  df[df$y == yy[y],'med'] <- qlx[2,]
  df[df$y == yy[y],'lq'] <- qlx[1,]
  df[df$y == yy[y],'uq'] <- qlx[3,]
  
}

pp <-  ggplot(df,aes(x=l,y=obs))+geom_point(colour='magenta',size=1)+facet_wrap(~y)+geom_line(aes(x=l,y=med),colour='blue')+geom_line(aes(x=l,y=uq),linetype='dashed',colour='blue')+geom_line(aes(x=l,y=lq),linetype='dashed',colour='blue')+xlab("length")+ylab("proportion")

# CPUE

df <- expand.grid(obs=NA,med=NA,lq=NA,uq=NA,y=yy)
iobs <- Ix/mean(Ix)
ihat <- t(matrix(recon$I[1,,1,1,1,],nrow=20,ncol=nits))
ihat <- apply(ihat,1,function(x){x <- x/mean(x)})
iq <- apply(ihat,1,quantile,c(0.025,0.5,0.975))

ymax <- max(c(max(iq[3,]),max(iobs)))
plot(yrs.abc,iobs,xlab='year',ylab='Index',type='p',pch=19,col='magenta',ylim=c(0,ymax))
lines(yrs.abc,iq[2,],col='blue',lwd=1.5)
lines(yrs.abc,iq[3,],col='blue',lwd=1.5,lty=2)
lines(yrs.abc,iq[1,],col='blue',lwd=1.5,lty=2)

# save it

save.image("filestore/ABC_example.rda")
