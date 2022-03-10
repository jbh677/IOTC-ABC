####################
# HCR utility code #
####################

# get the maturity and weight-at-age from length

get.wt.mat <- function(mula,mulbins,ml,ages,awl,bwl,lm50,numl,sdl) {

  wt <- mat <- rep(NA,nages)
  for(a in 1:nages) {
  
    ltmp <- seq(mula[a]*exp(-1.96*sdl),mula[a]*exp(1.96*sdl),length=50)
    wtmp <- awl*ltmp^bwl
    dw <- dnorm(log(wtmp),log(awl*mula[a]^bwl),sdl,log=FALSE)
    dw <- dw/sum(dw)
    wt[a] <- sum(wtmp*dw)
    mux <- lm50
    mtmp <- ltmp^numl/(mux^numl+ltmp^numl)
    #mtmp <-  1/(1+19^{-(ltmp-lm50)/(lm95-lm50)})
    dm <- dnorm(log(mtmp),log(ml[a]),sdl,log=FALSE)
    dm <- dm/sum(dm)
    mat[a] <- sum(mtmp*dm) 

  }
  
  # distribution of length given age: p(l | a)
  
  nbins <- length(mulbins)
  pla <- array(dim=c(nbins,nages))
  for(a in 1:nages) {
    
    dx <- dnorm(log(mulbins),log(mula[a]),sdl,FALSE)
    dx <- dx/sum(dx)
    pla[,a] <- dx
    
  }
 
  return(list(wt=wt,mat=mat,pla=pla))

}

# projection

projection <- function(hfix,Biol,Fleet,ny,nits,rngseed) {

  res <- initpop(Biol,Fleet)
  Ninit <- res$Ninit
  SPR <- res$SPR
  nages <- Biol$nages
  nf <- Fleet$nf
  self <- Fleet$self
  fshare <- Fleet$fshare
  initdel <- Biol$initdel
  B0 <- Biol$B0
  M <- Biol$M
  hh <- Biol$hh
  mat <- Biol$mat
  wt <- Biol$wt
  sigmar <- Biol$sigmar
  rho <- Biol$rho
  dm <- c(ny,nages,nf)
  N <- array(dim=c(nits,nages,ny))
  
  set.seed(rngseed)
  for(n in 1:nits) N[n,,] <- proj(Ninit,B0,hh,SPR,sigmar,rho,M,wt,mat,self,fshare,hfix,dm)

  # get SSB depletion

  xx <- apply(N,c(1,3),function(x,mat,wt){x <- sum(x*mat*wt)},mat,wt)
  dSSB <- xx[]/B0

  return(list('dSSB'=dSSB,'N'=N))

}

projmultiyearh <- function(hy,Biol,Fleet,ny,nits,rngseed) {

  res <- initpop(Biol,Fleet)
  Ninit <- res$Ninit
  SPR <- res$SPR
  nages <- Biol$nages
  nf <- Fleet$nf
  self <- Fleet$self
  fshare <- Fleet$fshare
  initdel <- Biol$initdel
  B0 <- Biol$B0
  M <- Biol$M
  hh <- Biol$hh
  mat <- Biol$mat
  wt <- Biol$wt
  sigmar <- Biol$sigmar
  rho <- Biol$rho
  dm <- c(ny,nages,nf)
  N <- array(dim=c(nits,nages,ny))
  SSB <- array(dim=c(nits,ny))
  C <- array(dim=c(nits,nf,ny))
  
  set.seed(rngseed)
  for(n in 1:nits) {
    
    rezz <- projhy(Ninit,B0,hh,SPR,sigmar,rho,M,wt,mat,self,fshare,hy,dm)
    N[n,,] <- rezz$N
    SSB[n,] <- rezz$SSB
    C[n,,] <- rezz$C
     
  }
  
  # get SSB depletion

  dSSB <- SSB[]/B0

  return(list('dSSB'=dSSB,'N'=N,'C'=C))

}

# projection DD

projectiondd <- function(hfix,Biol,Fleet,ny,nits,rngseed) {

  res <- initpop(Biol,Fleet)
  Ninit <- res$Ninit
  SPR <- res$SPR
  nages <- Biol$nages
  nf <- Fleet$nf
  self <- Fleet$self
  fshare <- Fleet$fshare
  initdel <- Biol$initdel
  B0 <- Biol$B0
  M <- Biol$M
  hh <- Biol$hh
  mat <- Biol$mat
  wt <- Biol$wt
  sigmar <- Biol$sigmar
  rho <- Biol$rho
  dm <- c(ny,nages,nf)
  N <- N0 <- array(dim=c(nits,nages,ny))
  
  set.seed(rngseed)
  for(n in 1:nits) {
    
    zzz <- projdd(Ninit,B0,hh,SPR,sigmar,rho,M,wt,mat,self,fshare,hfix,dm)
    N[n,,] <- zzz$N
    N0[n,,] <- zzz$N0 
     

  }
  
  # get dynamic SSB depletion

  xx <- apply(N,c(1,3),function(x,mat,wt){x <- sum(x*mat*wt)},mat,wt)
  xx0 <- apply(N0,c(1,3),function(x,mat,wt){x <- sum(x*mat*wt)},mat,wt) 
  dSSB <- xx/xx0

  return(list('dSSB'=dSSB,'N'=N,'N0'=N0))

}

# multi-year harvest rates

projddmultiyearh <- function(hy,Biol,Fleet,ny,nits,rngseed) {

  res <- initpop(Biol,Fleet)
  Ninit <- res$Ninit
  SPR <- res$SPR
  nages <- Biol$nages
  nf <- Fleet$nf
  self <- Fleet$self
  fshare <- Fleet$fshare
  initdel <- Biol$initdel
  B0 <- Biol$B0
  M <- Biol$M
  hh <- Biol$hh
  mat <- Biol$mat
  wt <- Biol$wt
  sigmar <- Biol$sigmar
  rho <- Biol$rho
  dm <- c(ny,nages,nf)
  N <- N0 <- array(dim=c(nits,nages,ny))
  
  set.seed(rngseed)
  for(n in 1:nits) {
    
    zzz <- projddhy(Ninit,B0,hh,SPR,sigmar,rho,M,wt,mat,self,fshare,hy,dm)
    N[n,,] <- zzz$N
    N0[n,,] <- zzz$N0 

  }
  
  # get dynamic SSB depletion

  xx <- apply(N,c(1,3),function(x,mat,wt){x <- sum(x*mat*wt)},mat,wt)
  xx0 <- apply(N0,c(1,3),function(x,mat,wt){x <- sum(x*mat*wt)},mat,wt) 
  dSSB <- xx/xx0

  return(list('dSSB'=dSSB,'N'=N,'N0'=N0))

}

# initial population

initpop <- function(Biol,Fleet) {

  nages <- Biol$nages
  nf <- Fleet$nf
  self <- Fleet$self
  fshare <- Fleet$fshare
  initdel <- Biol$initdel
  M <- Biol$M
  B0 <- Biol$B0
  hh <- Biol$hh
  m <- Biol$mat
  wt <- Biol$wt

  ntmp <- rep(NA,nages)

  init.fn <- function(hinit) {


    hf <- self
    for(f in 1:nf) hf[f,] <- hinit*fshare[f]*self[f,]
    h <- apply(hf,2,sum)
    ntmp[1] <- 1
    for(a in 2:nages) ntmp[a] <- ntmp[a-1] * exp(-M[a-1])
    ntmp[nages] <- ntmp[nages]/(1-exp(-M[nages]))
    SPR0 <-  sum(ntmp*wt*m)
    ntmp[1] <- 1
    for(a in 2:nages) ntmp[a] <- ntmp[a-1] * exp(-M[a-1]) * (1-h[a-1])
    ntmp[nages] <- ntmp[nages]/(1-exp(-M[nages])*(1-h[nages]))
    SPR <-  sum(ntmp*wt*m) 
    spratio <- SPR/SPR0
    dep <- (4*hh*spratio+hh-1)/(5*hh-1)
    return(dep-initdel)

  }

  xinit <- uniroot(init.fn,interval=c(0,0.4))$root 
  hf <- self
  for(f in 1:nf) hf[f,] <- xinit*fshare[f]*self[f,]
  h <- apply(hf,2,sum) 
  ntmp[1] <- 1
  for(a in 2:nages) ntmp[a] <- ntmp[a-1]*exp(-M[a-1])
  ntmp[nages] <- ntmp[nages]/(1-exp(-M[nages]))
  spr0 <- sum(ntmp*wt*m)

  # get exploited SPR

  ntmp[1] <- 1
  for(a in 2:nages) ntmp[a] <- ntmp[a-1]*exp(-M[a-1])*(1-h[a-1])
  ntmp[nages] <- ntmp[nages]/(1-exp(-M[nages])*(1-h[nages]))
  sprf <- sum(ntmp*wt*m)
 
  # SPR ratio

  spr.ratio <- sprf/spr0

  # stock-recruit parameters (R = alp S (1+ bet S))

  alp <- 4*hh/(spr0*(1-hh))
  bet <- (5*hh-1)/(B0*(1-hh))

  # get exploited eqm mean spawners & recruitment

  Sbar <- max((alp*sprf-1)/bet,0)
  Rbar <- Sbar/sprf
  Ninit <- Rbar * ntmp

  return(list(Ninit=Ninit,SPR=spr0))

}

# MSY calculator

get.MSY <- function(Biol,Fleet,msytype='deterministic',hmax=0.4) {
  
  if(msytype == 'stochastic') {
    
  }
  
  if(msytype == 'deterministic') {
   
    nages <- Biol$nages
    nf <- Fleet$nf
    self <- Fleet$self
    fshare <- Fleet$fshare
    initdel <- Biol$initdel
    M <- Biol$M
    B0 <- Biol$B0
    hh <- Biol$hh
    m <- Biol$mat
    wt <- Biol$wt
    
    ntmp <- rep(NA,nages)
    
    msy.fn <- function(x) {
      
      hf <- self
      for(f in 1:nf) hf[f,] <- x*fshare[f]*self[f,]
      h <- apply(hf,2,sum)
      ntmp[1] <- 1
      for(a in 2:nages) ntmp[a] <- ntmp[a-1] * exp(-M[a-1])
      ntmp[nages] <- ntmp[nages]/(1-exp(-M[nages]))
      SPR0 <-  sum(ntmp*wt*m)
      ntmp[1] <- 1
      for(a in 2:nages) ntmp[a] <- ntmp[a-1] * exp(-M[a-1]) * (1-h[a-1])
      ntmp[nages] <- ntmp[nages]/(1-exp(-M[nages])*(1-h[nages]))
      SPR <-  sum(ntmp*wt*m) 
      spratio <- SPR/SPR0
      alp <- 4*hh/(SPR0*(1-hh))
      bet <- (5*hh-1)/(B0*(1-hh))
      Sbar <- max((alp*SPR-1)/bet,0)
      Bratio <<- Sbar/B0
      Rbar <- Sbar/SPR
      N <- Rbar*ntmp
      C <- sum(h*wt*N) # catch in weight
      
      return(C)
      
    } 
    
    res.msy <- optimise(msy.fn,interval=c(0,hmax),maximum=TRUE)
    
    return(list(Hmsy=res.msy$maximum,Cmsy=res.msy$objective,Bmsy=Bratio))
  }
}

# generate CPUE indices

get.CPUE <- function(Biol,FLeet,N,obsCPUE) {
 
  f <- obsCPUE$f
  q <- obsCPUE$q
  sdx <- sqrt(log(1+obsCPUE$cpuecv^2))
  ysub <- as.character(obsCPUE$yrs.cpue)
  Nx <- N[,,ysub]
  nits <- dim(Nx)[1] 
  I <- matrix(nrow=nits,ncol=length(ysub))
  colnames(I) <- ysub
  selx <- Fleet$self[f,]
  wt <- Biol$wt
  
  if(obsCPUE$unit == 'numbers') Ihat <- apply(Nx,c(1,3),function(x,q,selx){x <- sum(q*selx*x)},q,selx)
  if(obsCPUE$unit == 'weight') Ihat <- apply(Nx,c(1,3),function(x,q,wt,selx){x <- sum(q*wt*selx*x)},q,wt,selx)

  Iobs <- Ihat[]*rlnorm(dim(Ihat),-sdx^2/2,sdx)
  
  return(Iobs)
  
}
    
# generate catch length composition

get.CAL <- function(Biol,FLeet,N,obsCAL) {
  
  f <- obsCAL$f
  neff <- obsCAL$neff
  lbins <- obsCAL$lbins
  nbins <- length(lbins)-1
  mulbins <- 0.5*(lbins[-1]+lbins[-(nbins+1)])
  mula <- Biol$mula
  sdla <- Biol$sdla
  nages <- Biol$nages
  ysub <- as.character(obsCAL$yrs.cal)
  Nx <- N[,,ysub]
  nits <- dim(Nx)[1] 
  CAL <- array(dim=c(nits,nbins,length(ysub)))
  dimnames(CAL)[[3]] <- ysub
  selx <- Fleet$self[f,]
  pla <- Biol$pla
  
  CAA <- aperm(apply(Nx,c(1,3),function(x,selx){x <- x*selx},selx),c(2,1,3))
  for(n in 1:nits)
    for(y in 1:length(ysub))
      for(l in 1:nbins) CAL[n,l,y] <- sum(CAA[n,,y]*pla[l,])
  pyl <- aperm(apply(CAL,c(1,3),function(x){x <- x/sum(x)}),c(2,1,3))
  
  # resample using multinomial
  
  lfdat <-  aperm(apply(pyl,c(1,3),function(x){x <- as.vector(rmultinom(1,neff,x))}),c(2,1,3))
    
  return(lfdat)
  
}

# ABC MCMC sampler - algorithm D from Wilkinson (2013)

abc.mcmc <- function(thetaold,zeta,parx,rwsdvec,Cx,LFx,Ix) {
  
  thetanew <- thetaold
  npar <- length(thetanew)
  
  # parameter location
  
  if(parx == 'B0hinit') idx <- c(1,2) 
  if(parx == 'recdevs') idx <- 3:(npar-3)
  if(parx == 'selectivity') idx <- (npar-2):npar
  
  # new proposal
  
  thetanew[idx] <- rnorm(length(idx),thetaold[idx],rwsdvec[idx])
  
  # estimate MSY variables for new proposal parameters
  
  biol <- Biol
  fleet <- Fleet
  biol$B0 <- exp(thetanew[1])
  biol$hh <- zeta[1]
  biol$M[] <- zeta[2]
  pla <- biol$pla
  nages <- biol$nages
  selx <- rep(NA,nages)
  selpars <- exp(thetanew[(npar-2):npar])
  for(a in 1:nages) selx[a] <- ifelse(a <= selpars[1],2^{-(a-selpars[1])^2/selpars[2]^2},2^{-(a-selpars[1])^2/selpars[3]^2})
  fleet$self[1,] <- selx
  MSYvars <- get.MSY(biol,fleet,'deterministic',0.5)
  Bmsy <- MSYvars$Bmsy*biol$B0
  Hmsy <- MSYvars$Hmsy
  Cmsy <- MSYvars$Cmsy
  
  # extract biological and fleet variables for population model
  
  nages <- biol$nages
  nf <- fleet$nf
  self <- fleet$self
  fshare <- fleet$fshare
  B0 <- biol$B0
  hinit <- ilogit(thetanew[2])
  M <- biol$M
  hh <- biol$hh
  mat <- biol$mat
  wt <- biol$wt
  recdevs <- thetanew[3:(npar-3)]
  sigmar <- zeta[3]
  dm <- c(length(recdevs),nages,nf,dim(pla)[1])
  fidx <- 1
  
  # get popn trajectory/predicted data given parameters
  
  rezz <- abcpdyn(B0,hh,M,hinit,fshare,mat,wt,recdevs,sigmar,pla,self,dm,Cx,fidx) 
  b0 <- exp(thetanew[1])
  dssb <- rezz$SSB/b0
  dbmsy <- rezz$SSB/Bmsy
  dhmsy <- rezz$hy/Hmsy
  names(dssb) <- names(dbmsy) <- names(dhmsy) <- as.character(yrs.abc)
  mudssb <- mean(log(dssb[as.character(yrs.bdep)]))
  mudbmsy <- mean(log(dbmsy[as.character(yrs.bmsy)]))
  mudhmsy <- mean(log(dhmsy[as.character(yrs.hmsy)]))
  
  ###################################
  # priors & discrepancy statistics #
  ###################################
  
  # priors:
  # 1. zeta: hh, M, sigmar
  # 2. theta: estimated parameters (B0,hinit,recdevs,selpars)
  # 3. status: bmsy and hmsy ratios and depletion
  
  zetaprior <- sum(dnorm(c(log((zeta[1]-0.2)/(1-zeta[1])),log(zeta[2])),c(mueta,muM),c(sdeta,sdM),T))+dgamma(x=1/zeta[3]^2,shape=gamr,rate=psir,log=T)
  mutheta <- c(mulnB0,muhinit,rep(0,npar-5),muselpars)
  sdtheta <- c(sdlnB0,sdhinit,rep(zeta[3],npar-5),sdselpars)
  thetaprior <- sum(dnorm(thetanew,mutheta,sdtheta,T))
  statusprior <- sum(dnorm(c(mudssb,mudbmsy,mudhmsy),c(mubdep,mubmsy,muhmsy),c(sdbdep,sdbmsy,sdhmsy),T))
  prior <- zetaprior+thetaprior+statusprior
  
  # CPUE data: piecewise log-linear trends + Welch's t-distro likelihood
  
  ysplit <- as.character(cpue.idx)
  nsplit <- length(ysplit)+1
  y1 <- as.character(ystart)
  y2 <- as.character(yend+1) 
  yaug <- c(y1,ysplit,y2)
  ihat <- log(rezz$I)
  iobs <- log(Ix)
  names(iobs) <- names(ihat) <- as.character(yrs.abc)
  zzz <- rep(NA,nsplit)
  for(k in 1:nsplit) {
    
    ysub <- as.character(as.integer(yaug[k]):(as.integer(yaug[k+1])-1))
    yobs <- unname(iobs[ysub])
    yhat <- unname(ihat[ysub])
    tt <- 1:length(ysub)
    ntt <- length(tt)
    lamobs <- sum((tt-mean(tt))*(yobs-mean(yobs)))/sum((tt-mean(tt))^2)
    sdlamobs <- sqrt(var(yobs-lamobs*tt)/sum((tt-mean(tt))^2))
    lamhat <- sum((tt-mean(tt))*(yhat-mean(yhat)))/sum((tt-mean(tt))^2)
    sdlamhat <- sqrt(var(yhat-lamhat*tt)/sum((tt-mean(tt))^2))
    
    # Welch's t-test likelihood - claaaaaassssic
    
    tstat <- sqrt(ntt)*(lamobs-lamhat)/(sqrt(sdlamobs^2+sdlamhat^2))
    dft <- ntt^2*(ntt-1)*((sdlamobs^2+sdlamhat^2)/ntt)^2/(sdlamobs^4+sdlamhat^4)
    zzz[k] <- dt(tstat,dft,log=TRUE)
    
  }
  
  disc.cpue <- sum(zzz)
  
  # LF data
  
  lfhat <- rezz$LF
  lfobs <- LFx
  disc.lf <- 0
  for(k in 1:dim(lfhat)[2]) {
    
    P <- lfobs[,k]
    Pz <- P
    Pz[Pz==0] <- 1e-6
    Q <- lfhat[,k]
    disc.lf <- disc.lf-sum(P*log(Pz/Q))
      
  }
  
  disc <- disc.cpue+disc.lf
  
  # penalty on overall harvest rate: max of 0.9
  
  hvx <- 0.9-rezz$hy
  hvx[hvx > 0] <- 0
  pen <- -1e+5*sum(hvx^2)
  
  # overall discrepancy function, priors, & penalties
  
  piret <- disc+prior+pen
  
  return(list(par=thetanew,objf=piret,prior=prior,disc=disc,pen=pen))
  
}

# ABC MCMCv2 sampler - algorithm D from Wilkinson (2013)

abc.mcmc2 <- function(thetaold,zeta,parx,rwsdvec,Cx,LFx,lam.lf=1,Ix,sd.cpue) {
  
  thetanew <- thetaold
  npar <- length(thetanew)
  
  # parameter location
  
  if(parx == 'B0hinit') idx <- c(1,2) 
  if(parx == 'recdevs') idx <- 3:(npar-3)
  if(parx == 'selectivity') idx <- (npar-2):npar
  
  # new proposal
  
  thetanew[idx] <- rnorm(length(idx),thetaold[idx],rwsdvec[idx])
  
  # estimate MSY variables for new proposal parameters
  
  biol <- Biol
  fleet <- Fleet
  biol$B0 <- exp(thetanew[1])
  biol$hh <- zeta[1]
  biol$M[] <- zeta[2]
  pla <- biol$pla
  nages <- biol$nages
  selx <- rep(NA,nages)
  selpars <- exp(thetanew[(npar-2):npar])
  for(a in 1:nages) selx[a] <- ifelse(a <= selpars[1],2^{-(a-selpars[1])^2/selpars[2]^2},2^{-(a-selpars[1])^2/selpars[3]^2})
  fleet$self[1,] <- selx
  MSYvars <- get.MSY(biol,fleet,'deterministic',0.8)
  Bmsy <- MSYvars$Bmsy*biol$B0
  Hmsy <- MSYvars$Hmsy
  Cmsy <- MSYvars$Cmsy
  
  # extract biological and fleet variables for population model
  
  nages <- biol$nages
  nf <- fleet$nf
  self <- fleet$self
  fshare <- fleet$fshare
  B0 <- biol$B0
  hinit <- ilogit(thetanew[2])
  M <- biol$M
  hh <- biol$hh
  mat <- biol$mat
  wt <- biol$wt
  recdevs <- thetanew[3:(npar-3)]
  sigmar <- zeta[3]
  dm <- c(length(recdevs),nages,nf,dim(pla)[1])
  fidx <- 1
  
  # get popn trajectory/predicted data given parameters
  
  rezz <- abcpdyn(B0,hh,M,hinit,fshare,mat,wt,recdevs,sigmar,pla,self,dm,Cx,fidx) 
  dssb <- rezz$SSB/B0
  dbmsy <- rezz$SSB/Bmsy
  dhmsy <- rezz$hy/Hmsy
  names(dssb) <- names(dbmsy) <- names(dhmsy) <- as.character(yrs.abc)
  mudssb <- mean(log(dssb[as.character(yrs.bdep)]))
  mudbmsy <- mean(log(dbmsy[as.character(yrs.bmsy)]))
  mudhmsy <- mean(log(dhmsy[as.character(yrs.hmsy)]))
  
  ###################################
  # priors & discrepancy statistics #
  ###################################
  
  # priors:
  # 1. zeta: hh, M, sigmar
  # 2. theta: estimated parameters (B0,hinit,recdevs,selpars)
  # 3. status: bmsy and hmsy ratios and depletion
  
  zetaprior <- sum(dnorm(c(log((zeta[1]-0.2)/(1-zeta[1])),log(zeta[2])),c(mueta,muM),c(sdeta,sdM),T))+dgamma(x=1/zeta[3]^2,shape=gamr,rate=psir,log=T)
  mutheta <- c(mulnB0,muhinit,rep(0,npar-5),muselpars)
  sdtheta <- c(sdlnB0,sdhinit,rep(zeta[3],npar-5),sdselpars)
  thetaprior <- sum(dnorm(thetanew,mutheta,sdtheta,T))
  statusprior <- sum(dnorm(c(mudssb,mudbmsy,mudhmsy),c(mubdep,mubmsy,muhmsy),c(sdbdep,sdbmsy,sdhmsy),T))
  prior <- zetaprior+thetaprior+statusprior
  
  # CPUE data: classic log-normal likelihood
  
  iobs <- log(Ix)-mean(log(Ix))
  ihat <- log(rezz$I)-mean(log(rezz$I))
  disc.cpue <- sum(dnorm(iobs,ihat,sd.cpue,T))
  
  # LF data
  
  lfhat <- rezz$LF
  lfobs <- LFx
  disc.lf <- 0
  for(k in 1:dim(lfhat)[2]) {
    
    P <- lfobs[,k]
    Pz <- P
    Pz[Pz==0] <- 1e-6
    Q <- lfhat[,k]
    disc.lf <- disc.lf-sum(P*log(Pz/Q))
    
  }
  
  disc <- disc.cpue+lam.lf*disc.lf
  
  # penalty on overall harvest rate: max of 0.8
  
  hvx <- 0.8-rezz$hy
  hvx[hvx > 0] <- 0
  pen <- -1e+6*sum(hvx^2)
  
  # overall discrepancy function, priors, & penalties
  
  piret <- disc+prior+pen
  
  return(list(par=thetanew,objf=piret,prior=prior,disc=disc,pen=pen))
  
}

# ABC MCMCv3 sampler - algorithm D from Wilkinson (2013)

abc.mcmc3 <- function(thetaold,parx,rwsdvec,Cx,LFx,lam.lf=1,Ix,sd.cpue) {
  
  thetanew <- thetaold
  npar <- length(thetanew)
  
  # parameter location
  
  if(parx == 'B0hinit') idx <- c(1,2) 
  if(parx == 'recdevs') idx <- 3:(npar-6)
  if(parx == 'selectivity') idx <- (npar-5):(npar-3)
  if(parx == 'steepM') idx <- (npar-2):(npar-1)
  
  # new proposal
  
  thetanew[idx] <- rnorm(length(idx),thetaold[idx],rwsdvec[idx])
  
  # estimate MSY variables for new proposal parameters
  
  biol <- Biol
  fleet <- Fleet
  biol$B0 <- exp(thetanew[1])
  biol$hh <- (exp(thetanew[npar-2])+0.2)/(1+exp(thetanew[npar-2]))
  biol$M[] <- exp(theta[npar-1])
  pla <- biol$pla
  nages <- biol$nages
  selx <- rep(NA,nages)
  selpars <- exp(thetanew[(npar-5):(npar-3)])
  for(a in 1:nages) selx[a] <- ifelse(a <= selpars[1],2^{-(a-selpars[1])^2/selpars[2]^2},2^{-(a-selpars[1])^2/selpars[3]^2})
  fleet$self[1,] <- selx
  MSYvars <- get.MSY(biol,fleet,'deterministic',0.8)
  Bmsy <- MSYvars$Bmsy*biol$B0
  Hmsy <- MSYvars$Hmsy
  Cmsy <- MSYvars$Cmsy
  
  # extract biological and fleet variables for population model
  
  nages <- biol$nages
  nf <- fleet$nf
  self <- fleet$self
  fshare <- fleet$fshare
  B0 <- biol$B0
  hinit <- ilogit(thetanew[2])
  M <- biol$M
  hh <- biol$hh
  mat <- biol$mat
  wt <- biol$wt
  recdevs <- thetanew[3:(npar-6)]
  sigmar <- thetanew[npar]
  dm <- c(length(recdevs),nages,nf,dim(pla)[1])
  fidx <- 1
  
  # get popn trajectory/predicted data given parameters
  
  rezz <- abcpdyn(B0,hh,M,hinit,fshare,mat,wt,recdevs,sigmar,pla,self,dm,Cx,fidx) 
  dssb <- rezz$SSB/B0
  dbmsy <- rezz$SSB/Bmsy
  dhmsy <- rezz$hy/Hmsy
  names(dssb) <- names(dbmsy) <- names(dhmsy) <- as.character(yrs.abc)
  mudssb <- mean(log(dssb[as.character(yrs.bdep)]))
  mudbmsy <- mean(log(dbmsy[as.character(yrs.bmsy)]))
  mudhmsy <- mean(log(dhmsy[as.character(yrs.hmsy)]))
  
  ###################################
  # priors & discrepancy statistics #
  ###################################
  
  # priors:
  # 1. theta: estimated parameters (B0,hinit,recdevs,selpars,hh,M)
  # 2. status: bmsy and hmsy ratios and depletion
  
  mutheta <- c(mulnB0,muhinit,rep(0,npar-8),muselpars,mueta,muM)
  sdtheta <- c(sdlnB0,sdhinit,rep(sigmar,npar-8),sdselpars,sdeta,sdM)
  thetaprior <- sum(dnorm(thetanew[-npar],mutheta,sdtheta,T))
  statusprior <- sum(dnorm(c(mudssb,mudbmsy,mudhmsy),c(mubdep,mubmsy,muhmsy),c(sdbdep,sdbmsy,sdhmsy),T))
  prior <- thetaprior+statusprior
  
  # CPUE data: classic log-normal likelihood
  
  iobs <- log(Ix)-mean(log(Ix))
  ihat <- log(rezz$I)-mean(log(rezz$I))
  disc.cpue <- sum(dnorm(iobs,ihat,sd.cpue,T))
  
  # LF data
  
  lfhat <- rezz$LF
  lfobs <- LFx
  disc.lf <- 0
  for(k in 1:dim(lfhat)[2]) {
    
    P <- lfobs[,k]
    Pz <- P
    Pz[Pz==0] <- 1e-6
    Q <- lfhat[,k]
    disc.lf <- disc.lf-sum(P*log(Pz/Q))
    
  }
  
  disc <- disc.cpue+lam.lf*disc.lf
  
  # penalty on overall harvest rate: max of 0.8
  
  hvx <- 0.8-rezz$hy
  hvx[hvx > 0] <- 0
  pen <- -1e+6*sum(hvx^2)
  
  # overall discrepancy function, priors, & penalties
  
  piret <- disc+prior+pen
  
  return(list(par=thetanew,objf=piret,prior=prior,disc=disc,pen=pen))
  
}

# reconstruct the popn dynamics

reconstruct <- function(Biol,Fleet,theta.mcmc,zeta.mcmc,Cx,yrs.abc) {

  # objectos
  
  biol <- Biol
  fleet <- Fleet
  pla <- biol$pla

  # dims  

  nits <- dim(theta.mcmc)[1]
  npar <- dim(theta.mcmc)[2]
  ny <- dim(Cx)[2]
  nages <- Biol$nages
  nbins <- dim(pla)[1]
  
  # return objs
  
  yy <- as.character(yrs.abc)
  SSB <- I <- H <- FLQuant(dim=c(1,ny,1,1,1,nits))
  N <- FLQuant(quant='age',dim=c(nages,ny,1,1,1,nits))
  LF <- FLQuant(quant='length',dim=c(nbins,ny,1,1,1,nits))
  Bmsy <- Cmsy <- Hmsy <- rep(NA,nits)
  
  dimnames(SSB)[[2]] <- dimnames(I)[[2]] <- dimnames(H)[[2]] <- dimnames(N)[[2]] <- dimnames(LF)[[2]] <- yy
  
  for(n in 1:nits) {
  
    biol$hh <- zeta.mcmc[n,1]
    biol$M[] <- zeta.mcmc[n,2]
    biol$B0 <- exp(theta.mcmc[n,1])
    selx <- rep(NA,nages)
    selpars <- exp(theta.mcmc[n,(npar-2):npar])
    for(a in 1:nages) selx[a] <- ifelse(a <= selpars[1],2^{-(a-selpars[1])^2/selpars[2]^2},2^{-(a-selpars[1])^2/selpars[3]^2})
    fleet$self[1,] <- selx
    MSYvars <- get.MSY(biol,fleet,'deterministic',0.8)
    Bmsy[n] <- MSYvars$Bmsy*biol$B0
    Hmsy[n] <- MSYvars$Hmsy
    Cmsy[n] <- MSYvars$Cmsy
    
    # extract biological and fleet variables for population model
    
    nages <- biol$nages
    nf <- fleet$nf
    self <- fleet$self
    fshare <- fleet$fshare
    B0 <- biol$B0
    hinit <- ilogit(theta.mcmc[n,2])
    M <- biol$M
    hh <- biol$hh
    mat <- biol$mat
    wt <- biol$wt
    recdevs <- theta.mcmc[n,3:(npar-3)]
    sigmar <- zeta.mcmc[n,3]
    dm <- c(length(recdevs),nages,nf,dim(pla)[1])
    fidx <- 1
    
    # get popn trajectory/predicted data given parameters
    
    rezz <- abcpdyn(B0,hh,M,hinit,fshare,mat,wt,recdevs,sigmar,pla,self,dm,Cx,fidx) 
    SSB[1,,1,1,1,n] <- rezz$SSB
    H[1,,1,1,1,n] <- rezz$hy
    I[1,,1,1,1,n] <- rezz$I
    LF[,,1,1,1,n] <- rezz$LF
    N[,,1,1,1,n] <- rezz$N
    
  }  
  
  return(list(N=N,SSB=SSB,H=H,I=I,LF=LF,Bmsy=Bmsy,Hmsy=Hmsy,Cmsy=Cmsy))
  
}

# reconstruction v2

reconstructv2 <- function(Biol,Fleet,theta.mcmc,Cx,yrs.abc) {
  
  # objectos
  
  biol <- Biol
  fleet <- Fleet
  pla <- biol$pla
  
  # dims  
  
  nits <- dim(theta.mcmc)[1]
  npar <- dim(theta.mcmc)[2]
  ny <- dim(Cx)[2]
  nages <- Biol$nages
  nbins <- dim(pla)[1]
  
  # return objs
  
  yy <- as.character(yrs.abc)
  SSB <- I <- H <- FLQuant(dim=c(1,ny,1,1,1,nits))
  N <- FLQuant(quant='age',dim=c(nages,ny,1,1,1,nits))
  LF <- FLQuant(quant='length',dim=c(nbins,ny,1,1,1,nits))
  Bmsy <- Cmsy <- Hmsy <- rep(NA,nits)
  
  dimnames(SSB)[[2]] <- dimnames(I)[[2]] <- dimnames(H)[[2]] <- dimnames(N)[[2]] <- dimnames(LF)[[2]] <- yy
  
  for(n in 1:nits) {
    
    biol$hh <- (exp(theta.mcmc[n,26])+0.2)/(1+exp(theta.mcmc[n,26]))
    biol$M[] <- exp(theta.mcmc[n,27])
    biol$B0 <- exp(theta.mcmc[n,1])
    selx <- rep(NA,nages)
    selpars <- exp(theta.mcmc[n,(npar-5):(npar-3)])
    for(a in 1:nages) selx[a] <- ifelse(a <= selpars[1],2^{-(a-selpars[1])^2/selpars[2]^2},2^{-(a-selpars[1])^2/selpars[3]^2})
    fleet$self[1,] <- selx
    MSYvars <- get.MSY(biol,fleet,'deterministic',0.8)
    Bmsy[n] <- MSYvars$Bmsy*biol$B0
    Hmsy[n] <- MSYvars$Hmsy
    Cmsy[n] <- MSYvars$Cmsy
    
    # extract biological and fleet variables for population model
    
    nages <- biol$nages
    nf <- fleet$nf
    self <- fleet$self
    fshare <- fleet$fshare
    B0 <- biol$B0
    hinit <- ilogit(theta.mcmc[n,2])
    M <- biol$M
    hh <- biol$hh
    mat <- biol$mat
    wt <- biol$wt
    recdevs <- theta.mcmc[n,3:(npar-6)]
    sigmar <- theta.mcmc[n,npar]
    dm <- c(length(recdevs),nages,nf,dim(pla)[1])
    fidx <- 1
    
    # get popn trajectory/predicted data given parameters
    
    rezz <- abcpdyn(B0,hh,M,hinit,fshare,mat,wt,recdevs,sigmar,pla,self,dm,Cx,fidx) 
    SSB[1,,1,1,1,n] <- rezz$SSB
    H[1,,1,1,1,n] <- rezz$hy
    I[1,,1,1,1,n] <- rezz$I
    LF[,,1,1,1,n] <- rezz$LF
    N[,,1,1,1,n] <- rezz$N
    
  }  
  
  return(list(N=N,SSB=SSB,H=H,I=I,LF=LF,Bmsy=Bmsy,Hmsy=Hmsy,Cmsy=Cmsy))
  
}
