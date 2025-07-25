library(RTMB)
# remotes::install_github("fishfollower/SAM/stockassessment")
library(stockassessment)
#fit<-fitfromweb("NEA_sei_21_v5Reca")
load("fit.RData")
dat<-list()
dat$logobs <- fit$data$logobs
dat$aux <- fit$data$aux
dat$idx1 <- fit$data$idx1
dat$idx2 <- fit$data$idx2
dat$minYear <- min(fit$data$years)
dat$minAge <- min(fit$data$minAgePerFleet)
dat$fleetTypes <- fit$data$fleetTypes
dat$sampleTimes <- fit$data$sampleTimes
dat$year <- fit$data$years
dat$age <-  min(fit$data$minAgePerFleet):max(fit$data$maxAgePerFleet)
dat$M <- fit$data$natMor
dat$SW <- fit$data$stockMeanWeight
dat$MO <- fit$data$propMat
dat$PF <- fit$data$propF
dat$PM <- fit$data$propM

dat$srmode <- 0
dat$fcormode <- 2
dat$keyF <- fit$conf$keyLogFsta[1,]
dat$keyQ <- fit$conf$keyLogFpar
dat$keySd <- fit$conf$keyVarObs
dat$keySd[dat$keySd<(-.1)] <- NA
dat$covType <- c(0,1,2) #as.integer(fit$conf$obsCorStruct)-1
dat$keyIGAR <- fit$conf$keyCorObs
dat$keyIGAR[fit$conf$keyCorObs==-1] <- NA
dat$keyIGAR[is.na(fit$conf$keyCorObs)] <- -1
dat$keyIGAR[2, 1:4]<-0
dat$noParUS <- sapply(1:length(dat$fleetTypes),
                      function(f){
                        A <- sum(!is.na(dat$keySd[f,]))
                        ifelse(dat$covType[f]==2, (A*A-A)/2, 0)
                      })
# AR type flag
dat$ar1code <- -1 # -2 (options: -2=no AR1, 0=centered,-1=non-centered)


#dat$logobs[c(7,9,13)+100]<-NA

par <- list()
par$logsdR <-0
par$logsdS <- 0
par$logsdF <- numeric(max(dat$keyF)+1)
par$rickerpar <- if(dat$srmode==1){c(1,1)}else{numeric(0)}
par$transRhoF <- if(dat$fcormode==0){numeric(0)}else{0.1}
par$bhpar <- if(dat$srmode==2){c(1,1)}else{numeric(0)}
par$logQ <- numeric(max(dat$keyQ, na.rm=TRUE)+1)
par$logsd <- numeric(max(dat$keySd, na.rm=TRUE)+1)
par$logIGARdist <- numeric(max(dat$keyIGAR, na.rm=TRUE)+1)
par$parUS <- numeric(sum(dat$noParUS))
par$logN <- matrix(0, nrow=length(dat$year), ncol=length(dat$age))
par$logF <- matrix(0, nrow=length(dat$year), ncol=max(dat$keyF)+1)
par$missing <- numeric(sum(is.na(dat$logobs)))
par$tPhi <- if(dat$ar1code==-2){numeric(0)}else{(1)} # AR phi for recruitment

##############

itrans <- function(x){
  2/(1 + exp(-2 * x)) - 1;
}

ssbFUN <- function(logN, logFF, M, SW, MO, PF, PM){
  nrow <- nrow(logN)
  ncol <- ncol(logN)
  ret <- numeric(nrow)
  for(y in 1:nrow){
    for(a in 1:ncol){
      ret[y] = ret[y]+SW[y,a]*MO[y,a]*exp(logN[y,a])*exp(-PF[y,a]*exp(logFF[y,a])-PM[y,a]*M[y,a])
    }
  }
  return(ret);
}

jnll<-function(par){
  getAll(par, dat)

  logobs <- OBS(logobs)

  nobs <- length(logobs)
  nrow <- nrow(M)
  ncol <- ncol(M)

  sdR <- exp(logsdR)
  sdS <- exp(logsdS)
  sdF <- exp(logsdF)
  sd <- exp(logsd)

  phi <- 2*plogis(tPhi)-1

  logobs <- logobs+numeric(1) ## hack to make advector
  logobs[is.na(logobs)] <- missing  ##patch missing

  logFF <- logF[,keyF+1] ## expand F


  ssb <- ssbFUN(logN,logFF,M,SW,MO,PF,PM)

  jnll <- 0

  # AR1 init conditions
  # centered + non-centered
  if(ar1code%in%c(0,-1)) {
    jnll <- jnll -dnorm(logN[1,1],0, sdR/sqrt(1-phi^2),log=TRUE)
  }

  for(y in 2:nrow){

    thisSSB <- ifelse((y-minAge-1)>(-.5),ssb[y-minAge],ssb[1])

    if(srmode==0){
      predN <- logN[y-1,1]
    }
    if(srmode==1){
      predN <- rickerpar[1]+log(thisSSB)-exp(rickerpar[2])*thisSSB
    }
    if(srmode==2){
      predN <- bhpar[1]+log(thisSSB)-log(1.0+exp(bhpar[2])*thisSSB)
    }
    # no AR1
    if(ar1code==-2){
      jnll <- jnll - dnorm(logN[y,1],predN,sdR,TRUE)
    }
    # centered
    if(ar1code==0){
      jnll <- jnll - dnorm(logN[y,1],phi*predN,sdR,TRUE)
    }
    # non-centered:
    # x[t] = p * x[t - 1] + sigma * z[t];
    # z[t]=(x[t]-p*x[t-1])/sigma
    if(ar1code==-1){
      jnll <- jnll - dnorm((logN[y,1]-phi*logN[y-1,1])/sdR,0,1,TRUE)
    }

  }
  logrec <- as.vector(logN[,1]) + numeric(1)
  # remaining N matrix
  for(y in 2:nrow){
    for(a in 2:ncol){
      predN <- logN[y-1,a-1]-exp(logFF[y-1,a-1])-M[y-1,a-1]
      if(a==ncol){
        predN <- log(exp(predN)+exp(logN[y-1,a]-exp(logFF[y-1,a])-M[y-1,a]))
      }
      jnll <- jnll - dnorm(logN[y,a],predN,sdS,TRUE)
    }
  }

  ############################################# F part
  SigmaF <- matrix(0,ncol(logF),ncol(logF))

  if(fcormode==0){
    diag(SigmaF) <- sdF*sdF
  }

  if(fcormode==1){
    diag(SigmaF) <- sdF*sdF
    rhoF <- itrans(transRhoF[1])
    for(i in 2:ncol(logF)){
      for(j in 1:(i-1)){
        SigmaF[i,j] <- rhoF*sdF[i]*sdF[j]
        SigmaF[j,i] <- SigmaF[i,j]
      }
    }
  }

  if(fcormode==2){
    diag(SigmaF) <- sdF*sdF
    rhoF <- itrans(transRhoF[1])
    for(i in 2:ncol(logF)){
      for(j in 1:(i-1)){
        SigmaF[i,j] <- sdF[i]*sdF[j]*(rhoF^(i-j))
        SigmaF[j,i] <- SigmaF[i,j]
      }
    }
  }

  for(y in 2:nrow){
    jnll <- jnll - dmvnorm(logF[y,],logF[y-1,],SigmaF,log=TRUE)
  }

  logPred <- numeric(nobs)
  for(i in 1:nobs){
    y <- aux[i,1]-minYear+1
    f <- aux[i,2]
    a <- aux[i,3]-minAge+1
    Z <- exp(logFF[y,a])+M[y,a]
    if(fleetTypes[f]==0){
      logPred[i] <- logN[y,a]-log(Z)+log(1-exp(-Z))+logFF[y,a]
    }
    if(fleetTypes[f]==2){
      logPred[i] <- logQ[keyQ[f,a]+1]+logN[y,a]-Z*sampleTimes[f]
    }
  }

  Svec <- list()
  for(f in 1:nrow(idx1)){
    thisdim <- sum(!is.na(keySd)[f,])
    S <- matrix(0,thisdim,thisdim)
    if(covType[f]==0){
      diag(S) <- sd[na.omit(keySd[f,])+1]^2
    }
    if(covType[f]==1){
      dist <- numeric(thisdim)
      d=2;
      for(a in 1:ncol(keyIGAR)){
        if(!is.na(keyIGAR[f,a])){
          dist[d] <- dist[d-1]+exp(logIGARdist[keyIGAR[f,a]+1])
          d <- d+1
        }
      }
      sdvec <- sd[na.omit(keySd[f,])+1]
      for(i in 1:nrow(S)){
        for(j in 1:(i-1)){
          S[i,j] <- sdvec[i]*sdvec[j]*(0.5^(dist[i]-dist[j]));
          S[j,i] <- S[i,j];
        }
        S[i,i] <- sdvec[i]^2;
      }
    }
    if(covType[f]==2){
      sdvec <- numeric(thisdim);
      d <- 1;
      for(a in 1:ncol(keySd)){
        if(!is.na(keySd[f,a])){
          sdvec[d] <- sd[keySd[f,a]+1]
          d <- d+1
        }
      }
      from <- 1
      ii <- 1
      while(ii<f){from=from+noParUS[ii]; ii<-ii+1}
      thispar <- parUS[from:(from+noParUS[f]-1)]
      U <- diag(thisdim)
      U[upper.tri(U)] <- thispar
      R <- cov2cor(t(U)%*%(U))
      D <- diag(sdvec)
      S <- D%*%R%*%D
    }
    Svec[[f]] <- S;
  }

  for(f in 1:nrow(idx1)){
    for(y in 1:ncol(idx1)){
      if(!is.na(idx1[f,y])){
        idx <- (idx1[f,y]+1):(idx2[f,y]+1)
        jnll <- jnll -  dmvnorm(logobs[idx],logPred[idx],Svec[[f]],log=TRUE)
      }
    }
  }
  REPORT(logPred)
  ADREPORT(ssb)
  ADREPORT(logrec)
  jnll
}

obj <- MakeADFun(jnll, par, random=c("logN", "logF", "missing"), map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=FALSE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=1000, iter.max=1000))
opt$objective

stockassessment::ssbplot(fit)
sdr<-sdreport(obj)
sdr
pl<-as.list(sdr,report=TRUE, "Est")

plsd<-as.list(sdr,report=TRUE, "Std")
lines(dat$year, pl$ssb, lwd=3, col="darkred")
lines(dat$year, pl$ssb-2*plsd$ssb, lwd=3, col="darkred", lty="dotted")
lines(dat$year, pl$ssb+2*plsd$ssb, lwd=3, col="darkred", lty="dotted")

library(ggplot2);library(dplyr)
pl$logrec
data.frame(year=dat$year,
           rec=exp(pl$logrec),
           up=exp(pl$logrec+2*plsd$logrec),
           lw=exp(pl$logrec-2*plsd$logrec)) %>%
  ggplot(aes(x=year, y=rec))+
           geom_bar(stat="identity")+
  geom_errorbar(aes(ymax=up,ymin=lw),col="black")

par <- as.list(sdr, what = "Est")
parsd <- as.list(sdr, what = "Std")
2*plogis(par$tPhi)-1



plot(dat$year, pl$ssb, lwd=3, col="darkred", ylim=c(min(pl$ssb-2*plsd$ssb-100), max(pl$ssb+2*plsd$ssb+100)))
lines(dat$year, pl$ssb-2*plsd$ssb, lwd=3, col="darkred", lty="dotted")
lines(dat$year, pl$ssb+2*plsd$ssb, lwd=3, col="darkred", lty="dotted")


odat <- dat # save a copy of the orig data set

# let's try to simulate some logobs based on phi=0.55, tphi~=1.238
2*plogis(1.238)-1


sdr<-sdreport(obj)
sdr
pl<-as.list(sdr,"Est")
pl$tPhi <- 1.238 # phi=0.55, tphi~=1.238
objcon

docondsim <- function(){
  dat <- list()  # Local dat
  sdr <- sdreport(obj)
  pl <- as.list(sdr, "Est")
  objcon <- MakeADFun(f, pl, silent=TRUE) # Simulate data conditional on MLE
  dat$y <<- objcon$simulate()$y  # need to update the dat object in the global env to make it accessible to MakeADFun
  objsim <- MakeADFun(jnll, par, random = "gamma", silent=TRUE)
  fitsim <- nlminb(objsim$par, objsim$fn, objsim$gr)
  as.list(sdreport(objsim),"Est")$gamma
}
sim <- replicate(100, docondsim())
dat <- odat  # Restore original data

