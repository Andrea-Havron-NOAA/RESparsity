load("simplefsa/fsa.RData") # gets "dat"
library(RTMB)
par_manual <- list(
  logN1Y=rep(0,nrow(dat$M)),
  x=rep(0,ncol(dat$M)-1),
  logFY=rep(0,ncol(dat$M)),
  logFA=rep(0,nrow(dat$M)),
  logSdCatch=0, 
  logQ=rep(0,length(unique(dat$age[dat$fleet==2]))),
  logSdSurvey=0,
  logSdR=0,
  tphiR=0,
  logMuR=0
)

nll_manual<-function(par){
  getAll(par, dat)
    
  na <- max(age)-min(age)+1
  ny <- max(year)-min(year)+1

  ## setup F
  F <- exp(outer(logFA,logFY,"+"))

  ## setup N

 
  ans<-0
  logN1A<-numeric(length(x))
  sdAR<-exp(logSdR)
  phiAR<-2*plogis(tphiR)-1
  ans <- ans -sum(dnorm(x,0,1,log=TRUE))
  ans <- ans-dnorm(logN1Y[1],logMuR, sdAR, log=TRUE)
  logN1A[1]<- phiAR*logN1Y[1]+logMuR*(1-phiAR)+sdAR*sqrt(1-phiAR^2)*x[1]
  for(i in 2:length(logN1A)){    
    logN1A[i] <- phiAR*logN1A[i-1]+logMuR*(1-phiAR)+sdAR*sqrt(1-phiAR^2)*x[i] 
  }
  
  logN <- matrix(0, nrow=na, ncol=ny)
  logN[,1] <- logN1Y  
  for(y in 2:ny){
    logN[1,y] <- logN1A[y-1]
    for(a in 2:na){
      logN[a,y] <- logN[a-1,y-1]-F[a-1,y-1]-M[a-1,y-1]
    }
  } 

  # Match to observations
  logObs <- log(obs)
  logPred <- numeric(length(logObs))
  sdvec <- numeric(length(logObs))
  for(i in 1:length(logObs)){
    a <- age[i]-min(age)+1
    y <- year[i]-min(year)+1
    if(fleet[i]==1){
      logPred[i] <- log(F[a,y])-log(F[a,y]+M[a,y])+log(1-exp(-F[a,y]-M[a,y]))+logN[a,y]
      sdvec[i] <- exp(logSdCatch)
    }else{
      logPred[i] <- logQ[a]-(F[a,y]+M[a,y])*surveyTime+logN[a,y]
      sdvec[i] <- exp(logSdSurvey)
    }    
  }
  
  ans <- ans-sum(dnorm(logObs,logPred,sdvec,TRUE))

  logssb <- log(apply(exp(logN)*stockMeanWeight*propMature,2,sum))
  ADREPORT(logssb)
  return(ans)
}

obj <- MakeADFun(nll_manual, par_manual, map=list(logFA=factor(c(1:4,NA,NA,NA))), silent=TRUE, random="x")

opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(iter.max=1000,eval.max=1000))
sdrep <- sdreport(obj)
pl <- as.list(sdrep, "Est", report=TRUE)
plsd <- as.list(sdrep, "Std", report=TRUE)


#yr<-sort(unique(dat$year))
#plot(yr, exp(pl$logssb), type="l", lwd=5, col="red", ylim=c(0,550000), xlab="Year", ylab="SSB")
#lines(yr, exp(pl$logssb-2*plsd$logssb), type="l", lwd=1, col="red")
#lines(yr, exp(pl$logssb+2*plsd$logssb), type="l", lwd=1, col="red")

Matrix:::image(obj$env$spHess()) 

