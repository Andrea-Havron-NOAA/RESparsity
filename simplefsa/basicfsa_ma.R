load("fsa.RData") # gets "dat"
library(RTMB)
par <- list(
  logN1Y=rep(0,nrow(dat$M)),
  eps=rep(0,ncol(dat$M)),
  logFY=rep(0,ncol(dat$M)),
  logFA=rep(0,nrow(dat$M)),
  logSdCatch=0,
  logQ=rep(0,length(unique(dat$age[dat$fleet==2]))),
  logSdSurvey=0,
  logSdR=0,
  tthetaR=0,
  logMuR=0
)

nll<-function(par){
  getAll(par, dat)

  na <- max(age)-min(age)+1
  ny <- max(year)-min(year)+1

  ## setup F
  F <- exp(outer(logFA,logFY,"+"))

  ## setup N
  ans <- -sum(dnorm(eps,0,sd=exp(logSdR), log=TRUE))
  logN1A <- numeric(length(eps)-1)
  theta <- 2*plogis(tthetaR)-1
  for(i in 1:length(logN1A)){
    logN1A[i] <- logMuR+eps[i+1]+theta*eps[i]
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

  ans <- ans -sum(dnorm(logObs,logPred,sdvec,TRUE))

  logssb <- log(apply(exp(logN)*stockMeanWeight*propMature,2,sum))

  ADREPORT(logssb)
  return(ans)
}

obj <- MakeADFun(nll, par, map=list(logFA=factor(c(1:4,NA,NA,NA))), silent=TRUE, random="eps")

opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(iter.max=1000,eval.max=1000))
sdrep <- sdreport(obj)
pl <- as.list(sdrep, "Est", report=TRUE)
plsd <- as.list(sdrep, "Std", report=TRUE)

yr<-sort(unique(dat$year))
plot(yr, exp(pl$logssb), type="l", lwd=5, col="red", ylim=c(0,550000), xlab="Year", ylab="SSB")
lines(yr, exp(pl$logssb-2*plsd$logssb), type="l", lwd=1, col="red")
lines(yr, exp(pl$logssb+2*plsd$logssb), type="l", lwd=1, col="red")


#################################################################
# EW added this to do the Stan sampling
remotes::install_github("kaskr/tmbstan/tmbstan")
library(tmbstan)
library(rstan)
library(broom.mixed)
library(ggplot2)
library(gridExtra)

# set up cores to sample in parallel
#cores <- parallel::detectCores()-1
#options(mc.cores = cores)
# fit the model
stan_fit <- tmbstan(obj, chains=4,
                    iter=500, # 2500 burn in
                    control=list(adapt_delta=0.99))

# diagnostics look pretty good
pars <- summary(stan_fit)$summary
which(pars[,"Rhat"] > 1.15)

# extract tidied parameters for plotting
tidied_pars <- tidy(stan_fit)

post <- rstan::extract(stan_fit)

df <- data.frame(
  logSdSurvey = post$logSdSurvey,
  logSdR      = post$logSdR,
  tphiR       = post$tthetaR,
  logMuR      = post$logMuR
)

post <- rstan::extract(stan_fit)

p1 <- ggplot(data.frame(x = post$logSdSurvey), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40) +
  geom_density() + ggtitle("logSdSurvey")

p2 <- ggplot(data.frame(x = post$logSdR), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40) +
  geom_density() + ggtitle("logSdR")

p3 <- ggplot(data.frame(x = post$tthetaR), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40) +
  geom_density() + ggtitle("tthetaR")

p4 <- ggplot(data.frame(x = post$logMuR), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40) +
  geom_density() + ggtitle("logMuR")

grid.arrange(p1, p2, p3, p4, ncol = 2)
