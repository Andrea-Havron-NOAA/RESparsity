library(RTMB)
# remotes::install_github("fishfollower/SAM/stockassessment")
library(stockassessment)
library(tictoc)
#fit<-fitfromweb("NEA_sei_21_v5Reca")
load("fit.RData")
dat <- list()
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
dat$covType <- c(0,1,2)
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
dat$ar1code <- -1  # 0 = centered, -1 = non-centered
dat$est_rec_intercept <- FALSE # Flag to estimate recruitment intercept
mulogobs <- dat$logobs

# Parameter section
par <- list()
par$logsdR <- 0
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
par$tPhi <- 1  # AR phi for recruitment
par$z <- rep(0, length(dat$year))  # Standard normal innovations
par$rec_intercept <- 0 # mapped off if not estimated

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

jnll <- function(par){
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

  ssb <- numeric(nrow)
  ssb[1] <- ssbFUN(logN[1,,drop=FALSE],logFF[1,,drop=FALSE],M[1,,drop=FALSE],
                   SW[1,,drop=FALSE],MO[1,,drop=FALSE],PF[1,,drop=FALSE],
                   PM[1,,drop=FALSE])

  jnll <- 0

  # NON-CENTERED PARAMETERIZATION
  if(ar1code == -1){
    # likelihood for standard normal innovations
    jnll <- jnll - dnorm(z[1], 0, sqrt(sdR*sdR/(1-phi*phi)), log=TRUE)
    jnll <- jnll - sum(dnorm(z[-1], 0, sdR, log=TRUE))

    # transform innovations (z) to get recruitment deviations
    rec_dev <- numeric(nrow)
    rec_dev[1] <- z[1]  # initial condition
    logN[1,1] <- rec_intercept + rec_dev[1]
    predN <- numeric((ncol-1)*(nrow-1))
    estN <- numeric((ncol-1)*(nrow-1))
    indN <- 0
    # 3. Set recruitment as random walk + transformed deviations
    for(y in 2:nrow){
      rec_dev[y] <- phi * rec_dev[y-1] + z[y]
      if(srmode==0){
        # Random walk with optional intercept
        predR <- rec_intercept
        logN[y,1] <- predR + rec_dev[y]
      }else{
        thisSSB <- ifelse((y-minAge-1)>(-.5), ssb[y-minAge], ssb[1])
        if(srmode==1){
          # Ricker: predicted recruitment + AR(1) deviation + intercept
          predR <- rickerpar[1] + log(thisSSB) - exp(rickerpar[2])*thisSSB
          logN[y,1] <- predR + rec_dev[y]
        }
        if(srmode==2){
          # Beverton-Holt: predicted recruitment + AR(1) deviation + intercept
          predR <- bhpar[1] + log(thisSSB) - log(1.0 + exp(bhpar[2])*thisSSB)
          logN[y,1] <- predR + rec_dev[y]
        }
      }


      for(a in 2:ncol){
        indN <- indN + 1
        predN[indN] <- logN[y-1,a-1] - exp(logFF[y-1,a-1]) - M[y-1,a-1]
        if(a==ncol){
          predN[indN] <- log(exp(predN[indN]) + exp(logN[y-1,a] - exp(logFF[y-1,a]) - M[y-1,a]))
        }
        estN[indN] <- logN[y,a]
      }
      ssb[y] <- ssbFUN(logN[y,,drop=FALSE],logFF[y,,drop=FALSE],M[y,,drop=FALSE],
                       SW[y,,drop=FALSE],MO[y,,drop=FALSE],PF[y,,drop=FALSE],
                       PM[y,,drop=FALSE])
    }
    jnll <- jnll - sum(dnorm(estN, predN, sdS, TRUE))
  }

  # CENTERED PARAMETERIZATION
  if(ar1code == 0){
    # AR1 init conditions
    jnll <- jnll - dnorm(logN[1,1], rec_intercept, sdR/sqrt(1-phi^2), log=TRUE)

    predN <- numeric((ncol-1)*(nrow-1))
    estN <- numeric((ncol-1)*(nrow-1))
    predR <- numeric((nrow-1))
    indN <- 0
    indR <- 0
    for(y in 2:nrow){
      indR <- indR + 1
      if(srmode==0){
        predR[indR] <- rec_intercept
      }else{
        thisSSB <- ifelse((y-minAge-1)>(-.5), ssb[y-minAge], ssb[1])
        if(srmode==1){
          # Ricker: predicted recruitment + AR(1) deviation + intercept
          predR[indR] <- rickerpar[1] + log(thisSSB) - exp(rickerpar[2])*thisSSB
        }
        if(srmode==2){
          # Beverton-Holt: predicted recruitment + AR(1) deviation + intercept
          predR[indR] <- bhpar[1] + log(thisSSB) - log(1.0 + exp(bhpar[2])*thisSSB)
        }
      }

      for(a in 2:ncol){
        indN <- indN + 1
        predN[indN] <- logN[y-1,a-1] - exp(logFF[y-1,a-1]) - M[y-1,a-1]
        if(a==ncol){
          predN[indN] <- log(exp(predN[indN]) + exp(logN[y-1,a] - exp(logFF[y-1,a]) - M[y-1,a]))
        }
        estN[indN] <- logN[y,a]
      }

      ssb[y] <- ssbFUN(logN[y,,drop=FALSE],logFF[y,,drop=FALSE],M[y,,drop=FALSE],
                       SW[y,,drop=FALSE],MO[y,,drop=FALSE],PF[y,,drop=FALSE],
                       PM[y,,drop=FALSE])

    }
    jnll <- jnll - sum(dnorm(estN, predN, sdS, TRUE))
    jnll <- jnll - sum(dnorm(logN[-1,1], predR + phi*(logN[-nrow,1] - predR), sdR, TRUE))
  }

  # F part
  SigmaF <- matrix(0, ncol(logF), ncol(logF))

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
    jnll <- jnll - dmvnorm(logF[y,], logF[y-1,], SigmaF, log=TRUE)
  }

  logPred <- numeric(nobs)
  for(i in 1:nobs){
    y <- aux[i,1] - minYear + 1
    f <- aux[i,2]
    a <- aux[i,3] - minAge + 1
    Z <- exp(logFF[y,a]) + M[y,a]
    if(fleetTypes[f]==0){
      logPred[i] <- logN[y,a] - log(Z) + log(1-exp(-Z)) + logFF[y,a]
    }
    if(fleetTypes[f]==2){
      logPred[i] <- logQ[keyQ[f,a]+1] + logN[y,a] - Z*sampleTimes[f]
    }
  }

  Svec <- list()
  for(f in 1:nrow(idx1)){
    thisdim <- sum(!is.na(keySd)[f,])
    S <- matrix(0, thisdim, thisdim)
    if(covType[f]==0){
      diag(S) <- sd[na.omit(keySd[f,])+1]^2
    }
    if(covType[f]==1){
      dist <- numeric(thisdim)
      d=2;
      for(a in 1:ncol(keyIGAR)){
        if(!is.na(keyIGAR[f,a])){
          dist[d] <- dist[d-1] + exp(logIGARdist[keyIGAR[f,a]+1])
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
        jnll <- jnll - dmvnorm(logobs[idx], logPred[idx], Svec[[f]], log=TRUE)
      }
    }
  }

  REPORT(logPred)
  ADREPORT(ssb)
  if(ar1code == -1) {
    ADREPORT(rec_dev)  # recruitment deviations
  }
  ADREPORT(logN[,1])  # report recruitment
  ADREPORT(rec_intercept)  # report intercept
  jnll
}

###################################################################################
# This was getting overwhelming, so adding this map function to map off things we don't estimate
###################################################################################
create_map_list <- function(est_rec_intercept = FALSE, ar1code = 0) {
  map_list <- list(
    logsdF = as.factor(rep(0, length(par$logsdF)))
  )

  # Map off rec_intercept if not estimated
  if (!est_rec_intercept) {
    map_list$rec_intercept <- as.factor(NA)
  }

  # For non-centered version, map off rec column of logN
  if (ar1code == -1) {
    nyr <- length(dat$year)
    nage <- length(dat$age)
    logN_map <- matrix(1:(nyr * nage), nyr, nage)
    logN_map[,1] <- NA  # Map off recruitment
    map_list$logN <- as.factor(logN_map)
  }

  return(map_list)
}


run_speed_test <- TRUE
###################################################################################
# Test centered version without intercept
###################################################################################
dat$ar1code <- 0
dat$est_rec_intercept <- FALSE
par_centered <- par
par_centered$z <- NULL  # remove z for centered version

obj_centered_no_int <- MakeADFun(jnll, par_centered,
                                 random=c("logN", "logF", "missing"),
                                 map=create_map_list(est_rec_intercept = FALSE, ar1code = 0),
                                 silent=TRUE)

opt_centered_no_int <- nlminb(obj_centered_no_int$par, obj_centered_no_int$fn, obj_centered_no_int$gr,
                              control=list(eval.max=1000, iter.max=1000))
sdr_centered_no_int <- sdreport(obj_centered_no_int)


if(run_speed_test) { # 348.225 # NV 67.09
  tic()
  for(i in 1:25) {
    set.seed(i)
    dat$logobs <- mulogobs + rnorm(length(dat$logobs), 0, sd = 0.03)
    obj_centered_no_int <- MakeADFun(jnll, par_centered,
                                     random=c("logN", "logF", "missing"),
                                     map=create_map_list(est_rec_intercept = FALSE, ar1code = 0),
                                     silent=TRUE)
    opt_centered_no_int <- nlminb(obj_centered_no_int$par, obj_centered_no_int$fn, obj_centered_no_int$gr,
                                  control=list(eval.max=1000, iter.max=1000))
  }
  toc()
}

run_speed_test <- FALSE
###################################################################################
# test centered version with intercept
###################################################################################
dat$est_rec_intercept <- TRUE

obj_centered_int <- MakeADFun(jnll, par_centered,
                              random=c("logN", "logF", "missing"),
                              map=create_map_list(est_rec_intercept = TRUE, ar1code = 0),
                              silent=TRUE)

opt_centered_int <- nlminb(obj_centered_int$par, obj_centered_int$fn, obj_centered_int$gr,
                           control=list(eval.max=1000, iter.max=1000))
sdr_centered_int <- sdreport(obj_centered_int)


if(run_speed_test) { # 312.867 sec # NV vec 18.66 (5 its no sdreport) # NV vec 21.72 (5 its with no sdreport)
  tic()
  for(i in 1:5) {
    set.seed(i)
    dat$logobs <- mulogobs + rnorm(length(dat$logobs), 0, sd = 0.03)
    obj_centered_int <- MakeADFun(jnll, par_centered,
                                  random=c("logN", "logF", "missing"),
                                  map=create_map_list(est_rec_intercept = TRUE, ar1code = 0),
                                  silent=TRUE)

    opt_centered_int <- nlminb(obj_centered_int$par, obj_centered_int$fn, obj_centered_int$gr,
                               control=list(eval.max=1000, iter.max=1000))
    sdr_centered_int <- sdreport(obj_centered_int)
  }
  toc()
}
###################################################################################
# test non-centered version without intercept
###################################################################################
dat$ar1code <- -1
dat$est_rec_intercept <- FALSE

obj_noncentered_no_int <- MakeADFun(jnll, par,
                                    random=c("z", "logN", "logF", "missing"),
                                    map=create_map_list(est_rec_intercept = FALSE, ar1code = -1),
                                    silent=TRUE)

opt_noncentered_no_int <- nlminb(obj_noncentered_no_int$par, obj_noncentered_no_int$fn, obj_noncentered_no_int$gr,
                                 control=list(eval.max=1000, iter.max=1000))
sdr_noncentered_noint <- sdreport(obj_noncentered_no_int)

if(run_speed_test) { # 167.642
  tic()
  for(i in 1:25) {
    set.seed(i)
    dat$logobs <- mulogobs + rnorm(length(dat$logobs), 0, sd = 0.03)
    obj_noncentered_no_int <- MakeADFun(jnll, par,
                                        random=c("z", "logN", "logF", "missing"),
                                        map=create_map_list(est_rec_intercept = FALSE, ar1code = -1),
                                        silent=TRUE)

    opt_noncentered_no_int <- nlminb(obj_noncentered_no_int$par, obj_noncentered_no_int$fn, obj_noncentered_no_int$gr,
                                     control=list(eval.max=1000, iter.max=1000))
  }
  toc()
}
###################################################################################
# test non-centered version with intercept
###################################################################################
dat$est_rec_intercept <- TRUE

obj_noncentered_int <- MakeADFun(jnll, par,
                                 random=c("z", "logN", "logF", "missing"),
                                 map=create_map_list(est_rec_intercept = TRUE, ar1code = -1),
                                 silent=TRUE)

opt_noncentered_int <- nlminb(obj_noncentered_int$par, obj_noncentered_int$fn, obj_noncentered_int$gr,
                              control=list(eval.max=1000, iter.max=1000))
sdr_noncentered_int <- sdreport(obj_noncentered_int)

if(run_speed_test) { # 192.298 sec # NV vec 44.75 (5 its no sdreport) # NV vec 51.78 (5 its with sdreport)
  tic()
  for(i in 1:5) {
    set.seed(i)
    dat$logobs <- mulogobs + rnorm(length(dat$logobs), 0, sd = 0.03)
    obj_noncentered_int <- MakeADFun(jnll, par,
                                     random=c("z", "logN", "logF", "missing"),
                                     map=create_map_list(est_rec_intercept = TRUE, ar1code = -1),
                                     silent=TRUE)

    opt_noncentered_int <- nlminb(obj_noncentered_int$par, obj_noncentered_int$fn, obj_noncentered_int$gr,
                                  control=list(eval.max=1000, iter.max=1000))

    sdr_noncentered_int <- sdreport(obj_noncentered_int)

  }
  toc()
}
###################################################################################
# We have two pairs of results that are identical. Models with the intercept have the same
# values, and models without the intercepts have the same values.
#
# Test that all fixed effects parameters are the same -- this includes phi
###################################################################################
all((abs(sdr_centered_int$par.fixed - sdr_noncentered_int$par.fixed) < 0.001))
all((abs(diag(sdr_centered_int$cov.fixed) - diag(sdr_noncentered_int$cov.fixed)) < 0.001))

all((abs(sdr_centered_no_int$par.fixed - sdr_noncentered_noint$par.fixed) < 0.001))
all((abs(diag(sdr_centered_no_int$cov.fixed) - diag(sdr_noncentered_noint$cov.fixed)) < 0.001))

###################################################################################
# Test that all recruitment is the same
###################################################################################
logN_vals <- sdr_centered_no_int$par.random[which(names(sdr_centered_no_int$par.random) == "logN")]
logN_matrix <- matrix(logN_vals, nrow = nyr, ncol = nage)

# For the centered model, this is the first column of logN
recruitment_centered <- logN_matrix[,1]
# For non - centered model, this is a little trickier
# xxtract z innovations
z <- sdr_noncentered_noint$par.random[names(sdr_noncentered_noint$par.random) == "z"]
phi <- 2*plogis(sdr_noncentered_noint$par.fixed["tPhi"]) - 1
# Transform z to rec_dev (the AR(1) process)
nyr <- length(recruitment_centered)
rec_dev <- numeric(nyr)
rec_dev[1] <- z[1]  # initial condition
for(i in 2:nyr) {
  rec_dev[i] <- phi * rec_dev[i-1] + z[i]
}
# For srmode=0, logN[,1] should equal rec_dev
recruitment_noncentered <- rec_dev
# Test whether these equal
all((abs(recruitment_noncentered - recruitment_centered) < 0.001))



# This is the same test for the versions of the model with the intercept
logN_vals <- sdr_centered_int$par.random[which(names(sdr_centered_int$par.random) == "logN")]
logN_matrix <- matrix(logN_vals, nrow = nyr, ncol = nage)

# For the centered model, this is the first column of logN
recruitment_centered <- logN_matrix[,1]
# For non - centered model, this is a little trickier
# xxtract z innovations
z <- sdr_noncentered_int$par.random[names(sdr_noncentered_int$par.random) == "z"]
phi <- 2*plogis(sdr_noncentered_int$par.fixed["tPhi"]) - 1
# Transform z to rec_dev (the AR(1) process)
nyr <- length(recruitment_centered)
rec_dev <- numeric(nyr)
rec_dev[1] <- z[1]  # initial condition
for(i in 2:nyr) {
  rec_dev[i] <- phi * rec_dev[i-1] + z[i]
}
# For srmode=0, logN[,1] should equal rec_dev
rec_intercept <- sdr_noncentered_int$par.fixed["rec_intercept"]
recruitment_noncentered <- rec_intercept + rec_dev
# Test whether these equal
all((abs(recruitment_noncentered - recruitment_centered) < 0.001))
