library(RTMB)
N<-150
sd <- 1
theta <- .9
sdo <- .5
eps <- rnorm(N+1,sd=sd)
x <- numeric(N)
for(i in 1:N){
  x[i] <- eps[i+1]+theta*eps[i]
}                  
y <- x+rnorm(length(x), sd=sdo)

dat <- list(y=y, code=0)
f<-function(par){
  getAll(dat,par)
  theta <- 2*plogis(tTheta)-1
  timeSteps <-length(y)
  sd <- exp(logSigma)
  sdo <- exp(logSigmaObs)
  nll <- 0

  if(code==-1){ # here x is epsilon
    nll <- nll -sum(dnorm(x,0,sd=sd, log=TRUE))
    lam <- numeric(timeSteps)
    for(i in 1:timeSteps){
      lam[i] <- x[i+1]+theta*x[i]
    }                  
    nll <- nll - sum(dnorm(y, lam, sdo, TRUE))
  }

  if(code==0){ # here x is the process
    eps <- numeric(timeSteps + 1) 
    eps[1] <- eps0 
    for (i in 1:timeSteps){
      eps[i+1] <- x[i] - theta*eps[i]
    }
    nll <- nll - sum(dnorm(eps, 0, sd, log = TRUE))
    nll <- nll - sum(dnorm(y, x, sdo, TRUE))
  }
  
  if(code==2){
    S <- diag(timeSteps)*(1+theta*theta)*sd*sd;
    S[abs(row(S)-col(S))==1] <- sd*sd*theta
    nll <- nll - dmvnorm(x,0,S,log=TRUE)
    nll <- nll - sum(dnorm(y, x, sdo, TRUE))
  }

  if(code==3){ # first attemp to speed up the dmvnorm - can likely be improved 
    d <- rep((1+theta*theta)*sd*sd,timeSteps)
    s <- rep(sd*sd*theta, timeSteps-1)
    dL <- numeric(timeSteps)
    sL <- numeric(timeSteps-1)
    dL[1] <- sqrt(d[1])
    for (i in 2:timeSteps){
      sL[i-1] <- s[i-1] / dL[i-1]   
      dL[i] <- sqrt(d[i] - sL[i-1]^2)
    }
    logDet <- 2*sum(log(dL))
    r <- numeric(timeSteps)
    r[1] <- (x[1] - dL[1] * r[1]) / dL[1]
    for (i in 2:timeSteps) {
      r[i] <- (x[i] - sL[i-1] * r[i-1]) / dL[i]
    } 
    nll <- nll + 0.5 * (timeSteps * log(2*pi) + logDet + sum(r^2))
    nll <- nll - sum(dnorm(y, x, sdo, TRUE))
  }

  
  REPORT(theta)
  REPORT(sd)
  REPORT(sdo)
  nll
}

pdf("ma1hess%03d.pdf", width=6, height=6, onefile=FALSE)
dat$code <- -1
par <- list(logSigma=.1, tTheta=3, logSigmaObs=log(0.5), x=rep(0,length(x)+1))
obj <- MakeADFun(f,par, silent=TRUE, random="x", map=list(logSigmaObs=as.factor(NA)))
cat("eps :                   ", paste0(round(system.time(opt<-nlminb(obj$par,obj$fn,obj$gr)),5)[3], "s, par = "))
cat(c(opt$obj, opt$par), "\n")
plot(Matrix::image(obj$env$spHess(random=TRUE), main="Epsilons"))
dat$code <- 0
par <- list(logSigma=.1, tTheta=3, logSigmaObs=log(0.5), x=rep(0,length(x)), eps0=0)
obj <- MakeADFun(f,par, silent=TRUE, random=c("x", "eps0"), map=list(logSigmaObs=as.factor(NA)))
cat("seq :                   ", paste0(round(system.time(opt<-nlminb(obj$par,obj$fn,obj$gr)),5)[3], "s, par = "))
cat(c(opt$obj, opt$par), "\n")
plot(Matrix::image(obj$env$spHess(random=TRUE), main="Process and eps0"))
# works, but very slow
#dat$code <- 2
#par <- list(logSigma=.1, tTheta=3, logSigmaObs=log(0.5), x=rep(0,length(x)))
#obj <- MakeADFun(f,par, silent=TRUE, random="x", map=list(logSigmaObs=as.factor(NA)))
#cat("dmvnorm:                ", paste0(round(system.time(opt<-nlminb(obj$par,obj$fn,obj$gr)),5)[3], "s, par = "))
#cat(c(opt$obj, opt$par), "\n")
#plot(Matrix::image(obj$env$spHess(random=TRUE), main="Process using dmvnorm"))
dat$code <- 3
par <- list(logSigma=.1, tTheta=3, logSigmaObs=log(0.5), x=rep(0,length(x)))
obj <- MakeADFun(f,par, silent=TRUE, random="x", map=list(logSigmaObs=as.factor(NA)))
cat("dmvnorm-sparse:         ", paste0(round(system.time(opt<-nlminb(obj$par,obj$fn,obj$gr)),5)[3], "s, par = "))
cat(c(opt$obj, opt$par), "\n")
plot(Matrix::image(obj$env$spHess(random=TRUE), main="Process"))

dev.off()
