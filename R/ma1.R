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

dat <- list(y=y, code=0, predict=50)

f<-function(par){
  getAll(dat,par)
  theta <- 2*plogis(tTheta)-1
  timeSteps <-length(y)
  sd <- exp(logSigma)
  sdo <- exp(logSigmaObs)
  io <- 1:length(y)
  nll <- 0

  if(code==-1){ # here x is epsilon
    nll <- nll -sum(dnorm(x,0,sd=sd, log=TRUE))
    lam <- numeric(timeSteps+predict)
    for(i in 1:(timeSteps+predict)){
      lam[i] <- x[i+1]+theta*x[i]
    }                  
    nll <- nll - sum(dnorm(y, lam[io], sdo, TRUE))
  }

  if(code==0){ # here x is the process
    eps <- numeric(timeSteps + predict + 1) 
    eps[1] <- eps0 
    for (i in 1:(timeSteps+predict)){
      eps[i+1] <- x[i] - theta*eps[i]
    }
    nll <- nll - sum(dnorm(eps, 0, sd, log = TRUE))
    nll <- nll - sum(dnorm(y, x[io], sdo, TRUE))
  }
  
  if(code==2){
    S <- diag(timeSteps+predict)*(1+theta*theta)*sd*sd;
    S[abs(row(S)-col(S))==1] <- sd*sd*theta
    nll <- nll - dmvnorm(x,0,S,log=TRUE)
    nll <- nll - sum(dnorm(y, x[io], sdo, TRUE))
  }

  if(code==3){ # first attemp to speed up the dmvnorm - can likely be improved
    #browser()  
    d <- rep((1+theta*theta)*sd*sd,timeSteps+predict)
    s <- rep(sd*sd*theta, timeSteps+predict-1)
    dL <- numeric(timeSteps+predict)
    sL <- numeric(timeSteps+predict-1)
    dL[1] <- sqrt(d[1])
    for (i in 2:(timeSteps+predict)){
      sL[i-1] <- s[i-1] / dL[i-1]   
      dL[i] <- sqrt(d[i] - sL[i-1]^2)
    }
    logDet <- 2*sum(log(dL))
    r <- numeric(timeSteps+predict)
    r[1] <- (x[1] - dL[1] * r[1]) / dL[1]
    for (i in 2:(timeSteps+predict)) {
      r[i] <- (x[i] - sL[i-1] * r[i-1]) / dL[i]
    } 
    nll <- nll + 0.5 * ((timeSteps+predict) * log(2*pi) + logDet + sum(r^2))
    nll <- nll - sum(dnorm(y, x[io], sdo, TRUE))
  }

  
  REPORT(theta)
  REPORT(sd)
  REPORT(sdo)
  nll
}

pdf("ma1hess%03d.pdf", width=6, height=6, onefile=FALSE)
dat$code <- -1
par <- list(logSigma=.1, tTheta=3, logSigmaObs=log(0.5), x=rep(0,length(x)+dat$predict+1))
obj <- MakeADFun(f,par, silent=TRUE, random="x", map=list(logSigmaObs=as.factor(NA)))
cat("eps :                   ", paste0(round(system.time(opt<-nlminb(obj$par,obj$fn,obj$gr)),5)[3], "s, par = "))
cat(c(opt$obj, opt$par), "\n")
plot(Matrix::image(obj$env$spHess(random=TRUE), main="Epsilons"))
# works, but very slow
dat$code <- 2
par <- list(logSigma=.1, tTheta=3, logSigmaObs=log(0.5), x=rep(0,length(x)+dat$predict))
obj <- MakeADFun(f,par, silent=TRUE, random="x", map=list(logSigmaObs=as.factor(NA)))
cat("dmvnorm:                ", paste0(round(system.time(opt<-nlminb(obj$par,obj$fn,obj$gr)),5)[3], "s, par = "))
cat(c(opt$obj, opt$par), "\n")
plot(Matrix::image(obj$env$spHess(random=TRUE), main="Process using dmvnorm"))
# optimized code sometimes unstable when doing predictions
#dat$code <- 3
#par <- list(logSigma=-1, tTheta=-1, logSigmaObs=log(0.5), x=rep(0,length(x)+dat$predict))
#obj <- MakeADFun(f,par, silent=TRUE, random="x", map=list(logSigmaObs=as.factor(NA)))
#cat("dmvnorm-sparse:         ", paste0(round(system.time(opt<-nlminb(obj$par,obj$fn,obj$gr)),5)[3], "s, par = "))
#cat(c(opt$obj, opt$par), "\n")
#plot(Matrix::image(obj$env$spHess(random=TRUE), main="Process"))
dat$code <- 0
par <- list(logSigma=.1, tTheta=3, logSigmaObs=log(0.5), x=rep(0,length(x)+dat$predict), eps0=0)
obj <- MakeADFun(f,par, silent=TRUE, random=c("x", "eps0"), map=list(logSigmaObs=as.factor(NA)))
cat("seq :                   ", paste0(round(system.time(opt<-nlminb(obj$par,obj$fn,obj$gr)),5)[3], "s, par = "))
cat(c(opt$obj, opt$par), "\n")
plot(Matrix::image(obj$env$spHess(random=TRUE), main="Process and eps0"))
dev.off()

pdf("predictma1.pdf", width=12, height=6)
  sdr <- sdreport(obj)
  pl <- as.list(sdr, "Est")
  plsd <- as.list(sdr, "Std")
  plot(y, xlim=c(1,length(pl$x)), ylim=c(min(pl$x-2*plsd$x), max(pl$x+2*plsd$x)))
  lines(pl$x, lwd=2, col="navy")
  lines(pl$x-2*plsd$x, lwd=2, lty="dashed", col="navy")
  lines(pl$x+2*plsd$x, lwd=2, lty="dashed", col="navy")
dev.off()
