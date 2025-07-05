library(RTMB)
N<-500
x <- numeric(N)
sd <- 1
sdo <- .5
phi <- .9
x[1] <- rnorm(1, sd=sqrt(sd*sd/(1-phi*phi)))
for(i in 2:N){
  x[i] <- phi*x[i-1]+rnorm(1,sd=sd)
}                  
y <- x+rnorm(length(x), sd=sdo)

dat <- list(y=y, code=0)
par <- list(logSigma=0, tPhi=1, logSigmaObs=log(0.1), x=rep(0,length(x)))

f<-function(par){
  getAll(dat,par)
  phi <- 2*plogis(tPhi)-1
  timeSteps <-length(x)
  sd <- exp(logSigma)
  sdo <- exp(logSigmaObs)
  nll<-0

  if(code==-1){  
    nll <- nll -dnorm(x[1],0,sqrt(sd*sd/(1-phi*phi)),log=TRUE)
    nll <- nll -sum(dnorm(x[-1],0,sd,log=TRUE))
    lam <- numeric(length(x))
    lam[1] <- 0+x[1]
    for(i in 2:timeSteps){    
      lam[i] <- phi*lam[i-1]+x[i] 
    }
    nll <- nll - sum(dnorm(y, lam, sdo, TRUE))
  }
  
  if(code==0){
    nll <- nll -dnorm(x[1],0,sqrt(sd*sd/(1-phi*phi)),log=TRUE)
    for(i in 2:timeSteps){    
      nll <- nll -dnorm(x[i],phi*x[i-1],sd,log=TRUE)
    }
    nll <- nll - sum(dnorm(y, x, sdo, TRUE))
  }
  
  if(code==1){
    nll <- nll - dautoreg(x,phi=phi,scale=sqrt(sd*sd/(1-phi*phi)), log=TRUE)
    nll <- nll - sum(dnorm(y, x, sdo, TRUE))  
  }

  if(code==2){
    t<-1:length(x)  
    S <- sd^2/(1-phi^2)*phi^abs(outer(t,t,"-"))  
    nll <- nll - dmvnorm(x,0,S,log=TRUE)
    nll <- nll - sum(dnorm(y, x, sdo, TRUE))
  }

  if(code==3){
    Q<-Matrix::spMatrix(timeSteps,timeSteps)
    diag(Q)<-c(1,rep(1+phi^2, timeSteps-2),1)
    Q[abs(row(Q)-col(Q))==1] <- -phi
    nll <- nll - dgmrf(x, Q=Q, scale=sd, log=TRUE)
    nll <- nll - sum(dnorm(y, x, sdo, TRUE))
  }
  ADREPORT(phi)
  nll
}

pdf("ar1hess%03d.pdf", width=6, height=6, onefile=FALSE)
dat$code = -1
obj <- MakeADFun(f,par, silent=TRUE, random="x")
cat("Deviance formulation:   ", paste0(round(system.time(opt<-nlminb(obj$par,obj$fn,obj$gr)),5)[3], "s, par = "))
cat(c(opt$obj, opt$par), "\n")
plot(Matrix::image(obj$env$spHess(random=TRUE), main="Deviances"))
dat$code=0
obj <- MakeADFun(f,par, silent=TRUE, random="x")
cat("Sequential:             ", paste0(round(system.time(opt<-nlminb(obj$par,obj$fn,obj$gr)),5)[3], "s, par = "))
cat(c(opt$obj, opt$par), "\n")
plot(Matrix::image(obj$env$spHess(random=TRUE), main="Sequential"))
dat$code=1
obj <- MakeADFun(f,par, silent=TRUE, random="x")
cat("dautoreg:               ", paste0(round(system.time(opt<-nlminb(obj$par,obj$fn,obj$gr)),5)[3], "s, par = "))
cat(c(opt$obj, opt$par), "\n")
plot(Matrix::image(obj$env$spHess(random=TRUE), main="dautoreg"))
# works but really slow so turned off
#dat$code=2
#obj <- MakeADFun(f,par, silent=TRUE, random="x")
#cat("dmvnorm:                ", paste0(round(system.time(opt<-nlminb(obj$par,obj$fn,obj$gr)),5)[3], "s, par = "))
#cat(c(opt$obj, opt$par), "\n")
#plot(Matrix::image(obj$env$spHess(random=TRUE), main="dmvnorm"))
dat$code=3
obj <- MakeADFun(f,par, silent=TRUE, random="x")
cat("dgmrf:                  ", paste0(round(system.time(opt<-nlminb(obj$par,obj$fn,obj$gr)),5)[3], "s, par = "))
cat(c(opt$obj, opt$par), "\n")
plot(Matrix::image(obj$env$spHess(random=TRUE), main="dgmrf"))
dev.off()
