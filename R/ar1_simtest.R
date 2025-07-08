library(RTMB)
library(tictoc)
library(ggplot2)
library(viridis)

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

grid_pars <- expand.grid(seed = 1:100, N = c(30, 50, 100),
                         sd = c(0.1, 0.5, 1),
                         sdo = c(0.1, 0.5, 1), phi = c(0, 0.5, 0.9),
                         all_eq = 0)

for(ii in 1:nrow(grid_pars)) {

set.seed(grid_pars$seed[ii])

N<-grid_pars$N[ii]
x <- numeric(N)
sd <- grid_pars$sd[ii]
sdo <- grid_pars$sdo[ii]
phi <- grid_pars$phi[ii]
# init condition
x[1] <- rnorm(1, sd=sqrt(sd*sd/(1-phi*phi)))
# simulate ss model
for(i in 2:N){
  x[i] <- phi*x[i-1]+rnorm(1,sd=sd)
}
y <- x+rnorm(length(x), sd=sdo)

dat <- list(y=y, code=0)
par <- list(logSigma=runif(1,0.5,1.5), tPhi=1,
            logSigmaObs=log(runif(1,0.5,1.5)), x=rep(0,length(x)))

dat$code = -1
f1 <- MakeADFun(f,par, silent=TRUE, random="x")

dat$code=0
f2 <- MakeADFun(f,par, silent=TRUE, random="x")

dat$code=1
f3 <- MakeADFun(f,par, silent=TRUE, random="x")

dat$code=3
f4 <- MakeADFun(f,par, silent=TRUE, random="x")
# Test whether parameters are identical
all_eq <- FALSE
if(all.equal(f1$par, f2$par) && all.equal(f1$par, f3$par) && all.equal(f1$par, f4$par)) all_eq <- TRUE
grid_pars$all_eq[ii] <- as.numeric(all_eq)
}

# Result: all 8100 parameter sets above are identical
grid_pars <- expand.grid(seed = 1:500, N = c(100, 500, 1000),
                         sd = c(1),
                         sdo = c(0.1), phi = c(0.9),
                         time = NA)
grid_pars_1 <- grid_pars
grid_pars_1$model <- "Centered"
grid_pars_2 <- grid_pars
grid_pars_2$model <- "Non-centered"
grid_pars_3 <- grid_pars
grid_pars_3$model <- "dautoreg"
grid_pars_4 <- grid_pars
grid_pars_4$model <- "Sparse matrix"

for(ii in 1:nrow(grid_pars)) {

  set.seed(grid_pars$seed[ii])

  N<-grid_pars$N[ii]
  x <- numeric(N)
  sd <- grid_pars$sd[ii]
  sdo <- grid_pars$sdo[ii]
  phi <- grid_pars$phi[ii]
  # init condition
  x[1] <- rnorm(1, sd=sqrt(sd*sd/(1-phi*phi)))
  # simulate ss model
  for(i in 2:N){
    x[i] <- phi*x[i-1]+rnorm(1,sd=sd)
  }
  y <- x+rnorm(length(x), sd=sdo)

  dat <- list(y=y, code=0)
  par <- list(logSigma=runif(1,0.5,1.5), tPhi=1,
              logSigmaObs=log(runif(1,0.5,1.5)), x=rep(0,length(x)))

  dat$code = -1
  tic()
  f1 <- MakeADFun(f,par, silent=TRUE, random="x")
  timr = toc()
  grid_pars_1$time[ii] <- timr$toc - timr$tic

  dat$code=0
  tic()
  f2 <- MakeADFun(f,par, silent=TRUE, random="x")
  timr = toc()
  grid_pars_2$time[ii] <- timr$toc - timr$tic

  dat$code=1
  tic()
  f3 <- MakeADFun(f,par, silent=TRUE, random="x")
  timr = toc()
  grid_pars_3$time[ii] <- timr$toc - timr$tic

  dat$code=3
  tic()
  f4 <- MakeADFun(f,par, silent=TRUE, random="x")
  timr = toc()
  grid_pars_4$time[ii] <- timr$toc - timr$tic
}

pars <- rbind(grid_pars_1, grid_pars_2, grid_pars_3, grid_pars_4)

dplyr::group_by(pars, model, N) |>
  dplyr::summarise(mean_time = mean(time), sd_time = sd(time)) |>
  ggplot(aes(N, mean_time, color=model)) +
  geom_line() +
  scale_color_viridis_d(option="magma", begin=0.2, end=0.8) +
  theme_bw() +
  xlab("Time series length") +
  ylab("Mean time (s)")
ggsave("figures/ar1_simtest.png", width=6, height=5)

