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

  if(code==-2){
    nll <- nll - dnorm(x[1],0,sqrt(sd*sd/(1-phi*phi)),log=TRUE)
    nll <- nll - sum(dnorm(x[-1],phi*x[-(length(x))],sd,log=TRUE))
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
                         all_eq = 0, time=NA, nll=NA, iters=0,
                         evals_F=0, evals_G=0, logSigma=NA,
                         tPhi=NA, logSigmaObs=NA, converged=NA)

grid_pars_1 <- grid_pars
grid_pars_1$model <- "Centered"
grid_pars_2 <- grid_pars
grid_pars_2$model <- "Non-centered"
grid_pars_3 <- grid_pars
grid_pars_3$model <- "dautoreg"
grid_pars_4 <- grid_pars
grid_pars_4$model <- "Sparse matrix"
grid_pars_5 <- grid_pars
grid_pars_5$model <- "Alt-Non-centered"

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
opt1 <- nlminb(f1$par, f1$fn, f1$gr)
timr = toc()
grid_pars_1$time[ii] <- timr$toc - timr$tic
grid_pars_1$converged[ii] <- opt1$convergence
grid_pars_1$nll[ii] <- opt1$objective
grid_pars_1$iters[ii] <- opt1$iterations
grid_pars_1$evals_F[ii] <- opt1$evaluations[["function"]]
grid_pars_1$evals_G[ii] <- opt1$evaluations[["gradient"]]
grid_pars_1$logSigma[ii] <- opt1$par[["logSigma"]]
grid_pars_1$tPhi[ii] <- opt1$par[["tPhi"]]
grid_pars_1$logSigmaObs[ii] <- opt1$par[["logSigmaObs"]]


dat$code=0
tic()
f2 <- MakeADFun(f,par, silent=TRUE, random="x")
opt2 <- nlminb(f2$par, f2$fn, f2$gr)
timr = toc()
grid_pars_2$time[ii] <- timr$toc - timr$tic
grid_pars_2$converged[ii] <- opt2$convergence
grid_pars_2$nll[ii] <- opt2$objective
grid_pars_2$iters[ii] <- opt2$iterations
grid_pars_2$evals_F[ii] <- opt2$evaluations[["function"]]
grid_pars_2$evals_G[ii] <- opt2$evaluations[["gradient"]]
grid_pars_2$logSigma[ii] <- opt2$par[["logSigma"]]
grid_pars_2$tPhi[ii] <- opt2$par[["tPhi"]]
grid_pars_2$logSigmaObs[ii] <- opt2$par[["logSigmaObs"]]

dat$code=1
tic()
f3 <- MakeADFun(f,par, silent=TRUE, random="x")
opt3 <- nlminb(f3$par, f3$fn, f3$gr)
timr = toc()
grid_pars_3$time[ii] <- timr$toc - timr$tic
grid_pars_3$converged[ii] <- opt3$convergence
grid_pars_3$nll[ii] <- opt3$objective
grid_pars_3$iters[ii] <- opt3$iterations
grid_pars_3$evals_F[ii] <- opt3$evaluations[["function"]]
grid_pars_3$evals_G[ii] <- opt3$evaluations[["gradient"]]
grid_pars_3$logSigma[ii] <- opt3$par[["logSigma"]]
grid_pars_3$tPhi[ii] <- opt3$par[["tPhi"]]
grid_pars_3$logSigmaObs[ii] <- opt3$par[["logSigmaObs"]]

dat$code=3
tic()
f4 <- MakeADFun(f,par, silent=TRUE, random="x")
opt4 <- nlminb(f4$par, f4$fn, f4$gr)
timr = toc()
grid_pars_4$time[ii] <- timr$toc - timr$tic
grid_pars_4$converged[ii] <- opt4$convergence
grid_pars_4$nll[ii] <- opt4$objective
grid_pars_4$iters[ii] <- opt4$iterations
grid_pars_4$evals_F[ii] <- opt4$evaluations[["function"]]
grid_pars_4$evals_G[ii] <- opt4$evaluations[["gradient"]]
grid_pars_4$logSigma[ii] <- opt4$par[["logSigma"]]
grid_pars_4$tPhi[ii] <- opt4$par[["tPhi"]]
grid_pars_4$logSigmaObs[ii] <- opt4$par[["logSigmaObs"]]

dat$code=-2
tic()
f5 <- MakeADFun(f,par, silent=TRUE, random="x")
opt5 <- nlminb(f5$par, f5$fn, f5$gr)
timr = toc()
grid_pars_5$time[ii] <- timr$toc - timr$tic
grid_pars_5$converged[ii] <- opt5$convergence
grid_pars_5$nll[ii] <- opt5$objective
grid_pars_5$iters[ii] <- opt5$iterations
grid_pars_5$evals_F[ii] <- opt5$evaluations[["function"]]
grid_pars_5$evals_G[ii] <- opt5$evaluations[["gradient"]]
grid_pars_5$logSigma[ii] <- opt5$par[["logSigma"]]
grid_pars_5$tPhi[ii] <- opt5$par[["tPhi"]]
grid_pars_5$logSigmaObs[ii] <- opt5$par[["logSigmaObs"]]

# Test whether parameters are identical
all_eq <- FALSE
if(all.equal(opt1$par, opt2$par)==TRUE && all.equal(opt1$par, opt3$par)==TRUE && all.equal(opt1$par, opt4$par)==TRUE && all.equal(opt1$par, opt5$par)==TRUE) all_eq <- TRUE
grid_pars$all_eq[ii] <- as.numeric(all_eq)
}

pars <- rbind(grid_pars_1, grid_pars_2, grid_pars_3, grid_pars_4, grid_pars_5)

dplyr::group_by(pars, model, N) |>
  dplyr::summarise(mean_time = mean(time), sd_time = sd(time)) |>
  ggplot(aes(N, mean_time, color=model)) +
  geom_line() +
  scale_color_viridis_d(option="magma", begin=0.2, end=0.8) +
  theme_bw() +
  xlab("Time series length") +
  ylab("Mean time (s)")
ggsave("figures/ar1_simtest_opt.png", width=6, height=5)

dplyr::group_by(pars, model, N) |>
  dplyr::summarise(mean_iter = mean(iters), sd_time = sd(iters)) |>
  ggplot(aes(N, mean_iter, color=model)) +
  geom_line() +
  scale_color_viridis_d(option="magma", begin=0.2, end=0.8) +
  theme_bw() +
  xlab("Time series length") +
  ylab("Mean iterations (n)")
ggsave("figures/ar1_simtest_iters.png", width=6, height=5)

dplyr::group_by(pars, model, N) |>
  dplyr::summarise(mean_evals_F = mean(evals_F), sd_time = sd(evals_F)) |>
  ggplot(aes(N, mean_evals_F, color=model)) +
  geom_line() +
  scale_color_viridis_d(option="magma", begin=0.2, end=0.8) +
  theme_bw() +
  xlab("Time series length") +
  ylab("Mean function evaluations (n)")
ggsave("figures/ar1_simtest_evals_F.png", width=6, height=5)


dplyr::group_by(pars, model, N) |>
  dplyr::summarise(mean_evals_G = mean(evals_G), sd_time = sd(evals_G)) |>
  ggplot(aes(N, mean_evals_G, color=model)) +
  geom_line() +
  scale_color_viridis_d(option="magma", begin=0.2, end=0.8) +
  theme_bw() +
  xlab("Time series length") +
  ylab("Mean gradient evaluations (n)")
ggsave("figures/ar1_simtest_evals_G.png", width=6, height=5)

