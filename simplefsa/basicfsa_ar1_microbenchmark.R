library(microbenchmark)
library(lineprof)
library(dplyr)

source("simplefsa/basicfsa_ar_good.R")
source("simplefsa/basicfsa_ar_bad.R")

results_intercept <- microbenchmark(
  sparse_intercept = {
    obj_sparse <- MakeADFun(nll_sparse, par_sparse,
                        map=list(logFA=factor(c(1:4,NA,NA,NA))),
                        silent=TRUE,
                        random="logN1A")
    opt_sparse <- nlminb(obj_sparse$par, obj_sparse$fn, obj_sparse$gr,
                  control=list(iter.max=1000,eval.max=1000))
  },
  manual_intercept = {
    obj_manual <- MakeADFun(nll_manual, par_manual,
                        map=list(logFA=factor(c(1:4,NA,NA,NA))),
                        silent=TRUE,
                        random="x")
    opt_manual <- nlminb(obj_manual$par, obj_manual$fn, obj_manual$gr,
                  control=list(iter.max=1000,eval.max=1000))
  },
  times = 100
)
results_intercept |> summary() |> dplyr::select(expr, median)



print(results_intercept)

results_nointercept <- microbenchmark(
  sparse_intercept = {

    obj <- MakeADFun(nll_sparse, par_sparse,
                     map=list(logFA=factor(c(1:4,NA,NA,NA)),
                              logMuR = factor(NA)),
                     silent=TRUE,
                     random="logN1A")

    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control=list(iter.max=1000,eval.max=1000))

  },

  manual_intercept = {

    obj <- MakeADFun(nll_manual, par_manual,
                     map=list(logFA=factor(c(1:4,NA,NA,NA)),
                              logMuR = factor(NA)),
                     silent=TRUE,
                     random="z")

    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control=list(iter.max=1000,eval.max=1000))


  },
  times = 100
)


print(results_nointercept)  |> summary() |> dplyr::select(expr, median)
