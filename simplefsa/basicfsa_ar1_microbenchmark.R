library(bench)
library(lineprof)
library(dplyr)

source("simplefsa/basicfsa_ar_good.R")
source("simplefsa/basicfsa_ar_bad.R")

# Manual warm-up
for (i in 1:10) {
  obj_sparse <- MakeADFun(nll_sparse, par_sparse,
                          map = list(logFA = factor(c(1:4, NA, NA, NA))),
                          silent = TRUE,
                          random = "logN1A")
  opt_sparse <- nlminb(obj_sparse$par, obj_sparse$fn, obj_sparse$gr,
                       control = list(iter.max = 1000, eval.max = 1000))
  obj_manual <- MakeADFun(nll_manual, par_manual,
                          map = list(logFA = factor(c(1:4, NA, NA, NA))),
                          silent = TRUE,
                          random = "x")
  opt_manual <- nlminb(obj_manual$par, obj_manual$fn, obj_manual$gr,
                       control = list(iter.max = 1000, eval.max = 1000))
}


results_intercept <- 
bench::mark(
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
  iterations = 100, check = FALSE, relative = TRUE, filter_gc = FALSE
)
results_intercept |> print(width = Inf)



results_nointercept <- bench::mark(
  sparse_no_intercept = {

    obj <- MakeADFun(nll_sparse, par_sparse,
                     map=list(logFA=factor(c(1:4,NA,NA,NA)),
                              logMuR = factor(NA)),
                     silent=TRUE,
                     random="logN1A")

    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control=list(iter.max=1000,eval.max=1000))

  },

  manual_no_intercept = {

    obj <- MakeADFun(nll_manual, par_manual,
                     map=list(logFA=factor(c(1:4,NA,NA,NA)),
                              logMuR = factor(NA)),
                     silent=TRUE,
                     random="x")

    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control=list(iter.max=1000,eval.max=1000))


  },
  iterations = 100, check = FALSE, relative = TRUE
)


results_nointercept |> print(width = Inf)
