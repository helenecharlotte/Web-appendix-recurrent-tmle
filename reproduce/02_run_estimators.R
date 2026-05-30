### 02_run_estimators.R --- 
#----------------------------------------------------------------------
### Code:


library(hdnom)
library(MASS)
library(data.table)
library(survival)
library(prodlim)
library(zoo)
library(nleqslv)
library(foreach)
library(doParallel)
library(xtable)
library(stringr)
library(glmnet)
library(Matrix)
library(rlist)
library(digest)

source("./R/tmle.estimation.fun.R")
source("./R/sim.data.recurrent.R")
source("./R/lebesgue.loss.fun.R")
source("./R/cv.fun.R")     
source("./R/basis.fun.R")
source("./R/fit.hal.R")
source("./R/predict.hal.R")

## simulation data from setting 2 in the paper: 

sim.dt <- sim.data.outer(n = 200, sim.setting = "7A", cens.percentage = "low", seed = 200)
sim.dt 


# The general version of the TMLE (based on Cox models for initial estimation) is fitted as follows:

cox.tmle <- tmle.est.fun(sim.dt, 
                         fit.type1 = list(model = "Surv(tstart, tstop, delta == 1)~A+L1+L2+L3+Y.dummy", 
                                          fit = "cox"),
                         fit.type2 = list(model = "Surv(tstart, tstop, delta == 2)~A+L1+L2+L3+Y.dummy",
                                          fit = "cox"),
                         fit.type0 = list(model = "Surv(tstart, tstop, delta == 0)~A+L1+L2+L3+Y.dummy",
                                          fit = "cox"),
                         fit.treatment = "A~L1+L2+L3")
print(cox.tmle)

# The TMLE based on HAL for initial estimation as follows:

hal.tmle <- tmle.est.fun(sim.dt, 
                         fit.type1 = list(model = "Surv(tstart, tstop, delta == 1)~A+L1+L2+L3+Y.1", 
                                          fit = "hal"),
                         fit.type2 = list(model = "Surv(tstart, tstop, delta == 2)~A+L1+L2+L3+Y.1",
                                          fit = "hal"),
                         fit.type0 = list(model = "Surv(tstart, tstop, delta == 0)~A+L1+L2+L3+Y.1",
                                          fit = "hal"),
                         cut.time.varying = 8, 
                         fit.treatment = "A~L1+L2+L3")
print(hal.tmle)

######################################################################
### 02_run_estimators.R ends here
