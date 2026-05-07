### 01_simulate_data.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: May  7 2026 (13:47) 
## Version: 
## Last-Updated: May  7 2026 (14:07) 
##           By: Helene
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
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

## setting 1 in the paper: 

sim.dt <- sim.data.outer(n = 200, sim.setting = "4A", cens.percentage = "low", seed = 200)
sim.dt 

## setting 2 in the paper: 

sim.dt <- sim.data.outer(n = 200, sim.setting = "7A", cens.percentage = "low", seed = 200)
sim.dt 

######################################################################
### 01_simulate_data.R ends here
