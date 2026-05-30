### 03_real_data_analysis.R --- 
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

# Load and prepare data: 

library(frailtypack)
data(readmission)

dt.readm <- data.table(readmission)

dt.readm[event == 1 & death == 1, death := 0]

dt.readm[, time := t.stop]
dt.readm[, delta := event]
dt.readm[death == 1, delta := 2]

dt.readm[, max.time := max(time), by = "id"]

#-- treatment variable:
dt.readm[, chemo := 1*(chemo == "Treated")]

#-- use only baseline value of charlson:
dt.readm[, charlson := charlson[1], by = "id"]

dt.readm[, table(event, death)]

## NB: this is not the same call that was used to get the result in the main paper;
## for example, in the main paper, CV was repeated several times to robustify against seed dependence;
## here we don't do that, so that the code will run faster.

set.seed(202002)

est.hal.1 <- 
  tmle.est.fun(dt.readm, tau = 5*365, 
               fit.type1 = list(model = "Surv(tstart, tstop, delta == 1)~chemo+sex+dukes+charlson+Y.1",
                                fit = "hal"),
               fit.type2 = list(model = "Surv(tstart, tstop, delta == 2)~chemo+sex+dukes+charlson+Y.1",
                                fit = "hal"),
               fit.type0 = list(model = "Surv(tstart, tstop, delta == 0)~chemo+sex+dukes+charlson+Y.1",
                                fit = "hal"),
               fit.treatment = "chemo~sex+dukes+charlson",
               cut.time = 35,
               cut.time.varying = 8,
               cut.two.way = 1,
               min.no.of.ones = 0.05,
               event.dependent.cv = TRUE,
               reduce.seed.dependence = FALSE,
               V = 10,
               return.eic = TRUE,
               lambda.cvs = c((9:1)/10,(9:2)/10^2, seq(1/10^2, 1/10^5, length = 100)))

est.hal.0 <- 
  tmle.est.fun(dt.readm, tau = 5*365,
               intervention.A = 0,                           
               fit.type1 = list(model = "Surv(tstart, tstop, delta == 1)~chemo+sex+dukes+charlson+Y.1",
                                fit = "hal"),
               fit.type2 = list(model = "Surv(tstart, tstop, delta == 2)~chemo+sex+dukes+charlson+Y.1",
                                fit = "hal"),
               fit.type0 = list(model = "Surv(tstart, tstop, delta == 0)~chemo+sex+dukes+charlson+Y.1",
                                fit = "hal"),
               fit.treatment = "chemo~sex+dukes+charlson",
               cut.time = 35,
               cut.time.varying = 8,
               cut.two.way = 1,
               min.no.of.ones = 0.05,
               event.dependent.cv = TRUE,
               reduce.seed.dependence = FALSE,
               V = 10,
               return.eic = TRUE,
               lambda.cvs = c((9:1)/10,(9:2)/10^2, seq(1/10^2, 1/10^5, length = 100)))

est.hal.1[[1]]["tmle.est"]
est.hal.1[[1]]["tmle.se"]

est.hal.0[[1]]["tmle.est"]
est.hal.0[[1]]["tmle.se"]

est.hal.1[[1]]["tmle.est"] + 1.96*c(-1,1)*est.hal.1[[1]]["tmle.se"]
est.hal.0[[1]]["tmle.est"] + 1.96*c(-1,1)*est.hal.0[[1]]["tmle.se"]

(est.diff <- est.hal.1[[1]][["tmle.est"]] - est.hal.0[[1]][["tmle.est"]])
(se.diff.eic <- sqrt(mean((est.hal.1$eic - est.hal.0$eic)^2)/length(unique(dt.readm[["id"]]))))
est.diff + c(-1,1)*1.96*se.diff.eic 

est.diff*length(unique(dt.readm[["id"]]))
(est.diff + c(-1,1)*1.96*se.diff.eic)*length(unique(dt.readm[["id"]]))

######################################################################
### 03_real_data_analysis.R ends here
