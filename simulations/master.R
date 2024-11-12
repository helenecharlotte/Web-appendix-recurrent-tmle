### master.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Nov 12 2024 (15:29) 
## Version: 
## Last-Updated: Nov 12 2024 (15:31) 
##           By: Helene
##     Update #: 2
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
library(SuperLearner)
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

#-------------------------------------------------------------------------------------------#
## source relevant scripts
#-------------------------------------------------------------------------------------------#

source("./R/tmle.estimation.fun.R")
source("./R/sim.data.recurrent.R")
source("./R/lebesgue.loss.fun.R")
source("./R/cv.fun.R")     
source("./R/basis.fun.R")
source("./R/fit.hal.R")
source("./R/predict.hal.R")

#-------------------------------------------------------------------------------------------#
## specify parameters for simulation setting
#-------------------------------------------------------------------------------------------#

n <- 500
M <- 250

verbose <- TRUE

get.truth <- FALSE
get.uncensored.truth <- FALSE
get.cens.fraction <- FALSE

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){
    no.cores <- 30
} else {
    no.cores <- 3
}



if (get.truth) {
    for (intervention.A in c(1)) {
        for (sim.setting in c("1A")) {
            true.psi <- sim.data.outer(sim.setting = sim.setting, n = 1e4, intervention.A = intervention.A, rep.true = 50)
            saveRDS(true.psi,
                    file=paste0("./output/",
                                "save-true-psi-A",
                                intervention.A,
                                "-recurrent", 
                                "-sim.setting.", sim.setting,
                                ".rds"))
        }
    }
}

if (get.uncensored.truth) {
    for (sim.setting in c("1A")) {
        true.psi <- sim.data.outer(sim.setting = sim.setting, n = 1e4, censoring = FALSE, rep.true = 50)
        saveRDS(true.psi,
                file=paste0("./output/",
                            "save-true-psi-no-cens",
                            "-recurrent", 
                            "-sim.setting.", sim.setting,
                            ".rds"))
    }
}

if (get.cens.fraction) {
    for (sim.setting in c("1A")) {
        for (cens.percentage in c("low", "high")) {
            print(paste0("sim.setting ", sim.setting, ": ",
                         cens.fraction <- sim.data.outer(sim.setting = sim.setting,
                                                         cens.percentage = cens.percentage,
                                                         n = 1e5,
                                                         get.cens.fraction = TRUE)*100, "%"))
            saveRDS(cens.fraction,
                    file=paste0("./output/",
                                "save-cens-fraction",
                                "-recurrent", 
                                "-sim.setting.", sim.setting,
                                "-cens.percentage.", cens.percentage,
                                ".rds"))
        }
    }
}


    for (sim.setting in c("1Ax")) {

        #--- amount of censoring ---#
        for (cens.percentage in c("low")) { #, "high"

            print("-----------------------------------------------------------------------------------")
            print("-----------------------------------------------------------------------------------")
            print(paste0("sim.setting = ", sim.setting))
            print(paste0("cens.percentage = ", cens.percentage))
            print("-----------------------------------------------------------------------------------")
        
            for (misspecify.list in list(
                                        c(misspecify.T = FALSE, misspecify.C = FALSE)
                                    )) {

                misspecify.T <- misspecify.list["misspecify.T"][[1]]
                misspecify.C <- misspecify.list["misspecify.C"][[1]]
                use.hal <- misspecify.list["use.hal"][[1]]

                registerDoParallel(no.cores)

                est.list <- foreach(m=1:M, .errorhandling="pass"
                                    ) %dopar% {

                                        print(paste0("m = ", m))
                                        set.seed(100+m)
                                        dt <- sim.data.outer(n = n, sim.setting = sim.setting)
                                        print(summary(dt))

                                        #--- use HAL ---#
                                        if (!is.na(use.hal)) {

                                            return(list(
                                                tmle.est = tmle.est.fun(dt, intervention.A = intervention.A,
                                                                        verbose = verbose,
                                                                        cv.glmnet = FALSE,
                                                                        fit.type1 = list(model = "Surv(tstart, tstop, delta == 1)~A+L2+L1+L3+Y.dummy+Y.1",
                                                                                         fit = "hal"),
                                                                        fit.type2 = list(model = "Surv(tstart, tstop, delta == 2)~A+L2+L1+L3+Y.dummy",
                                                                                         fit = "hal"),
                                                                        fit.type0 = list(model = "Surv(tstart, tstop, delta == 0)~A+L2+L1+L3+Y.dummy",
                                                                                         fit = "hal"),
                                                                        cut.time.varying = 4)
                                            ))
                                        }
                                    
                                        #--- use Cox ---#

                                        if (substr(sim.setting, 1, 1) == "1") {
                                            fit.type1.model <- "Surv(tstart, tstop, delta == 1)~A+L2+L1.squared+L3+Y.dummy"
                                            fit.type2.model <- "Surv(tstart, tstop, delta == 2)~A+L2+L1+L3+Y.dummy"
                                            fit.type0.model <- "Surv(tstart, tstop, delta == 0)~A+L2+L1+L3+Y.dummy"
                                        } else if (substr(sim.setting, 1, 1) == "2") {
                                            fit.type1.model <- "Surv(tstart, tstop, delta == 1)~A+L2+L1.squared+L3+Y.dummy"
                                            fit.type2.model <- "Surv(tstart, tstop, delta == 2)~A+L2+L1+L3+Y.dummy"
                                            fit.type0.model <- "Surv(tstart, tstop, delta == 0)~A+L2+L1+L3"
                                        } else if (substr(sim.setting, 1, 1) == "3") {
                                            fit.type1.model <- "Surv(tstart, tstop, delta == 1)~A+L2+L1.squared+L3"
                                            fit.type2.model <- "Surv(tstart, tstop, delta == 2)~A+L2+L1+L3"
                                            fit.type0.model <- "Surv(tstart, tstop, delta == 0)~A+L2+L1+L3"
                                        }

                                        if (misspecify.T) {
                                            if (!substr(sim.setting, 2, 2) == "A") {
                                                next
                                            }
                                            fit.type1.model <- gsub("L1.squared", "L1", fit.type1.model)
                                        } else if (misspecify.C) {
                                            if (!substr(sim.setting, 1, 1) == "1") {
                                                next
                                            } else {
                                                fit.type0.model <- gsub("\\+Y.dummy", "", fit.type0.model)
                                            }
                                        }
                                    
                                        return(list(
                                            tmle.est = tmle.est.fun(dt, intervention.A = intervention.A,
                                                                    verbose = verbose,
                                                                    fit.type1 = fit.type1.model,
                                                                    fit.type2 = fit.type2.model,
                                                                    fit.type0 = fit.type0.model),
                                            baseline.tmle = tmle.est.fun(dt, intervention.A = intervention.A,
                                                                         verbose = verbose,
                                                                         baseline.tmle = TRUE,
                                                                         fit.type1 = gsub("\\+Y.dummy", "", fit.type1.model),
                                                                         fit.type2 = gsub("\\+Y.dummy", "", fit.type2.model),
                                                                         fit.type0 = gsub("\\+Y.dummy", "", fit.type0.model)),
                                            standard.np = tmle.est.fun(dt, intervention.A = intervention.A,
                                                                       verbose = verbose,
                                                                       standard.np = TRUE),
                                            naive1 = tmle.est.fun(dt, intervention.A = intervention.A,
                                                                  verbose = verbose,
                                                                  naive.gcomp = TRUE,
                                                                  fit.type1 = gsub("\\+Y.dummy", "", fit.type1.model)),
                                            naive2 = tmle.est.fun(dt, intervention.A = intervention.A,
                                                                  verbose = verbose,
                                                                  naive.gcomp = TRUE,
                                                                  fit.type1 = paste0(gsub("\\+Y.dummy", "", fit.type1.model), "+Y.dummy"))
                                        ))
                                    }


                stopImplicitCluster()

                saveRDS(est.list,
                        file=paste0("./output/",
                                    "save-est-list-A",
                                    intervention.A,
                                    "-recurrent",
                                    "-n", n, "-M", M,
                                    "-sim.setting.", sim.setting,
                                    "-cens.percentage.", cens.percentage,
                                    ifelse(!is.na(use.hal), "-useHAL", ""),
                                    ifelse(misspecify.T, "-misspecify.T", ""),
                                    ifelse(misspecify.C, "-misspecify.C", ""),
                                    ".rds"))

                if (M <= 5) print(est.list)
                
            }
        }
    }
}



######################################################################
### master.R ends here
