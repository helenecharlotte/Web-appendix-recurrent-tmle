### master.R --- 
#----------------------------------------------------------------------

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
library(parallel)
library(xtable)
library(stringr)
library(glmnet)
library(Matrix)
library(digest)
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
## get true values, amount of censoring, etc
#-------------------------------------------------------------------------------------------#


for (intervention.A in c(1,0)) {
    for (sim.setting in c("8A")) { #c("1A", "1B", "2A", "2B", "3A")
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

for (sim.setting in c("8A")) { #c("1A", "1B", "2A", "2B", "3A")
    for (cens.percentage in c("low", "high", "10%", "30%")) {
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

for (sim.setting in c("8A")) { #c("1A", "1B", "2A", "2B", "3A")
    true.psi <- sim.data.outer(sim.setting = sim.setting, n = 1e4, censoring = FALSE, rep.true = 50)
    saveRDS(true.psi,
            file=paste0("./output/",
                        "save-true-psi-no-cens",
                        "-recurrent", 
                        "-sim.setting.", sim.setting,
                        ".rds"))
}

#-------------------------------------------------------------------------------------------#
## specify parameters for simulation setting
#-------------------------------------------------------------------------------------------#

n <- 500
M <- 1000

verbose <- TRUE

#-------------------------------------------------------------------------------------------#
## loop over simulation repetitions
#-------------------------------------------------------------------------------------------#

for (intervention.A in c(1)) {

    #--- simulation setting ---#
    for (sim.list in list(
                         list(sim.setting = "4A", ## <-- setting 1
                              misspecify.list = list(
                                  c(misspecify.T = FALSE, misspecify.C = FALSE),
                                  c(misspecify.T = TRUE, misspecify.C = FALSE),
                                  c(misspecify.T = 2, misspecify.C = 2),
                                  c(use.hal = "no-interaction-3")
                              )),
                         list(sim.setting = "7A", ## <-- setting 2
                              misspecify.list = list(
                                  c(misspecify.T = FALSE, misspecify.C = FALSE),
                                  c(misspecify.T = 2, misspecify.C = TRUE),
                                  c(misspecify.T = FALSE, misspecify.C = TRUE),
                                  c(use.hal = "no-interaction-8")
                              ))
                     )) {
 
        sim.setting <- sim.list[["sim.setting"]]
        misspecify.list <- sim.list[["misspecify.list"]]

        #--- amount of censoring ---# 
        for (cens.percentage in c("low", "high", "10%", "30%")) { #

            print("-----------------------------------------------------------------------------------")
            print("-----------------------------------------------------------------------------------")
            print(paste0("sim.setting = ", sim.setting))
            print(paste0("cens.percentage = ", cens.percentage))
            print("-----------------------------------------------------------------------------------")
        
            for (misspecify in misspecify.list) {

                misspecify.T <- misspecify["misspecify.T"][[1]]
                misspecify.C <- misspecify["misspecify.C"][[1]]
                use.hal <- misspecify["use.hal"][[1]]
                btmle <- misspecify["btmle"][[1]]

                if (!is.na(use.hal)) {
                    misspecify.T <- misspecify.C <- FALSE

                    if (length(grep("no-interaction", use.hal)) > 0) {
                        cut.two.way <- 0
                        cut.time.treatment <- 0
                        cut.time.covar <- 0
                        cut.time.time.varying.covar <- 0
                    } else {
                        cut.two.way <- 10
                        cut.time.treatment <- 10
                        cut.time.covar <- 10
                        cut.time.time.varying.covar <- 10
                    }
                        
                    if ((num <- gregexpr('[0-9]+',use.hal)[[1]][1])>0) {
                        cut.time.varying <- as.numeric(substr(use.hal, num, nchar(use.hal)))                                          
                    } else {
                        cut.time.varying <- 6
                    }
                    
                    cut.time <- misspecify["cut.time"][[1]]
                    if (is.na(cut.time)) {
                        cut.time <- 35
                    } else {
                        use.hal <- paste0(use.hal, "-cut-time-", cut.time)
                        cut.time <- as.numeric(cut.time)
                    }                                                

                    cut.one.way <- misspecify["cut.one.way"][[1]]
                    if (is.na(cut.one.way)) {
                        cut.one.way <- 35
                    } else {
                        use.hal <- paste0(use.hal, "-cut-one-way-", cut.one.way)
                        cut.one.way <- as.numeric(cut.one.way)
                    }

                    use.exponential <- misspecify["use.exponential"][[1]]
                    if (is.na(use.exponential)) {
                        use.exponential <- FALSE
                    } else {
                        use.exponential <- (use.exponential == "TRUE")
                        use.hal <- paste0(use.hal, "-use-exponential-", use.exponential)
                    }

                    cv.hal.fit <- misspecify["cv.hal.fit"][[1]]
                    if (is.na(cv.hal.fit)) {
                        cv.hal.fit <- FALSE
                    } else {
                        cv.hal.fit <- (cv.hal.fit == "TRUE")
                        use.hal <- paste0(use.hal, "-cv.hal.fit-", cv.hal.fit)
                    }

                    
                    V <- misspecify["V"][[1]]
                    if (is.na(V)) {
                        V <- 10
                    } else {
                        V <- as.numeric(V)
                        use.hal <- paste0(use.hal, "-V", V)
                    }

                    reduce.seed.dependence <- misspecify["reduce.seed.dependence"][[1]]
                    if (is.na(reduce.seed.dependence)) {
                        reduce.seed.dependence <- FALSE
                    } else {
                        reduce.seed.dependence <- as.numeric(reduce.seed.dependence)
                        use.hal <- paste0(use.hal, "-reduce.seed.dependence", reduce.seed.dependence)
                    }
                    
                    event.dependent.cv <- misspecify["event.dependent.cv"][[1]]
                    if (is.na(event.dependent.cv)) {
                        event.dependent.cv <- FALSE
                    } else {
                        event.dependent.cv <- (event.dependent.cv == "TRUE")
                        use.hal <- paste0(use.hal, "-event.dependent.cv", event.dependent.cv)
                    }
                }

                parallelize.Z <- ifelse(detectCores()>15, 5, 3)
                no.cores <- floor(min(detectCores()/parallelize.Z, ifelse(n <= 500, 22,
                                                                          ifelse(!is.na(use.hal), 10, 15))))

                track_m <- function(m) paste0("./output/",
                                              paste0("m=",m,
                                                     ", sim.setting=", sim.setting,
                                                     ", cens.percentage=", cens.percentage,
                                                     ", use.hal=", use.hal,
                                                     ", misspecify.T=", misspecify.T, 
                                                     ", misspecify.C=", misspecify.C, 
                                                     ", btmle=", btmle), ".txt")

                tw <- ceiling(M / 25)

                registerDoParallel(no.cores)

                est.list <- foreach(m=1:M, .errorhandling="pass"
                                    ) %dopar%  {

                                        print(paste0("m = ", m))
                                        set.seed(100+m)
                                        dt <- sim.data.outer(n = n, sim.setting = sim.setting, cens.percentage = cens.percentage)
                                        print(summary(dt))

                                        if (M >= tw & m %% floor(M / tw) == 0) {
                                            write.csv(m, file=track_m(m))
                                            if (file.exists(track_m(m-floor(M / tw)))) {
                                                file.remove(track_m(m-floor(M / tw)))
                                            }
                                        }

                                        #--- use HAL ---#
                                        if (!is.na(use.hal)) {

                                            if (!is.na(btmle)) {
                                                if (btmle) {
                                                    return(list(
                                                        baseline.tmle = tmle.est.fun(dt, intervention.A = intervention.A,
                                                                                     verbose = verbose,
                                                                                     cv.glmnet = FALSE,
                                                                                     fit.type1 = list(model = "Surv(tstart, tstop, delta == 1)~A+L2+L1+L3",
                                                                                                      fit = "hal"),
                                                                                     fit.type2 = list(model = "Surv(tstart, tstop, delta == 2)~A+L2+L1+L3",
                                                                                                      fit = "hal"),
                                                                                     fit.type0 = list(model = "Surv(tstart, tstop, delta == 0)~A+L2+L1+L3",
                                                                                                      fit = "cox"),
                                                                                     cut.time.varying = 4,
                                                                                     parallelize.Z = parallelize.Z,
                                                                                     baseline.tmle = TRUE,
                                                                                     cut.one.way = 30,
                                                                                     #output.lambda.cvs = TRUE,
                                                                                     cut.two.way = 0,
                                                                                     cut.time = 35)
                                                    ))
                                                }
                                            } else {

                                                return(list(
                                                    tmle.est = tmle.est.fun(dt, intervention.A = intervention.A,
                                                                            verbose = verbose,
                                                                            cv.glmnet = FALSE,
                                                                            fit.type1 = list(model = "Surv(tstart, tstop, delta == 1)~A+L2+L1+L3+Y.1",
                                                                                             fit = "hal"),
                                                                            fit.type2 = list(model = "Surv(tstart, tstop, delta == 2)~A+L2+L1+L3+Y.1",
                                                                                             fit = "hal"),
                                                                            fit.type0 = list(model = "Surv(tstart, tstop, delta == 0)~A+L2+L1+L3+Y.1",
                                                                                             fit = "hal"),
                                                                            parallelize.Z = parallelize.Z,
                                                                            use.exponential = use.exponential,
                                                                            cv.hal.fit = cv.hal.fit,
                                                                            V = V,
                                                                            event.dependent.cv = event.dependent.cv,
                                                                            reduce.seed.dependence = reduce.seed.dependence,
                                                                            cut.time.varying = cut.time.varying,
                                                                            cut.one.way = cut.one.way,
                                                                            cut.two.way = cut.two.way,
                                                                            cut.time = cut.time)
                                                ))
                                            }
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
                                        } else if (substr(sim.setting, 1, 1) == "5") {
                                            fit.type1.model <- "Surv(tstart, tstop, delta == 1)~A+L2+L1.squared+L3+Y.dummy+Y.dummy:L2"
                                            fit.type2.model <- "Surv(tstart, tstop, delta == 2)~A+L2+L1+L3+Y.dummy"
                                            fit.type0.model <- "Surv(tstart, tstop, delta == 0)~A+L2+L1+L3+Y.dummy+Y.dummy:L2"
                                        } else if (substr(sim.setting, 1, 1) %in% c("6", "7", "8")) {
                                            fit.type1.model <- "Surv(tstart, tstop, delta == 1)~A+L2+L1.squared+L3+Y.dummy"
                                            fit.type2.model <- "Surv(tstart, tstop, delta == 2)~A+L2+L1+L3+Y.dummy"
                                            if (((substr(sim.setting, 2, 2) == "B" & substr(sim.setting, 1, 1) %in% c("6")) | (substr(sim.setting, 2, 2) == "A" & substr(sim.setting, 1, 1) %in% c("7", "8"))) & !misspecify.C) {
                                                fit.type0.model <- "Surv(tstart, tstop, delta == 0)~A+L2+L1+L3+Y.dummy"
                                            } else {
                                                fit.type0.model <- "Surv(tstart, tstop, delta == 0)~A+L2+L1+L3"
                                            }
                                        } else if (substr(sim.setting, 1, 1) == "3") {
                                            fit.type1.model <- "Surv(tstart, tstop, delta == 1)~A+L2+L1.squared+L3"
                                            fit.type2.model <- "Surv(tstart, tstop, delta == 2)~A+L2+L1+L3"
                                            fit.type0.model <- "Surv(tstart, tstop, delta == 0)~A+L2+L1+L3"
                                        } else if (substr(sim.setting, 1, 1) == "4") {
                                            fit.type1.model <- "Surv(tstart, tstop, delta == 1)~A+L2+L1.squared+L3+Y.dummy+Y.dummy.3"
                                            fit.type2.model <- "Surv(tstart, tstop, delta == 2)~A+L2+L1+L3+Y.dummy+Y.dummy.3"
                                            fit.type0.model <- "Surv(tstart, tstop, delta == 0)~A+L2+L1+L3+Y.dummy+Y.dummy.3"
                                        }

                                        if (misspecify.T == 2) {
                                            if (substr(sim.setting, 1, 1) == "4") {
                                                fit.type1.model <- gsub("\\+Y.dummy.3", "", fit.type1.model)
                                            } else if (substr(sim.setting, 1, 1) == "5") {
                                                fit.type1.model <- gsub("\\+Y.dummy:L2", "", fit.type1.model)
                                            } else {
                                                fit.type1.model <- gsub("\\+Y.dummy", "", fit.type1.model)
                                            }
                                        } else if (misspecify.T) {
                                            if (!substr(sim.setting, 2, 2) == "A" & substr(sim.setting, 1, 1) %in% c("1", "2", "3")) {
                                                next
                                            }
                                            fit.type1.model <- gsub("L1.squared", "L1", fit.type1.model)
                                        }

                                        if (misspecify.C == 2) {
                                            if (substr(sim.setting, 1, 1) == "5") {
                                                fit.type0.model <- gsub("\\+Y.dummy:L2", "", fit.type0.model)
                                            } else {
                                                fit.type0.model <- gsub("\\+Y.dummy.3", "", fit.type0.model)
                                            }
                                        } else if (misspecify.C) {
                                            if (substr(sim.setting, 1, 1) %in% c("2", "3")) {
                                                next
                                            } else {
                                                fit.type0.model <- gsub("\\+Y.dummy", "", fit.type0.model)
                                            }
                                        }

                                        message("---------------------------------")
                                        print(fit.type0.model)
                                        message("---------------------------------")
                                        print(fit.type1.model)
                                        message("---------------------------------")
                                        print(fit.type2.model)
                                        message("---------------------------------")

                                        return(list(
                                            tmle.est = tmle.est.fun(dt, intervention.A = intervention.A,
                                                                    verbose = verbose,
                                                                    fit.type1 = list(model = fit.type1.model, fit = "cox"),
                                                                    fit.type2 = list(model = fit.type2.model, fit = "cox"),
                                                                    fit.type0 = list(model = fit.type0.model, fit = "cox"),
                                                                    parallelize.Z = parallelize.Z),
                                            baseline.tmle = tmle.est.fun(dt, intervention.A = intervention.A,
                                                                         verbose = verbose,
                                                                         baseline.tmle = TRUE,
                                                                         fit.type1 = list(model = gsub("\\+Y.dummy", "", gsub("\\+Y.dummy:L2", "", gsub("\\+Y.dummy.3", "", fit.type1.model))), fit = "cox"),
                                                                         fit.type2 = list(model = gsub("\\+Y.dummy", "", gsub("\\+Y.dummy:L2", "", gsub("\\+Y.dummy.3", "", fit.type2.model))), fit = "cox"),
                                                                         fit.type0 = list(model = gsub("\\+Y.dummy", "", gsub("\\+Y.dummy:L2", "", gsub("\\+Y.dummy.3", "", fit.type0.model))), fit = "cox")),
                                            standard.np = tmle.est.fun(dt, intervention.A = intervention.A,
                                                                       verbose = verbose,
                                                                       standard.np = TRUE),
                                            naive1 = tmle.est.fun(dt, intervention.A = intervention.A,
                                                                  verbose = verbose,
                                                                  naive.gcomp = TRUE,
                                                                  fit.type1 = list(model = gsub("\\+Y.dummy", "", gsub("\\+Y.dummy:L2", "", gsub("\\+Y.dummy.3", "", fit.type1.model))), fit = "cox")),
                                            naive2 = tmle.est.fun(dt, intervention.A = intervention.A,
                                                                  verbose = verbose,
                                                                  naive.gcomp = TRUE,
                                                                  fit.type1 = list(model = gsub("\\+Y.dummy", "", paste0(gsub("\\+Y.dummy:L2", "", gsub("\\+Y.dummy.3", "", fit.type1.model))), "+Y.dummy"), fit = "cox"))
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
                                    ifelse(!is.na(use.hal), ifelse(use.hal >= 2, paste0("-useHAL", use.hal), "-useHAL"), ""),
                                    ifelse(!is.na(btmle), ifelse(btmle, "-btmle", ""), ""),
                                    ifelse(misspecify.T, ifelse(misspecify.T == 2, "-misspecify2.T", "-misspecify.T"), ""),
                                    ifelse(misspecify.C, ifelse(misspecify.C == 2, "-misspecify2.C", "-misspecify.C"), ""),
                                    ".rds"))

                if (file.exists(track_m(M))) {
                    file.remove(track_m(M))
                }
                

                if (M <= 5) print(est.list)
                
            }
        }
    }
}





######################################################################
### master.R ends here
