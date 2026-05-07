### summarize.output.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: May  7 2026 (13:50) 
## Version: 
## Last-Updated: May  7 2026 (14:05) 
##           By: Helene
##     Update #: 16
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

summarize.output <- function(sim.setting = "1A",
                           cens.percentage = "low",
                           intervention.A = 1,
                           n = 500, M = 100,
                           use.hal = NA,
                           btmle = NA,
                           misspecify.T = FALSE, misspecify.C = FALSE,
                           baseline.tmle = FALSE,
                           standard.np = FALSE,
                           naive.gcomp = FALSE,
                           browse = FALSE, firstM = NULL,
                           plot.output = FALSE,
                           print.oracle.coverage = TRUE,
                           file.path = "",
                           test.seed.dependence = FALSE) {

    if (!is.na(btmle)) {
        if (btmle) baseline.tmle <- TRUE
    }

    true.psi <- readRDS(file=paste0("./", file.path, "output/",
                                    "save-true-psi-A",
                                    intervention.A,
                                    "-recurrent", 
                                    "-sim.setting.", sim.setting,
                                    ".rds"))

    print(paste0("sim setting ", sim.setting, ": true psi = ",
                 true.psi["E1", "psi0"]))

    mse.fun <- function(x) mean((x - true.psi["E1", "psi0"])^2) #

    cens.fraction <- readRDS(file=paste0("./", file.path, "saved_output/",
                                         "save-cens-fraction",
                                         "-recurrent", 
                                         "-sim.setting.", sim.setting,
                                         "-cens.percentage.", cens.percentage,
                                         ".rds"))

    print(paste0("sim setting ", sim.setting, ": cens fraction = ",
                 cens.fraction, "%"))

    if (test.seed.dependence) {
        est.list <- readRDS(file=paste0("./", file.path, "saved_output/",
                                        "test-seed-dependence-est-list-A",
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
        lapply(est.list, function(est) {
            print(sapply(est, function(seed) seed["tmle.est"]))
            print(sapply(est, function(seed) seed["tmle.se"]))
        })
        return(lapply(est.list, function(est) c(mean = mean(sapply(est, function(seed) seed["tmle.est"])),
                                                sd = sd(sapply(est, function(seed) seed["tmle.est"])),
                                                mean.se = mean(sapply(est, function(seed) seed["tmle.se"])),
                                                sd.se = sd(sapply(est, function(seed) seed["tmle.se"])))))
        
    } else{
        est.list <- readRDS(file=paste0("./", file.path, "saved_output/",
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
    }

    if (length(firstM)>0) {
        est.list <- est.list[1:firstM]
    }

    print(file.info(paste0("./", file.path, "saved_output/",
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
                           ".rds"))$mtime)

    if ("positivity.issues" %in% names(est.list[[1]][[1]])) {
        print(paste0("sim setting ", sim.setting, ": positivity issues = ",
                     round(100*mean(as.numeric(unlist(lapply(est.list, function(x) x[[1]]["positivity.issues"]))), na.rm = TRUE), 4), "%"))
        (pos.jj <- (1:length(est.list))[unlist(lapply(est.list, function(x) x[[1]]["positivity.issues"])) == 1])
    }

    if (browse) browser()

    if (mean(unlist(lapply(est.list, function(x) length(x) == 0))) > 0) {
        print(paste0("output NULL in ", 100*mean(unlist(lapply(est.list, function(x) length(x) == 0))),
                     " % of the repetitions"))
        (1:length(est.list))[unlist(lapply(est.list, function(x) length(x) == 0))]
    }

    if (mean(unlist(lapply(est.list, function(x) is.character(x[[1]])))) > 0) {
        print(paste0("output ERROR in ", 100*mean(unlist(lapply(est.list, function(x) is.character(x[[1]])))),
                     " % of the repetitions"))
        (1:length(est.list))[unlist(lapply(est.list, function(x) is.character(x[[1]])))]        
    }

    if (plot.output) {
        if (length(est.list[[1]]) == 1) {
            #if (length(pos.jj)>0) par(mfrow = c(1,2))
            hist(unlist(lapply(est.list, function(x) x[["tmle.est"]]["tmle.est"])), breaks = 10)
            #if (length(pos.jj)>0)
            #    hist(unlist(lapply(est.list[-pos.jj], function(x) x[["tmle.est"]]["tmle.est"])))
        } else {
            #if (length(pos.jj)>0) par(mfrow = c(1,2))
            hist(unlist(lapply(est.list, function(x) x[[1]]["tmle.est"])), breaks = 10)
            #if (length(pos.jj)>0)
            #    hist(unlist(lapply(est.list[-pos.jj], function(x) x[[1]]["tmle.est"])))            
        }

    }
    
    if (standard.np) {
        if (print.oracle.coverage) {
            tmle.est <- unlist(lapply(est.list, function(x) x$standard.np["np.est"]))
            tmle.se <- unlist(lapply(est.list, function(x) x$standard.np["np.se"]))
            tmle.sd <- sd(unlist(lapply(est.list, function(x) x$standard.np["np.est"])))
            print(paste0("oracle coverage = ",
                         round(mean(tmle.est-1.96*tmle.sd <= true.psi["E1", "psi0"] &
                                    tmle.est+1.96*tmle.sd >= true.psi["E1", "psi0"]), 4)))
        }
        return(results <- data.table(sim.setting = as.character(sim.setting),
                                     which = "",
                                     estimator = c("standard np"),
                                     t(cbind(snp = c(mean = round(mean(unlist(lapply(est.list, function(x) x$standard.np[1]))[!is.na(unlist(lapply(est.list, function(x) x$standard.np["np.se"])))]), 4),
                                                     bias = round(mean(unlist(lapply(est.list, function(x) x$standard.np[1]))[!is.na(unlist(lapply(est.list, function(x) x$standard.np["np.se"])))]) - true.psi["E1", "psi0"], 4),
                                                     mse = round(mse.fun(unlist(lapply(est.list, function(x) x$standard.np[1]))[!is.na(unlist(lapply(est.list, function(x) x$standard.np["np.se"])))]), 4), #"-",                                   
                                                     sd = round(sd(unlist(lapply(est.list, function(x) x$standard.np[1]))[!is.na(unlist(lapply(est.list, function(x) x$standard.np["np.se"])))]), 4), #"-",
                                                     se = round(mean(unlist(lapply(est.list, function(x) x$standard.np["np.se"]))[!is.na(unlist(lapply(est.list, function(x) x$standard.np["np.se"])))]), 4),
                                                     cov = round(mean((unlist(lapply(est.list, function(x) x$standard.np["np.est"]))-1.96*unlist(lapply(est.list, function(x) x$standard.np["np.se"])) <= true.psi["E1", "psi0"] &
                                                                       unlist(lapply(est.list, function(x) x$standard.np["np.est"]))+1.96*unlist(lapply(est.list, function(x) x$standard.np["np.se"])) >= true.psi["E1", "psi0"])[!is.na(unlist(lapply(est.list, function(x) x$standard.np["np.se"])))]), 4),
                                                     oracle.cov = round(mean(tmle.est-1.96*tmle.sd <= true.psi["E1", "psi0"] &
                                                                             tmle.est+1.96*tmle.sd >= true.psi["E1", "psi0"]), 4)))
                                       )))
    } else if (naive.gcomp) {
        return(results <- data.table(sim.setting = as.character(sim.setting),
                                     which = ifelse(misspecify.T & misspecify.C, "both missp",
                                             ifelse(misspecify.T, "missp",
                                             ifelse(misspecify.C, "cens missp", "correct"))),
                                     estimator = c("naive g-comp 1", "naive g-comp 2"),
                                     t(cbind(naive1 = c(mean = round(mean(unlist(lapply(est.list, function(x) x$naive1))), 4),
                                                        bias = round(mean(unlist(lapply(est.list, function(x) x$naive1))) - true.psi["E1", "psi0"], 4),
                                                        mse = round(mse.fun(unlist(lapply(est.list, function(x) x$naive1))), 4), #"-",                                   
                                                        sd = round(sd(unlist(lapply(est.list, function(x) x$naive1))), 4), #"-",
                                                        se = "-",
                                                        cov = "-",
                                                        oracle.cov = "-"),
                                             naive2 = c(mean = round(mean(unlist(lapply(est.list, function(x) x$naive2))), 4),
                                                        bias = round(mean(unlist(lapply(est.list, function(x) x$naive2))) - true.psi["E1", "psi0"], 4),
                                                        mse = round(mse.fun(unlist(lapply(est.list, function(x) x$naive2))), 4), #"-",
                                                        sd = round(sd(unlist(lapply(est.list, function(x) x$naive2))), 4), #"-",
                                                        se = "-",
                                                        cov = "-",
                                                        oracle.cov = "-"))
                                       )))
    } else if (baseline.tmle) {

        if (!is.na(btmle)) {
            if (names(est.list[[1]]) == "tmle.est") {
                for (mm in 1:length(est.list)) {
                    names(est.list[[mm]]) <- "baseline.tmle"
                }
            }
        }
        
        if (print.oracle.coverage) {
            tmle.est <- unlist(lapply(est.list, function(x) x[["baseline.tmle"]]["tmle.est"]))
            tmle.se <- unlist(lapply(est.list, function(x) x[["baseline.tmle"]]["tmle.se"]))
            tmle.sd <- sd(unlist(lapply(est.list, function(x) x[["baseline.tmle"]]["tmle.est"])))
            print(paste0("oracle coverage = ",
                         round(mean(tmle.est-1.96*tmle.sd <= true.psi["E1", "psi0"] &
                                    tmle.est+1.96*tmle.sd >= true.psi["E1", "psi0"]), 4)))
        }

        return(results <- data.table(sim.setting = as.character(sim.setting),
                                     which = ifelse(misspecify.T & misspecify.C, "both missp",
                                             ifelse(misspecify.T, "missp",
                                             ifelse(misspecify.C, "cens missp", "correct"))),
                                     estimator = c(ifelse(!is.na(use.hal), "hal-baseline-tmle", "baseline-tmle")),
                                     t(cbind(tmle.est = c(mean = round(mean(unlist(lapply(est.list, function(x) x[["baseline.tmle"]]["tmle.est"]))), 4),
                                                          bias = round(mean(unlist(lapply(est.list, function(x) x[["baseline.tmle"]]["tmle.est"]))) - true.psi["E1", "psi0"], 4),
                                                          mse = round(mse.fun(unlist(lapply(est.list, function(x) x[["baseline.tmle"]]["tmle.est"]))), 4),
                                                          sd = round(sd(unlist(lapply(est.list, function(x) x[["baseline.tmle"]]["tmle.est"]))), 4),
                                                          se = round(mean(unlist(lapply(est.list, function(x) x[["baseline.tmle"]]["tmle.se"]))), 4),
                                                          cov = round(mean(unlist(lapply(est.list, function(x) x[["baseline.tmle"]]["tmle.est"]))-1.96*unlist(lapply(est.list, function(x) x[["baseline.tmle"]]["tmle.se"])) <= true.psi["E1", "psi0"] &
                                                                           unlist(lapply(est.list, function(x) x[["baseline.tmle"]]["tmle.est"]))+1.96*unlist(lapply(est.list, function(x) x[["baseline.tmle"]]["tmle.se"])) >= true.psi["E1", "psi0"]), 4),
                                                          oracle.cov = round(mean(tmle.est-1.96*tmle.sd <= true.psi["E1", "psi0"] &
                                                                                  tmle.est+1.96*tmle.sd >= true.psi["E1", "psi0"]), 4))
                                             ))))
    } else {
        
        if (print.oracle.coverage) {
            tmle.est <- unlist(lapply(est.list, function(x) x[["tmle.est"]]["tmle.est"]))
            tmle.se <- unlist(lapply(est.list, function(x) x[["tmle.est"]]["tmle.se"]))
            tmle.sd <- sd(unlist(lapply(est.list, function(x) x[["tmle.est"]]["tmle.est"])))
            print(paste0("oracle coverage = ",
                         round(mean(tmle.est-1.96*tmle.sd <= true.psi["E1", "psi0"] &
                                    tmle.est+1.96*tmle.sd >= true.psi["E1", "psi0"]), 4)))
        }

        if (length(est.list[[1]]) == 1) {
            return(results <- data.table(sim.setting = as.character(sim.setting),
                                         which = ifelse(!is.na(use.hal), "", ifelse(misspecify.T & misspecify.C, "both missp",
                                                                             ifelse(misspecify.T, "missp",
                                                                             ifelse(misspecify.C, "cens missp", "correct")))),
                                         estimator = c(ifelse(!is.na(use.hal), c("hal-tmle"), c("tmle")),
                                                       ifelse(!is.na(use.hal), c("hal g-comp"), c("g-comp"))),
                                         t(cbind(tmle.est = c(mean = round(mean(as.numeric(unlist(lapply(est.list, function(x) x[[1]]["tmle.est"]))), na.rm = TRUE), 4),
                                                              bias = round(mean(as.numeric(unlist(lapply(est.list, function(x) x[[1]]["tmle.est"]))), na.rm = TRUE) - true.psi["E1", "psi0"], 4),
                                                              mse = round(mse.fun(na.omit(as.numeric(unlist(lapply(est.list, function(x) x[[1]]["tmle.est"]))))), 4),
                                                              sd = round(sd(as.numeric(unlist(lapply(est.list, function(x) x[[1]]["tmle.est"]))), na.rm = TRUE), 4),
                                                              se = round(mean(as.numeric(unlist(lapply(est.list, function(x) x[[1]]["tmle.se"]))), na.rm = TRUE), 4),
                                                              cov = round(mean(as.numeric(unlist(lapply(est.list, function(x) x[[1]]["tmle.est"])))-1.96*as.numeric(unlist(lapply(est.list, function(x) x[[1]]["tmle.se"]))) <= true.psi["E1", "psi0"] &
                                                                               as.numeric(unlist(lapply(est.list, function(x) x[[1]]["tmle.est"])))+1.96*as.numeric(unlist(lapply(est.list, function(x) x[[1]]["tmle.se"]))) >= true.psi["E1", "psi0"], na.rm = TRUE), 4),
                                                              oracle.cov = round(mean(tmle.est-1.96*tmle.sd <= true.psi["E1", "psi0"] &
                                                                                      tmle.est+1.96*tmle.sd >= true.psi["E1", "psi0"]), 4)),
                                                 g.est = c(mean = round(mean(as.numeric(unlist(lapply(est.list, function(x) x[[1]]["g.est"]))), na.rm = TRUE), 4),
                                                           bias = round(mean(as.numeric(unlist(lapply(est.list, function(x) x[[1]]["g.est"]))), na.rm = TRUE) - true.psi["E1", "psi0"], 4),
                                                           mse = round(mse.fun(na.omit(as.numeric(unlist(lapply(est.list, function(x) x[[1]]["g.est"]))))), 4) ,#"-",
                                                           sd = round(sd(as.numeric(unlist(lapply(est.list, function(x) x[[1]]["g.est"]))), na.rm = TRUE), 4) ,#"-",
                                                           se = "-",
                                                           cov = "-",
                                                           oracle.cov = "-")
                                                 ))))
        } else {
            return(results <- data.table(sim.setting = as.character(sim.setting),
                                         which = ifelse(!is.na(use.hal), "", ifelse(misspecify.T & misspecify.C, "both missp",
                                                                             ifelse(misspecify.T, "missp",
                                                                             ifelse(misspecify.C, "cens missp", "correct")))),
                                         estimator = c(ifelse(!is.na(use.hal), c("hal-tmle"), c("tmle")),
                                                       ifelse(!is.na(use.hal), c("hal g-comp"), c("g-comp"))),
                                         t(cbind(tmle.est = c(mean = round(mean(unlist(lapply(est.list, function(x) x[["tmle.est"]]["tmle.est"]))), 4),
                                                              bias = round(mean(unlist(lapply(est.list, function(x) x[["tmle.est"]]["tmle.est"]))) - true.psi["E1", "psi0"], 4),
                                                              mse = round(mse.fun(unlist(lapply(est.list, function(x) x[["tmle.est"]]["tmle.est"]))), 4),
                                                              sd = round(sd(unlist(lapply(est.list, function(x) x[["tmle.est"]]["tmle.est"]))), 4),
                                                              se = round(mean(unlist(lapply(est.list, function(x) x[["tmle.est"]]["tmle.se"]))), 4),
                                                              cov = round(mean(unlist(lapply(est.list, function(x) x[["tmle.est"]]["tmle.est"]))-1.96*unlist(lapply(est.list, function(x) x[["tmle.est"]]["tmle.se"])) <= true.psi["E1", "psi0"] &
                                                                               unlist(lapply(est.list, function(x) x[["tmle.est"]]["tmle.est"]))+1.96*unlist(lapply(est.list, function(x) x[["tmle.est"]]["tmle.se"])) >= true.psi["E1", "psi0"]), 4),
                                                              oracle.cov = round(mean(tmle.est-1.96*tmle.sd <= true.psi["E1", "psi0"] &
                                                                                      tmle.est+1.96*tmle.sd >= true.psi["E1", "psi0"]), 4)),
                                                 g.est = c(mean = round(mean(unlist(lapply(est.list, function(x) x[["tmle.est"]]["g.est"]))), 4),
                                                           bias = round(mean(unlist(lapply(est.list, function(x) x[["tmle.est"]]["g.est"]))) - true.psi["E1", "psi0"], 4),
                                                           mse = round(mse.fun(unlist(lapply(est.list, function(x) x[["tmle.est"]]["g.est"]))), 4) ,#"-",
                                                           sd = round(sd(unlist(lapply(est.list, function(x) x[["tmle.est"]]["g.est"]))), 4) ,#"-",
                                                           se = "-",
                                                           cov = "-",
                                                           oracle.cov = "-")
                                                 ))))
        }
    }    
}

######################################################################
### summarize.output.R ends here
