### predict.hal.R --- 
#----------------------------------------------------------------------


predict.hal <- function(seed = NULL, #13349,
                        fit.hals, ## output from fit.hal
                        pseudo.dt, ## tmp.hal,
                        delta.var = "delta",
                        time.var = "time",
                        treatment = "Aobs", ## the *observed* treatment variable
                        treatment.prediction = "A", ## the *counterfactual* treatment variable
                        parallelize.predict = 1,
                        cv.fit = FALSE,
                        verbose = FALSE, verbose2 = FALSE,
                        browse0 = FALSE,
                        browse1 = FALSE,
                        browse2 = FALSE,
                        browse3 = FALSE
                        ) {

    if (length(seed)>0) set.seed(seed)

    if (browse0) browser()

    #--------------------------------
    #--- detect dependence on recurrent event process;

    history.Y <- depend.Y <- depend.Y.all <- unique(unlist(lapply(fit.hals, function(fit.hal) {
        indicator.names <- coef(fit.hal[["hal.fit"]], s=fit.hal[["lambda.cv"]])@Dimnames[[1]]
        nonzero.vars <- coef(fit.hal[["hal.fit"]], s=fit.hal[["lambda.cv"]])@i
        return(indicator.names[nonzero.vars+1][(substr(indicator.names[nonzero.vars+1], 1, 1) == "Y") |
                                               (indicator.names[nonzero.vars+1] %in% grep(":Y", indicator.names[nonzero.vars+1], value = TRUE))])
    })))

    history.main.Y <- unique(c(history.Y[history.Y %in% grep("Y.1", history.Y[!(history.Y %in% grep(":", history.Y, value = TRUE))], value = TRUE)],
                               sapply(strsplit(history.Y, ":"), function(vvar) {
                                   vvar[substr(vvar, 1, 1) == "Y"]
                               })))

    if (verbose2) {
        print("--------------------------------------------")
        print(history.main.Y)
    }

    if (length(depend.Y)>0) {
        max.Y <- max(as.numeric(gsub("Y.1 >= ", "", gsub("TRUE|FALSE", "", history.main.Y))))
        history.Y <- depend.Y <- depend.Y.all <- c(paste0("Y.1 >= ", 1:max.Y, "TRUE"),
                                                   history.Y[history.Y%in%grep(":", history.Y, value = TRUE)])
    } else{
        max.Y <- 1
    }
    
    if (verbose2) {
        print("--------------------------------------------")
        print(history.Y)
        print("--------------------------------------------")
    }
    
    if (length(depend.Y)>0) {
        
        depend.Y <- unique(unlist(sapply(depend.Y.all, function(yvar) {
            depend.tmp <- strsplit(yvar, ":")[[1]]
            return(depend.tmp[substr(depend.tmp, 1, 1) == "Y"])
        })))

        depend.Y <- c(
            grep("Y.dummy", depend.Y, value = TRUE),
            grep("Y.1", depend.Y, value = TRUE)[order(sapply(unlist(gsub("TRUE", "", grep("Y.1", depend.Y, value = TRUE))), function(yvar) {
                as.numeric(strsplit(yvar, ">=")[[1]][2])
            }))],
            grep("Y.time", depend.Y, value = TRUE)[order(sapply(unlist(gsub("TRUE", "", grep("Y.time", depend.Y, value = TRUE))), function(yvar) {
                as.numeric(strsplit(yvar, ">=")[[1]][2])
            }))])       

        depend.Y[substr(depend.Y, 1, 6) == "Y.time"] <- paste0("Y.dummy >= 1TRUE:", depend.Y[substr(depend.Y, 1, 6) == "Y.time"])

        depend.Y.list <- lapply(depend.Y, function(yvar) {
            tmp.Y <- unique(c(yvar, depend.Y.all[depend.Y.all %in% grep(yvar, depend.Y.all, value = TRUE) &
                                                 !(depend.Y.all %in% depend.Y)]))
            if (yvar == "Y.dummy >= 1TRUE" & length(grep("Y.dummy >= 1TRUE:Y.time", tmp.Y)) > 0) {
                return(tmp.Y[-grep("Y.dummy >= 1TRUE:Y.time", tmp.Y)])
            } else {
                return(tmp.Y)
            }
        })
        
        history.Y <- do.call("rbind", lapply(do.call(list.expand, replicate(length(depend.Y), 0:1, simplify = FALSE)), function(jj) unlist(jj)))
        #-- find "lower triangle" for Y.1;
        for (jj in (1:length(depend.Y))[substr(depend.Y, 1, 3) == "Y.1"]) {
            if (jj == grep("Y.1", depend.Y)[1]) {
                if (length(grep("Y.dummy", depend.Y)) > 0)
                    history.Y <- history.Y[!(history.Y[, grep("Y.dummy", depend.Y)[1]] == 0 & history.Y[, jj] == 1), ]
            } else {
                history.Y <- history.Y[!(history.Y[, jj-1] == 0 & history.Y[, jj] == 1), ]
            }
        }
        #-- find "lower triangle" for Y.time;
        for (jj in (1:length(depend.Y))[grep("Y.time", depend.Y)]) {
            if (jj == grep("Y.time", depend.Y)[1]) {
                if (length(grep("Y.dummy", depend.Y)) > 0)
                    history.Y <- history.Y[!(history.Y[, grep("Y.dummy", depend.Y)[1]] == 0 & history.Y[, jj] == 1), ]
            } else {
                history.Y <- history.Y[!(history.Y[, jj-1] == 0 & history.Y[, jj] == 1), ]
            }
        }
        colnames(history.Y) <- depend.Y
    }
    
    #--------------------------------
    #--- prediction part;

 
    if (browse1) browser()

    by.vars1 <- c("id", treatment.prediction)
    by.vars2 <- c("id", "observed.Y", treatment.prediction)

    hal.vars.list <- lapply(1:length(fit.hals), function(kk) {
        fit.hals[[kk]][["hal.fit"]]$beta@Dimnames[[1]]
    })

    hal.vars <- unique(c(unlist(hal.vars.list)))

    hal.vars.static <- sapply(hal.vars[!(hal.vars %in% grep("Y.|:Y", hal.vars, value = TRUE))],
                              function(xstatic) paste0("(", gsub(":", "):(", xstatic), ")"))
    hal.vars.dynamic <- sapply(hal.vars[hal.vars %in% grep("Y.|:Y", hal.vars, value = TRUE)],
                               function(xstatic) paste0("(", gsub(":", "):(", xstatic), ")"))

    X.hal.static <- Matrix(
        model.matrix(formula(paste0(delta.var, "~-1+", paste(gsub("FALSE|TRUE", "", hal.vars.static),
                                                             collapse = "+"))), 
                     data=pseudo.dt[,
                     (treatment) := get(treatment.prediction)]), sparse=TRUE) 

    ## pseudo.dt[, risk.time := diff(c(get(time.var), max.time[.N])), by = by.vars1]

    if (browse2) browser()

    if (length(history.Y) == 0) {
        history.Y <- matrix(c(0), nrow = 1)
    }
    
    for (ii in 1:nrow(history.Y)) {
        pseudo.dt[Y.1 == rowSums(history.Y)[ii], which.Y := ii]
    }
        
    pseudo.dt <- cbind(pseudo.dt,
                       do.call("cbind",
                               mclapply(
                                   X = 1:nrow(history.Y),
                                   mc.cores = min(detectCores()-1, parallelize.predict),#parallelize.cve
                                   FUN = function(ii) {

                                       pseudo.ii <- copy(pseudo.dt)[,
                                           (treatment) := get(treatment.prediction)][, which.Y := ii][, Y.1 := rowSums(history.Y)[ii]]

                                       if (length(hal.vars.dynamic)>0) {
                                           X.ii <- cbind(X.hal.static, Matrix(
                                                                           model.matrix(formula(paste0(delta.var, "~-1+", paste(gsub("FALSE|TRUE", "", hal.vars.dynamic), collapse = "+"))), 
                                                                                        data=pseudo.ii), sparse=TRUE))
                                       } else {
                                           X.ii <- X.hal.static
                                       }

                                       colnames.missing <- colnames(X.ii)[!colnames(X.ii) %in% hal.vars]

                                       colnames(X.ii)[colnames(X.ii) %in% colnames.missing] <- sub(
                                           "^([^:]+):([^:]+)$",
                                           "\\2:\\1",
                                           colnames(X.ii)[colnames(X.ii) %in% colnames.missing]
                                       )
               
                                       out.vars <- rep(NA, length(fit.hals))

                                       if (browse3) browser()

                                       for (kk in 1:length(fit.hals)) {
        
                                           delta.value <- fit.hals[[kk]][["delta.value"]]

                                           if (cv.fit) {
                                               for (jjj in 1:length(fit.hals[[kk]]$train.fit)) {
                                                   test.set <- fit.hals[[kk]]$train.fit[[jjj]]$test.set
                                                   train.fit <- fit.hals[[kk]]$train.fit[[jjj]]$train.fit 
                                                   pseudo.ii[id %in% test.set, (paste0("fit.lambda", delta.value)) := exp(predict(train.fit, X.ii[pseudo.ii[["id"]] %in% test.set, hal.vars.list[[kk]]],
                                                                                                                  newoffset=0, s=fit.hals[[kk]][["lambda.cv"]]))]
                                               }
                                           } else {
                                               pseudo.ii[, (paste0("fit.lambda", delta.value)) := exp(predict(fit.hals[[kk]][["hal.fit"]], X.ii[, hal.vars.list[[kk]]],
                                                                                                              newoffset=0, s=fit.hals[[kk]][["lambda.cv"]]))]
                                           }
                                           pseudo.ii[, (paste0("fit.dLambda", delta.value)) := get(paste0("fit.lambda", delta.value))*risk.time]
                                           pseudo.ii[, (paste0("fit.Lambda", delta.value)) := cumsum(get(paste0("fit.dLambda", delta.value))), by = by.vars2]

                                           pseudo.ii[, (paste0("P", delta.value)) := get(paste0("fit.dLambda", delta.value))]
                                           pseudo.ii[, (paste0("surv", delta.value)) := exp(-get(paste0("fit.Lambda", delta.value)))]
                                           pseudo.ii[, (paste0("surv", delta.value, ".1")) := get(paste0("surv", delta.value)), by = "id"]

                                           setnames(pseudo.ii, paste0("P", delta.value), paste0("P", delta.value, ".Y", ii))

                                           out.vars[kk] <- paste0("P", delta.value, ".Y", ii)
                                       }

                                       return(pseudo.ii[, out.vars, with = FALSE])
                                   })))

    if (browse2) browser()

    return(pseudo.dt)
}

######################################################################
### predict.hal.R ends here
