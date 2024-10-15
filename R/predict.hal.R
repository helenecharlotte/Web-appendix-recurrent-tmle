### predict.hal.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Oct 15 2024 (09:34) 
## Version: 
## Last-Updated: Oct 15 2024 (09:34) 
##           By: Helene
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

predict.hal <- function(seed = 13349,
                        fit.hals, ## output from fit.hal
                        pseudo.dt, ## tmp3,
                        delta.var = "delta",
                        delta.value = 1,
                        time.var = "time",
                        covars,
                        treatment = "Aobs", ## the *observed* treatment variable
                        treatment.prediction = "A", ## the *counterfactual* treatment variable
                        cut.one.way = 10,
                        cut.time.varying = 5,
                        Y.time.grid = NULL,
                        cut.time = 15,
                        two.way = cbind(var1="", var2=""),
                        cut.two.way = 5,
                        cut.time.treatment = NULL,
                        cut.time.covar = NULL,
                        penalize.treatment = FALSE,
                        penalize.time = FALSE,
                        cv.glmnet = FALSE,
                        V = 5,
                        lambda.cvs = c(sapply(1:5, function(jjj) (9:1)/(10^jjj))),
                        verbose = FALSE,
                        browse0 = FALSE,
                        browse1 = FALSE,
                        browse2 = FALSE,
                        browse3 = FALSE
                        ) {

    set.seed(seed)

    if (browse0) browser()

    #--------------------------------
    #--- detect dependence on recurrent event process;

    # - following should incorporate interactions ----------

    history.Y <- depend.Y <- depend.Y.all <- unique(unlist(lapply(fit.hals, function(fit.hal) {
        indicator.names <- coef(fit.hal[["hal.fit"]], s=fit.hal[["lambda.cv"]])@Dimnames[[1]]
        nonzero.vars <- coef(fit.hal[["hal.fit"]], s=fit.hal[["lambda.cv"]])@i
        return(indicator.names[nonzero.vars+1][(substr(indicator.names[nonzero.vars+1], 1, 1) == "Y") |
                                               (indicator.names[nonzero.vars+1] %in% grep(":Y", indicator.names[nonzero.vars+1], value = TRUE))])
    })))
    
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

    #pseudo.dt[, check.order := 1:nrow(pseudo.dt)]

    X.hal.a <- basis.fun(pseudo.dt = copy(pseudo.dt)[pseudo.dt[["observed.Y"]] == 1][, (treatment) := get(treatment.prediction)],
                         delta.var = delta.var,
                         delta.value = delta.value,
                         covars = covars,
                         treatment = treatment,
                         cut.one.way = cut.one.way,
                         cut.time.varying = cut.time.varying,
                         Y.time.grid = Y.time.grid,
                         cut.time = cut.time,
                         two.way = two.way,
                         cut.two.way = cut.two.way,
                         cut.time.treatment = cut.time.treatment,
                         cut.time.covar = cut.time.covar,
                         predict = TRUE)

    #pseudo.dt[, grid.period := X.hal.a$hal.pseudo.dt[["grid.period"]]]
    #pseudo.dt[, check.order1 := X.hal.a$hal.pseudo.dt[["check.order"]]]

    if (browse1) browser()
    
    X.hal.a <- X.hal.a$X

    if (length(history.Y)>0) {

        # --- NB - should think interactions into this?

        pseudo.dt <- do.call("rbind", lapply(1:nrow(history.Y), function(ii) {
            pseudo.ii <- copy(pseudo.dt)[pseudo.dt[["observed.Y"]] == 1][, (treatment) := get(treatment.prediction)][, which.Y := ii]
            for (jj in 1:ncol(history.Y)) {
                #pseudo.ii[, (colnames(history.Y)[jj]) := history.Y[ii,jj]]
                for (yvar in depend.Y.list[[jj]]) {
                    if (yvar == colnames(history.Y)[jj]) {
                        pseudo.ii[, (yvar) := history.Y[ii,jj]]
                    } else { #--- this to handle interactions
                        #print(yvar)
                        #print(strsplit(yvar, ":")[[1]][substr(strsplit(yvar, ":")[[1]], 1, 1) != "Y"])
                        pseudo.ii[, (yvar) := history.Y[ii,jj]*X.hal.a[, strsplit(yvar, ":")[[1]][substr(strsplit(yvar, ":")[[1]], 1, 1) != "Y"]]]
                    }
                }
            }
            return(pseudo.ii)
        }))

        X.hal.a <- do.call("rbind", lapply(1:nrow(history.Y), function(ii) {
            X.ii <- copy(X.hal.a)
            observed.Y <- rep(1, length = nrow(X.hal.a))
            for (jj in 1:ncol(history.Y)) {
                observed.Y <- observed.Y*(X.ii[, colnames(history.Y)[jj]] == history.Y[ii,jj])
                for (yvar in depend.Y.list[[jj]]) {
                    if (yvar == colnames(history.Y)[jj]) {
                        X.ii[, yvar] <- history.Y[ii,jj]
                    } else {
                        X.ii[, yvar] <- history.Y[ii,jj]*X.ii[, strsplit(yvar, ":")[[1]][substr(strsplit(yvar, ":")[[1]], 1, 1) != "Y"]]
                    }
                }
            }
            return(cbind(X.ii, observed.Y))
        }))

        pseudo.dt[, observed.Y := X.hal.a[, "observed.Y"]]

        if (FALSE) {
            pseudo.dt[observed.Y == 1]
            pseudo.dt[observed.Y == 1, table(get("Y.dummy >= 1TRUE"), Y.dummy)]
            pseudo.dt[observed.Y == 1, table(get("Y.dummy >= 1TRUE"), get("Y.dummy >= 1TRUE:Y.time >= 903TRUE"), Y.dummy, Y.time >= 903)]
            pseudo.dt[observed.Y == 0, table(get("Y.dummy >= 1TRUE"), get("Y.dummy >= 1TRUE:Y.time >= 903TRUE"), Y.dummy, Y.time >= 903)]
            pseudo.dt[, table(get("Y.dummy >= 1TRUE"), get("Y.dummy >= 1TRUE:Y.time >= 903TRUE"))]
            pseudo.dt[, table(Y.dummy, Y.time >= 903)]
            pseudo.dt[observed.Y == 0, table(Y.dummy, Y.time >= 903)]
            pseudo.dt[observed.Y == 1, table(Y.dummy, Y.time >= 903)]
            pseudo.dt[observed.Y == 1, table(Y.dummy, Y.dummy & Y.time >= 903)]
            #pseudo.dt[pseudo.dt[["observed.Y"]] == 1 & pseudo.dt[["Y.dummy >= 1TRUE"]]]
        }
        
        X.hal.a <- X.hal.a[, colnames(X.hal.a) != "observed.Y"]
        
        by.vars1 <- c("id", colnames(history.Y), treatment.prediction)
        by.vars2 <- c("id", "observed.Y", treatment.prediction)

    } else {
        by.vars1 <- by.vars2 <- c("id", treatment.prediction)
    }

    pseudo.dt[, risk.time := diff(c(0,get(time.var))), by = by.vars1]
    
    if (browse2) browser()
    
    for (kk in 1:length(fit.hals)) {
        
        delta.value <- fit.hals[[kk]][["delta.value"]]
        pseudo.dt[, (paste0("fit.lambda", delta.value)) := exp(predict(fit.hals[[kk]][["hal.fit"]], X.hal.a,
                                              newoffset=0, s=fit.hals[[kk]][["lambda.cv"]]))]
        pseudo.dt[, (paste0("fit.dLambda", delta.value)) := get(paste0("fit.lambda", delta.value))*risk.time]
        pseudo.dt[, (paste0("fit.Lambda", delta.value)) := cumsum(get(paste0("fit.dLambda", delta.value))), by = by.vars2]

        pseudo.dt[, (paste0("P", delta.value)) := get(paste0("fit.dLambda", delta.value))]
        pseudo.dt[, (paste0("surv", delta.value)) := exp(-get(paste0("fit.Lambda", delta.value)))]
        pseudo.dt[, (paste0("surv", delta.value, ".1")) := get(paste0("surv", delta.value)), by = "id"]
    }

    return(pseudo.dt)
}

######################################################################
### predict.hal.R ends here
