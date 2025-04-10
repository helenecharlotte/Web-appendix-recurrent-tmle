### tmle.estimation.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Oct 14 2024 (15:01) 
## Version: 
## Last-Updated: Apr  1 2025 (13:26) 
##           By: Helene
##     Update #: 509
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

tmle.est.fun <- function(dt,
                         intervention.A = 1,
                         tau = 1.2,
                         fit.type1 = list(model = "Surv(tstart, tstop, delta == 1)~A+L2+L1.squared+L3",
                                          fit = "cox"),
                         fit.type2 = list(model = "Surv(tstart, tstop, delta == 2)~A+L2+L1+L3",
                                          fit = "cox"),
                         fit.type0 = list(model = "Surv(tstart, tstop, delta == 0)~A+L2+L1+L3",
                                          fit = "cox"),
                         fit.treatment = list(model = "A~L1+L2+L3",
                                              fit = "glm"),
                         verbose = FALSE,
                         browse.hal = FALSE,
                         browse2 = FALSE, 
                         #--- simpler estimators: 
                         baseline.tmle = FALSE,
                         naive.gcomp = FALSE, 
                         standard.np = FALSE,
                         max.iter = FALSE,
                         #--- HAL parameters:
                         lambda.cvs = c(sapply(1:5, function(jjj) (9:1)/(10^jjj))),
                         cut.one.way = 15,
                         cut.time.varying = 5,
                         Y.time.grid = NULL,
                         cut.time = 35,
                         two.way = cbind(var1="", var2=""),
                         cut.two.way = 5,
                         cut.time.treatment = 0,
                         cut.time.covar = 0,
                         penalize.time = FALSE,
                         reduce.seed.dependence = FALSE, 
                         stime = 0,
                         scovar = 0,
                         tmle.criterion = 1,
                         V = 10, 
                         cv.glmnet = FALSE,
                         seed.hal = NULL,
                         output.lambda.cvs = FALSE
                         ) {

    n <- length(unique(dt[["id"]]))

    #-- if no "fit" is specified, then specify it as cox: 
    if (!is.list(fit.type1)) {
        fit.type1 <- list(model = fit.type1, fit = "cox")
    }

    if (!is.list(fit.type2)) {
        fit.type2 <- list(model = fit.type2, fit = "cox")
    }

    if (!is.list(fit.type0)) {
        fit.type0 <- list(model = fit.type0, fit = "cox")
    }

    #-- if no "fit" is specified for treatment model, then specify it as glm: 
    if (!is.list(fit.treatment)) {
        fit.treatment <- list(model = fit.treatment, fit = "glm")
    }
    
    #-- covariates/predictors extracted from models: 
    varnames <- unique(unlist(lapply(list(fit.type1[["model"]], fit.type2[["model"]], fit.type0[["model"]], fit.treatment[["model"]]), function(fit.type) {
        strsplit(strsplit(fit.type, "~")[[1]][2], "\\+|\\*")[[1]]
    })))

    if ("Y.dummy." %in% substr(varnames,1,nchar("Y.dummy."))) {
        Y.dummy.max <- max(as.numeric(gsub("Y.dummy.", "", varnames[substr(varnames,1,nchar("Y.dummy.")) == "Y.dummy."])))
    } else {
        Y.dummy.max <- NULL
    }

    varnames <- varnames[varnames != "1"]
    varnames <- varnames[!(varnames %in% grep(">=", varnames, value = TRUE))]

    #-- name of treatment variable: 
    Aname <- strsplit(fit.treatment[["model"]], "~")[[1]][1] 
   
    if (typeof(dt[["time"]]) != "double") {
        warning("NB: the time variable is not numeric - will be converted")
        dt[["time"]] <- as.numeric(dt[["time"]])
    }

    if ("L1" %in% names(dt)) dt[, L1.squared := L1^2]
    
    dt[, tstart := c(0, time[-.N]), by = "id"]
    dt[, tstop := time]

    #-- dummy variables for history of recurrent event process:
    dt[, Y := cumsum(1*(delta == 1)), by = "id"]
    dt[, Y.1 := c(0, Y[-.N]), by = "id"]
    dt[, Y.dummy := (Y.1>0)]
    dt[, Y2.dummy := (Y.1 >= 2)]
    dt[, Y3.dummy := (Y.1 >= 3)]
    if (length(Y.dummy.max)>0) {
        for (kY in 1:Y.dummy.max) {
            dt[, (paste0("Y.dummy.", kY)) := (Y.1 >= kY)]
        }
    }
    dt[, Y.time := time-tstart, by = "id"]

    #--------------------------------
    #-- "G-part"; for clever weight estimation:
    # (we start with cox models, if HAL is specified this is fitted later)
    fit.cox0 <- coxph(as.formula(fit.type0[["model"]]), data = dt, 
                      control = coxph.control(timefix = FALSE))
    if (verbose) print(fit.cox0)
    dt[, idN := 1:.N, by = "id"]
    if (fit.treatment[["fit"]] == "glm") {
        fitA <- glm(as.formula(fit.treatment[["model"]]), data=dt[idN == 1], family=binomial)
        if (verbose) print(summary(fitA))
    } else {
        print("NB: need to incorporate other estimations methods than glm for treatment")
    }
    dt[, probA := predict(fitA, newdata = dt, type = "response")]

    #--------------------------------
    #-- "Q-part":
    # (we start with cox models, if HAL is specified this is fitted later)
    fit.cox1 <- coxph(as.formula(fit.type1[["model"]]), data = dt, 
                      control = coxph.control(timefix = FALSE))
    fit.cox2 <- coxph(as.formula(fit.type2[["model"]]), data = dt, 
                      control = coxph.control(timefix = FALSE))
    if (verbose) print(fit.cox1)
    if (verbose) print(fit.cox2)

    #-- with dependence on Y.time>t0: 
    if ("t0" %in% names(fit.type1)) {
        t0 <- fit.type1[["t0"]]
        dt[, Y.time.dummy := 0]
        dt2 <- rbind(copy(dt)[Y.1 >= 1 & tstop-tstart>t0][, tstop := tstart+t0][, delta := 0],
                     copy(dt)[tstop-tstart>t0 & Y.1 >= 1, Y.time.dummy := 1][tstop-tstart>t0 & Y.1 >= 1, tstart := tstart+t0])[order(id, tstart)]
        dt2[, time := tstop]
        fit.cox1 <- coxph(as.formula(fit.type1[["model"]]), data = dt2)
        if (verbose) print(fit.cox1)
        fit.cox2 <- coxph(as.formula(fit.type2[["model"]]), data = dt2)
        if (verbose) print(fit.cox2)
    } else {
        dt2 <- copy(dt)
    }

    #-- get baseline intensities: 
    tmp.type1 <- suppressWarnings(setDT(basehaz(fit.cox1, centered=TRUE)))[, dhazard1 := c(hazard[1],diff(hazard))][, hazard1 := hazard][, -"hazard", with = FALSE]
    tmp.type2 <- suppressWarnings(setDT(basehaz(fit.cox2, centered=TRUE)))[, dhazard2 := c(hazard[1],diff(hazard))][, hazard2 := hazard][, -"hazard", with = FALSE]
    tmp.type0 <- suppressWarnings(setDT(basehaz(fit.cox0, centered=TRUE)))[, dhazard0 := c(hazard[1],diff(hazard))][, hazard0 := hazard][, -"hazard", with = FALSE]

    #-- get all unique times; 
    unique.times <- sort(unique(dt2[["time"]]))
    unique.times <- unique.times[unique.times <= tau]
    all.times <- data.table(expand.grid(time = unique.times,
                                        id = unique(dt[["id"]])))

    #-- collect data with all time-points;
    tmp.inner <- merge(dt2[, time.obs := time], all.times, by = c("id", "time"), all = TRUE)[order(id, time.obs)]

    for (varname in varnames) {
        tmp.inner[, (varname) := na.locf(get(varname)), by = "id"]
    }

    #--------------------------------
    #-- expanded dataset to work with for TMLE:
    # (should rename tmp3:))
    tmp3 <- merge(merge(merge(tmp.inner[order(id, time)][is.na(delta), delta := 0][, -c("tstop"), with = FALSE],
                              tmp.type1[dhazard1>0], by = "time", all = TRUE),
                        tmp.type2[dhazard2>0], by = "time", all = TRUE),
                  tmp.type0[dhazard0>0], by = "time", all = TRUE)[order(id,time)][!is.na(id)]

    #-- get last observed time for each line: 
    tmp3[delta == 1, time.obs2 := time.obs]
    tmp3[, time.obs.last := nafill(time.obs2, "locf"), by = "id"]
    tmp3[is.na(time.obs.last), time.obs.last := 0]
    
    tmp3[, time.obs := nafill(time.obs, "nocb"), by = "id"]
    tmp3[is.na(time.obs), time.obs := -Inf]

    tmp3[is.na(dhazard0), dhazard0 := 0]
    tmp3[is.na(dhazard1), dhazard1 := 0]
    tmp3[is.na(dhazard2), dhazard2 := 0]

    if ("dhazard1.t0.0" %in% names(tmp3)) {
        tmp3[is.na(dhazard1.t0.0), dhazard1.t0.0 := 0]
        tmp3[is.na(dhazard1.t0.1), dhazard1.t0.1 := 0]
    }

    #-- dummy variables for history of recurrent event process:
    tmp3[, Y := cumsum(1*(delta == 1)), by = "id"]
    tmp3[, Y.1 := c(0, Y[-.N]), by = "id"]
    tmp3[, Y.dummy := (Y.1>0)]
    tmp3[, Y2.dummy := (Y.1 >= 2)]
    tmp3[, Y3.dummy := (Y.1 >= 3)]
    if (length(Y.dummy.max)>0) {
        tmp3[, Y.dummy.index := findInterval(Y.1, 0:Y.dummy.max)]
        for (kY in 1:Y.dummy.max) {
            tmp3[, (paste0("Y.dummy.", kY)) := (Y.1 >= kY)]
        }
    }
    tmp3[, Y.time := time-c(0, time.obs.last[-.N]), by = "id"]
    tmp3[Y.dummy == 0, Y.time := 0]
    tmp3[, final.time := max(time.obs), by = "id"]
    tmp3[time > final.time, Y.time := time-final.time]
    tmp3 <- tmp3[time <= tau] #- only include data up to time tau!
    tmp3[, before.tau := 1*(time <= tau)]
    tmp3[, time.diff := c(time[-1], tau) - time, by = "id"]
    tmp3[Y.dummy == 0, Y.time := 0] #- only if there was a jump in past

    if ("t0" %in% names(fit.type1)) {
        tmp3[, Y.time.dummy := (Y.1 >= 1)*(Y.time > t0+0.00000001)]
    }

    #-- remame treatment variable to "Aobs" in tmp3 data: 
    setnames(tmp3, Aname, "Aobs")
    #-- the treatment variable is then set to the interventional level: 
    tmp3[, (Aname) := intervention.A]

    #--------------------------------
    #-- standard nonparametric estimator:
    if (standard.np) {
        if (sessionInfo()$otherPkgs$prodlim$Version == "2019.11.13") {
            print("NB: remove this hack at some point; problem with prodlim versions")
            dt.km <- copy(dt)[, N := .N, by = "id"][idN == N]
            #-- kaplan-meier estimator for survival;
            km.mod <- "Hist(time, delta==2)~A"
            km.fit <- summary(prodlim(formula(km.mod), data=dt.km),
                              times=unique.times[unique.times<=tau], asMatrix=TRUE)$table
            #-- kaplan-meier estimator for censoring;
            km.cens.mod <- "Hist(time, delta==0)~A"
            km.cens.fit <- summary(prodlim(formula(km.cens.mod), data=dt.km),
                                   times=unique.times[unique.times<=tau], asMatrix=TRUE)$table
            #-- nelson-aalen estimator for intensity of recurrent events;
            na.mod <- "Hist(time, delta)~A"
            na.fit <- summary(prodlim(formula(na.mod), data=dt),
                              cause=1,times=unique.times[unique.times<=tau], asMatrix=TRUE)$table
            dt.non.fit <- merge(data.table(A = as.numeric(gsub("A\\=", "", km.fit[, "X"])),
                                           time = as.numeric(km.fit[, "time"]),
                                           n.risk.km = as.numeric(km.fit[, "n.risk"]),
                                           n.event.km = as.numeric(km.fit[, "n.event"]),
                                           S.km = as.numeric(km.fit[, "surv"])),
                                data.table(A = as.numeric(gsub("A\\=", "", na.fit[, "X"])),
                                           time = as.numeric(na.fit[, "time"]),
                                           n.risk.na = as.numeric(na.fit[, "n.risk"]),
                                           n.event.na = as.numeric(na.fit[, "n.event"])),
                                by = c("A", "time"))[, S.km.1 := c(1, S.km[-.N]), by = c("A")]
            dt.cens.non.fit <- merge(dt.non.fit,
                                     data.table(A = as.numeric(gsub("A\\=", "", km.cens.fit[, "X"])),
                                                time = as.numeric(km.cens.fit[, "time"]),
                                                n.risk.censkm = as.numeric(km.cens.fit[, "n.risk"]),
                                                n.event.cens.km = as.numeric(km.cens.fit[, "n.event"]),
                                                S.cens.km = as.numeric(km.cens.fit[, "surv"])),
                                     by = c("A", "time"))[, S.cens.km.1 := c(1, S.cens.km[-.N]), by = c("A")]
            setnames(dt.non.fit, "A", Aname)
            setnames(dt.cens.non.fit, "A", Aname)
        } else {
            dt.km <- copy(dt)[, N := .N, by = "id"][idN == N]
            #-- kaplan-meier estimator for survival;
            km.mod <- paste0("Hist(time, delta==2)~", Aname)
            km.fit <- data.table(
                summary(prodlim(formula(km.mod), data=dt.km),
                        times=unique.times[unique.times<=tau])[, c(Aname, "time", "n.risk", "n.event", "surv"), with = FALSE])
            setnames(km.fit, c("n.risk", "n.event", "surv"),  c("n.risk.km", "n.event.km", "S.km"))
            #-- kaplan-meier estimator for censoring;
            km.cens.mod <- paste0("Hist(time, delta==0)~", Aname)
            km.cens.fit <- data.table(
                summary(prodlim(formula(km.cens.mod), data=dt.km),
                        times=unique.times[unique.times<=tau])[, c(Aname, "time", "n.risk", "n.event", "surv"), with = FALSE])
            setnames(km.cens.fit, c("n.risk", "n.event", "surv"),  c("n.risk.cens.km", "n.event.cens.km", "S.cens.km"))
            #-- nelson-aalen estimator for intensity of recurrent events;
            na.mod <- paste0("Hist(time, delta)~", Aname)
            na.fit <- data.table(
                summary(prodlim(formula(na.mod), data=dt),
                        cause=1, times=unique.times[unique.times<=tau])[, c(Aname, "time", "n.risk", "n.event"), with = FALSE])
            setnames(na.fit, c("n.risk", "n.event"),  c("n.risk.na", "n.event.na"))
            dt.non.fit <- merge(km.fit, na.fit,
                                by = c(Aname, "time"))[, S.km.1 := c(1, S.km[-.N]), by = Aname]
            dt.cens.non.fit <- merge(dt.non.fit, km.cens.fit,
                                     by = c(Aname, "time"))[, S.cens.km.1 := c(1, S.cens.km[-.N]), by = Aname]
        }
        dt.cens.non.fit[, target.tau := sum(S.km.1*n.event.na/n.risk.km), by = Aname]
        dt.cens.non.fit[, target.t := cumsum(S.km.1*n.event.na/n.risk.km), by = Aname]
        tmp3 <- merge(tmp3,  dt.cens.non.fit, by = c(Aname, "time"))
        tmp3[, final.time := max(time.obs), by = "id"]
        dt.non.fit.eic <- tmp3[time <= tau, sum((time <= final.time & get(Aname) == Aobs)/S.cens.km.1*(
            ((delta == 1) - n.event.na/n.risk.km) -
            (target.tau - target.t)/S.km*((delta == 2) - n.event.km/n.risk.km)
        )), by = c("id", "Aobs")][, eic := V1]
        n.intervention <- dt[get(Aname) == intervention.A, length(unique(id))]
        dt.non.fit.eic[Aobs == intervention.A, sqrt(mean(eic^2/n.intervention))]
        return(c(np.est = dt.non.fit[, sum(S.km.1*n.event.na/n.risk.km), by = Aname][get(Aname) == intervention.A][,2][[1]],
                 np.se = dt.non.fit.eic[Aobs == intervention.A, sqrt(mean(eic^2/n.intervention))]))
    }

    #--------------------------------
    #-- compute needed quantities for (initial) estimation and targeting
    
    tmp3[, exp1 := exp(predict(fit.cox1, newdata=tmp3, type="lp"))]
    tmp3[, surv1 := exp(-cumsum(exp1*dhazard1)), by = "id"]

    tmp3[, exp2 := exp(predict(fit.cox2, newdata=tmp3, type="lp"))]
    tmp3[, surv2 := exp(-cumsum(exp2*dhazard2)), by = "id"]
    tmp3[, surv2.1 := c(1, surv2[-.N]), by = "id"]

    tmp3[, P1 := dhazard1*exp1]
    tmp3[, P2 := dhazard2*exp2]

    tmp3[, exp0 := exp(predict(fit.cox0, newdata=tmp3, type="lp"))]
    tmp3[, surv0 := exp(-cumsum(exp0*dhazard0)), by = "id"]
    tmp3[, surv0.1 := c(1, surv0[-.N]), by = "id"]

    tmp3[, surv0.cox := surv0]
    tmp3[, surv0.1.cox := surv0.1]

    #--------------------------------
    #-- naive g-comp estimator:
    
    if (naive.gcomp) {
        return(mean(tmp3[time <= tau, sum(P1*surv2.1), by = "id"][[2]]))
    }

    #--------------------------------
    #-- HAL estimation:

    any.hal <- unlist(lapply(list(fit.type1, fit.type2, fit.type0), function(each) length(grep("hal", each[["fit"]]))>0))

    if (any(any.hal)) {

        covars <- varnames[varnames != Aname]
        if (length(grep("\\.squared", covars))>0) {
            covars <- covars[!(covars %in% grep("\\.squared", covars, value = TRUE))]
        }

        #-- interactions: 
        if (cut.two.way > 0) {
            time.varying.covars <-
                covars[sapply(covars, function(covar) {
                    (nrow(tmp3[, unique(get(covar)), by = "id"])>n)
                })]
            two.way <- t(combn(c(covars, "Aobs"), 2))
            colnames(two.way) <- c("var1", "var2")
            two.way <- two.way[!(two.way[,1] %in% time.varying.covars & two.way[,2] %in% time.varying.covars),]
        }

        tmp.hal <- copy(tmp3)[, observed.Y := 1]

        #-- fit HAL: 

        fit.hals <- lapply(c(1,2,0)[any.hal], function(delta.value) {
            fit.hal(
                pseudo.dt = tmp.hal,
                delta.var = "delta",
                delta.value = delta.value,
                time.var = "time",
                covars = covars,
                treatment = "Aobs",
                treatment.prediction = Aname,
                lambda.cvs = lambda.cvs,
                cut.one.way = cut.one.way,
                cut.time.varying = cut.time.varying,
                Y.time.grid = Y.time.grid,
                cut.time = cut.time,
                two.way = two.way,
                cut.two.way = cut.two.way,
                cut.time.treatment = cut.time.treatment,
                cut.time.covar = cut.time.covar,
                penalize.time = penalize.time,
                stime = stime,
                scovar = scovar,
                reduce.seed.dependence = reduce.seed.dependence,
                V = V,
                cv.glmnet = cv.glmnet,
                verbose = verbose,
                seed = seed.hal)
        })


        if (output.lambda.cvs) {
            lambda.cvs <- lapply(fit.hals, function(fh) c(delta = fh$delta.value, lambda.cv = fh$lambda.cv))
            lambda.cvs.names <- sapply(lambda.cvs, function(lc) {
                paste0("delta = ", lc["delta"])
            })
            lambda.cvs <- sapply(lambda.cvs, function(lc) {
                lc["lambda.cv"]
            })
            names(lambda.cvs) <- lambda.cvs.names
            hal.coefs <- lapply(fit.hals, function(fh) {
                delta <- fh$delta.value
                get.coef <- coef(fh$hal.fit, s = fh$lambda)
                coef.dt <- data.table(non.zero = get.coef@Dimnames[[1]][get.coef@i+1],
                                      coef = get.coef@x)[non.zero %in% c("Aobs", "Y.dummy >= 1TRUE")]
                coefs <- coef.dt[["coef"]]
                names(coefs) <- paste0(coef.dt[["non.zero"]], "(d = ", delta, ")")
                return(coefs)
            })
        }
        
        #-- predict in expanded data: 
        tmp.hal <- predict.hal(
            fit.hals = fit.hals,
            pseudo.dt = tmp.hal,
            delta.var = "delta",
            time.var = "time",
            covars = covars,
            treatment = "Aobs",
            treatment.prediction = Aname,
            cut.one.way = cut.one.way,
            cut.time.varying = cut.time.varying,
            Y.time.grid = Y.time.grid,
            cut.time = cut.time,
            two.way = two.way,
            cut.two.way = cut.two.way,
            cut.time.treatment = cut.time.treatment,
            cut.time.covar = cut.time.covar,
            stime = stime,
            scovar = scovar,
            verbose = verbose,
            seed = seed.hal)

        if (browse.hal) browser()

    }

    #--------------------------------    
    #-- prediction part is over, so remove data after tau:
    #-- NB CHECK
    
    tmp3 <- tmp3[time <= tau]
    if (any(any.hal)) tmp.hal <- tmp.hal[time <= tau]
    unique.times <- unique.times[unique.times <= tau]
    
    #--------------------------------    
    #-- to handle dependence on jump in the past:
    
    if (!baseline.tmle) { #-- only relevant for the general tmle:

        if (length(Y.dummy.max)>0) { #-- to handle effect of Y.1  (no of jumps in the past)
          
            index.j <- 1:(Y.dummy.max+1)
            index.j1 <- sapply(index.j+1, function(ij) min(ij, max(as.numeric(index.j))))

            for (kY in index.j) {

                tmp3.Y.kY <- copy(tmp3)

                for (kY2 in index.j) {
                    if (kY2 < kY) {
                        tmp3.Y.kY[, (paste0("Y.dummy.", kY2)) := 1]
                    } else {
                        tmp3.Y.kY[, (paste0("Y.dummy.", kY2)) := 0]
                    }
                    if (kY2 == 1) {
                        tmp3.Y.kY[, Y.dummy := get(paste0("Y.dummy.", kY2))]
                    }
                }

                tmp3[, (paste0("exp1.Y", kY)) := exp(predict(fit.cox1, newdata=tmp3.Y.kY, type="lp"))]
                tmp3[, (paste0("exp2.Y", kY)) := exp(predict(fit.cox2, newdata=tmp3.Y.kY, type="lp"))]                

                tmp3[, (paste0("P1.Y", kY)) := dhazard1*get(paste0("exp1.Y", kY))]
                tmp3[, (paste0("P2.Y", kY)) := dhazard2*get(paste0("exp2.Y", kY))]
                
            }
            
        } else {  #-- to handle only effect of Y.dummy (jump or not in the past)
        
            tmp3[, Y1 := (Y >= 1)*Y+(Y == 0)*1]
            tmp3.Y1 <- copy(tmp3)
            tmp3.Y0 <- copy(tmp3)
            tmp3.Y1[, Y.dummy := 1]
            tmp3.Y0[, Y.dummy := 0]

            tmp3[, exp1.Y1 := exp(predict(fit.cox1, newdata=tmp3.Y1, type="lp"))]
            tmp3[, exp1.Y0 := exp(predict(fit.cox1, newdata=tmp3.Y0, type="lp"))]

            tmp3[, exp2.Y1 := exp(predict(fit.cox2, newdata=tmp3.Y1, type="lp"))]
            tmp3[, exp2.Y0 := exp(predict(fit.cox2, newdata=tmp3.Y0, type="lp"))]

            tmp3[, P1.Y1 := dhazard1*exp1.Y1]
            tmp3[, P1.Y0 := dhazard1*exp1.Y0]

            tmp3[, P2.Y1 := dhazard2*exp2.Y1]
            tmp3[, P2.Y0 := dhazard2*exp2.Y0]

            if (browse2) {
                tmp3[, P1.Y1.cox := dhazard1*exp1.Y1]
                tmp3[, P1.Y0.cox := dhazard1*exp1.Y0]

                tmp3[, P2.Y1.cox := dhazard2*exp2.Y1]
                tmp3[, P2.Y0.cox := dhazard2*exp2.Y0]
            }

        }
        
    }
    
    #--------------------------------    
    #-- update estimators if fitted with HAL:

    if (any(any.hal)) {

        if (baseline.tmle) { #-- for the simple baseline version of tmle: 
            
            if ("Y.dummy >= 1TRUE" %in% names(tmp.hal) | "Y.1 >= 1TRUE" %in% names(tmp.hal)) {
                print("NB: do not use time-dependent variables with the baseline tmle")
            }

            for (kk in 1:length(fit.hals)) {
                delta.value <- fit.hals[[kk]][["delta.value"]]
                tmp3 <- merge(tmp3[, !(names(tmp3) %in% c(paste0("P", delta.value), paste0("surv", delta.value), paste0("surv", delta.value, ".1"))), with = FALSE],
                              tmp.hal[, names(tmp.hal) %in% c("time", "id", paste0("P", delta.value), paste0("surv", delta.value), paste0("surv", delta.value, ".1")), with = FALSE], by = c("id", "time"))[order(id, time)]
            }

        } else { #-- for the general version of tmle:

            #-- NB: when using HAL, we here go from long to wide format
            #-- NB: remember that there is a faster version which works with "Y.dummy" dependence alone

            if (length(tmp.hal[["Y.dummy >= 1TRUE"]])>0) {
                tmp.hal[["Y.dummy"]] <- tmp.hal[["Y.dummy >= 1TRUE"]]
            }

            Y.dummy.vars <- grep(">=", names(tmp.hal)[substr(names(tmp.hal), 1, 3)  %in% c("Y.d", "Y.1")], value = TRUE)

            if (length(grep("Y.1", Y.dummy.vars))>0) {
                
                Y.dummy.max <- max(as.numeric(gsub("TRUE|Y.1 >= ", "", grep("Y.1", Y.dummy.vars, value = TRUE))))

                index.value <- as.numeric(gsub("Y.dummy >= |Y.1 >= |TRUE", "", Y.dummy.vars))
                
                tmp3[, Y.dummy.index := findInterval(Y.1, 0:Y.dummy.max)]
                
                index.j <- 1:(Y.dummy.max+1)
                index.j1 <- sapply(index.j+1, function(ij) min(ij, max(as.numeric(index.j))))

                for (kY in index.j) {

                    kY.vector <- apply(cbind(
                        do.call("cbind", lapply((1:length(index.value))[index.value %in% index.j[index.j<kY]], function(kY.below) {
                            tmp.hal[[Y.dummy.vars[kY.below]]] == 1
                        })),
                        do.call("cbind", lapply((1:length(index.value))[index.value %in% index.j[index.j >= kY]], function(kY.above) {
                            tmp.hal[[Y.dummy.vars[kY.above]]] == 0
                        }))), 1, prod)
                    
                    tmp.hal[, (paste0("kY.", kY)) := kY.vector]

                }
            }

            for (kk in 1:length(fit.hals)) {

                delta.value <- fit.hals[[kk]][["delta.value"]]

                tmp3 <- merge(tmp3[, !(names(tmp3) %in% c(paste0("P", delta.value), paste0("surv", delta.value), paste0("surv", delta.value, ".1"))), with = FALSE],
                              tmp.hal[observed.Y == 1, names(tmp.hal) %in% c("time", "id", paste0("P", delta.value), paste0("surv", delta.value), paste0("surv", delta.value, ".1")), with = FALSE], by = c("id", "time"))

                if (length(grep("Y.1", Y.dummy.vars))>0) { #-- to handle dependence on Y.1:
                    
                    for (kY in index.j) {
                        
                        tmp.hal.kY <- tmp.hal[get(paste0("kY.", kY)) == 1]
                        tmp.hal.kY[, (paste0("P", delta.value, ".Y", kY)) := get(paste0("P", delta.value))]

                        tmp3 <- merge(tmp3[, !(names(tmp3) %in% paste0("P", delta.value, ".Y", kY)), with = FALSE],
                                      tmp.hal.kY[, c(paste0("P", delta.value, ".Y", kY), "id", "time"), with = FALSE], by = c("id", "time"))

                    }
                    
                } else {

                    if (length(Y.dummy.vars)>0) {
                        
                        tmp.Y1 <- tmp.hal[Y.dummy == 1, c("id", "time",
                                                          paste0("P", delta.value)), with = FALSE]

                        setnames(tmp.Y1, paste0("P", delta.value), paste0("P", delta.value, ".Y1"))
            
                        tmp.Y0 <- tmp.hal[Y.dummy == 0, c("id", "time",
                                                          paste0("P", delta.value)), with = FALSE]

                        setnames(tmp.Y0, paste0("P", delta.value), paste0("P", delta.value, ".Y0"))
                        
                    } else { ##--- when there is no dependence on past of recurrent event process. (this can be improved/simplified)

                        tmp.Y1 <- tmp.hal[, c("id", "time",
                                              paste0("P", delta.value)), with = FALSE]

                        setnames(tmp.Y1, paste0("P", delta.value), paste0("P", delta.value, ".Y1"))
            
                        tmp.Y0 <- tmp.hal[, c("id", "time",
                                              paste0("P", delta.value)), with = FALSE]

                        setnames(tmp.Y0, paste0("P", delta.value), paste0("P", delta.value, ".Y0"))

                    }

                    tmp3 <- merge(tmp3[, !(names(tmp3) %in% paste0("P", delta.value, c(".Y1", ".Y0"))), with = FALSE],
                                  merge(tmp.Y1, tmp.Y0, by = c("id", "time")), by = c("id", "time"))

                }

            }
        }
    }

    if (browse2) browser()
    
    #--------------------------------    
    #-- compute clever weights:
    tmp3[, C := cumsum(1*(time == time.obs & delta == 0 & time %in% dt[["time"]])), by = "id"]
    tmp3[, C.1 := c(0, C[-.N]), by = "id"] 
    tmp3[, probA := predict(fitA, newdata = tmp3, type = "response")]
    tmp3[, clever.weight := (Aobs == get(Aname))/((probA^(Aobs == 1)*(1-probA)^(Aobs == 0)))*(C.1 == 0)/surv0.1]
    
    tmp3[, final.time := max(time.obs), by = "id"]

    if (verbose) {
        print("--------------------------------------------")
        print("clever weights:")
        print(tmp3[time <= final.time, summary(clever.weight)])
        print(tmp3[time <= final.time, summary(surv0)])
        print(tmp3[clever.weight>0 & time <= final.time, summary(clever.weight)])
        print("--------------------------------------------")
    }

    #--------------------------------    
    #-- collect info on positivity issues:
    
    if (tmp3[clever.weight>0 & time <= final.time, max(clever.weight)]>100) {
        positivity.issues <- tmp3[clever.weight>0 & time <= final.time & clever.weight > 100, length(unique(id))]
        print("NB: positivity issues! weights will be truncated")
        tmp3[clever.weight>0 & time <= final.time & clever.weight > 100, clever.weight := 100]
    } else {
        positivity.issues <- 0
    }

    #--------------------------------    
    #-- simple "baseline" version of TMLE:
    
    if (baseline.tmle) {

        for (iter in 1:max.iter) {

            if (verbose) print(paste0("tmle iter = ", iter))

            tmp3[, Z.tau := sum(P1*surv2.1), by = "id"]

            #-- g.est: 
            if (iter == 1) g.est <- mean(tmp3[, sum(P1*surv2.1), by = "id"][[2]])

            #-- clever covariates: 
            tmp3[, clever.D := -c(rev(cumsum(rev(P1*surv2.1))))/surv2, by = "id"]
            tmp3[, clever.Y := 1, by = "id"]

            #-- current estimator for target parameter: 
            target.est <- mean(tmp3[, sum(P1*surv2.1), by = "id"][[2]])

            #-- efficient influence curve:         
            eic <- tmp3[time <= final.time, sum(clever.weight*clever.Y*((delta == 1) - P1)) +
                                            sum(clever.weight*clever.D*((delta == 2) - P2)) +
                                            Z.tau[1]-target.est,
                        by = "id"][[2]]
            if (iter == 1) target.se <- sqrt(mean(eic^2/n))

            #-- check if solved well enough: 
            if (verbose) print(paste0("eic equation solved at = ", abs(mean(eic))))
            if (abs(mean(eic)) <= target.se/(log(n))) break(print(paste0("finished after ", iter, " iterations")))

            #-- otherwise update: 
            target.fun.Y <- function(eps) {
                mean(tmp3[time <= final.time, sum(clever.weight*(clever.Y)*((delta == 1) - P1*exp(eps))), by = "id"][[2]])}
            target.fun.D <- function(eps) {
                mean(tmp3[time <= final.time, sum(clever.weight*(clever.D)*((delta == 2) - P2*exp(eps))), by = "id"][[2]])}

            eps.Y <- nleqslv(0.00, target.fun.Y)$x
            eps.D <- nleqslv(0.00, target.fun.D)$x
            
            if (verbose) print(paste0("eps.Y = ", eps.Y))
            if (verbose) print(paste0("eps.D = ", eps.D))

            tmp3[, P1 := P1*exp(eps.Y)]
            tmp3[, P2 := P2*exp(eps.D)]

            tmp3[, surv2 := exp(-cumsum(P2)), by = "id"] 
            tmp3[, surv2.1 := c(1, surv2[-.N]), by = "id"]
        }

        return(c(g.est = g.est, tmle.est = target.est, tmle.se = target.se, positivity.issues = positivity.issues))
    }

    #--------------------------------    
    #-- "general" version of TMLE: 

    for (iter in 1:max.iter) {

         if (verbose) print(paste0("tmle iter = ", iter))

         tmp3[time == rev(unique.times)[1], Z := P1 + Y.1]

         tmp3[time == rev(unique.times)[1], clever.Z.D.0 := Z] 
         tmp3[time == rev(unique.times)[1], clever.Z.D.1 := Y.1]

         if (length(Y.dummy.max)>0) {
             for (kY in index.j) {
                 tmp3[time == rev(unique.times)[1], (paste0("Z.Y", kY)) := get(paste0("P1.Y", kY)) + Y.1]
             }
         } else {
             tmp3[time == rev(unique.times)[1], Z.Y1 := P1.Y1 + Y.1]
             tmp3[time == rev(unique.times)[1], Z.Y0 := P1.Y0 + Y.1]
         }

        tmp3[time == rev(unique.times)[1], clever.Z.Y.0 := Y.1] 
        tmp3[time == rev(unique.times)[1], clever.Z.Y.1 := Y.1+1]

         tmp3[, Z.next := c(Z[-1], Z[.N]), by = "id"]
         if (length(Y.dummy.max)>0) {
             for (kY in index.j) {
                 tmp3[, (paste0("Z.Y", kY, ".next")) := c(get(paste0("Z.Y", kY))[-1], Y[.N]), by = "id"]
             }
         } else {
             tmp3[, Z.Y1.next := c(Z.Y1[-1], Y[.N]), by = "id"]
             tmp3[, Z.Y0.next := c(Z.Y0[-1], 0), by = "id"]
         }

        ii.count <- 0

         for (ii in length(unique.times):2) {

            ii.count <- ii.count+1

            if (verbose) if (round(ii.count/length(length(unique.times):2)*100) %% 10 == 0 & round((ii.count-1)/length(length(unique.times):2)*100) != round(ii.count/length(length(unique.times):2)*100))
                             print(paste0(round(ii.count/length(length(unique.times):2)*100), "% done with current iterations"))

             if (length(Y.dummy.max)>0) {

                 for (kY in index.j) {

                     #-- compute clever covariates for Y:
                     tmp3[time == unique.times[ii-1] & Y.dummy.index == kY,
                          clever.Z.Y.0 := (get(paste0("Z.Y", kY, ".next")) - (delta == 1))*(1-P2)+Y.1*P2]
                     tmp3[time == unique.times[ii-1] & Y.dummy.index == kY,
                          clever.Z.Y.1 := (get(paste0("Z.Y", index.j1[index.j == kY], ".next")) + 1 - (delta == 1))*(1-P2)+Y.1*P2]

                     tmp3[time == unique.times[ii-1], 
                     (paste0("Z.Y", kY)) := (get(paste0("Z.Y", kY, ".next")) - (delta == 1))*(1-get(paste0("P1.Y", kY))) + 
                         (get(paste0("Z.Y", index.j1[index.j == kY], ".next")) - (delta == 1) + 1)*get(paste0("P1.Y", kY))]
                     
                     #-- compute clever covariates for D:
                     tmp3[time == unique.times[ii-1] & Y.dummy.index == kY,
                          clever.Z.D.0 := get(paste0("Z.Y", kY))] 
                     tmp3[time == unique.times[ii-1] & Y.dummy.index == kY,
                          clever.Z.D.1 := Y.1]
                
                     tmp3[time == unique.times[ii-1],
                     (paste0("Z.Y", kY)) := get(paste0("Z.Y", kY))*(1-get(paste0("P2.Y", kY)))+Y.1*get(paste0("P2.Y", kY))]

                     tmp3[time == unique.times[ii-1] & Y.dummy.index == kY,
                          Z := get(paste0("Z.Y", kY))]

                     tmp3[, Z.next := c(Z[-1], 1), by = "id"]
                     tmp3[, (paste0("Z.Y", kY, ".next")) := c(get(paste0("Z.Y", kY))[-1], Y[.N]), by = "id"]

                 }
                 
             } else {
                              
                 #-- compute clever covariates for Y:
                 tmp3[time == unique.times[ii-1], clever.Z.Y.0 := ((Z.Y0.next - (delta == 1))*(Y.1 == 0) + (Z.Y1.next - (delta == 1))*(Y.1 >= 1))*(1-P2)+Y.1*P2]
                 tmp3[time == unique.times[ii-1], clever.Z.Y.1 := (Z.Y1.next + 1 - (delta == 1))*(1-P2)+Y.1*P2]

                 tmp3[time == unique.times[ii-1], Z.Y0 := (Z.Y0.next - (delta == 1))*(1-P1.Y0)+(Z.Y1.next + 1 - (delta == 1))*P1.Y0] 
                 tmp3[time == unique.times[ii-1], Z.Y1 := (Z.Y1.next - (delta == 1))*(1-P1.Y1)+(Z.Y1.next + 1 - (delta == 1))*P1.Y1] 
                 
                 #-- compute clever covariates for D:
                 tmp3[time == unique.times[ii-1], clever.Z.D.0 := Z.Y1*(Y.1 >= 1) + Z.Y0*(Y.1 == 0)] 
                 tmp3[time == unique.times[ii-1], clever.Z.D.1 := Y.1]
                
                 tmp3[time == unique.times[ii-1], Z.Y1 := (Z.Y1*(1-P2.Y1)+Y.1*P2.Y1)]
                 tmp3[time == unique.times[ii-1], Z.Y0 := (Z.Y0*(1-P2.Y0)+Y.1*P2.Y0)]

                 tmp3[time == unique.times[ii-1], Z := Z.Y1*(Y.1 >= 1) + Z.Y0*(Y.1 == 0)]

                 tmp3[, Z.next := c(Z[-1], 1), by = "id"]
                 tmp3[, Z.Y1.next := c(Z.Y1[-1], 1), by = "id"]
                 tmp3[, Z.Y0.next := c(Z.Y0[-1], 1), by = "id"]
                 
             }      
         }
        
        #-- g.est: 
        if (iter == 1) g.est <- tmp3[time == unique.times[1], mean(Z)]
        
        #-- current estimator for target parameter:  
        target.est <- tmp3[time == unique.times[1], mean(Z)]

        #-- efficient influence curve:         
        eic <- tmp3[time <= final.time, sum(clever.weight*(clever.Z.Y.1-clever.Z.Y.0)*((delta == 1) - P1)) +
                                        sum(clever.weight*(clever.Z.D.1-clever.Z.D.0)*((delta == 2) - P2)) +
                                        Z[1]-target.est,
                    by = "id"][[2]]
        if (iter == 1) target.se <- sqrt(mean(eic^2/n))

        #-- check if solved well enough: 
        if (verbose) print(paste0("eic equation solved at = ", abs(mean(eic))))
        if (tmle.criterion == 1) {
            if (abs(mean(eic)) <= target.se/(log(n))) break(print(paste0("finished after ", iter, " iterations")))
        } else {
            if (abs(mean(eic)) <= target.se/(sqrt(n)*log(n))) break(print(paste0("finished after ", iter, " iterations")))
        }

         #-- otherwise update: 
         target.fun.Y <- function(eps) {
             mean(tmp3[time <= final.time, sum(clever.weight*(clever.Z.Y.1-clever.Z.Y.0)*((delta == 1) - P1*exp(eps))), by = "id"][[2]])}
         target.fun.D <- function(eps) {
             mean(tmp3[time <= final.time, sum(clever.weight*(clever.Z.D.1-clever.Z.D.0)*((delta == 2) - P2*exp(eps))), by = "id"][[2]])}

         eps.Y <- nleqslv(0.00, target.fun.Y)$x
         eps.D <- nleqslv(0.00, target.fun.D)$x
         
         if (verbose) print(paste0("eps.Y = ", eps.Y))
         if (verbose) print(paste0("eps.D = ", eps.D))
         
         tmp3[, P1 := P1*exp(eps.Y)]
         tmp3[, P2 := P2*exp(eps.D)]

         if (length(Y.dummy.max)>0) {

             for (kY in index.j) {
                 tmp3[, (paste0("P1.Y", kY)) := get(paste0("P1.Y", kY))*exp(eps.Y)]
                 tmp3[, (paste0("P2.Y", kY)) := get(paste0("P2.Y", kY))*exp(eps.D)]
             }

         } else {
             
             tmp3[, P1.Y1 := P1.Y1*exp(eps.Y)]
             tmp3[, P1.Y0 := P1.Y0*exp(eps.Y)]

             tmp3[, P2.Y1 := P2.Y1*exp(eps.D)]
             tmp3[, P2.Y0 := P2.Y0*exp(eps.D)]
             
         }
        
    }

    if (output.lambda.cvs) {
        return(c(g.est = g.est, tmle.est = target.est, tmle.se = target.se, positivity.issues = positivity.issues, lambda.cvs = lambda.cvs, hal.coefs = hal.coefs))
    } else {
        return(c(g.est = g.est, tmle.est = target.est, tmle.se = target.se, positivity.issues = positivity.issues))
    }
    
}


######################################################################
### tmle.estimation.fun.R ends here

