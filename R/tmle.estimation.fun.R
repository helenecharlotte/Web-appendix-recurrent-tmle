### tmle.estimation.fun.R --- 
#----------------------------------------------------------------------

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
                         use.exponential = FALSE,
                         verbose = FALSE, verbose2 = FALSE,
                         browse.hal = FALSE,
                         browse2 = FALSE,
                         browse3 = FALSE,
                         #--- simpler estimators: 
                         baseline.tmle = FALSE,
                         naive.gcomp = FALSE, 
                         standard.np = FALSE,
                         max.iter = 100,
                         #--- HAL parameters:
                         lambda.cvs = c(sapply(1:5, function(jjj) (9:1)/(10^jjj))),
                         min.no.of.ones = 0.01,
                         cut.one.way = 15,
                         cut.time.varying = 5,
                         cut.time = 35,
                         cut.two.way = 5,
                         penalize.time = FALSE,
                         reduce.seed.dependence = FALSE,
                         event.dependent.cv = FALSE,
                         cv.hal.fit = FALSE,
                         parallelize.Z = 1,
                         parallelize.cve = parallelize.Z,
                         parallelize.predict = parallelize.Z,
                         tmle.criterion = 1,
                         V = 10, 
                         cv.glmnet = FALSE,
                         seed.hal = NULL,
                         output.lambda.cvs = FALSE,
                         return.eic = FALSE
                         ) {

    dt <- copy(dt)

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

    if (length(interaction.names <- grep(":", varnames, value = TRUE))>0) {
        varnames <- unique(c(varnames[!varnames %in% interaction.names],
                             str_split(interaction.names, ":")[[1]]))
    }

    if ("Y.dummy." %in% substr(varnames,1,nchar("Y.dummy."))) {
        Y.dummy.max <- max(as.numeric(gsub("Y.dummy.", "", varnames[substr(varnames,1,nchar("Y.dummy.")) == "Y.dummy."])))
    } else {
        if ("Y.dummy" %in% substr(varnames,1,nchar("Y.dummy"))) {
            Y.dummy.max <- 1
            print(Y.dummy.max)
        } else {
            Y.dummy.max <- NULL
        }
    }

    varnames <- varnames[varnames != "1"]
    varnames <- varnames[!(varnames %in% grep(">=", varnames, value = TRUE))]

    #-- name of treatment variable: 
    Aname <- strsplit(fit.treatment[["model"]], "~")[[1]][1]
  
    if (typeof(dt[["time"]]) != "double") {
        warning("NB: the time variable is not numeric - will be converted")
        dt[, time := as.numeric(time)]
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

        if ("Y.dummy" %in% covars) covars[covars == "Y.dummy"] <- "Y.1"
        
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

        tmp.hal <- copy(tmp.hal)

        tmp.hal[, tstart := c(0, time[-.N]), by = "id"]
        tmp.hal[, final.time := max(time.obs), by = "id"]

        time.varying.covars <- covars[covars %in% "Y.1"]
        baseline.covars <- covars[!covars %in% "Y.1"]

        tmp.hal[, (baseline.covars):=lapply(.SD, function(x) {
            if (is.character(x)) {
                return(as.numeric(as.factor(x)))
            } else if (!is.numeric(x)) {
                return(as.numeric(x))
            } else return(x)
        }), .SDcols=baseline.covars]

        grid.times <- unique.times[floor(seq(1, length(unique.times), length=cut.time+1))[-(cut.time+1)]]

        grid.covars <- lapply(baseline.covars, function(covar) {
            sort(unique(tmp.hal[idN == 1, covar, with = FALSE][floor(seq(1, n, length=cut.one.way))])[[1]])
        })
        names(grid.covars) <- baseline.covars
     
        tmp.hal[, grid.period:=as.numeric(cut(tstart, c(0, grid.times, Inf), include.lowest=TRUE, right=FALSE))]
        tmp.hal[, grid.time:=c(0, grid.times, Inf)[grid.period]]

        by.vars <- c("id", "grid.period")
        if (length(time.varying.covars)>0) by.vars <- c(by.vars, time.varying.covars)

        for (delta.value in c(1,2,0)[any.hal]) {
            if (delta.value == 0) {
                tmp.hal[, (paste0("ddd.", delta.value)) := sum(time==final.time & delta==delta.value), by = by.vars]
            } else {
                tmp.hal[, (paste0("ddd.", delta.value)) := sum((time<=final.time)*(delta==delta.value)), by = by.vars]
            }

            ## tmp.hal[, risk.time := sum(time - tstart), by = by.vars]
            tmp.hal[time <= final.time, risk.time := sum(time - tstart), by = by.vars]
        }

        tmp.hal.reduced <- unique(tmp.hal[time<=final.time], by=by.vars)

        ## tmp.hal.reduced[grid.time == 0, time:=0]
        ## tmp.hal.reduced[, risk.time := diff(c(time, final.time[.N])), by = "id"]

        hal.formula.main <- paste0("delta ~ -1 + Aobs + ",
                                   paste0("(grid.time >=", grid.times, ")", collapse = "+"),
                                   " + ", paste0(sapply(baseline.covars, function(covar)
                                       paste0("(", covar, ">=", grid.covars[[covar]], ")", collapse = "+")), collapse = "+"),
                                   ifelse(length(time.varying.covars)>0, paste0(" + ", sapply(time.varying.covars, function(covar)
                                       paste0("(", covar, ">=", 1:cut.time.varying, ")", collapse = "+"))), ""))

        if (!(cut.two.way == 0)) {
                       
            screening.covars.interaction <- lapply(1:length(grid.covars), function(cc) paste0("(", names(grid.covars)[cc], " >= ", min(grid.covars[[cc]][grid.covars[[cc]] >= median(grid.covars[[cc]])]), ")"))
            screening.times.interaction <- paste0("(grid.time >= ", median(grid.times), ")")

            screening.two.way <- do.call("rbind", combn(c(screening.covars.interaction,
                                                          screening.times.interaction,
                                                          "Aobs",
                                                          "(Y.1>= 1)"), 2, simplify = FALSE))

            hal.formula.main.without.two.way <- hal.formula.main

            hal.formula.main <- paste0(hal.formula.main,
                                       " + ", paste0(screening.two.way[,1], ":", screening.two.way[,2], collapse = "+"))
        }

        X <- Matrix(
            model.matrix(formula(hal.formula.main),
                         data = tmp.hal.reduced), sparse = TRUE)
        
        col_ones <- Matrix::colMeans(X != 0)
        X <- X[, col_ones >= min.no.of.ones]
        
        x.vector <- hash_sparse_rows_dgC(X)
        tmp.hal.reduced[, x:=x.vector]

        if (browse.hal) browser()

        if (length(seed.hal) == 0) seed.hal <- sample(1e8, 1)

        ## print(paste0("seed.hal = ", seed.hal))

        fit.hals <- lapply(c(1,2,0)[any.hal], function(delta.value) {
           
            first.hal <- fit.hal(
                hal.pseudo.dt = tmp.hal.reduced, 
                X.hal = X, 
                delta.var = "delta",
                delta.value = delta.value,
                time.var = "time",
                treatment = "Aobs",
                lambda.cvs = lambda.cvs,
                penalize.time = !(cut.two.way == 0),
                penalize.treatment = !(cut.two.way == 0),
                event.dependent.cv = event.dependent.cv,
                reduce.seed.dependence = reduce.seed.dependence,
                V = V,
                parallelize.cve = parallelize.cve,
                cv.glmnet = cv.glmnet,
                verbose = verbose,
                seed = seed.hal)

            if (cut.two.way == 0) {

                return(first.hal)

            } else {

                indicator.names <- coef(first.hal[["hal.fit"]], s=first.hal[["lambda.cv"]])@Dimnames[[1]]
                nonzero.vars <- coef(first.hal[["hal.fit"]], s=first.hal[["lambda.cv"]])@i
                nonzero.coefs <- coef(first.hal[["hal.fit"]], s=first.hal[["lambda.cv"]])@x
                selected.bases <- indicator.names[nonzero.vars+1]
                selected.bases <- selected.bases[rev(order(abs(nonzero.coefs)))][1:min(length(selected.bases), 50)]
                selected.bases <- selected.bases[!selected.bases %in% "(Intercept)" &
                                                 !selected.bases %in% grep("FALSE", selected.bases, value = TRUE)]
                selected.bases <- gsub("TRUE", "", selected.bases)

                two.way <- lapply(str_split(selected.two.way.bases <- grep(":", selected.bases, value = TRUE), ":"),
                       function(str_interaction) {
                           c(str_split(str_interaction[1], " >= ")[[1]][1],
                             str_split(str_interaction[2], " >= ")[[1]][1])
                       })

                selected.one.way.bases <- selected.bases[!selected.bases %in% selected.two.way.bases]

                if (length(two.way)>0) {
                    two.way <- do.call("rbind", lapply(two.way, function(two.ways) {
                        if (length(grep(two.ways[1], selected.one.way.bases))>0 & length(grep(two.ways[2], selected.one.way.bases))>0)  
                            do.call("rbind", lapply(grep(two.ways[1], selected.one.way.bases, value = TRUE),
                                                    function(var1) 
                                                        cbind(var1, grep(two.ways[2], selected.one.way.bases, value = TRUE))))
                    }))
                }

                if (length(two.way)>0) {
                    two.way.bases <- unique(c(
                        paste0("(", gsub(":", "):(", selected.two.way.bases), ")"),
                        apply(two.way, 1,
                              function(row) paste0("(", row[1], "):(", row[2], ")"))))
                }
                
                hal.formula.two.way <-
                    paste0(hal.formula.main.without.two.way,
                           ifelse(length(two.way)>0, paste0(" + ", paste0(two.way.bases, collapse = "+")), ""))

                X.two.way <- Matrix(
                    model.matrix(formula(hal.formula.two.way),
                                 data = tmp.hal.reduced), sparse = TRUE)
            
                col_ones <- Matrix::colMeans(X.two.way != 0)
                X.two.way <- X.two.way[, col_ones >= min.no.of.ones]

                return(fit.hal(
                    hal.pseudo.dt = tmp.hal.reduced, 
                    X.hal = X.two.way, 
                    delta.var = "delta",
                    delta.value = delta.value,
                    time.var = "time",
                    treatment = "Aobs",
                    lambda.cvs = lambda.cvs,
                    penalize.time = penalize.time,
                    reduce.seed.dependence = reduce.seed.dependence,
                    event.dependent.cv = event.dependent.cv,
                    V = V,
                    parallelize.cve = parallelize.cve,
                    cv.glmnet = cv.glmnet,
                    verbose = verbose,
                    seed = seed.hal))
            }
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

        if (browse.hal) browser()

        tmp.hal[, risk.time := time-tstart]

        tmp.hal <- predict.hal(
            fit.hals = fit.hals,
            pseudo.dt = tmp.hal,
            delta.var = "delta",
            time.var = "time",
            treatment = "Aobs",
            treatment.prediction = Aname,
            parallelize.predict = parallelize.predict,
            cv.fit = cv.hal.fit,
            verbose = verbose,
            verbose2 = verbose2,
            seed = seed.hal)
    }

    #--------------------------------    
    #-- prediction part is over
    
    tmp3 <- tmp3[time <= tau]
    if (any(any.hal)) tmp.hal <- tmp.hal[time <= tau]
    
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
            
            if ("P1.Y" %in% substr(names(tmp.hal), 1, 4) | "Y.dummy >= 1TRUE" %in% names(tmp.hal) | "Y.1 >= 1TRUE" %in% names(tmp.hal)) {
                print("NB: do not use time-dependent variables with the baseline tmle")
            }

            for (kk in 1:length(fit.hals)) {
                delta.value <- fit.hals[[kk]][["delta.value"]]
                tmp3 <- merge(tmp3[, !(names(tmp3) %in% c(paste0("P", delta.value), paste0("surv", delta.value), paste0("surv", delta.value, ".1"))), with = FALSE],
                              tmp.hal[, names(tmp.hal) %in% c("time", "id", paste0("P", delta.value), paste0("surv", delta.value), paste0("surv", delta.value, ".1")), with = FALSE], by = c("id", "time"))[order(id, time)]
            }

        } else { #-- for the general version of tmle:

            index.j <- as.numeric(gsub("P1.Y", "", names(tmp.hal)[substr(names(tmp.hal),1,4) == "P1.Y"]))

            for (kk in 1:length(fit.hals)) {

                delta.value <- fit.hals[[kk]][["delta.value"]]

                for (kY in index.j) {
                    tmp.hal[Y.1 >= kY-1, (paste0("P", delta.value)) := get(paste0("P", delta.value, ".Y", kY))]
                    ##tmp.hal[which.Y == kY, (paste0("P", delta.value)) := get(paste0("P", delta.value, ".Y", kY))]
                }

                tmp.hal[, (paste0("surv", delta.value)) := exp(-cumsum(get((paste0("P", delta.value))))), by = "id"]
                #tmp.hal[, (paste0("surv", delta.value, ".1")) := exp(-cumsum(get((paste0("P", delta.value))))), by = "id"]
                tmp.hal[, (paste0("surv", delta.value, ".1")) := c(1,get(paste0("surv", delta.value))[-.N]), by = "id"]
                #tmp.hal[, (paste0("surv", delta.value, ".1.test2")) := get(paste0("surv", delta.value, ".1")), by = "id"]

                tmp3 <- merge(tmp3[, !(names(tmp3) %in% c(paste0("P", delta.value), paste0("P", delta.value, ".Y", index.j), paste0("surv", delta.value), paste0("surv", delta.value, ".1"))), with = FALSE],
                              tmp.hal[, names(tmp.hal) %in% c("time", "id", paste0("P", delta.value), paste0("P", delta.value, ".Y", index.j), paste0("surv", delta.value), paste0("surv", delta.value, ".1")), with = FALSE], by = c("id", "time"))
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
    # browser()
    # Number of Y-history states
    ### browser()
    if (any(any.hal)) {
        index.j <- sort(as.numeric(gsub("P1.Y", "", names(tmp.hal)[substr(names(tmp.hal),1,4) == "P1.Y"])))
    } else {
        index.j <- sort(as.numeric(gsub("P1.Y", "", names(tmp3)[substr(names(tmp3),1,4) == "P1.Y"])))
    }

    if (min(index.j) == 0) {
        index.j <- index.j+1
        for (kY in rev(index.j)) {
            setnames(tmp3, paste0("P1.Y", kY-1), paste0("P1.Y", kY))
            setnames(tmp3, paste0("P2.Y", kY-1), paste0("P2.Y", kY))
 
        }
    }

    if (use.exponential) {

        for (kY in index.j) {
            tmp3[, (paste0("P1.Y", kY)) := 1-exp(-tmp3[[paste0("P1.Y", kY)]])]
            tmp3[, (paste0("P2.Y", kY)) := 1-exp(-tmp3[[paste0("P2.Y", kY)]])]
        }

        tmp3[, P1 := 1-exp(-tmp3[["P1"]])]
        tmp3[, P2 := 1-exp(-tmp3[["P2"]])]
                           
    }

    if (verbose) print(index.j)
    if (verbose2) print(summary(tmp3))

    for (kY in index.j) {
        tmp3[Y.1 >= kY-1, Y.dummy.index := kY]
    }

    K <- max(index.j)
    states <- 1:K
    
    if (browse3) browser()

    ## browser()

    for (iter in 1:max.iter) {

        if (verbose) print(paste0("tmle iter = ", iter))

        setkey(tmp3, id, time)
        dt_list <- split(tmp3, by="id", keep.by=TRUE)

        t2 <- system.time({
            dt_list <- mclapply(
                dt_list,
                compute_Z_and_clever_per_id,
                states = states,
                mc.cores = min(detectCores()-1, parallelize.Z)
            )
        })

        test.error <- try(tmp3 <- rbindlist(dt_list))

        if (any(class(test.error) == "try-error")) {
            t2 <- system.time({
                dt_list <- mclapply(
                    dt_list,
                    compute_Z_and_clever_per_id,
                    states = states,
                    mc.cores = 1
                )
            })
        }

        if (verbose) print(t2)
        
        #-- current estimator for target parameter:  
        target.est <- mean(tmp3[, Z[1], by = "id"][[2]])

        #-- g.est: 
        if (iter == 1) g.est <- target.est

        #-- efficient influence curve:         
        eic <- tmp3[time <= final.time, sum(clever.weight*(clever.Z.Y.1-clever.Z.Y.0)*((delta == 1) - P1)) +
                                        sum(clever.weight*(clever.Z.D.1-clever.Z.D.0)*((delta == 2) - P2)) +
                                        Z[1]-target.est,
                    by = "id"][[2]]
        if (iter == 1) {
            target.se <- sqrt(mean(eic^2/n))
            if (return.eic) {
                eic.out <- eic
            }
        }

        #-- check if solved well enough: 
        if (verbose) print(paste0("eic equation solved at = ", abs(mean(eic))))
        if (tmle.criterion == 1) {
            if (abs(mean(eic)) <= target.se/(log(n))) break(print(paste0("finished after ", iter, " iterations")))
        } else {
            if (abs(mean(eic)) <= target.se/(sqrt(n)*log(n))) break(print(paste0("finished after ", iter, " iterations")))
        }

        t3 <- system.time({

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

            for (kY in index.j) {
                tmp3[, (paste0("P1.Y", kY)) := get(paste0("P1.Y", kY))*exp(eps.Y)]
                tmp3[, (paste0("P2.Y", kY)) := get(paste0("P2.Y", kY))*exp(eps.D)]
            }
        })

        if (verbose) print(t3)
        
    }

    if (return.eic) {
        if (output.lambda.cvs) {
            return(list(c(g.est = g.est, tmle.est = target.est, tmle.se = target.se, positivity.issues = positivity.issues, lambda.cvs = lambda.cvs, hal.coefs = hal.coefs),
                        eic = eic.out))
        } else {
            return(list(c(g.est = g.est, tmle.est = target.est, tmle.se = target.se, positivity.issues = positivity.issues),
                        eic = eic.out))
        }
    } else {
        if (output.lambda.cvs) {
            return(c(g.est = g.est, tmle.est = target.est, tmle.se = target.se, positivity.issues = positivity.issues, lambda.cvs = lambda.cvs, hal.coefs = hal.coefs))
        } else {
            return(c(g.est = g.est, tmle.est = target.est, tmle.se = target.se, positivity.issues = positivity.issues))
        }
    }
    
}


library(data.table)

# -----------------------------------------------------------------------------
# Function: compute_Z_and_clever_per_id
# Purpose : For one individual, run backward matrix recursion over K states
# -----------------------------------------------------------------------------

compute_Z_and_clever_per_id <- function(dt_id,
                                        states,
                                        P1.col.pref = "P1.Y",
                                        P2.col.pref = "P2.Y",
                                        ycount.col = "Y.1",
                                        delta.col = "delta",
                                        state.idx.col = "Y.dummy.index") {


    K <- length(states)
    Tn <- nrow(dt_id)

    P1_cols <- paste0(P1.col.pref, states)    # "P1.Y1","P1.Y2",...
    P2_cols <- paste0(P2.col.pref, states)

    P1_mat <- as.matrix(dt_id[, ..P1_cols])   # Tn x K
    P2_mat <- as.matrix(dt_id[, ..P2_cols])   # Tn x K

    # containers for Z vectors per time (each is length K)
    Z_by_time <- vector("list", Tn)          # will store numeric vectors length K
    Znext_by_time <- vector("list", Tn)      # store Z_next used for clever covariates (optional)
    # output columns per row (these are scalars per row)
    Z_row <- numeric(Tn)
    clever_Y0_row <- numeric(Tn)
    clever_Y1_row <- numeric(Tn)
    clever_D0_row <- numeric(Tn)
    clever_D1_row <- numeric(Tn)   

    # Terminal time initialization (last row) 
    last_idx <- Tn
    P1_last <- P1_mat[last_idx, ]   # vector length K
    Z_T <- as.numeric(P1_last)
    Z_by_time[[last_idx]] <- Z_T
    Znext_by_time[[last_idx]] <- Z_T
  
    # Backwards recursion
    # Loop tt from Tn-1 down to 1; for each tt we compute Z_t (vector length K)
    for (tt in (Tn - 1):1) {

        # Next-step vector
        Z_next <- Z_by_time[[tt + 1]]        # length K
        Z_next_j1 <- c(Z_next[-1], Z_next[K])
        
        # Extract P1 and P2 for this time (vectors length K)
        P1_vec <- as.numeric(P1_mat[tt, ])
        P2_vec <- as.numeric(P2_mat[tt, ])
        
        ## Z_t <- (1-P1_vec)*(1-P2_vec) * Z_next + P1_vec*(1-P2_vec) * Z_next_j1 + (1-P2_vec)*P1_vec
        Z_t <- (1-P1_vec-P2_vec) * Z_next + P1_vec * Z_next_j1 + P1_vec
        
        Z_by_time[[tt]] <- Z_t
        Znext_by_time[[tt]] <- Z_next
    }

    for (tt in seq_len(Tn)) {
        state_pos <- as.integer(dt_id[[state.idx.col]][tt])  

        # Current Z for the observed state
        Z_row[tt] <- Z_by_time[[tt]][state_pos]

        # clever.Z.D.1 is 0
        clever_D1_row[tt] <- 0

        # clever.Z.D.0 := Z at the state
        clever_D0_row[tt] <- Z_row[tt]

        # clever.Z.Y.* come from Z_next
        Znext_tt <- Znext_by_time[[tt]]
        # Z.Y0: use Z_next for same state
        clever_Y0_row[tt] <- Znext_tt[state_pos]
        # Z.Y1: use Z_next of next state (cap at last)
        next_state_pos <- min(state_pos + 1, K)
        clever_Y1_row[tt] <- (Znext_tt[next_state_pos] + 1) 
    }

    out_dt <- cbind(
        data.table(id = dt_id$id,
                   time = dt_id$time,
                   final.time = dt_id$final.time,
                   clever.weight = dt_id$clever.weight,
                   delta = dt_id$delta,
                   P1 = dt_id$P1,
                   P2 = dt_id$P2,
                   Y.dummy.index = dt_id$Y.dummy.index,
                   Z = Z_row,
                   clever.Z.Y.0 = clever_Y0_row,
                   clever.Z.Y.1 = clever_Y1_row,
                   clever.Z.D.0 = clever_D0_row,
                   clever.Z.D.1 = clever_D1_row),
        P1_mat, P2_mat)[]

    # Keep original row order
    setorder(out_dt, time)
    return(out_dt)
}



######################################################################
### tmle.estimation.fun.R ends here

library(digest)

hash_sparse_rows_dgC <- function(M) {
    stopifnot(inherits(M, "dgCMatrix"))

    p  <- M@p
    i  <- M@i
    x  <- M@x
    nr <- M@Dim[1]

    # storage: list of integer/values per row
    idx_list <- vector("list", nr)
    val_list <- vector("list", nr)

    # loop over columns, add entries to each row
    for (col in seq_len(M@Dim[2])) {
        start <- p[col] + 1
        end   <- p[col + 1]

        if (end >= start) {
            rows <- i[start:end] + 1        # row indices (1-based)
            vals <- x[start:end]            # values

            # append col index + value to each row's list
            for (k in seq_along(rows)) {
                r <- rows[k]
                idx_list[[r]] <- c(idx_list[[r]], col)
                val_list[[r]] <- c(val_list[[r]], vals[k])
            }
        }
    }

    # hash each row using only the sparse pattern
    out <- character(nr)
    for (r in seq_len(nr)) {
        out[r] <- digest::digest(
                              list(idx_list[[r]], val_list[[r]]),
                              algo = "xxhash64",
                              serialize = TRUE
                          )
    }

    out
}

