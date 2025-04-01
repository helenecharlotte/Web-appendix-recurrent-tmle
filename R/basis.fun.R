### basis.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Oct 15 2024 (10:09) 
## Version: 
## Last-Updated: Mar 27 2025 (18:23) 
##           By: Helene
##     Update #: 131
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

basis.fun <- function(pseudo.dt, ## tmp3,
                      delta.var = "delta",
                      delta.value = 1,
                      time.var = "time",
                      covars,
                      treatment = "Aobs",
                      cut.one.way = 10,
                      cut.time.varying = 5,
                      Y.time.grid = NULL,
                      cut.time = 15,
                      two.way = cbind(var1="", var2=""),
                      cut.two.way = 5,
                      cut.time.treatment = NULL,
                      cut.time.covar = NULL,
                      predict = FALSE,
                      stime = 0,
                      scovar = 0) {

    n <- nrow(unique(pseudo.dt, by = "id"))

    time.varying.covars <-
        covars[sapply(covars, function(covar) {
            (nrow(pseudo.dt[, unique(get(covar)), by = "id"])>n)
        })]

    pseudo.dt[, (covars):=lapply(.SD, function(x) {
        if (is.character(x)) {
            return(as.numeric(as.factor(x)))
        } else if (!is.numeric(x)) {
            return(as.numeric(x))
        } else return(x)
    }), .SDcols=covars]

    covars <- covars[!covars %in% time.varying.covars]

    if (length(time.var)>0) {

        pseudo.dt[, final.time := max(time.obs), by = "id"]
        grid.times <- c(0, indicator.basis.fun(pseudo.dt[get(time.var)<=final.time], time.var, cut.time, return.grid=TRUE), Inf)
        pseudo.dt[, grid.period:=as.numeric(cut(get(time.var), grid.times, include.lowest=TRUE, right=FALSE))]

        by.vars <- c("id", "grid.period")
        if (length(time.varying.covars)>0) by.vars <- c(by.vars, time.varying.covars)

        if (delta.value == 0) {
            pseudo.dt[, ddd := sum(get(time.var)==final.time & get(delta.var)==delta.value), by = by.vars]
        } else {
            pseudo.dt[, ddd := sum((get(time.var)<=final.time)*(get(delta.var)==delta.value)), by = by.vars]
        }

        #-- NB: following solution for predict option is *slow*;
        if (!predict) hal.pseudo.dt <- unique(pseudo.dt[get(time.var)<=final.time], by=by.vars) else hal.pseudo.dt <- copy(pseudo.dt)
        hal.pseudo.dt[, grid.time:=grid.times[grid.period]]
        hal.pseudo.dt[grid.time == 0, (time.var):=0]

        hal.pseudo.dt[, risk.time := diff(c(get(time.var), final.time[.N])), by = "id"]

        if (FALSE) {
            hal.pseudo.dt[id == 3, c("time", "grid.time", "risk.time", "Y.dummy", "ddd")] #, "test.risk.time"
            hal.pseudo.dt[id == 3, sum(ddd)]
            dt1[id == 3, sum(delta == 1)]
        }
        
    }

    options(na.action="na.pass")

    X <- Matrix(
        model.matrix(formula(paste0(
            delta.var, "~-1",
            ifelse(length(treatment)>0, paste0("+", treatment), ""),
            ifelse(length(time.var)>0, paste0("+",
                                              paste0(indicator.basis.fun(hal.pseudo.dt, "grid.time", cut.time), collapse="+")), ""),
            ifelse(length(time.var)>0 & cut.time.treatment>0,
                   paste0("+", paste0(paste0(indicator.basis.fun(hal.pseudo.dt, "grid.time", cut.time.treatment), ":", treatment, "+"), collapse="")), ""),
            ifelse(length(time.var)>0 & cut.time.covar>0,
                   paste0("+", paste0(sapply(
                                   covars, function(covar) paste0(paste0(indicator.basis.fun(hal.pseudo.dt, "grid.time", cut.time.covar), ":", indicator.basis.fun(hal.pseudo.dt, covar, cut.time.covar), "+"), collapse="")), collapse="+"), ""), ""),
            ifelse(length(covars)>0 & cut.one.way>2,
                   paste0("+", paste0(sapply(covars, function(covar) paste0(indicator.basis.fun(hal.pseudo.dt, covar, cut.one.way), collapse="+")), collapse="+")),
                   ""),
            ifelse(length(time.varying.covars)>0 & cut.time.varying>2,
                   paste0("+", paste0(sapply(time.varying.covars, function(covar) {
                       paste0(paste0(ifelse(covar == "Y.time" & TRUE, "(Y.dummy>=1):", ""), indicator.basis.fun(pseudo.dt, covar, cut.time.varying, xgrid=Y.time.grid, delta.var = delta.var)), collapse="+")
                   }), collapse="+")),
                   ""),
            ifelse(two.way[1,1]!="" & cut.two.way>2,
                   paste0(apply(two.way, 1, function(row) {
                       if (row[1]!=row[2]) {
                           if (any(row==treatment)) {
                               if (row[1] == treatment) {
                                   if (row[2] %in% time.varying.covars) {
                                       paste0("+", paste0(apply(expand.grid(treatment,
                                                                            paste0(ifelse(row[2] == "Y.time" & TRUE, "(Y.dummy>=1):", ""), indicator.basis.fun(pseudo.dt, row[2], cut.time.varying, xgrid=Y.time.grid, delta.var = delta.var))), 1,
                                                                function(row2) paste0(row2, collapse=":")), collapse="+"))
                          
                                   } else {
                                       paste0("+", paste0(apply(expand.grid(treatment,
                                                                            na.omit(indicator.basis.fun(hal.pseudo.dt, row[2], cut.one.way)[seq(1, cut.one.way-2, length = cut.two.way-2)])), 1,
                                                                function(row2) paste0(row2, collapse=":")), collapse="+"))
                                   }
                               } else {
                                   if (row[1] %in% time.varying.covars) {
                                       paste0("+", paste0(apply(expand.grid(paste0(ifelse(row[1] == "Y.time" & TRUE, "(Y.dummy>=1):", ""), indicator.basis.fun(pseudo.dt, row[1], cut.time.varying, xgrid=Y.time.grid, delta.var = delta.var)),
                                                                            treatment), 1,
                                                                function(row2) paste0(row2, collapse=":")), collapse="+"))
                                   } else {
                                       paste0("+", paste0(apply(expand.grid(na.omit(indicator.basis.fun(hal.pseudo.dt, row[1], cut.one.way)[seq(1, cut.one.way-2, length = cut.two.way-2)]),
                                                                            treatment), 1,
                                                                function(row2) paste0(row2, collapse=":")), collapse="+"))
                                   }
                               }
                           } else if (any(row %in% time.varying.covars)) { 
                               if (row[1] %in% time.varying.covars) {
                                   paste0("+", paste0(apply(expand.grid(paste0(ifelse(row[1] == "Y.time" & TRUE, "(Y.dummy>=1):", ""), indicator.basis.fun(pseudo.dt, row[1], cut.time.varying, xgrid=Y.time.grid, delta.var = delta.var)),
                                                                        na.omit(indicator.basis.fun(hal.pseudo.dt, row[2], cut.one.way)[seq(1, cut.one.way-2, length = cut.two.way-2)])), 1,
                                                            function(row2) paste0(row2, collapse=":")), collapse="+"))
                               } else {
                                   paste0("+", paste0(apply(expand.grid(na.omit(indicator.basis.fun(hal.pseudo.dt, row[1], cut.one.way)[seq(1, cut.one.way-2, length = cut.two.way-2)]),
                                                                        paste0(ifelse(row[2] == "Y.time" & TRUE, "(Y.dummy>=1):", ""), indicator.basis.fun(pseudo.dt, row[2], cut.time.varying, xgrid=Y.time.grid, delta.var = delta.var))), 1,
                                                            function(row2) paste0(row2, collapse=":")), collapse="+"))
                               }
                           } else {
                               paste0("+", paste0(apply(expand.grid(na.omit(indicator.basis.fun(hal.pseudo.dt, row[1], cut.one.way)[seq(1, cut.one.way-2, length = cut.two.way-2)]),
                                                                    na.omit(indicator.basis.fun(hal.pseudo.dt, row[2], cut.one.way)[seq(1, cut.one.way-2, length = cut.two.way-2)])), 1,
                                                        function(row2) paste0(row2, collapse=":")), collapse="+"))
                           }
                       } else return("")
                   }), collapse=""), "")
        )), 
        data=hal.pseudo.dt), sparse=FALSE)

    if (length(grep("Y.dummy >= 1FALSE:Y.time", colnames(X))) > 0) {
        remove.X <- grep("Y.dummy >= 1FALSE:Y.time", colnames(X))
        X <- X[, -remove.X]
    }
    
    if (length(time.var)>0) {
        if (!predict) {
            x.vector <- apply(X, 1, function(x) paste0(x, collapse=","))
            hal.pseudo.dt[, x:=x.vector]
        }
        return(list(X=X, hal.pseudo.dt=hal.pseudo.dt))
    } 

}

indicator.basis.fun <- function(pseudo.dt, xvar, xcut, xgrid=NULL, type="obs", return.grid=FALSE, delta.var=NULL) {
    if (length(xgrid)>0 & xvar == "Y.time") {
    } else if (type=="obs") {
        if (xvar == "Y.1") { # fix to remove heavy tail: 
            xvar.values <- 1:xcut
        } else if (FALSE & xvar == "Y.time") { # 
            xvar.values <- pseudo.dt[before.tau == 1 & get(delta.var) == 1 & Y.dummy >= 1, sort(get(xvar))]
        } else {
            xvar.values <- pseudo.dt[before.tau == 1, sort(unique(get(xvar)))]
        }
        xvar.pick <- floor(seq(1, length(xvar.values), length=min(xcut, length(xvar.values))))
        xgrid <- xvar.values[xvar.pick][-c(1,xcut)]
    } else {
        xgrid <- round(seq(pseudo.dt[before.tau == 1, min(get(xvar))],
                           pseudo.dt[before.tau == 1, max(get(xvar))],
                           length=xcut)[-c(1,xcut)], 2)
    }
    if (return.grid) {
        return(c(unique(xgrid)))
    } else {
        return(paste0("(", xvar, ">=", unique(xgrid), ")"))
    }
}

######################################################################
### basis.fun.R ends here
