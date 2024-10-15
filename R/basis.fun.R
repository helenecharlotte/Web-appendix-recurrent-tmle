### basis.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Oct 15 2024 (10:09) 
## Version: 
## Last-Updated: Oct 15 2024 (10:15) 
##           By: Helene
##     Update #: 12
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
                      predict = FALSE) {

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

        if (delta.value == 0) {
            pseudo.dt[, ddd := sum(get(time.var)==final.time & get(delta.var)==delta.value), by = c("id", "grid.period")]
        } else {
            pseudo.dt[, ddd := sum((get(time.var)<=final.time)*get(delta.var)==delta.value), by = c("id", "grid.period")]
        }
        pseudo.dt[, risk.time := diff(c(0,get(time.var))), by = "id"]
        pseudo.dt[, risk.time := sum((get(time.var)<=final.time)*risk.time), by = c("id", "grid.period")]

        #-- NB: following solution for predict option is *slow*;
        if (!predict) hal.pseudo.dt <- unique(pseudo.dt[get(time.var)<=final.time], by=c("id", "grid.period")) else hal.pseudo.dt <- copy(pseudo.dt)
        hal.pseudo.dt[, grid.time:=grid.times[grid.period]]
        #-- NB: could consider if there is a better way to allow time-varying covariates change where
        #-- they actually do +++ consider if this step makes good sense;
        if (!predict) for (covar in time.varying.covars) hal.pseudo.dt[, (covar):=get(covar)[1], by = c("id", "grid.period")]
        hal.pseudo.dt[, (time.var):=grid.time]
        
    }

    options(na.action="na.pass")
    
    X <- Matrix(
        model.matrix(formula(paste0(
            delta.var, "~-1",
            ifelse(length(treatment)>0, paste0("+", treatment), ""),
            ifelse(length(time.var)>0, paste0("+",
                                              paste0(indicator.basis.fun(hal.pseudo.dt, "grid.period", cut.time), collapse="+")), ""),
            ifelse(length(time.var)>0 & cut.time.treatment>0,
                   paste0("+", paste0(paste0(indicator.basis.fun(hal.pseudo.dt, "grid.period", cut.time.treatment), ":", treatment, "+"), collapse="")), ""),
            ifelse(length(time.var)>0 & cut.time.covar>0,
                   paste0("+", paste0(sapply(
                                   covars, function(covar) paste0(paste0(indicator.basis.fun(hal.pseudo.dt, "grid.period", cut.time.covar), ":", indicator.basis.fun(hal.pseudo.dt, covar, cut.time.covar), "+"), collapse="")), collapse="+"), ""), ""),
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

indicator.basis.fun <- function(pseudo.dt, xvar, xcut, xgrid=NULL, type="obs", return.grid=FALSE, delta.var=NULL, seed=220) {
    set.seed(seed)
    if (length(xgrid)>0 & xvar == "Y.time") {
    } else if (type=="obs") {
        if (xvar == "Y.1") { # fix to remove heavy tail: 
            xvar.values <- 1:xcut
        } else if (FALSE & xvar == "Y.time") { # 
            xvar.values <- pseudo.dt[get(delta.var) == 1 & Y.dummy >= 1, sort(get(xvar))]
        } else {
            xvar.values <- pseudo.dt[, sort(unique(get(xvar)))]
        }
        xvar.pick <- floor(seq(1, length(xvar.values), length=min(xcut, length(xvar.values))))
        xgrid <- xvar.values[xvar.pick][-c(1,xcut)]
    } else {
        xgrid <- round(seq(pseudo.dt[, min(get(xvar))],
                           pseudo.dt[, max(get(xvar))],
                           length=xcut)[-c(1,xcut)], 2)
    }
    if (return.grid) return(c(unique(xgrid))) else return(paste0("(", xvar, ">=", unique(xgrid), ")"))
}

######################################################################
### basis.fun.R ends here
