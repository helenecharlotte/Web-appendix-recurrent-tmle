### fit.hal.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Oct 15 2024 (09:34) 
## Version: 
## Last-Updated: Mar 27 2025 (18:23) 
##           By: Helene
##     Update #: 111
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

fit.hal <- function(seed = NULL,#13349,
                    pseudo.dt, 
                    delta.var = "delta",
                    delta.value = 1,
                    time.var = "time",
                    covars,
                    treatment = "Aobs", #- the *observed* treatment variable
                    treatment.prediction = "A", #- the *counterfactual* treatment variable
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
                    stime = 0,
                    scovar = 0,
                    cv.glmnet = FALSE,
                    V = 5,
                    reduce.seed.dependence = FALSE,
                    lambda.cvs = c(sapply(1:5, function(jjj) (9:1)/(10^jjj))),
                    verbose = FALSE
                    ) {

    if (length(seed)>0) set.seed(seed)

    if (length(unique(pseudo.dt[[treatment.prediction]]))>1) {
        pseudo.subset <- pseudo.dt[[treatment]] == pseudo.dt[[treatment.prediction]]
    } else {
        pseudo.subset <- rep(TRUE, nrow(pseudo.dt))
    }

    if (!"observed.Y" %in% names(pseudo.dt)) {
        pseudo.dt[, observed.Y := 1]
    }

    #-- basis matrix for hal: 
    X.hal <- basis.fun(pseudo.dt = pseudo.dt[pseudo.subset & pseudo.dt[["observed.Y"]] == 1],
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
                       stime = stime,
                       scovar = scovar)

    #-- compute risk time and number of events for each grid period:
    hal.pseudo.dt <- X.hal$hal.pseudo.dt
    X.hal <- X.hal$X
    
    hal.pseudo.dt[, D:=sum(ddd), by="x"]
    hal.pseudo.dt[, RT:=sum(risk.time), by="x"]

    #-- penalty.factor; penalize treatment and time indicators, or not:
    penalty.factor <- rep(1, ncol(X.hal))

    if (!penalize.treatment & length(grep(treatment, colnames(X.hal)))>0) { #-- treatment:
        not.penalize <- grep(treatment, colnames(X.hal), value=TRUE)
        if (length(grep(":", not.penalize))>0) not.penalize <- not.penalize[-grep(":", not.penalize)]
        penalty.factor[colnames(X.hal) %in% not.penalize] <- 0
    }

    if (!penalize.time & length(grep("grid.time", colnames(X.hal)))>0) {  #-- time:
        not.penalize <- grep("grid.time", colnames(X.hal), value=TRUE)
        not.penalize <- not.penalize[!(not.penalize %in% grep("I\\(", colnames(X.hal), value=TRUE))]
        if (length(grep(":", not.penalize))>0) not.penalize <- not.penalize[-grep(":", not.penalize)]
        penalty.factor[colnames(X.hal) %in% not.penalize] <- 0
    }    
   
    #-- data as input to glmnet:
    tmp.dt <- unique(hal.pseudo.dt[, c("x", "D", "RT")])
    Y2.hal <- tmp.dt[RT>0, D]
    offset2 <- tmp.dt[RT>0, log(RT)]
    X2.hal <- unique.matrix(X.hal)[tmp.dt$RT>0,]

    #-- fit glmnet: 
    if (cv.glmnet) {
        
        hal.fit <- cv.glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                             offset=offset2,
                             family="poisson",
                             penalty.factor=penalty.factor,
                             maxit=10000)

        lambda.cv <- hal.fit$lambda.min
        
    } else {

        if (length(lambda.cvs) == 0) {
           
            hal.fit <- glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                              offset=offset2,
                              family="poisson",
                              penalty.factor=penalty.factor,
                              maxit=10000)

            lambda.cvs <- hal.fit$lambda

            cve.hal <- cv.fun(loss.fun = lebesgue.loss.fun,
                              hal.pseudo.dt = hal.pseudo.dt,
                              X = X.hal,
                              penalty.factor = penalty.factor,
                              V = V,
                              lambda.cvs = lambda.cvs,
                              delta.var = delta.var,
                              delta.value = delta.value,
                              seed = if(is.null(seed)) NULL else {seed+5843})
        
            lambda.cv <- cve.hal$min$lambda.cv
            
        } else {

            if (reduce.seed.dependence>1) {
                cve.hal.list <- list()
                for (mm in 1:reduce.seed.dependence) {
                    cve.hal.list[[mm]] <- cv.fun(loss.fun = lebesgue.loss.fun,
                                                 hal.pseudo.dt = hal.pseudo.dt,
                                                 X = X.hal,
                                                 penalty.factor = penalty.factor,
                                                 V = V,
                                                 lambda.cvs = lambda.cvs,
                                                 delta.var = delta.var,
                                                 delta.value = delta.value,
                                                 seed = if(is.null(seed)) NULL else {seed+5843})

                    
                }
                cve <- rowMeans(do.call("cbind", lapply(cve.hal.list, function(xlist) xlist$all[, "cve"])))
                lambda.cv <- cve.hal.list[[1]]$all[, "lambda"][cve == min(cve)]
            } else {
                cve.hal <- cv.fun(loss.fun = lebesgue.loss.fun,
                                  hal.pseudo.dt = hal.pseudo.dt,
                                  X = X.hal,
                                  penalty.factor = penalty.factor,
                                  V = V,
                                  lambda.cvs = lambda.cvs,
                                  delta.var = delta.var,
                                  delta.value = delta.value,
                                  seed = if(is.null(seed)) NULL else {seed+5843})
                lambda.cv <- cve.hal$min$lambda.cv
            }

            xxx <- 0 
            tryCatch(
                hal.fit <- glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                                  offset=offset2,
                                  lambda=lambda.cv, 
                                  family="poisson",
                                  penalty.factor=penalty.factor,
                                  maxit=10000),
                warning = function(w) { xxx <<- xxx + 1 }
            )

            if (xxx == 1) {
                print("NB: needed more iterations for glmnet to converge")
                tryCatch(
                    hal.fit <- glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                                      offset=offset2,
                                      lambda=lambda.cv, 
                                      family="poisson",
                                      penalty.factor=penalty.factor,
                                      maxit=100000),
                    warning = function(w) { xxx <<- xxx + 1 }
                )            
            }

            if (xxx == 2) {
                print("NB: try glmnet with all lambdas")
                tryCatch(
                    hal.fit <- glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                                      offset=offset2,
                                      lambda=lambda.cvs, 
                                      family="poisson",
                                      penalty.factor=penalty.factor,
                                      maxit=100000),
                    warning = function(w) { xxx <<- xxx + 1 }
                )            
            }

            if (xxx == 3) {
            
                print("NB: glmnet lack of fit not fixed; try cv.glmnet")
            
                hal.fit <- cv.glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                                     offset=offset2,
                                     family="poisson",
                                     penalty.factor=penalty.factor,
                                     maxit=10000)
            
                lambda.cv <- hal.fit$lambda.min
            
            } 

        }

    }

    if (verbose) {
        print("--------------------------------------------")
        print(paste0(delta.var, " = ", delta.value))
        print("--------------------------------------------")
        print(coef(hal.fit, s=lambda.cv))
    }

    return(list(delta.value = delta.value, hal.fit = hal.fit, lambda.cv = lambda.cv))

}

######################################################################
### fit.hal.R ends here
