### fit.hal.R --- 
#----------------------------------------------------------------------

fit.hal <- function(seed = NULL,#13349,
                    hal.pseudo.dt,
                    X.hal, 
                    delta.var = "delta",
                    delta.value = 1,
                    time.var = "time",
                    treatment = "Aobs", #- the *observed* treatment variable
                    penalize.treatment = FALSE,
                    penalize.time = FALSE,
                    cv.glmnet = FALSE,
                    V = 5,
                    browse = FALSE, 
                    parallelize.cve = 1,
                    event.dependent.cv = FALSE,
                    reduce.seed.dependence = FALSE,
                    lambda.cvs = c(sapply(1:5, function(jjj) (9:1)/(10^jjj))),
                    verbose = FALSE
                    ) {

    if (length(seed)>0) set.seed(seed)

    hal.pseudo.dt[, ddd:=hal.pseudo.dt[[paste0("ddd.", delta.value)]]]
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
                              event.dependent.cv = event.dependent.cv,
                              verbose = verbose,
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
                                                 event.dependent.cv = event.dependent.cv,
                                                 verbose = verbose,
                                                 seed = if(is.null(seed)) NULL else {seed+5843+mm*20})

                }
                cve <- rowMeans(do.call("cbind", lapply(cve.hal.list, function(xlist) xlist$all[, "cve"])))
                lambda.cv <- cve.hal.list[[1]]$all[, "lambda"][cve == min(cve)]
            } else {
                cve.hal <- cv.fun(loss.fun = lebesgue.loss.fun,
                                  hal.pseudo.dt = hal.pseudo.dt,
                                  X = X.hal,
                                  penalty.factor = penalty.factor,
                                  V = V,
                                  parallelize.cve = parallelize.cve,
                                  lambda.cvs = lambda.cvs,
                                  delta.var = delta.var,
                                  delta.value = delta.value,
                                  event.dependent.cv = event.dependent.cv,
                                  verbose = verbose,
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
        message("--------------------------------------------")
        message(paste0(delta.var, " = ", delta.value))
        message("--------------------------------------------")
        print(coef(hal.fit, s=lambda.cv))

        if (FALSE) {
            print(coef(cve.hal$train.fit[[1]]$train.fit, s=lambda.cv))
        }
    }
    
    if (reduce.seed.dependence>1) {
        return(list(delta.value = delta.value, hal.fit = hal.fit, lambda.cv = lambda.cv,
                    train.fit = cve.hal.list[[1]]$train.fit))
    } else {
        return(list(delta.value = delta.value, hal.fit = hal.fit, lambda.cv = lambda.cv,
                    train.fit = cve.hal$train.fit))
    }

}

######################################################################
### fit.hal.R ends here
