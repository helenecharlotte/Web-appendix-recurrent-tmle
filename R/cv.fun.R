### cv.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Oct 15 2024 (10:12) 
## Version: 
## Last-Updated: Mar 20 2025 (09:54) 
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

cv.fun <- function(loss.fun = lebesgue.loss.fun,
                   hal.pseudo.dt,
                   V = 5,
                   seed = NULL,
                   X,
                   time.var = "time",
                   penalty.factor = rep(1, ncol(X)),
                   delta.var = "delta",
                   delta.value = 1,
                   verbose = FALSE,
                   browse = FALSE, 
                   lambda.cvs = c(sapply(1:5, function(jjj) (9:1)/(10^jjj)))) {

    if (length(seed)>0) set.seed(seed)

    if (browse) browser()
    
    #-- split sample: 
    n <- length(hal.pseudo.dt[, unique(id)])
    cv <- sample(rep(1:V, length = n), size = n)
    cv.split <- lapply(1:V, function(vv) (1:n)[cv == vv])

    cve.list <- lapply(1:V, function(vv) {

        test.set <- cv.split[[vv]]
        train.set <- hal.pseudo.dt[, id][!hal.pseudo.dt[, id] %in% test.set]
        dt.train <- hal.pseudo.dt[id %in% train.set]

        dt.train[, D:=sum(ddd), by="x"]
        dt.train[, RT:=sum(risk.time), by="x"]

        tmp.dt <- unique(dt.train[, c("x", "D", "RT")])
        Y.train <- tmp.dt[RT>0, D]
        offset.train <- tmp.dt[RT>0, log(RT)]
        X.train <- unique.matrix(X[hal.pseudo.dt$id %in% train.set,])[tmp.dt$RT>0,]

        if (length(lambda.cvs)>0) {

            train.fit <- try(glmnet(x=as.matrix(X.train), y=Y.train,
                                    offset=offset.train,
                                    family="poisson",
                                    maxit=100000,
                                    penalty.factor=penalty.factor,
                                    lambda=lambda.cvs), silent=TRUE)

            if (class(train.fit)[1]=="try-error") {
                train.fit <- glmnet(x=as.matrix(X.train), y=Y.train,
                                    offset=offset.train,
                                    family="poisson",
                                    maxit=100000,
                                    penalty.factor=penalty.factor,
                                    lambda=lambda.cvs)
            } else if (all(train.fit$lambda==Inf)) {
                train.fit <- glmnet(x=as.matrix(X.train), y=Y.train,
                                    offset=offset.train,
                                    family="poisson",
                                    maxit=100000,
                                    penalty.factor=penalty.factor,
                                    lambda=lambda.cvs)
            }
        }
        
        return(sapply(lambda.cvs, function(lambda.cv) {
            return(loss.fun(train.fit = train.fit,
                            hal.pseudo.dt = hal.pseudo.dt,
                            test.set = test.set,
                            X = X,
                            time.var = time.var, 
                            lambda.cv = lambda.cv,
                            delta.var = delta.var,
                            delta.value = delta.value))
        }))
        
    })

    cve <- unlist(lapply(1:length(lambda.cvs), function(mm) {
        sum(unlist(lapply(cve.list, function(out) out[[mm]])))
    }))
    
    lambda.cv <- min(lambda.cvs[cve==min(cve)])
    
    return(list(min=list(lambda.cv=lambda.cv,
                         cve=unique(cve[cve==min(cve)])),
                all=cbind(lambda=lambda.cvs, cve=cve)))
}

######################################################################
### cv.fun.R ends here
