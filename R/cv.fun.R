### cv.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Oct 15 2024 (10:12) 
## Version: 
## Last-Updated: Dec 10 2025 (10:28) 
##           By: Helene
##     Update #: 147
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
                   parallelize.cve = 1,
                   event.dependent.cv = FALSE,
                   verbose = FALSE,
                   browse = FALSE, 
                   lambda.cvs = c(sapply(1:5, function(jjj) (9:1)/(10^jjj)))) {

    if (length(seed)>0) set.seed(seed)

    ## print(paste0("seed = ", seed))

    if (browse) browser()

    #-- split sample: 
    
    if (event.dependent.cv) {
        
        hal.pseudo.dt[, Ny := max(Y.1), by = "id"][Ny >= 3, Ny := 3]
        hal.pseudo.dt[, Nc := sum(ddd.0), by = "id"]
        hal.pseudo.dt[, Nd := sum(ddd.2), by = "id"]

        sample.events <- unique(hal.pseudo.dt, by = c("id", "Ny", "Nc", "Nd"))
        unique.sample.events <- unique(sample.events, by = c("Ny", "Nc", "Nd"))[, id.unique := 1:.N]
        sample.events <- merge(sample.events, unique.sample.events[, c("Ny", "Nc", "Nd", "id.unique"), with = FALSE], by = c("Ny", "Nc", "Nd"))

        unique(sample.events[, count := .N, by = c("Ny", "Nc", "Nd")], by = c("Ny", "Nc", "Nd"))[, c("Ny", "Nc", "Nd", "id.unique", "count"), with = FALSE]

        sample.events[, random.row := sample(1:nrow(sample.events), nrow(sample.events))]
        sample.events <- sample.events[order(id.unique, random.row)]

        sampV <- rep(sample(1:V, V), length = nrow(sample.events)) ##unlist(lapply(1:ceiling(nrow(sample.events)/V), function(vv) sample(1:V, V)))[1:nrow(sample.events)]
        sample.events[, cv := sampV]
       
        ## for (jj in sample.events[, unique(id.unique)]) {
        ##     sample.events[id.unique == jj, cv := sample(rep(sample(1:V, V), length = nrow(sample.events[id.unique == jj])), size = nrow(sample.events[id.unique == jj]))]
        ## }

        if (verbose&FALSE) {
            print(sample.events[, table(cv)])
            for (vv in 1:V) {
                print(sample.events[cv == vv, table(id.unique)])
                print(sample.events[cv == vv, sum(Nc)])
                print(sample.events[cv == vv, sum(Ny)])
                print(sample.events[cv == vv, sum(Nd)])
            }
        }
        
        cv.split <- lapply(1:V, function(vv) sample.events[cv == vv, id])

    } else {

        n <- length(hal.pseudo.dt[, unique(id)])
        cv <- sample(rep(1:V, length = n), size = n)
        cv.split <- lapply(1:V, function(vv) (1:n)[cv == vv])

    }

    cve.list <- parallel::mclapply(
                              X = 1:V,
                              mc.cores = min(detectCores()-1, parallelize.cve),
                              FUN = function(vv) {

    #cve.list <- lapply(1:V, function(vv) {

                                  test.set <- cv.split[[vv]]
                                  train.set <- hal.pseudo.dt[, id][!(hal.pseudo.dt[, id] %in% test.set)]
                                  dt.train <- hal.pseudo.dt[id %in% train.set]

                                  if (FALSE) {
                                      print("--------------")
                                      print(paste0("v = ", vv))
                                      print(paste0("delta = 1: ", dt.train[, sum(ddd.1)]))
                                      print(paste0("delta = 2: ", dt.train[, sum(ddd.2)]))
                                      print(paste0("delta = 0: ", dt.train[, sum(ddd.0)]))
                                      print(mean(dt.train[, final.time[1], by = "id"][[2]]))
                                      print(table(dt.train[, max(Y.1), by = "id"][[2]]))
                                      print("--------------")
                                  }
                                  
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

                                  return(list(
                                      cve = sapply(lambda.cvs, function(lambda.cv) {
                                          return(loss.fun(train.fit = train.fit,
                                                          hal.pseudo.dt = hal.pseudo.dt,
                                                          test.set = test.set,
                                                          X = X,
                                                          time.var = time.var, 
                                                          lambda.cv = lambda.cv,
                                                          delta.var = delta.var,
                                                          delta.value = delta.value))
                                      }),
                                      train.fit = train.fit,
                                      test.set = test.set))
        
                              })

    test.error <- try(cve <- unlist(lapply(1:length(lambda.cvs), function(mm) {
        sum(unlist(lapply(cve.list, function(out) {
            if (is.numeric(out$cve)) {
                return(out$cve[[mm]])
            } else {
                return(0)
            }
        })))
    })))

    if (!is.numeric(test.error)) {
        message("Error in CV step")
        cve <- 0
    }
    
    lambda.cv <- min(lambda.cvs[cve==min(cve)])

    return(list(min=list(lambda.cv=lambda.cv,
                         cve=unique(cve[cve==min(cve)])),
                all=cbind(lambda=lambda.cvs, cve=cve),
                train.fit = lapply(cve.list, function(out) {
                    return(list(train.fit = out$train.fit,
                                test.set = out$test.set))
                })))
}

######################################################################
### cv.fun.R ends here
