### lebesgue.loss.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Oct 15 2024 (10:13) 
## Version: 
## Last-Updated: Oct 15 2024 (10:13) 
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

lebesgue.loss.fun <- function(train.fit,
                              hal.pseudo.dt,
                              test.set,
                              X,
                              lambda.cv,
                              delta.var = "delta",
                              delta.value = 1) {

    tmp <- copy(hal.pseudo.dt)

    tmp[id %in% test.set, fit.lambda:=exp(predict(train.fit, X[tmp$id %in% test.set,],
                                                  newoffset=0, s=lambda.cv))]
    tmp[id %in% test.set, fit.dLambda:=fit.lambda*risk.time]
    tmp[id %in% test.set, fit.Lambda:=cumsum(fit.dLambda), by="id"]

    tmp[id %in% test.set, dN:=1*ddd]

    return(tmp[id %in% test.set, -sum(log(fit.lambda)*dN - fit.Lambda)])
}

######################################################################
### lebesgue.loss.fun.R ends here
