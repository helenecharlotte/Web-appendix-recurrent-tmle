### sim.data.recurrent.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Oct  2 2024 (14:47) 
## Version: 
## Last-Updated: Apr  1 2025 (07:32) 
##           By: Helene
##     Update #: 150
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

sim.data.outer <- function(n = 200,
                           intervention.A = NULL, censoring = TRUE,
                           sim.setting = "1A",
                           cens.percentage = "low",
                           tau = 1.2,
                           rep.true = 10,
                           get.cens.fraction = FALSE,
                           verbose = FALSE,
                           seed = NULL) {

    if (length(intervention.A) == 0) {
        if (!censoring) {
            if (length(seed)>0) set.seed(seed)
            out.true <- lapply(1:rep.true, function(ii) {
                sim.data.fun(n = n, tau = tau,
                             intervention.A = intervention.A,
                             censoring = censoring,
                             sim.setting = sim.setting,
                             cens.percentage = cens.percentage,
                             verbose = verbose)
            })
            out.true.1 <- out.true[[1]]
            print(cumsum(sapply(out.true, function(x) x[1, 2]))/(1:rep.true))
            plot(cumsum(sapply(out.true, function(x) x[1, 2]))/(1:rep.true))
            out.true.1[1, 2] <- mean(sapply(out.true, function(x) x[1, 2]))
            out.true.1[2, 2] <- mean(sapply(out.true, function(x) x[2, 2]))
            return(out.true.1)
        } else if (get.cens.fraction) {  
            sim.dt <- sim.data.fun(n = n, tau = tau,
                                   intervention.A = intervention.A,
                                   randomize.A = TRUE,
                                   sim.setting = sim.setting,
                                   cens.percentage = cens.percentage,
                                   verbose = verbose,
                                   seed = seed)
            sim.dt[, end.time := max(time), by = "id"]
            return(sim.dt[time == end.time, mean(delta == 0 & time <= tau)])
        } else {
            return(sim.data.fun(n = n, tau = tau,
                                intervention.A = intervention.A,
                                randomize.A = TRUE,
                                sim.setting = sim.setting,
                                cens.percentage = cens.percentage,
                                verbose = verbose,
                                seed = seed))
        }
    } else {
        if (length(seed)>0) set.seed(seed)
        out.true <- lapply(1:rep.true, function(ii) {
            sim.data.fun(n = n, tau = tau,
                         intervention.A = intervention.A,
                         sim.setting = sim.setting,
                         cens.percentage = cens.percentage,
                         verbose = verbose)
        })
        out.true.1 <- out.true[[1]]
        print(cumsum(sapply(out.true, function(x) x[1, 2]))/(1:rep.true))
        plot(cumsum(sapply(out.true, function(x) x[1, 2]))/(1:rep.true))
        out.true.1[1, 2] <- mean(sapply(out.true, function(x) x[1, 2]))
        out.true.1[2, 2] <- mean(sapply(out.true, function(x) x[2, 2]))
        return(out.true.1)
    }

}


sim.data.fun <- function(n = 200,
                         randomize.A = TRUE, 
                         sim.setting = "1A",
                         cens.percentage = "low",
                         intervention.A = NULL,
                         censoring = TRUE,
                         loop.max = 100,
                         tau = 1.2,
                         endoffollowup = 3,
                         rep.true = 10,
                         get.cens.fraction = FALSE,
                         seed = sample(1e9, 1),
                         verbose = FALSE) {

    if (length(intervention.A)>0) {
        censoring <- FALSE
    }

    #-- weibull distribution parameters:
    eta <- 0.7
    nu <- 1.7

    #-- effects on recurrent event process:
    betaT.A <- -0.8
    betaT.L1 <- 1.2
    betaT.k <- 0

    betaT.k3 <- 0
    t0 <- 0
    betaT.time <- 0
    betaT.time.t0 <- 0
    betaT.A.t0 <- 0

    #-- effects on terminal process:
    betaT2.A <- -0.4
    betaT2.L1 <- 0.7
    betaT2.k <- 0

    #-- effects on censoring process: 
    betaC.L3 <- 0
    betaC.L1 <- 1.4
    betaC.A <- 0
    betaC.k <- 0

    betaC.k3 <- 0

    #-- interaction effects: 
    betaC.AxL1 <- 0
    betaC.L1xL3 <- 0
    betaT.AxL1 <- 0
    betaT.L1xL3 <- 0

    if (sim.setting == "1A") {

        if (cens.percentage == "low") {
            alpha.C <- -1.8
        } else {
            alpha.C <- -0.8
        }
        
        betaT.k <- 2.1
        betaT2.k <- 1.4
        betaC.k <- 1.8
        alpha.T <- 0.8
        alpha.T2 <- 0

    } else if (sim.setting == "1B") {

        if (cens.percentage == "low") {
            alpha.C <- -1.7
        } else {
            alpha.C <- -0.7
        }

        betaT.k <- 2.1
        betaT2.k <- 1.4
        betaC.k <- 1.8
        alpha.T <- 0.8
        alpha.T2 <- 0
        
        betaT.L1 <- betaT2.L1 <- 0
        betaC.L3 <- betaC.L1 <- betaC.A <- 0

    } else if (sim.setting == "1Ax") {

        if (cens.percentage == "low") {
            alpha.C <- -1.7
        } else {
            alpha.C <- -0.7
        }
        
        betaT.k <- 3.1
        betaT2.k <- 1.4
        betaC.k <- 1.8
        alpha.T <- 0.4
        alpha.T2 <- 0

    } else if (sim.setting == "1Bx") {

        if (cens.percentage == "low") {
            alpha.C <- -1.7
        } else {
            alpha.C <- -0.75
        }

        betaT.k <- 3.1
        betaT2.k <- 1.4
        betaC.k <- 1.8
        alpha.T <- 0.9
        alpha.T2 <- 0
        
        betaT.L1 <- betaT2.L1 <- 0
        betaC.L3 <- betaC.L1 <- betaC.A <- 0

    } else if (sim.setting == "2A") {

        if (cens.percentage == "low") {
            alpha.C <- -0.9
        } else {
            alpha.C <- 0
        }
        
        betaT.k <- 3.1
        betaT2.k <- 1.4
        alpha.T <- 0.4
        alpha.T2 <- 0

        betaC.L3 <- betaC.L1 <- betaC.A <- 0

    } else if (sim.setting == "2B") {

        if (cens.percentage == "low") {
            alpha.C <- -0.9
        } else if (cens.percentage == "extreme") {
            alpha.C <- 1
        } else {
            alpha.C <- 0
        }

        betaT.k <- 3.1
        betaT2.k <- 1.4
        alpha.T <- 0.9
        alpha.T2 <- 0
        
        betaT.L1 <- betaT2.L1 <- 0
        betaC.L3 <- betaC.L1 <- betaC.A <- 0

    } else if (sim.setting == "3A") {

        if (cens.percentage == "low") {
            alpha.C <- -1
        } else {
            alpha.C <- -0.2
        }

        betaT.k <- 0
        betaT2.k <- 0
        
        alpha.T <- 1.1
        alpha.T2 <- 0.4

        betaT.A <- -0.8
        betaT2.A <- -0.4
        
        betaC.k <- 0

    } else if (sim.setting == "3Ax") {

        if (cens.percentage == "low") {
            alpha.C <- -0.8
        } else {
            alpha.C <- 0.2
        }

        betaT.k <- 0
        betaT2.k <- 0
        
        alpha.T <- 1.8
        alpha.T2 <- 0.8

        betaT.A <- -0.8
        betaT2.A <- -0.4
        
        betaC.k <- 0

    } else if (sim.setting == "3A1") {

        if (cens.percentage == "low") {
            alpha.C <- -0.95
        } else {
            alpha.C <- -0.05
        }

        betaT.k <- 0
        betaT2.k <- 0
        
        alpha.T <- 1.1
        alpha.T2 <- 0.4

        betaT.A <- -0.8
        betaT2.A <- -0.4
        
        betaC.k <- 0
        betaC.L3 <- betaC.L1 <- betaC.A <- 0

    } else if (sim.setting == "3A2") {

        if (cens.percentage == "low") {
            alpha.C <- -2.3
        } else {
            alpha.C <- -1.2
        }

        betaT.k <- 0
        betaT2.k <- 0
        
        alpha.T <- 1.1
        alpha.T2 <- 0.4

        betaT.A <- -0.8
        betaT2.A <- -0.4
        
        betaC.k <- 1.8

    } else if (sim.setting == "4A") {

        if (cens.percentage == "low") {
            alpha.C <- -0.5
        } else {
            alpha.C <- 0.5
        }

        betaT.k <- -0.8
        betaT.k3 <- 2.1
        betaT2.k <- 1.4
        alpha.T <- 2.1
        alpha.T2 <- 0

        betaC.k3 <- 1.8
        betaC.k <- -0.8

    }

    if (verbose) {

        print(paste0("alpha.T = ", alpha.T))
        print(paste0("betaT.A = ", betaT.A))
        print(paste0("betaT.L1 = ", betaT.L1))
        print(paste0("betaT.k = ", betaT.k))

        #-- effects on terminal process:
        print(paste0("alpha.T2 = ", alpha.T2))
        print(paste0("betaT2.A = ", betaT2.A))
        print(paste0("betaT2.L1 = ", betaT2.L1))
        print(paste0("betaT2.k = ", betaT2.k))

        #-- effects on censoring process:
        print(paste0("alpha.C = ", alpha.C))
        print(paste0("betaC.L1 = ", betaC.L1))
        print(paste0("betaC.L3 = ", betaC.L3))
        print(paste0("betaC.A = ", betaC.A))
        print(paste0("betaC.k = ", betaC.k))

    }

    ##U <- runif(n, -1, 1)
    if (length(seed)>0) set.seed(seed)
    L1 <- runif(n, -1, 1)
    U <- runif(n, -1, 1)
    L2 <- runif(n, 0, 1)
    L3 <- runif(n, 0, 1)

    if (randomize.A) {
        A <- rbinom(n, 1, plogis(qlogis(0.5)))
    } else {
        A <- rbinom(n, 1, plogis(0.2+1.3*L1))
    }

    if (length(intervention.A)>0) A <- intervention.A

    #-- recurrent event intensity:
    
    phiT <- function(t, A, L1, L2, L3, betaT.A, betaT.time, k) {
        return(exp(alpha.T+A*betaT.A+L1^2*betaT.L1+betaT.k*(k>1)+
                   betaT.time*(k>1)+betaT.k3*(k>3)+
                   betaT.L1xL3*L1*L3+betaT.AxL1*A*(L1 >= 0.5)))
    }

    lambdaT <- function(t, A, L1, L2, L3, betaT.A, betaT.time, eta, nu, k) {
        return(phiT(t, A, L1, L2, L3, betaT.A, betaT.time, k)*eta*nu*t^{nu-1})
    }

    #-- censoring intensity:

    phiC <- function(t, A, L1, L2, L3, k) {
        return(exp(betaC.L3*L3 + betaC.L1*L1 + betaC.A*A + 0.1 + alpha.C +
                   betaC.k*(k>1) + betaC.k3*(k>3) +
                   betaC.AxL1*A*L1^2 + betaC.L1xL3*L1*L3))
    }

    lambdaC <- function(t, A, L1, L2, L3, eta, nu, k) {
        return(phiC(t, A, L1, L2, L3, k)*eta*nu*t^{nu-1})
    }

    #-- terminal process intensity:
    
    phiT2 <- function(t, A, L1, L2, L3, k) {
        return(exp(alpha.T2+betaT2.L1*L1+betaT2.A*A+betaT2.k*(k>1)))
    }

    lambdaT2 <- function(t, A, L1, L2, L3, eta, nu, k) {
        return(phiT2(t, A, L1, L2, L3, k)*eta*nu*t^{nu-1})
    }

    #-- add up contributions:
    
    if (censoring) {
        phi <- function(t, A, L1, L2, L3, betaT.A, betaT.time, k) phiT(t, A, L1, L2, L3, betaT.A, betaT.time, k) +
                                                                         phiC(t, A, L1, L2, L3, k) +
                                                                         phiT2(t, A, L1, L2, L3, k)
    } else {
        phi <- function(t, A, L1, L2, L3, betaT.A, betaT.time, k) phiT(t, A, L1, L2, L3, betaT.A, betaT.time, k) +
                                                                         phiT2(t, A, L1, L2, L3, k)
    }

    if (t0 > 0) { #-- if there is an effect of time since last event: 

        Lambda.inv <- function(u, t, A, L1, L2, L3, nu, eta, k, phi) {
            return( rowSums(cbind((k == 1 | u <= (eta*phi(t, A, L1, L2, L3, betaT.A, betaT.time, k))*( (t+t0)^{nu} - t^{nu} )) *
                                  (( (u + eta*phi(t, A, L1, L2, L3, betaT.A, betaT.time, k)*t^{nu}) /
                                     (eta*phi(t, A, L1, L2, L3, betaT.A, betaT.time, k)) )^{1/nu} - t),
            (k>1 & u > (eta*phi(t, A, L1, L2, L3, betaT.A, betaT.time, k))*( (t+t0)^{nu} - t^{nu} )) *
            (( (u - (eta*phi(t, A, L1, L2, L3, betaT.A, betaT.time, k))*(t+t0)^{nu} +
                eta*phi(t, A, L1, L2, L3, betaT.A.t0*betaT.A, betaT.time.t0, k)*(t+t0)^{nu} +
                eta*phi(t, A, L1, L2, L3, betaT.A, betaT.time, k)*t^{nu}) /
               (eta*phi(t, A, L1, L2, L3, betaT.A.t0*betaT.A, betaT.time.t0, k)) )^{1/nu} - t)), na.rm=TRUE) )
        }
        
    } else {
        
        Lambda.inv <- function(u, t, A, L1, L2, L3, nu, eta, k, phi) {
            return( ( (u + eta*phi(t, A, L1, L2, L3, betaT.A, betaT.time, k)*t^{nu}) /
                      (eta*phi(t, A, L1, L2, L3, betaT.A, betaT.time, k)) )^{1/nu} - t)
        }
        
    }

    
    #-- intialize event times:
    Tlist <- list(cbind(time=rep(0, n), delta=rep(0, n), id=1:n))

    #-- function to loop over for time-points: 
    if (verbose) set.seed(seed+120202)
    loop.fun <- function(k, Tprev) {

        #-- simulate event time: 
        Unif <- -log(runif(n))

        Tout <- Lambda.inv(Unif, Tprev, A, L1, L2, L3, nu, eta, k, phi) + Tprev

        if (censoring) {
            denom <- (lambdaT(Tout, A, L1, L2, L3, betaT.A, betaT.time.t0*(Tout-Tprev>t0), eta, nu, k) +
                      lambdaC(Tout, A, L1, L2, L3, eta, nu, k) +
                      lambdaT2(Tout, A, L1, L2, L3, eta, nu, k))
        } else {
            denom <- (lambdaT(Tout, A, L1, L2, L3, betaT.A, betaT.time.t0*(Tout-Tprev>t0), eta, nu, k) +
                      lambdaT2(Tout, A, L1, L2, L3, eta, nu, k))
        }
        
        probT <- lambdaT(Tout, A, L1, L2, L3, betaT.A, betaT.time.t0*(Tout-Tprev>t0), eta, nu, k) / (denom)
        probT2 <- lambdaT2(Tout, A, L1, L2, L3, eta, nu, k) / (denom)

        if (censoring) {
            probC <- lambdaC(Tout, A, L1, L2, L3, eta, nu, k) / (denom)
            which <- apply(cbind(probC, probT, probT2), 1, function(p) sample(0:2, size=1, prob=p))
        } else {
            which <- apply(cbind(probT, probT2), 1, function(p) sample(1:2, size=1, prob=p))
        }
        
        #-- return event time
        return(cbind(time=Tout, delta=which, id=1:n))
    }

    for (k in 1:loop.max) {
        Tlist[[k+1]] <- loop.fun(k, Tlist[[k]][,1])
    }

    dt <- data.table(do.call("rbind", Tlist))

    #-- order & throw away obs after end of followup: 
    setorder(dt, id, time, delta)
    dt <- dt[time<=endoffollowup]

    dt[time>0, Y := cumsum(delta == 2 | delta == 0), by = "id"]
    dt[time>0, idN:=1:.N, by=c("id", "Y")]
    dt <- dt[(idN==1 & Y == 1) | (Y == 0)][, -c("idN", "Y"), with=FALSE][time>0]

    baseline <- data.table(A=A, L1=L1, L2=L2, L3=L3, id=as.numeric(1:n))
    setkey(baseline, id); setkey(dt, id)
    dt <- merge(dt, baseline, by="id")

    if (length(intervention.A)>0 | !censoring) {
        return(do.call("rbind", lapply(tau, function(tau.kk) {
            out <- do.call("rbind", lapply(1:2, function(each) {
                c(tau=tau.kk,
                  psi0=mean(dt[, sum(time<=tau.kk & delta==each), by = "id"][[2]]))
            }))
            rownames(out) <- c("E1", "F2")
            return(out)
        })))
    } else {
        return(dt)
    }    
}


######################################################################
### sim.data.recurrent.R ends here
