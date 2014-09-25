## These are some of the core functions used in the analyses

## Some initial setup
library(RColorBrewer)
palette(brewer.pal(8,"Dark2"))
library("rootSolve")
library("deSolve")
library(xtable)
library(fields)

############################
## The model            ###
###########################

##' Single age class model for adult TB
##' This model has an explicit Tx compartment and a presymptomatic compartment
##' @param t
##' @param y
##' @param parms
##' @return
##' @author Andrew Azman
dxdt.TBHIV3 <- function(t,y,parms){
    with(as.list(c(parms,y)),{

        ac <- 1
        hivc <- 4
        tbc <- 9

        inds <- seq(1,tbc*hivc*ac+1,by=hivc*ac) #indices for state arrays below
        S <- array(y[1:(inds[2]-1)],dim=c(ac,4))
        Lf <- array(y[inds[2]:(inds[3]-1)],dim=c(ac,4))
        Ls <- array(y[inds[3]:(inds[4]-1)],dim=c(ac,4))
        Ps <- array(y[inds[4]:(inds[5]-1)],dim=c(ac,4))
        Asp <- array(y[inds[5]:(inds[6]-1)],dim=c(ac,4))
        Asn <- array(y[inds[6]:(inds[7]-1)],dim=c(ac,4))
        Aep <- array(y[inds[7]:(inds[8]-1)],dim=c(ac,4))
        Tx <-  array(y[inds[8]:(inds[9]-1)],dim=c(ac,4))
        Rtx <- array(y[inds[9]:(inds[10]-1)],dim=c(ac,4))

        N <- sum(S + Lf + Ls + Ps + Asp + Asn + Aep + Tx + Rtx)

        ## may want to add a real hiv force of infection here later
        foi <- as.numeric(Asp %*% c(beta.sp/N*rep(1,4)) +
                          Asn %*% c((beta.sp/N)*phi.sn) +
                          Ps %*% c((beta.sp/N)*phi.ps))

        theta.sp.c <- theta.sp + theta.spI
        theta.sn.c <- theta.sn + theta.snI
        theta.ep.c <- theta.ep + theta.epI

        dS <- dLf <- dLs <- dPs <- dAsp <- dAsn <- dAep <- dTx <- dRtx <- array(0,dim=c(ac,hivc))

        ##hiv uninfected susceptibles
        dS  <- S*(nu - foi - foi.hiv  - mu.hiv - chi.elg - chi.tx - delta) +
        c(0,(S*foi.hiv)[-hivc]) +
        c(0,(S*chi.elg)[-hivc]) +
        c(0,(S*chi.tx)[-hivc])

        ## keeping population size constant
        dS[1,1] <- dS[1,1] +
        Asp %*% mu.sp + Asn %*% mu.sn + Aep %*% mu.ep +    ## TB Deaths
        (S + Lf + Ls + Ps + Asp + Asn + Aep + Tx + Rtx) %*% delta +  ## Old Age
        (S + Lf + Ls + Ps + Asp + Asn + Aep + Tx + Rtx) %*% mu.hiv   ## HIV Deaths

        ## Latent fast
        dLf <-Lf*(nu - gamma.lf.ls - rho.lf -
                   foi.hiv - mu.hiv - chi.elg - chi.tx - delta) +
                   foi*(Ls*phi.l + Rtx*phi.l + S) +
                   c(0,(Lf*foi.hiv)[-hivc]) +
                   c(0,(Lf*chi.elg)[-hivc]) +
                   c(0,(Lf*chi.tx)[-hivc])

        ## Latent slow (remote infection)
        dLs <- Ls*(nu - foi*phi.l - rho.ls -
                   foi.hiv - mu.hiv - chi.elg - chi.tx - delta) +
                   Lf * gamma.lf.ls +
                   Rtx*gamma.rtx.ls +
                   c(0,(Ls*foi.hiv)[-hivc]) +
                   c(0,(Ls*chi.elg)[-hivc]) +
                   c(0,(Ls*chi.tx)[-hivc])

        ## Pre-symptomatic period
        dPs <- Ps*(nu - rho.ps - zeta.sn -
                   foi.hiv - mu.hiv - chi.elg - chi.tx - delta) +
                   Lf*rho.lf + Ls*rho.ls +
                   c(0,(Ps*foi.hiv)[-hivc]) +
                   c(0,(Ps*chi.elg)[-hivc]) +
                   c(0,(Ps*chi.tx)[-hivc])
        ## Smear Positive
        dAsp <- Asp*(nu - mu.sp - theta.sp.c - zeta.sp -
                     foi.hiv - mu.hiv - chi.elg - chi.tx - delta) +
                     (Ps*rho.ps + Rtx*rho.rel)*pi.sp*(1-pi.ep) +
                     c(0,(Asp*foi.hiv)[-hivc]) +
                     c(0,(Asp*chi.elg)[-hivc]) +
                     c(0,(Asp*chi.tx)[-hivc])

        dAsn <- Asn*(nu - mu.sn - theta.sn.c - zeta.sn -
                     foi.hiv - mu.hiv - chi.elg - chi.tx - delta) +
                     (Ps*rho.ps + Rtx*rho.rel)*(1-pi.sp)*(1-pi.ep) +
                     c(0,(Asn*foi.hiv)[-hivc]) +
                     c(0,(Asn*chi.elg)[-hivc]) +
                     c(0,(Asn*chi.tx)[-hivc])

        dAep <- Aep*(nu - mu.ep - theta.ep.c - zeta.ep -
                     foi.hiv - mu.hiv - chi.elg - chi.tx - delta) +
                     (Ps*rho.ps + Rtx*rho.rel)*pi.ep+
                     c(0,(Aep*foi.hiv)[-hivc]) +
                     c(0,(Aep*chi.elg)[-hivc]) +
                     c(0,(Aep*chi.tx)[-hivc])

        dTx <- Tx*(nu - gamma.tx.rtx -
                   foi.hiv - mu.hiv - chi.elg - chi.tx - delta) +
                   Asp*theta.sp.c +
                   Asn*theta.sn.c +
                   Aep*theta.ep.c +
                   c(0,(Tx*foi.hiv)[-hivc]) +
                   c(0,(Tx*chi.elg)[-hivc]) +
                   c(0,(Tx*chi.tx)[-hivc])

        dRtx <- Rtx*(nu - gamma.rtx.ls - rho.rel - foi*phi.l -
                     foi.hiv - mu.hiv - chi.elg - chi.tx - delta) +
                     Asp*zeta.sp + (Asn + Ps)*zeta.sn + Aep*zeta.ep +
                     Tx*(gamma.tx.rtx) +
                     c(0,(Rtx*foi.hiv)[-hivc]) +
                     c(0,(Rtx*chi.elg)[-hivc]) +
                     c(0,(Rtx*chi.tx)[-hivc])

        list(c(dS,dLf,dLs,dPs,dAsp,dAsn,dAep,dTx,dRtx))

    })}


##' take dxdt.TBHIV3 odes and appends some summary statistics to each time step
##' @param t
##' @param state
##' @param params
##' @return vector of state changes
##' @author Andrew Azman
dxdt.TBHIV.CI <- function(t,state,params){

    ## a little pre-processing
    ac <- 1 ## number of age classes
    hivc <- 4
    tbc <- 9

    inds <- seq(1,tbc*hivc*ac+1,by=hivc*ac) #indices for state arrays

    S <- array(state[1:(inds[2]-1)],dim=c(ac,4))
    Lf <- array(state[inds[2]:(inds[3]-1)],dim=c(ac,4))
    Ls <- array(state[inds[3]:(inds[4]-1)],dim=c(ac,4))
    Ps <- array(state[inds[4]:(inds[5]-1)],dim=c(ac,4))
    Asp <- array(state[inds[5]:(inds[6]-1)],dim=c(ac,4))
    Asn <- array(state[inds[6]:(inds[7]-1)],dim=c(ac,4))
    Aep <- array(state[inds[7]:(inds[8]-1)],dim=c(ac,4))
    Tx <-  array(state[inds[8]:(inds[9]-1)],dim=c(ac,4))
    Rtx <- array(state[inds[9]:(inds[10]-1)],dim=c(ac,4))


    with(as.list(c(state,params)),{
        ## rho.lf <- array(rho.lf,dim=c(ac,length(rho.lf)/2))
        ## rho.ls <- array(rho.ls,dim=c(ac,length(rho.ls)/2))
        ## rho.rel <- array(rho.rel,dim=c(ac,length(rho.rel)/2))

        dCI <- c((Lf * rho.lf) + (Ls * rho.ls) + (Rtx * rho.rel)) #1x4 number of new cases of each type
        #dCI <- c((Ps * rho.ps) + (Rtx * rho.rel)) #1x4 number of new cases of each type

        dCIall <- sum(dCI) #1x1 sum of all incident tb types
        ## tb deaths in each age class and hiv status 1x8
        dMtb <- c((Asp * mu.sp) + (Asn * mu.sn) + (Aep * mu.sn)) # number of new TB deaths

        ## cases detected of each type (this is what we will fit to)
        dN.Asp <- c(Asp * (theta.sp + theta.spI)) #1x4
        dN.Asn <- c(Asn * (theta.sn + theta.snI)) #1x4
        dN.Aep <- c(Aep * (theta.ep + theta.epI)) #1x4
        dReTx <-  c(Rtx * rho.rel)  #1x4

        c(dCI,dCIall,dMtb,dN.Asp,dN.Asn,dN.Aep,dReTx)
    }) -> dInds

    ##run TB model
    TBout <- dxdt.TBHIV3(t,state,params)

    ##Return the results
    rc <- list(c(TBout[[1]],dInds))
    return(rc)
}

##' take dxdt.TBHIV3 odes and appends some summary statistics to each time step
##' and allows beta to vary by a fixed amount per year
##' @param t
##' @param state
##' @param params
##' @return vector of state changes
##' @author Andrew Azman
dxdt.TBHIV.CI.var.beta <- function(t,state,params){

    ## a little pre-processing
    ac <- 1 ## number of age classes
    hivc <- 4
    tbc <- 9

    inds <- seq(1,tbc*hivc*ac+1,by=hivc*ac) #indices for state arrays

    S <- array(state[1:(inds[2]-1)],dim=c(ac,4))
    Lf <- array(state[inds[2]:(inds[3]-1)],dim=c(ac,4))
    Ls <- array(state[inds[3]:(inds[4]-1)],dim=c(ac,4))
    Ps <- array(state[inds[4]:(inds[5]-1)],dim=c(ac,4))
    Asp <- array(state[inds[5]:(inds[6]-1)],dim=c(ac,4))
    Asn <- array(state[inds[6]:(inds[7]-1)],dim=c(ac,4))
    Aep <- array(state[inds[7]:(inds[8]-1)],dim=c(ac,4))
    Tx <-  array(state[inds[8]:(inds[9]-1)],dim=c(ac,4))
    Rtx <- array(state[inds[9]:(inds[10]-1)],dim=c(ac,4))

    params$beta.sp <- params$beta.sp*exp(params$beta.delta*t)
    #cat(sprintf("beta.sp = %.2f, and beta.delta = %.3f \n",params$beta.sp[1],params$beta.delta[1]))


    with(as.list(c(state,params)),{
        ## rho.lf <- array(rho.lf,dim=c(ac,length(rho.lf)/2))
        ## rho.ls <- array(rho.ls,dim=c(ac,length(rho.ls)/2))
        ## rho.rel <- array(rho.rel,dim=c(ac,length(rho.rel)/2))

        dCI <- c((Lf * rho.lf) + (Ls * rho.ls) + (Rtx * rho.rel)) #1x4
        #dCI <- c((Ps * rho.ps) + (Rtx * rho.rel)) #1x4

        dCIall <- sum(dCI) #1x1
        ## tb deaths in each age class and hiv status 1x8
        dMtb <- c((Asp * mu.sp) + (Asn * mu.sn) + (Aep * mu.sn))

        ## cases detected of each type in formal sector (this is what we will fit to)
        dN.Asp <- c(Asp * (theta.sp + theta.spI)) #1x4
        dN.Asn <- c(Asn * (theta.sn + theta.snI)) #1x4
        dN.Aep <- c(Aep * (theta.ep + theta.epI)) #1x4
        dReTx <-  c(Rtx * rho.rel)  #1x4

        c(dCI,dCIall,dMtb,dN.Asp,dN.Asn,dN.Aep,dReTx)
    }) -> dInds


    ##run TB model
    TBout <- dxdt.TBHIV3(t,state,params)

    ##Return the results
    rc <- list(c(TBout[[1]],dInds))
    return(rc)
}

######################
## Helper functions ##
######################

##' Adds column names to output from ode
##' @param mod
##' @param time
##' @param ext
##' @param ac
##' @return
##' @author Andrew Azman
addColNames <- function(mod,time=T,ext=F,ac=1){
    ts <- c()
    if (time) ts <- "time"


        tmp  <-  c(ts,paste0("S",1:ac),
                   paste0("hS",1:ac),
                   paste0("aS",1:ac),
                   paste0("nS",1:ac),
                   paste0("Lf",1:ac),
                   paste0("hLf",1:ac),
                   paste0("aLf",1:ac),
                   paste0("nLf",1:ac),
                   paste0("Ls",1:ac),
                   paste0("hLs",1:ac),
                   paste0("aLs",1:ac),
                   paste0("nLs",1:ac),
                   paste0("Ps",1:ac),
                   paste0("hPs",1:ac),
                   paste0("aPs",1:ac),
                   paste0("nPs",1:ac),
                   paste0("Asp",1:ac),
                   paste0("hAsp",1:ac),
                   paste0("aAsp",1:ac),
                   paste0("nAsp",1:ac),
                   paste0("Asn",1:ac),
                   paste0("hAsn",1:ac),
                   paste0("aAsn",1:ac),
                   paste0("nAsn",1:ac),
                   paste0("Aep",1:ac),
                   paste0("hAep",1:ac),
                   paste0("aAep",1:ac),
                   paste0("nAep",1:ac),
                   paste0("Tx",1:ac),
                   paste0("hTx",1:ac),
                   paste0("aTx",1:ac),
                   paste0("nTx",1:ac),
                   paste0("Rtx",1:ac),
                   paste0("hRtx",1:ac),
                   paste0("aRtx",1:ac),
                   paste0("nRtx",1:ac))

    if (ext) {

        tmp <- c(tmp,paste0("CI",1:ac),paste0("hCI",1:ac),paste0("aCI",1:ac),
                 paste0("nCI",1:ac),"CIall",paste0("Mtb",1:ac),
                 paste0("hMtb",1:ac),paste0("aMtb",1:ac),paste0("nMtb",1:ac),
                 paste0("N.Asp",1:ac),paste0("hN.Asp",1:ac),paste0("aN.Asp",1:ac),
                 paste0("nN.Asp",1:ac),paste0("N.Asn",1:ac),paste0("hN.Asn",1:ac),
                 paste0("aN.Asn",1:ac),paste0("nN.Asn",1:ac),
                 paste0("N.Aep",1:ac),paste0("hN.Aep",1:ac),paste0("aN.Aep",1:ac),paste0("nN.Aep",1:ac),
                 paste0("ReTx",1:ac),paste0("hReTx",1:ac),paste0("aReTx",1:ac),paste0("nReTx",1:ac))
    }

    if (!is.null(nrow(mod))){
        colnames(mod) <- tmp
    } else {
        names(mod) <-  tmp
    }

    return(mod)
}

##' takes parameters from csv file
##' @param country name of country whose parameters we want (assumes in common form)
##' @param cols column numbers for data
##' @return list with each entry being a vector for that parameter
##' @author Andrew Azman
make.params <- function(country,cols=2:5){
    filename <- sprintf("Data/%s_params.csv",country)
    tmp <- read.csv(filename)
    params.block <- tmp[cols]
    rownames(params.block) <- tmp[,1]
    params.list <- do.call("list",as.data.frame(t(params.block)))
    return(params.list)
}

## runs TBHIV.CI model
##' @param params
##' @param initial.state
##' @param max.time
##' @param var.beta
##' @return output of lsoda or other ode solver
runTBHIVMod <- function(params,
                        initial.state,
                        max.time=1,
                        var.beta = FALSE
                        ){
    library(deSolve)

    times <- seq(0,max.time,by=0.1)
    ##print(params)
    if (var.beta){
        mod.out <- ode(initial.state,times,dxdt.TBHIV.CI.var.beta,params)
    } else {
        mod.out <- ode(initial.state,times,dxdt.TBHIV.CI,params)
    }
    return(mod.out)
}


##' takes a matrix with column names of model
##' and outputs just the columns needed for prevalence
##' @param run.mat
##' @param hiv.only
##' @return matrix of only columns of prev cases
##' @author Andrew Azman
getPrevCols <- function(run.mat,hiv.only=F){

    ## in case it is a vector
    if (is.null(nrow(run.mat))) run.mat <- t(as.matrix(run.mat))

    if (hiv.only){
        run.mat[,grep("(n|a|h)(A(sp|sn|ep)|Tx|Ps)1$",colnames(run.mat))]
    } else {
        run.mat[,grep("(n|a|h|^)(A(sp|sn|ep)|Tx|Ps)1$",colnames(run.mat))]
    }}

##' Objective function for fitting Incidence and CDR
##' @param params.fit
##' @param params
##' @param state
##' @param target.ci
##' @param target.cdr
##' @param target.prev.tb
##' @param plot.it
##' @param beta.or.theta - if we want to only fit one param ("beta" if we only want to fit beta, "theta" if we want to fit theta only)
##' @param weight.ci
##' @param weight.other
##' @return
##' @author Andrew Azman
incObFunc <- function(params.fit,
                      params,
                      state,
                      target.ci,
                      target.cdr,
                      target.prev.tb,
                      plot.it=FALSE,
                      beta.or.theta="",
                      weight.ci = 3,
                      weight.other=1
                      ){

    if (length(params.fit) == 1 && !missing(beta.or.theta)){
        if (beta.or.theta == "beta") params$beta.sp <- rep(params.fit,4)
        else if (beta.or.theta == "theta") params$theta.sp <- rep(params.fit,4)
        else stop("beta.or.theta is mispecficified")
    } else {
        params$beta.sp <- rep(params.fit[1],4)
        params$theta.sp <- rep(params.fit[2],4)
    }
    ##  cat(sprintf("fit.pars (post optim) = %f, %f \n",exp(params.fit)[1],exp(params.fit)[2]))

    ## assuming that the case detection rate of ep is same as sp
    ## and that sn is 0.75* sp
    ep.sn.mult <- 1
    params$theta.ep <- params$theta.sp*ep.sn.mult
    params$theta.sn <- params$theta.sp*ep.sn.mult

    tryCatch(
        RS <- runsteady(y=state[1:36],
                        fun=dxdt.TBHIV3,
                        parms=params,
                        verbose=F)
        ,
        error = function(e){
            ss.vals <- state
            cat(sprintf(e$message))
        }
        )

    if (attr(RS,"steady")){
        ss.vals <- c(RS$y,state[37:length(state)])
    } else {
        print("Couldn't reach steady state but proceeding to next set of paramters in optimization")
        ss.vals <- state
    }

    run <- runTBHIVMod(params,initial.state=ss.vals,max.time=1,var.beta=FALSE)
    run <- addColNames(run,ext=T,time=T)

    ci <- run[11,"CIall"] - run[1,"CIall"]

    if (!missing(target.prev.tb)){
        ## calc prevalance stats
        prev <- sum(getPrevCols(run)[11,])

        if (!missing(beta.or.theta) && beta.or.theta == "theta"){
            obj <-  ((prev/target.prev.tb) - 1)^2
            obj.no.trans <- 1 # value if there is no tranmission
        } else if (!missing(beta.or.theta) && beta.or.theta == "beta"){
            obj <- ((ci/target.ci) - 1)^2
            obj.no.trans <- 1
        } else {
            obj <- weight.ci*((ci/target.ci) - 1)^2 + weight.other*((prev/target.prev.tb) - 1)^2
            obj.no.trans <- 2
        }
        print(c(ci,target.ci=target.ci,prev=prev,target.prev=target.prev.tb))
    } else {
        cd <- (run[11,grep("N.Asp",colnames(run))] +
               run[11,grep("N.Asn",colnames(run))] +
               run[11,grep("N.Aep",colnames(run))]) -
                   (run[1,grep("N.Asp",colnames(run))] +
                    run[1,grep("N.Asn",colnames(run))] +
                    run[1,grep("N.Aep",colnames(run))])
        ## but we really want to fit to cases detected which is not implicitly a function of ci
        cd.num <- sum(cd)
        cdr <- (sum(cd)/ci)*100
        cd.num.target <- target.cdr*target.ci

        print(c(ci,target.ci=target.ci,cdr=cdr,target.cdr=100*target.cdr))

        if (!missing(beta.or.theta) && beta.or.theta == "theta"){
            obj <-  (cdr - target.cdr*100)^2
            obj.no.trans <- 1000000 # value if there is no tranmission
        } else if (!missing(beta.or.theta) && beta.or.theta == "beta"){
            print("beta")
            obj <- (ci - target.ci)^2
            obj.no.trans <- 1000000
        } else {
            obj <- weight.ci*((ci/target.ci) - 1)^2 + weight.other*((cd.num/cd.num.target) - 1)^2
            obj.no.trans <- 2
        }
    }
    print(c(params$beta.sp[1],params$theta.sp[1]))

    if (is.nan(obj) || obj == obj.no.trans) obj <- Inf #when we get no tranmission the ob func = 2
    cat(sprintf("objective func = %f \n",obj))

    if (plot.it){
        points(params$theta.sp[1],obj,col=2)
    }
    return(obj) # may think about scaling the objective function
}


##' For fitting incidence and % cases detected to thetea and beta
##' @param initial.state
##' @param params
##' @param target.ci
##' @param target.cdr
##' @return
##' @author Andrew Azman
fitIncCDR <- function(initial.state,
                      params,
                      target.ci,
                      target.cdr,
                      epsilon.cdr.inc.target=0.1
                      ){


    require("rootSolve")
    ## set all theta's to theta sp
    fit.pars <- c(params$beta.sp[1],params$theta.sp[1])
    print(fit.pars)

    ##fit each serperatley and iterate between em.
    epsilon.cdr.inc <- Inf
    while (epsilon.cdr.inc >= epsilon.cdr.inc.target){

        cur.beta <- params$beta.sp[1]
        cur.theta <- params$theta.sp[1]
        out.beta <- optim(fit.pars[1],
                          fn=incObFunc,
                          params=params,
                          state=initial.state,
                          target.ci=target.ci,
                          target.cdr=target.cdr,
                          beta.or.theta = "beta",
                          method="Brent",
                          lower=2,upper=100, #optimization is finicky! adjust lower bound
                          control=list(trace=T,abstol=1))
                                        #update beta

        params$beta.sp <- rep(out.beta$par,4)

                                        #update initial state

        out.theta <- optim(fit.pars[2],
                           fn=incObFunc,
                           params=params,
                           state=initial.state,
                           target.ci=target.ci,
                           target.cdr=target.cdr,
                           beta.or.theta = "theta",
                           method="Brent",
                           lower=0.1,
                           upper=2.5, #optimization is finicky! adjust lower bound
                           control=list(trace=T,abstol=1))

        ## update thetas
        ep.sn.mult <- 1 ## Assuming equal impcat on all tb types
        params$theta.sp <- rep(out.theta$par,4)
        params$theta.sn <- ep.sn.mult*rep(out.theta$par,4)
        params$theta.ep <- ep.sn.mult*rep(out.theta$par,4)

        ## now calculate the change
        epsilon.cdr.inc <- max(c(abs(cur.theta - out.theta$par)/cur.theta,abs(cur.beta - out.beta$par)/cur.beta))
    }

    ## start.state.min <- initial.state
    tryCatch(RS <- runsteady(y=initial.state,fun=dxdt.TBHIV.CI,parms=params,times=c(0,10000),verbose=F),
             error = function(e){
                 stop("Sorry can't reach steady state from optimized params")
             })

    ss.vals <- RS$y

    return(list(final.pars=params,ss=ss.vals))
}


##' Function to fit theta.sp and beta to TB preva and incidence
##' @param initial.state
##' @param params
##' @param target.ci
##' @param target.prev.tb
##' @return
##' @author Andrew Azman
fitIncPrev <- function(initial.state,
                       params,
                       target.ci,
                       target.prev.tb,
                       lowers=c(4,.1),
                       uppers=c(20,7)
                       ){

    require("rootSolve")
    ## set all theta's to theta sp
    fit.pars <- c(params$beta.sp[1],params$theta.sp[1])
    print(fit.pars)

    out <- optim(fit.pars,
                 fn=incObFunc,
                 params=params,
                 state=initial.state,
                 target.ci=target.ci,
                 target.prev.tb=target.prev.tb,
                 method="L-BFGS-B",
                 lower=lowers,upper=uppers, #optimization is finicky! adjust lower bound
                 control=list(trace=T,parscale=c(10,1),maxit=1000))


    final.pars <- params
    final.pars$beta.sp <- rep(out$par[1],4)
    final.pars$theta.sp <- rep(out$par[2],4)
    ep.sn.mult <- 1
    final.pars$theta.ep <-  final.pars$theta.sp*ep.sn.mult
    final.pars$theta.sn <-  final.pars$theta.sp*ep.sn.mult

    tryCatch(RS <- runsteady(y=initial.state,fun=dxdt.TBHIV.CI,parms=final.pars,times=c(0,10000),verbose=F),
             error = function(e){
                 stop("Sorry can't reach steady state from optimized params")
             })

    ss.vals <- RS$y

    return(list(final.pars=final.pars,ss=ss.vals))
}

## Runs intervention and control with a specfified increase in the detection rates
##' @param ss starting state for runs, should include the main states and claculated ones
##' @param params list of parameters to use in the simulations
##' @param time how long to run the models
##' @param int.theta.sp - increased rate of detection of sp TB
##' @param int.theta.sn - increased rate of detection of sn TB
##' @param int.theta.ep - increased rate of detection of ep TB
##' @return
runIntCont <- function(ss,
                       params,
                       time,
                       int.theta.sp,
                       int.theta.sn,
                       int.theta.ep,
                       var.beta=FALSE,
                       intervention.duration=time){

    ## make sure all the stats for the ss are set to zero
                                        #ss[37:length(ss)] <- 0

    cont <- runTBHIVMod(params,initial.state=ss,max.time=time,var.beta=var.beta)
    cont <- addColNames(cont,ext=T)

    params.int <- params
    params.int$theta.snI <- rep(int.theta.sn,4)
    params.int$theta.spI <- rep(int.theta.sp,4)
    params.int$theta.epI <- rep(int.theta.ep,4)

    ## first we will run the intervention
    int <- runTBHIVMod(params.int,initial.state=ss,max.time=intervention.duration,var.beta=var.beta)

    if (intervention.duration < time){
        int.part2 <- runTBHIVMod(params,initial.state=tail(int,1)[-1],max.time=time-intervention.duration,var.beta=var.beta)
        int <- rbind(int,int.part2[-1,])
        int[,1] <- seq(0,time,by=0.1)
    }

    int <- addColNames(int,ext=T)

    return(list(int=int,cont=cont))
}

#takes a run and plots incdience and cases detected
plotOut <- function(out,pop.adj=T,overlay=FALSE,legend=TRUE){
    if (pop.adj){
        limit <- grep("CI",colnames(out)) ##which is the first col of stats
        pa <- rowSums(out[,2:(limit-1)])/100000 ## pop.size / 100k
        pa <- pa[-1]  #since we are starting after 2008.0
    } else {
        pa <- rep(1,nrow(out)-1)
    }

    cd <- grep("N.",colnames(out))
    ## get cases detected per 100k (if adjusted)
    cases.detected <- (diff(rowSums(out[,cd]))/pa)*10
    times <- out[,1]
    ## get prevalence
    prev <- rowSums(getPrevCols(out))/c(1,pa)
    ##get incidence
    inc <- (diff(out[,"CI"])/pa)*10
    if (!overlay){
        plot(times,prev,col=1,type="l",ylim=c(0,700),lty=1,xlab="",ylab="Rate per 100k per year")
        lines(times[-1],inc,col=2,type="l",lty=1)
        lines(times[-1],cases.detected,col=3,type="l",lty=1)
    } else {
        lty <- 2
        lines(times,prev,col=1,type="l",lty=lty)
        lines(times[-1],inc,col=2,type="l",lty=lty)
        lines(times[-1],cases.detected,col=3,type="l",lty=lty)
    }

    if(legend & overlay){
        legend("topright",c("Prevalence, Intervention","Incidence, Intervention","Cases Detected, Intervention","Prevalence, No Intervention","Incidence, No Intervention","Cases Detected, No Intervention"),col=c(1:3,1:3),lty=c(rep(1,3),rep(2,3)),bty="n")
    } else if (legend){
        legend("topright",c("Prev","Inc","Detected"),col=1:3,lty=1,bty="n")
    }
}

##' Calculates HIV related summary statistics given model state
##' @param mod model state
##' @param full a flag for whether or not we are giving a full model output ot the function or not (or jsut a single line)
##' @return vector, prevalance and prop.on ARTs for both age classes
##' @author Andrew Azman
hIVStats <- function(mod,full=F){

    if(!is.null(nrow(mod)) && colnames(mod)[1] == "time") mod <- mod[,-1]

    if(is.null(nrow(mod)) && names(mod)[1] == "time") mod <- mod[-1]

    if(!is.null(nrow(mod))){
        #recover()
        ## assuming that the first CI column is the first one of cumulative statistics
        first.column.of.cum.stats <- grep("CI",colnames(mod))

        if (length(first.column.of.cum.stats) > 0){
            mod <- mod[,-c(first.column.of.cum.stats[1]:ncol(mod))]
        }

        prev.1 <-
        apply(mod[,grep("^[han]",colnames(mod))],1,sum)/
        rowSums(mod[,grep(".+1$",colnames(mod))])
        ## note the the labels for n and a are actually reveresed
        prop.on.art.1 <-
        rowSums(mod[,grep("^n.+1$",colnames(mod))])/
        rowSums(mod[,grep("(^a.+1$)|(^n.+1$)",colnames(mod))])
        ## only considering those eligible

        if (full) {
            return(list(prev.1,prop.on.art.1))
        } else {
            return(c(hiv.prev.1=tail(prev.1,1),prop.art.1=tail(prop.on.art.1,1)))
        }
    } else {

        ## assuming that the first CI column is the first one of cumulative statistics
        first.column.of.cum.stats <- grep("CI",names(mod))

        if (length(first.column.of.cum.stats) > 0){
            mod <- mod[first.column.of.cum.stats[1]:ncol(mod)]
        }

        prev.1 <- sum(mod[grep("^[han]",names(mod))])/sum(mod[grep(".+1$",names(mod))])
        ## note the the labels for n and a are actually reveresed
        prop.on.art.1 <- sum(mod[grep("^n.+1$",names(mod))])/
        sum(mod[grep("(^a.+1$)|(^n.+1$)",names(mod))])
        return(c(hiv.prev.1=prev.1,prop.art.1=prop.on.art.1))
    }
}


##' Takes parameters and model starting state, runs to steady state and estimates the sum of squared errors for HIV STAT output
##' @param fit.params
##' @param full.params
##' @param state
##' @param prev.1 true HIV prevalence for 1st age clas
##' @param prop.art.1 true propirtion of hiv eligible that are on ARTs (<15)
##' @return sum of squared errors for hiv.prev and prop.on.art for each age class (equally weighted and not scaled)
##' @author Andrew Azman
hIVObjective <- function(fit.params,
                         full.params,
                         state,
                         prev.1,
                         prop.art.1){

    full.params$chi.tx[3] <- fit.params[1]
#    full.params$chi.tx[2] <- fit.params[2]
    full.params$foi.hiv[1] <- fit.params[2]
    ## full.params$foi.hiv[2] <- fit.params[4]
    #print(fit.params)

    RS <- runsteady(y=state,fun=dxdt.TBHIV3,parms=full.params,verbose=F)

    tmp <- addColNames(RS$y,time=F)

    (stats <- hIVStats(tmp))

   # print(matrix(c(stats,prev.1,prop.art.1),nrow=2,byrow=T))
#    recover()
    sum((stats/c(prev.1,prop.art.1) - 1)^2)
}



##' Fits the chi.tx (rate of flow from eligble to ART) for each age class and foi.hiv (the constant rate of new hiv infections)
##' @param start.pars
##' @param params
##' @param start.state
##' @param prev.1
##' @param prop.art.1
##' @return final parameters of optimization routine
##' @author Andrew Azman
fitHIV <- function(params,
                   start.state,
                   prev.1,
                   prop.art.1){

    start.pars <- c(params$chi.tx[3],params$foi.hiv[1])
    fit <- optim(start.pars,
                 fn=hIVObjective,
                 full.params=params,
                 state=start.state,
                 prev.1=prev.1,
                 prop.art.1=prop.art.1,
                 method="L-BFGS-B",
                 lower=c(1e-5,1e-10),
                 upper=c(365,1),
                 control=list(parscale=c(1,.1)))
    fit
}


##' Gets percentage of people of each age for a given model output
##' @title
##' @param mod.out
##' @param classes number of age classes in the model
##' @return
getAgeDistribution <- function(mod.out,classes=2){
    ages <- c()
    for (i in 1:classes){
        ages[i] <- sum(mod.out[nrow(mod.out),grep(paste0(i,"$"),colnames(mod.out))])
    }
    ages/sum(ages)
}


##' Function takes a single year of data and returns some key TB related stats
##' prevalence , incidence, mortality
##' cases detected per year
##' percent of new TB infections that are HIV positive
##' @title
##' @param mod
##' @return
##' @author Andrew Azman
getTBStats <- function(mod,add.names=T,row.final,row.init){

    if (add.names) mod <- addColNames(mod,time=T,ext=T)

    if(missing(row.final) || missing(row.init)){
        row.final <- nrow(mod)
        row.init <- row.final - 10
    }

    ## overall TB mortality
    tb.mort <- sum(mod[row.final,grep("Mtb",colnames(mod))] -
                   mod[row.init,grep("Mtb",colnames(mod))])

    tb.hiv.mort <- sum(mod[row.final,grep("(a|h|n)Mtb",colnames(mod))] -
                       mod[row.init,grep("(a|h|n)Mtb",colnames(mod))])

    tb.prev <- sum(getPrevCols(mod)[row.final,])
    tb.hiv.prev <- sum(getPrevCols(mod,hiv.only=T)[row.final,])

    tb.inc <- mod[row.final,"CIall"] - mod[row.init,"CIall"]
    tb.hiv.inc <- sum(mod[row.final,grep("(a|h|n)CI",colnames(mod))] -
                      mod[row.init,grep("(a|h|n)CI",colnames(mod))])

    return(round(c(tb.mort.nohiv=tb.mort-tb.hiv.mort,
                   tb.hiv.mort=tb.hiv.mort,
                   tb.hiv.prev=tb.hiv.prev,
                   tb.prev=tb.prev,
                   tb.inc=tb.inc,tb.hiv.inc=tb.hiv.inc),1))
}


iterativeHIVTBFit <- function(start.state,
                              params.start,
                              target.ci=993,
                              target.cdr=0.69,
                              target.prev.tb = 768,
                              target.prev.hiv = 0.178,
                              target.art = 0.55,
                              epsilon.target=1e-2,
                              uppers.tb=c(20,4),
                              lowers.tb=c(5,.1)){

    ## initialize parameters
    epsilon <- Inf
    tmp.state <- start.state
    params.tmp <- params.start
    ## params.hiv.tmp <- params.hiv.start
    ## params.tb.tmp <- params.tb.start
    ## set up proposed parameter vector
    par.cur <- c(params.tmp$chi.tx[3],
                 params.tmp$foi.hiv[1],
                 params.tmp$beta.sp[1],
                 params.tmp$theta.sp[1])

    par.new <- rep(NA,4)

    while(epsilon > epsilon.target){

        hiv.fit.sa <- fitHIV(params.tmp,
                             tmp.state[1:36],
                             prev.1=target.prev.hiv,
                             prop.art.1=target.art)

        par.new[1] <- params.tmp$chi.tx[3] <- hiv.fit.sa$par[1]
        par.new[2] <- params.tmp$foi.hiv[1] <- hiv.fit.sa$par[2]

        if(!missing(target.prev.tb)){
            tb.fit.tmp <- fitIncPrev(initial.state=tmp.state,
                                     params=params.tmp,
                                     target.ci=target.ci,
                                     target.prev.tb=target.prev.tb,
                                     uppers=uppers.tb,lowers=lowers.tb)
        } else {
            tb.fit.tmp <- fitIncCDR(initial.state=tmp.state,
                                    params=params.tmp,
                                    target.ci=target.ci,
                                    target.cdr=target.cdr                                                                )
        }

        params.tmp$beta.sp <- tb.fit.tmp$final.pars$beta.sp
        params.tmp$theta.sp <- tb.fit.tmp$final.pars$theta.sp

        par.new[3] <- tb.fit.tmp$final.pars$beta.sp[1]
        par.new[4] <- tb.fit.tmp$final.pars$theta.sp[1]

        ## change if we alter relations hsip between theta.sp and the others
        params.tmp$theta.sn <- tb.fit.tmp$final.pars$theta.sp*1
        params.tmp$theta.ep <- tb.fit.tmp$final.pars$theta.sp*1

        epsilon <- max(abs(par.new - par.cur)/par.cur)
        par.cur <- par.new
        tmp.state <- tb.fit.tmp$ss

        cat(sprintf("Pct change in params from last optim is %f \n",epsilon))
    }

    list(params=params.tmp,
         state=tmp.state,
         epsilon=epsilon)
}


##' Takes output from runIntCont
##' @param out
##' @param times
##' @param costs
##' @param params
##' @param ...
##' @return
##' @author Andrew Azman
makeHorizonICERPlot <- function(out,times,costs,params,...){
    cols <- brewer.pal(6, name="Greens")
    cols <-colorRampPalette(cols, space = "Lab")
    colors<-cols(length(times)+3)
    plot(-100,-100,xlim=range(costs),ylim=c(0,600),xlab="",ylab="")
    sapply(1:length(times),function(horiz){
        lines(costs,sapply(1:length(costs),function(cost)
                     calcStats(out,eval.times=1:((horiz*10)+1),dtx.cost=cost,params=params,...)["ICER"]),col=horiz)
              #colors[horiz+2])
              })
}

##' Makes a levelplot of ICERs by cost and analystic time horizon
##' @param out output of runIntCont
##' @param times time horzozons
##' @param costs costs we want to evaluate it at
##' @param params parameters vector
##' @param xlabs
##' @param ylabs
##' @param ...
##' @return plot
##' @author Andrew Azman
makeLevelPlotICER <- function(out,times,costs,params,xlabs,ylabs,...){
    require(fields)
    cols <- brewer.pal(9, name="Greens")
    cols <-colorRampPalette(cols[-1], space = "Lab")
    grid <- expand.grid(times,costs)
    ICERS <- mapply(getICER,horiz=grid[,1],cost=grid[,2],MoreArgs= list(params=params,out=out,...))
    mat <- matrix(ICERS,nrow=length(times),ncol=length(costs))
#    layout(matrix(c(1,2),nrow=1),widths = c(.9,.1))
#    par(mar=c(2,2,2,2))
    par(mar=c(5,4.5,4,7))
    image(mat,col=cols(15),axes=F,xlab="Time Horizon (years)",ylab="Diagnosis Cost (USD)")
    axis(1,at=seq(0,1,length=length(xlabs)),labels=xlabs)
    axis(2,at=seq(0,1,length=length(ylabs)),labels=ylabs)
    image.plot(col=cols(15),zlim=range(ICERS),legend.only=T,horizontal=F,width=5)
}

##' Helper function
##' @param horiz
##' @param cost
##' @param params
##' @param out
##' @param fixed true if we are fixing
##' @param ...
##' @return
##' @author Andrew Azman
getICER <- function(horiz,cost,params,out,fixed,...){
    if (fixed){
        calcICERFixedCosts(out,eval.times=1:((horiz*10)+1),dtx.cost=cost,params=params,...)["ICER"]
    } else {
        calcICER(out,eval.times=1:((horiz*10)+1),dtx.cost=cost,params=params,...)["ICER"]
    }
}

##' objective function for fitting annual percent change in beta to change in CI
##' @title
##' @param beta.delta
##' @param params
##' @param ss
##' @param target.ci
##' @param years
##' @return
##' @author Andrew Azman
fitAnnualBetaDeltaObjFunc <- function(beta.delta,params,ss,target.ci,years){
    params$beta.delta <- rep(beta.delta,4)
    out <- runTBHIVMod(params,ss,years,T)
    ret <- (target.ci - getTBStats(out)[5])^2
    cat(sprintf("Target = %f, Current = %f \n",target.ci,getTBStats(out)[5]))
    ret
}

##' Fits annual pct change in beta
##' @param params
##' @param ss
##' @param target.ci
##' @param years
##' @return
##' @author Andrew Azman
fitAnnualBetaDelta <- function(params,
                               ss,
                               target.ci,
                               years){
    optim(params$beta.delta[1],
          fn=fitAnnualBetaDeltaObjFunc,
          ss=ss,params=params,target.ci=target.ci,years=years,
          method="Brent",lower=-10,upper=10,control=list(trace=T))
}

##' Returns data for a specific country for a specific year
##' @title
##' @return
##' @author Andrew Azman
getWHOStats <- function(target.country,years){
    dat <- read.csv("Data/TB_burden_countries_2012-12-10.csv")
    subset(dat,country == target.country & year %in% years)
}

##' Just to check that runsteady actually does what I hope it does
##' @param state
##' @param fun
##' @param params
##' @param check.every
##' @param var.beta
##' @return
##' @author Andrew Azman
runSteady <- function(state,fun,params,check.every=500,var.beta=FALSE){
    steady <- F
    while(!steady){
        tmp <- runTBHIVMod(params,state,check.every,var.beta=var.beta)
        if (abs(tail(tmp,10)[10] - tail(tmp,10)[1]) < 1){
            steady <- TRUE
        }
    }
    tail(tmp,1)[-1]
}


##' Fits increased theta to match a specific number increased cases detected in the first year
##' @param target.detection.increase number per 100k
##' @param duration
##' @param params
##' @param starting.state
##' @param ep.sn.muliplier
##' @param var.beta
##' @return
##' @author Andrew Azman
fitIncreasedDetectionRate <- function(target.detection.increase,
                                      duration,
                                      params,
                                      starting.state,
                                      ep.sn.multiplier,
                                      var.beta){

    optim(params$theta.spI[1]+.1,
          fn=fitIncreasedDetectionRateObjFunc,
          params=params,
          state=starting.state,
          duration=duration,
          ep.sn.multiplier=ep.sn.multiplier,
          target.detection.increase=target.detection.increase,
          var.beta=var.beta,method="Brent",lower=0,upper=10)
}

##' Objective function for fitting increased theta to increase in number of detected cases
##' @param theta.spI
##' @param params
##' @param state
##' @param duration
##' @param ep.sn.muliplier what percent of the sp rate increase shoudl be assigned to ep and sn?
##' @param var.beta
##' @param target.detection.increase
##' @return
##' @author Andrew Azman
fitIncreasedDetectionRateObjFunc <- function(theta.spI,
                                             params,
                                             state,
                                             duration,
                                             ep.sn.multiplier,
                                             var.beta,
                                             target.detection.increase){

    ## first run the model without an increased detection rate
    run.pre <- runTBHIVMod(params,state,duration,var.beta=var.beta)
    run.pre <- addColNames(run.pre,ext=T)
    last.time <- nrow(run.pre)

    ## now update the rates
    params$theta.spI <- rep(theta.spI,4)
    params$theta.snI <- rep(theta.spI,4)*ep.sn.multiplier
    params$theta.epI <- rep(theta.spI,4)*ep.sn.multiplier

    run.post <- runTBHIVMod(params,state,duration,var.beta=var.beta)
    run.post <- addColNames(run.post,ext=T)

                                        #how many additional cases are detected?
    cd.pre <- (run.pre[last.time,grep("N.Asp",colnames(run.pre))] +
               run.pre[last.time,grep("N.Asn",colnames(run.pre))] +
               run.pre[last.time,grep("N.Aep",colnames(run.pre))]) -
                   (run.pre[1,grep("N.Asp",colnames(run.pre))] +
                    run.pre[1,grep("N.Asn",colnames(run.pre))] +
                    run.pre[1,grep("N.Aep",colnames(run.pre))])

    cd.post <- (run.post[last.time,grep("N.Asp",colnames(run.post))] +
               run.post[last.time,grep("N.Asn",colnames(run.post))] +
               run.post[last.time,grep("N.Aep",colnames(run.post))]) -
                   (run.post[1,grep("N.Asp",colnames(run.post))] +
                    run.post[1,grep("N.Asn",colnames(run.post))] +
                    run.post[1,grep("N.Aep",colnames(run.post))])
#    cat(sprintf("pre = %.0f \n post = %.0f, \n increase = %.3f \n",sum(cd.pre),sum(cd.post),params$theta.spI[1]))
    ((sum(cd.post) - sum(cd.pre)) - target.detection.increase )^2
}

##' Calculates ICER for the output of intervention and counterfactual run
##' @param out output from runIntCont
##' @param eval.times - times to extract (in units of 1/10 year) and to analysze
##' @param dtx.cost - cost of finding cases in the first year (total - NOT per case)
##' @param tx.cost
##' @param tx.cost.mdr
##' @param tx.suc
##' @param tx.cost.partial
##' @param tx.cost.partial.mdr
##' @param discount
##' @param dis.wt.tx
##' @param dis.wt.tb
##' @param pct.mdr
##' @param params
calcICERFixedCosts <- function(out,
                               eval.times=1:11,
                               dtx.cost=20*100, #full cost in year 1
                               tx.cost=120,
                               tx.cost.mdr=120,
                               tx.suc=c(1),
                               tx.cost.partial=80,
                               tx.cost.partial.mdr=80,
                               discount=.03,
                               dis.wt.tx = c((0.331+0)/2,(0.399+0.221)/2,0.547,(0.399+0.053)/2), ## Weighted averages from solomon et al 2013
                               dis.wt.tb = c(0.331,0.399,0.547,0.399), ##using DB for AIDs only for HIV/TB from salomon et al 2013
                               pct.mdr = 0.023,
                               params){
    require(plyr)
    ## number of age classes (can put this as a param later)
    ac <- 1
    ## reduce to output for times of interest

    ## helper vectors
    types <- c("Asp","Asn","Aep")
    hivstatus <- c("","h","a","n")

    ## extract only the relavent evaluation times
    out <- lapply(out,function(x) x[eval.times,])
    ## get the times vector
    times <- out[[1]][,1][-1]

    ## get differences in stats over time
    diffs <- lapply(out,function(x) {
        diff(x)[,14:ncol(x)]
    })

    ## get unit costs through time
    ## NB: dtx.costs are total costs through time NOT per case
    ## taking Reimann integral here assuming step size of 0.1
    dtx.costs <- dtx.cost*exp(-times*discount)*0.1

    ## how many did we detect
    dtxs.int <- diffs$int[,grep("N.A(sp|sn|ep)",colnames(diffs$int))]
    dtxs.cont <- diffs$cont[,grep("N.A(sp|sn|ep)",colnames(diffs$cont))]

    ## get our costs discounted through time for DS and MDR TB
    tx.unitcosts <- tx.cost*exp(-times*discount)
    tx.part.unitcosts <- tx.cost.partial*exp(-times*discount)

    tx.unitcosts.mdr <- tx.cost.mdr*exp(-times*discount)
    tx.part.unitcosts.mdr <- tx.cost.partial.mdr*exp(-times*discount)

    ## Now we get the cost of full and partial treatment over time
    ## subtracting the control costs from the intervetion costs for case finding and then adding the diagnosis costs.
    ## for the treatment group
    txcost.int <- sum((rowSums(dtxs.int) - rowSums(dtxs.cont))*
                      (tx.suc*(tx.unitcosts*(1-pct.mdr) + tx.unitcosts.mdr*pct.mdr) +
                       (1-tx.suc)*(tx.part.unitcosts*(1-pct.mdr) + tx.part.unitcosts.mdr*pct.mdr)))
    ## Deaths

    death.cont <- diffs$cont[,grep("Mtb",colnames(diffs$cont))]
    death.int <- diffs$int[,grep("Mtb",colnames(diffs$int))]

    ## Years of life lost by hiv class
    ## taking a conservative approach where people
    ## can at most contribute horizon - time.step years
    YLL.cont <-  apply(death.cont,2,function(hiv.class) {
         hiv.class * (max(times) - times) * exp(-discount *times)
        ## hiv.class * which.max(times) - 1:length(times) * exp(-discount *times)

    })


    YLL.int <-  apply(death.int,2,function(hiv.class) {
        hiv.class * (max(times) - times) * exp(-discount *times)
       # hiv.class * which.max(times) - 1:length(times) * exp(-discount *times)

    })

    YLL.cont.minus.int <- YLL.cont - YLL.int

    ## from the model not accounting for deaths
    ## only considering symtomatic time not PS period
    ## NOTE: dis.dur.tx is not actually used anywhere anymore
    with(params,{
        dur.sp <- (theta.sp+theta.spI)*eta.sp+zeta.sp
        dur.sn <- (theta.sn+theta.snI)*eta.sn+zeta.sn
        dur.ep <- (theta.ep+theta.epI)*eta.ep+zeta.sn
        tmp <- 1/rbind(sp=dur.sp,sn=dur.sn,ep=dur.ep)
        colnames(tmp) <- c("","h","n","a")
        tmp
    }) -> dis.dur.tx

    with(params,{
        dur.sp <- theta.sp+zeta.sp
        dur.sn <- theta.sn+zeta.sn
        dur.ep <- theta.ep+zeta.sn
        tmp <- 1/rbind(sp=dur.sp,sn=dur.sn,ep=dur.ep)
        colnames(tmp) <- c("","h","n","a")
        tmp
    }) -> dis.dur.notx

    # taking mean treatment duration
    # assuming that all TB types have same duration of TX
    tx.dur <- 1/params$gamma.tx.rtx[1]

    ## Disability Years YLD = I * D * DW
    ## may need to split this by TB type
    prop.each.TB.type <- sapply(1:4,function(x) {
        with(params, c(pi.sp[x]*(1-pi.ep[x]),(1-pi.sp[x])*(1-pi.ep[x]),
                       pi.ep[x]))
    })

    colnames(prop.each.TB.type) <- c("","h","n","a")

    ## We consider prevalent cases are those contributing to YLD
    hiv.types <- c("^","h","a","n")

    ## list with each element as hiv type. matrices rows = time, columns = Asp, Asn, Aep
    prev.cases.cont <- sapply(1:4,function(x)
                                        (out[[2]])[,grep(paste0(hiv.types[x],"(Asp|Aep|Asn)"),
                                                                    colnames(out[[2]]))],simplify=F)
    prev.cases.int <- sapply(1:4,function(x)
                             (out[[1]])[,grep(paste0(hiv.types[x],"(Asp|Aep|Asn)"),
                                              colnames(out[[1]]))],simplify=F)
    prev.cases.cont.minus.int <- llply(1:4,function(x) prev.cases.cont[[x]] - prev.cases.int[[x]])

    ## these output lists of matrices. list by HIV, matrix columns by TB
    time.step <- 0.1

    YLD.notx.cont <- sapply(1:4,function(hiv){
        sum(sapply(1:3,function(tb) {
            prev.cases.cont[[hiv]][,tb] *
                time.step * dis.wt.tb[hiv]
        }))})

    YLD.notx.int <- sapply(1:4,function(hiv){
        sum(sapply(1:3,function(tb) {
            prev.cases.int[[hiv]][,tb] *
                time.step * dis.wt.tb[hiv]
        }))})

    YLD.notx.cont.minus.int <- YLD.notx.cont - YLD.notx.int

    #just getting them into a different form
    det.cases.int <- sapply(1:4,function(x){
        dtxs.int[,grep(paste0(hiv.types[x],"(N.Asp|N.Aep|N.Asn)"),
                       colnames(dtxs.int))]}
                            ,simplify=F)

    det.cases.cont <- sapply(1:4,function(x){
        dtxs.cont[,grep(paste0(hiv.types[x],"(N.Asp|N.Aep|N.Asn)"),
                        colnames(dtxs.cont))]}
                             ,simplify=F)

    ## NB: not discounting for time on treatment since it is SO short
    YLD.tx.int <-  sapply(1:4,function(hiv){
        sum(sapply(1:3,function(tb){
            det.cases.int[[hiv]][,tb] * pmin(0.5,max(times) - times) * dis.wt.tb[hiv]
        }))})

    YLD.tx.cont <-  sapply(1:4,function(hiv){
        sum(sapply(1:3,function(tb){
            det.cases.cont[[hiv]][,tb] * pmin(0.5,max(times) - times) * dis.wt.tb[hiv]
        }))})

    YLD.tx.cont.minus.int <- YLD.tx.cont - YLD.tx.int
    DALYs.int <- sum(YLL.int) + sum(YLD.notx.int) + sum(YLD.tx.int)
    DALYs.cont <- sum(YLL.cont) + sum(YLD.notx.cont) + sum(YLD.tx.cont)

    DALYs.averted <- sum(YLL.cont.minus.int) + sum(YLD.notx.cont.minus.int) + sum(YLD.tx.cont.minus.int)

    ret <- c(txcost.int=txcost.int,
             ICER=(txcost.int+sum(dtx.costs))/sum(DALYs.averted),
             DALYs.averted = DALYs.averted,
             DALYS.int = DALYs.int,
             DALYs.cont = DALYs.cont)

    return(ret)
}


##' @title
##' @param start.beta
##' @param state
##' @param params
##' @param target.ci
##' @return
##' @author Andrew Azman
ManualTuneBeta <- function(start.beta,state,params,target.ci){
    params$beta.sp <- start.beta
    RS <- runsteady(y=state[1:61],fun=dxdt.TBHIV.CI,parms=params,times=c(0,10000),verbose=F)
    run <- runTBHIVMod(params,initial.state=RS$y,max.time=1,var.beta=FALSE)
    getTBStats(run,add.names=T)
}


##' Makes fixed ICER plot
##' @param icer.min
##' @param icer.max
##' @param case.dt.dif fit to thi many numebr of cases detected int he first year
##' @param plot.params
##' @param start.state
##' @param tx.cost
##' @param tx.cost.partial
##' @param tx.cost.mdr
##' @param pct.mdr
##' @param tx.cost.partial.mdr
##' @param my.title
##' @param intcont.run if we give it this output from runIncCont we don't do it automatically
##' @param gdp per capita GDP
##' @param ICERS icers calaculated over the grid (calculated and returned in this function if not provided)
##' @param contours
##' @param xlab
##' @param ylab
##' @param leg legend?
##' @param ep.sn.multiplier
##' @param max.icer.cutoff
##' @param ...
##' @return list with ICERS and int.cont.run and makes plot
##' @author Andrew Azman
makeICERPlotFixed <- function(icer.min=0.00001,
                              icer.max=6000,
                              case.dt.dif,
                              plot.params  = india2011_params,
                              start.state = start.state.2011,
                              tx.cost = 81,
                              tx.cost.partial = tx.cost*.75,
                              tx.cost.mdr = 350,
                              pct.mdr = 0.023, # default for india
                              tx.cost.partial.mdr = tx.cost.mdr*.75,
                              my.title =  "",
                              intcont.run,
                              gdp,
                              ICERS,
                              contours,
                              xlab="",
                              ylab="",
                              leg=FALSE,
                              ep.sn.multiplier=1,
                              truncate.color=TRUE,
                              ...
                              ){

    ## fitting increased detetion rate that will give us X additional cases in the first year
    if (missing(intcont.run)){
        cat(sprintf("Fitting increased detection rate for %d case increase in year 1 \n",case.dt.dif))
        fit.tmp <- fitIncreasedDetectionRate(target.detection.increase = case.dt.dif,
                                             duration = 1,
                                             params = plot.params,
                                             starting.state = start.state,
                                             ep.sn.multiplier = ep.sn.multiplier,
                                             var.beta=FALSE)
        theta.reduction <- fit.tmp$par
        tmp  <- runIntCont(start.state,plot.params,10,
                           int.theta.sp= theta.reduction,
                           int.theta.sn = theta.reduction*ep.sn.multiplier,
                           int.theta.ep = theta.reduction*ep.sn.multiplier)
        plot.params$theta.spI <- rep(theta.reduction,4)
        plot.params$theta.snI <- rep(theta.reduction,4)*ep.sn.multiplier
        plot.params$theta.epI <- rep(theta.reduction,4)*ep.sn.multiplier

    } else {
        tmp <- intcont.run
    }

    times <- seq(1,10,by=.1)
    costs <- seq(50,5000,by=5)*case.dt.dif
    xlabs <- 1:10
    ylabs <- seq(50,5000,by=350)
    zlims <- c(0,log10(icer.max))
#    zlims <- c(icer.min,icer.max)
#    breaks <- seq(icer.min,icer.max,by=50)

    breaks <- seq(0,log10(icer.max),length=50)
    cols <- colorRampPalette(brewer.pal(9, name="Greens"))
    grid <- expand.grid(times,costs)

    # params for the mapply statement
    args.for.mapply <- list(params=plot.params,
                            out=tmp,
                            tx.cost=tx.cost,
                            tx.cost.partial=tx.cost.partial,
                            tx.cost.mdr=tx.cost.mdr,
                            tx.cost.partial.mdr=tx.cost.partial.mdr,
                            pct.mdr=pct.mdr,
                            fixed=TRUE)

    #only estiamte if we didn't supply ICERS
    if (missing(ICERS)) ICERS <- mapply(getICER,horiz=grid[,1],cost=grid[,2],MoreArgs=args.for.mapply)

    mat <- matrix(ICERS,nrow=length(times),ncol=length(costs))

                                        # if truncate color then we are going to set all values larger to icer.max
    if (truncate.color) {
        mat[which(mat > icer.max)] <- icer.max
        mat[which(mat < icer.min)] <- icer.min
    }

    ## par(mar=c(5,4.5,4,7))

    image(log10(mat),col=cols(length(breaks)-1),axes=F,xlab=xlab,ylab=ylab,zlim=zlims,breaks=breaks)

    if (!missing(gdp)){
        contour(log10(mat),levels=c(log10(0.0001),log10(gdp),log10(3*gdp)),col=addAlpha("black",.5),labcex=.5,lwd=1,lty=2,add=T,drawlabels=TRUE,method="edge",labels=c("cost saving","highly cost effective","cost effective"))
    }

    if (!missing(contours)){
        contour(log10(mat),
                levels=log10(contours), #[[1]]
                col=addAlpha("black",.5),
                labcex=.5,
                lwd=1,
                lty=2,
                labels=contours,
                add=TRUE)
                                        #,method="edge")
    }

    time.labs <- cbind(seq(0,1,length=length(times)),seq(1,10,length=length(times)))[seq(1,length(times),by=5),]
    axis(1,at=time.labs[,1],labels=time.labs[,2])
                                        #    axis(1,at=seq(0,1,length=length(xlabs)),labels=xlabs)
    costs.labs <- cbind(seq(0,1,length=length(costs)),costs/case.dt.dif)[seq(1,991,by=50),]
    axis(2,at=costs.labs[,1],labels=costs.labs[,2])

    #axis(2,at=seq(0,1,length=length(ylabs)),labels=ylabs)
    if (leg){
        legend.seq <- round(seq(min(zlims),max(zlims),length=5),0)
        image.plot(col=cols(length(breaks)-1),zlim=zlims,
                   ## breaks=seq(min(zlims),max(zlims),length=length(breaks)),
                   ## lab.breaks=round(10^seq(min(zlims),max(zlims),length=length(breaks)),0),
                   legend.only=T,horizontal=F,width=7,smallplot = c(.95,1,.05,.9),
                   axis.args=list(at=legend.seq, labels=10^legend.seq))

     }
    title(my.title)

    list("ICERS"=ICERS,"intcont.run"=tmp)
}

##' Plots Overview of Outputs from runIntCont
##' @param intcont
##' @param legend
##' @param by.TB - not implemented yet
##' @return
##' @author Andrew Azman
plotTBIncMort <- function(intcont,
                           legend=TRUE,
                           col1=1,
                           col2=2,
                          cd=FALSE,...){
    #CI, Prev # mortality, retx

    times <- intcont[[1]][,1]

    ci1 <- diff(intcont[[1]][,"CIall"])*10
    ci2 <- diff(intcont[[2]][,"CIall"])*10

    ## ##now prevalence
    ## prev1 <- rowSums(getPrevCols(intcont[[1]]))
    ## prev2 <- rowSums(getPrevCols(intcont[[2]]))

    ##mortality
    mort1 <- diff(rowSums(intcont[[1]][,grep("(n|a|h|^)(Mtb1)",colnames(intcont[[1]]))]))*10
    mort2 <- diff(rowSums(intcont[[2]][,grep("(n|a|h|^)(Mtb1)",colnames(intcont[[2]]))]))*10

    ## ##retreatment
    ## retx1 <- diff(rowSums(intcont[[1]][,grep("(n|a|h|^)(ReTx1)",colnames(intcont[[1]]))]))*10
    ## retx2 <- diff(rowSums(intcont[[2]][,grep("(n|a|h|^)(ReTx1)",colnames(intcont[[2]]))]))*10

    ## ## cases detected
    if (cd){
        cases.detected.1 <- diff(rowSums(intcont[[1]][,grep("N.A(sp|sn|ep)",colnames(intcont[[1]]))]))*10
        cases.detected.2 <- diff(rowSums(intcont[[2]][,grep("N.A(sp|sn|ep)",colnames(intcont[[2]]))]))*10
    }

    all.data.points <- c(ci1,ci2,mort1,mort2)
    if (cd) all.data.points <- c(all.data.points,cases.detected.1,cases.detected.2)

    plot(-100,-100,xlim=range(times),ylim=c(min(all.data.points),max(all.data.points)),xlab="",...)
    lines(times[-1],ci1,col=col1)
    lines(times[-1],ci2,col=col1,lty=2)

    ## lines(times,prev1,col=2)
    ## lines(times,prev2,col=2,lty=2)

    lines(times[-1],mort1,col=col2)
    lines(times[-1],mort2,col=col2,lty=2)

    ## lines(times[-1],retx1,col=4)
    ## lines(times[-1],retx2,col=4,lty=2)

    if (cd){
        lines(times[-1],cases.detected.1,col=5)
        lines(times[-1],cases.detected.2,col=5,lty=2)
    }

    if (legend){
        legend("topright",paste0(rep(c("CI","mort"),each=2),c(" - Interv."," - Baseline")),
               lty=rep(1:2,2),
               col=rep(1:2,each=2),
               bty="n")
    }
}

##' adds alpha to a set of colors
##' @title
##' @param COLORS
##' @param ALPHA
##' @return
addAlpha <- function(COLORS, ALPHA){
    if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
    RGB <- col2rgb(COLORS, alpha=TRUE)
    RGB[4,] <- round(RGB[4,]*ALPHA)
    NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
    return(NEW.COLORS)
}



plotCumTBIncMort <- function(intcont,
                             legend=TRUE,
                             col1=1,
                             col2=2,
                             diffs=FALSE,
                             poly=TRUE,
                             ...){

    times <- intcont[[1]][,1]

    ci1 <- intcont[[1]][,"CIall"]*10
    ci2 <- intcont[[2]][,"CIall"]*10

    ##mortality
    #mort1 <- rowSums(intcont[[1]][,grep("(n|a|h|^)(Mtb1)",colnames(intcont[[1]]))])*10
    #mort2 <- rowSums(intcont[[2]][,grep("(n|a|h|^)(Mtb1)",colnames(intcont[[2]]))])*10

    all.data.points <- c(ci1,ci2)#,mort1,mort2)

    if (diffs){
        plot(ci2 - ci1,col=col1,lty=6,...)
    } else {
        plot(-100,-100,xlim=range(times),ylim=c(0,max(all.data.points)),xlab="",...)
        lines(times,ci1,col=col1,lty=1)
        lines(times,ci2,col=col1,lty=2)
        if (poly){
            polygon(x=c(times,rev(times)),y=c(ci1,rev(ci2)),col=addAlpha(col1,.2),border=FALSE)
        }

    }

    ## lines(times,mort1,col=col2)
    ## lines(times,mort2,col=col2,lty=2)

}


##' Compares stats from model output to those from WHO
##' @param run
##' @param year
##' @param country
##' @return
compareStats <- function(run,year,country){
    tb.hiv.stats <- c(getTBStats(run),hIVStats(addColNames(run,ext=T)))
    who.stats <- getWHOStats(country,year)
    #mort

    who.stat.colnames <- c("e_mort_exc_tbhiv_100k","e_mort_exc_tbhiv_100k_lo","e_mort_exc_tbhiv_100k_hi",

                           "e_prev_100k","e_prev_100k_lo","e_prev_100k_hi",
                           "e_inc_100k","e_inc_100k_lo","e_inc_100k_hi",
                           "e_inc_tbhiv_100k","e_inc_tbhiv_100k_lo","e_inc_tbhiv_100k_hi")
    cbind(who.stats[who.stat.colnames])
}

##' Runs a short term ACF intervention then continues on for some years
##' @param country string with "india", "sa", or "china"
##' @param pct.incidence extra cases found in year one should be pct.incidence X incidence
##' @param int.dur total number of years we want to run the intervention
##' @param total.dur total number of years we want to run the smiluation
##' @param fits named (by country) list of fitted objects
##' @return intcont list for simulation
##' @author Andrew Azman
runNYearACF <- function(country,
                        pct.incidence,
                        case.dt.dif,
                        int.dur=2,
                        total.dur=10,
                        fits){
    #require(Hmisc)
    ## number of cases detecgted in year 1 proportional to incidence
    if (missing(case.dt.dif)){
        case.dt.dif <- c(round(getWHOStats("China",2011)[,"e_inc_100k"]*pct.incidence,0),
                         round(getWHOStats("India",2011)[,"e_inc_100k"]*pct.incidence,0),
                         round(getWHOStats("South Africa",2011)[,"e_inc_100k"]*pct.incidence,0))
    }

    case.dt.dif <- switch(country,
                          "india" = case.dt.dif[2],
                          "china" = case.dt.dif[1],
                          "sa" = case.dt.dif[3])

    fit.tmp <- fitIncreasedDetectionRate(target.detection.increase = case.dt.dif,
                                         duration = 1,
                                         params = fits[[country]]$params,
                                         starting.state = fits[[country]]$state,
                                         ep.sn.multiplier = 1,
                                         var.beta=FALSE)

    theta.reduction <- fit.tmp$par

    return(runIntCont(ss=fits[[country]]$state,
               params=fits[[country]]$params,
               time=total.dur,
               int.theta.sp=theta.reduction,
               int.theta.sn=theta.reduction*1,
               int.theta.ep=theta.reduction*1,
               intervention.duration = int.dur))
}



## Sens/Uncertainty Analyses Functions


##' Makes list of update functions for every param (or non-param) in the params list
##' used for running sesntivity analyses and dealing with dependent params
##' @return list of params suitable for use in the models
##' @author Andrew Azman
makeUpFuncs <- function(){
    up.funcs <- vector("list",length=148)

    up.funcs[[1]] <-
        update.func <- function(para,new.value) {
            para$beta.sp <- rep(new.value,4)
            para
        }
    up.funcs[[2]] <-
        update.func <- function(para,new.value) {
            para
        }
    up.funcs[[3]] <-
        update.func <- function(para,new.value) {
            para
    }

    up.funcs[[4]] <-
        update.func <- function(para,new.value) {
            para
        }
    up.funcs[[5]] <-
        update.func <- function(para,new.value) {
            para$phi.sn <- rep(new.value,4)
            para
        }
    up.funcs[[6]] <-
        update.func <- function(para,new.value) {
        para
    }
    up.funcs[[7]] <-
        update.func <- function(para,new.value) {
            para
        }
    up.funcs[[8]] <-
        update.func <- function(para,new.value) {
            para
        }
    up.funcs[[9]] <-
    update.func <- function(para,new.value) {
        para$phi.l[1] <- new.value
        para$phi.l[c(2,4)] <- new.value*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$phi.l[3]
        para
    }
    up.funcs[[10]] <-
        update.func <- function(para,new.value) {
            para
        }
    up.funcs[[11]] <-
        update.func <- function(para,new.value) {
            para$phi.l[3] <- new.value
            para$phi.l[c(2,4)] <- para$phi.l[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*new.value
            para
        }
    up.funcs[[12]] <-
        update.func <- function(para,new.value) {
            para
        }
    up.funcs[[13]] <-
        update.func <- function(para,new.value) {
            para$phi.ps[1:4] <- new.value
            para
        }
    up.funcs[[14]] <-
        update.func <- function(para,new.value) {
            para
    }
    up.funcs[[15]] <-
        update.func <- function(para,new.value) {
            para
        }
    up.funcs[[16]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[17]] <-
    update.func <- function(para,new.value) {
        para$gamma.lf.ls[1:4] <- new.value
        para
    }
up.funcs[[18]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[19]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[20]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[21]] <-
    update.func <- function(para,new.value) {
        para$gamma.rtx.ls[1:4] <- new.value
        para
    }
up.funcs[[22]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[23]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[24]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[25]] <-
    update.func <- function(para,new.value) {
        para$gamma.tx.rtx[1:4] <- new.value
        para
    }
up.funcs[[26]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[27]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[28]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[29]] <-
    update.func <- function(para,new.value) {
        para$rho.lf[1] <- new.value
        para$rho.lf[c(2,4)] <- new.value*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$rho.lf[3]
        para
    }
up.funcs[[30]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[31]] <-
    update.func <- function(para,new.value) {
        para$rho.lf[3] <- new.value
        para$rho.lf[c(2,4)] <- para$rho.lf[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$rho.lf[3]
        para
    }
up.funcs[[32]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[33]] <-
    update.func <- function(para,new.value) {
        para$rho.ls[1] <- new.value
        para$rho.ls[c(2,4)] <- new.value*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$rho.ls[3]
        para
    }
up.funcs[[34]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[35]] <-
    update.func <- function(para,new.value) {
        para$rho.ls[3] <- new.value
        para$rho.ls[c(2,4)] <- para$rho.ls[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$rho.ls[3]
        para
    }
up.funcs[[36]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[37]] <-
    update.func <- function(para,new.value) {
        para$rho.rel[1:4] <- new.value
        para
    }
up.funcs[[38]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[39]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[40]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[41]] <-
    update.func <- function(para,new.value) {
        para$rho.ps[1] <- new.value
        para$rho.ps[c(2,4)] <- para$rho.ps[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$rho.ps[3]
        para
    }
up.funcs[[42]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[43]] <-
    update.func <- function(para,new.value) {
        para$rho.ps[3] <- new.value
        para$rho.ps[c(2,4)] <- para$rho.ps[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$rho.ps[3]
        para
    }
up.funcs[[44]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[45]] <-
    update.func <- function(para,new.value) {
        para$pi.sp[1] <- new.value
        para$pi.sp[c(2,4)] <- para$pi.sp[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$pi.sp[3]
        para
    }
up.funcs[[46]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[47]] <-
    update.func <- function(para,new.value) {
        para$pi.sp[3] <- new.value
        para$pi.sp[c(2,4)] <- para$pi.sp[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$pi.sp[3]
        para
    }
up.funcs[[48]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[49]] <-
    update.func <- function(para,new.value) {
        para$pi.ep[1] <- new.value
        para$pi.ep[c(2,4)] <- para$pi.ep[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$pi.ep[3]
        para
    }
up.funcs[[50]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[51]] <-
    update.func <- function(para,new.value) {
        para$pi.ep[3] <- new.value
        para$pi.ep[c(2,4)] <- para$pi.ep[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$pi.ep[3]
        para
    }
up.funcs[[52]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[53]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[54]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[55]] <-
    update.func <- function(para,new.value) {
        para$mu.sp[3] <- new.value
        para$mu.sp[c(2,4)] <- para$mu.sp[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$mu.sp[3]
        para
    }
up.funcs[[56]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[57]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[58]] <-
    update.func <- function(para,new.value) {
        para
    }
up.funcs[[59]] <-
    update.func <- function(para,new.value) {
        para$mu.sn[3] <- new.value
        para$mu.sn[c(2,4)] <- para$mu.sn[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$mu.sn[3]
        para
    }
up.funcs[[60]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[61]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[62]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[63]] <-
    update.func <- function(para,new.value) {
        para$mu.ep[3] <- new.value
        para$mu.ep[c(2,4)] <- para$mu.ep[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$mu.ep[3]
        para
    }

up.funcs[[64]] <-
    update.func <- function(para,new.value) {
        para
    }

## zeta.sps
up.funcs[[65]] <-
    update.func <- function(para,new.value) {
        para$zeta.sp[1] <- new.value
        para$zeta.sp[c(2,4)] <- para$zeta.sp[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$zeta.sp[3]
        para$mu.sp[1] <- 1/3 - new.value
        para$mu.sp[c(2,4)] <- para$mu.sp[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$mu.sp[3]
        para
    }

up.funcs[[66]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[67]] <-
    update.func <- function(para,new.value) {
        para$zeta.sp[3] <- new.value
        para$zeta.sp[c(2,4)] <- para$zeta.sp[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$zeta.sp[3]
        para
    }

up.funcs[[68]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[69]] <-
    update.func <- function(para,new.value) {
        para$zeta.sn[1] <- new.value
        para$zeta.sn[c(2,4)] <- para$zeta.sn[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$zeta.sn[3]
        para$mu.sn[1] <- 1/3 - new.value
        para$mu.sn[c(2,4)] <- para$mu.sn[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$mu.sn[3]
        para
    }

up.funcs[[70]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[71]] <-
    update.func <- function(para,new.value) {
        para$zeta.sn[3] <- new.value
        para$zeta.sn[c(2,4)] <- para$zeta.sn[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$zeta.sn[3]
        para
    }

up.funcs[[72]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[73]] <-
    update.func <- function(para,new.value) {
        para$zeta.ep[1] <- new.value
        para$zeta.ep[c(2,4)] <- para$zeta.ep[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$zeta.ep[3]
        para$mu.ep[1] <- 1/3 - new.value
        para$mu.ep[c(2,4)] <- para$mu.ep[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$mu.ep[3]
        para
    }

up.funcs[[74]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[75]] <-
    update.func <- function(para,new.value) {
        para$zeta.ep[3] <- new.value
        para$zeta.ep[c(2,4)] <- para$zeta.ep[1]*para$`ART mulitplier`[1] + (1-para$`ART mulitplier`[1])*para$zeta.ep[3]
        para
    }

up.funcs[[76]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[77]] <-
    update.func <- function(para,new.value) {
        para$theta.sp[1:4] <- new.value
        para
    }

up.funcs[[78]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[79]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[80]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[81]] <-
    update.func <- function(para,new.value) {
        para$theta.sn[1:4] <- new.value
        para
    }

up.funcs[[82]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[83]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[84]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[85]] <-
    update.func <- function(para,new.value) {
        para$theta.ep[1:4] <- new.value
        para
    }

up.funcs[[86]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[87]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[88]] <-
    update.func <- function(para,new.value) {
        para
    }

for (i in 89:(89+(4*6)-1)){
    up.funcs[[i]] <-
        update.func <- function(para,new.value) {
            para
        }
}

up.funcs[[113]] <-
    update.func <- function(para,new.value) {
        para$foi.hiv[1] <- new.value
        para
    }

for (i in c(114:116,117,119,120,121,122,124,129:((129+4*4)-1),146:148)){
    up.funcs[[i]] <-
        update.func <- function(para,new.value) {
            para
        }
}

up.funcs[[118]] <-
    update.func <- function(para,new.value) {
        para$chi.elg[2] <- new.value
        para
    }

up.funcs[[123]] <-
    update.func <- function(para,new.value) {
        para$chi.tx[3] <- new.value
        para
    }

up.funcs[[125]] <-
    update.func <- function(para,new.value) {
        para
    }

up.funcs[[126]] <-
    update.func <- function(para,new.value) {
        para$mu.hiv[2] <- new.value
        para
    }

up.funcs[[127]] <-
    update.func <- function(para,new.value) {
        para$mu.hiv[3] <- new.value
        para
    }

up.funcs[[128]] <-
    update.func <- function(para,new.value) {
        para$mu.hiv[4] <- new.value
        para
    }

up.funcs[[145]] <-
    update.func <- function(para,new.value) {
        para$`ART mulitplier`[1:4] <- new.value
        para
    }
    up.funcs
}

##' helper function to generate array of parameters for sensitivty analyses
##' @param fits
##' @param country
##' @param p
##' @param seq.lengths
##' @param true.param.index
##' @return
genParamSeqs <- function(fits,country,
                         p=max.pct.change,
                         seq.lengths=num.points,
                         true.param.index=true.param.index){
    param.seq.array <- array(dim=c(seq.lengths,length(true.param.index)))
    for (i in seq_along(true.param.index)){
        orig.value <- c(t(do.call(rbind,fits[[country]]$params)))[true.param.index[i]]
        if (i %in% c(16:19,38)){ ## 38 is the ART multiplier
            param.seq.array[,i] <- seq(orig.value*p,min(orig.value*(1+p),1),length=seq.lengths)
        } else {
            param.seq.array[,i] <- seq(orig.value*p,orig.value*(1+p),length=seq.lengths)
        }
    }
    param.seq.array
}


##' For running on-way sensitivity analyses
##' @param country
##' @param fits
##' @param max.pct.change
##' @param num.points
##' @param cost.per.case
##' @param analytic.horizon
##' @param min.tx.costs
##' @param max.tx.costs
##' @param min.mdr.tx.costs
##' @param max.mdr.tx.costs
##' @return
##' @author Andrew Azman
runOneWaySens <- function(country,
                          fits,
                          max.pct.change,
                          num.points=5,
                          cost.per.case=2000,
                          analytic.horizon = 5,
                          min.tx.costs,
                          max.tx.costs,
                          min.mdr.tx.costs,
                          max.mdr.tx.costs
                          ){

    up.funcs <- makeUpFuncs()
    true.params <-1 - sapply(up.funcs,function(x) all.equal(c(do.call(rbind,x(fits[[country]]$params,-10))),c(do.call(rbind,fits[[country]]$params))) == TRUE)
    true.param.index <- which(true.params == 1)

    original.values <- c(t(do.call(rbind,fits[[country]]$params)))[true.param.index]
    seq.lengths <- num.points
    fits.orig <- fits
    out <- array(dim=c(seq.lengths,length(true.param.index)+2))
    ## 1. Let's first explore how the ICER for fixed cost per case detected in a single country varies by parameter
    param.array <- genParamSeqs(fits.orig,country,
                                p=max.pct.change,
                                seq.lengths = num.points,
                                true.param.index=true.param.index)

    ## get number of cases that will be detected
    pct.increase.in.yr1 <- 0.25
    cases.detected <- getIncreasedCasesDetected(TRUE,pct.increase.in.yr1)
    ## define ranges for parameters
    for (j in 1:ncol(param.array)){
        param.seq <- param.array[,j]
        for (i in seq_along(param.seq)){
            ## update param and any additional dependent params (e.g. HIV states)
            new.params <- up.funcs[[true.param.index[j]]](fits.orig[[country]]$params,param.seq[i])
            fits[[country]]$params <- new.params
            ## run 2 year ACF
            ## not we are not useing pct.incidence here as it is overridden by case.dt.fid
            run <- runNYearACF(country,pct.incidence = 0.15,case.dt.dif=cases.detected,int.dur = 2,total.dur = 10,fits=fits)
            ## Calculate and store ICER
            out[i,j] <- calcICERFixedCosts(out=run,
                                           eval.times = 1:(10*analytic.horizon + 1),
                                           dtx.cost=cases.detected[country]*cost.per.case,
                                           tx.suc=c(1),
                                           tx.cost = tx.cost.pc[country],
                                           tx.cost.partial = tx.cost.partial.pc[country],
                                           tx.cost.mdr = tx.cost.mdr.pc[country],
                                           pct.mdr= pct.mdr.pc[country],
                                           tx.cost.partial.mdr = tx.cost.partial.mdr[country],
                                           params=fits[[country]]$params)[2]
        }
        ## next paramater value
    }

    ## now for costs
    tx.costs <- seq(min.tx.costs,max.tx.costs,length=seq.lengths)
    mdr.tx.costs <- seq(min.mdr.tx.costs,max.mdr.tx.costs,length=seq.lengths)
    for (i in 1:seq.lengths){
        run <- runNYearACF(country,pct.incidence = 0.15,case.dt.dif=cases.detected,int.dur = 2,total.dur = 10,fits=fits.orig)
        ## Calculate and store ICER
        out[i,j+1] <- calcICERFixedCosts(out=run,
                                         eval.times = 1:(10*analytic.horizon + 1),
                                         dtx.cost=cases.detected[country]*cost.per.case,
                                         tx.suc=c(1),
                                         tx.cost = tx.costs[i],
                                         tx.cost.partial = tx.cost.partial.pc[country],
                                         tx.cost.mdr = tx.cost.mdr.pc[country],
                                         pct.mdr= pct.mdr.pc[country],
                                         tx.cost.partial.mdr = tx.cost.partial.mdr[country],
                                         params=fits[[country]]$params)[2]

        out[i,j+2] <- calcICERFixedCosts(out=run,
                                         eval.times = 1:(10*analytic.horizon + 1),
                                         dtx.cost=cases.detected[country]*cost.per.case,
                                         tx.suc=c(1),
                                         tx.cost = tx.cost.pc[country],
                                         tx.cost.partial = tx.cost.partial.pc[country],
                                         tx.cost.mdr = mdr.tx.costs[i],
                                         pct.mdr= pct.mdr.pc[country],
                                         tx.cost.partial.mdr = tx.cost.partial.mdr[country],
                                         params=fits[[country]]$params)[2]
    }

    param.array <- cbind(param.array,tx.costs,mdr.tx.costs)
    list(out,param.array)
}


##' @param sens.mat
##' @param param.array
##' @param country string of country
##' @param fits.orig
##' @param analytic.horizon
##' @param cost.per.case
##' @param lwd
##' @param top.n.params
##' @return pdf of tornado plot
##' @author Andrew Azman
makeTornadoPlot <- function(sens.mat,
                            param.array,
                            country,
                            fits.orig,
                            analytic.horizon,
                            cost.per.case,
                            lwd=10,
                            top.n.params=10){

    param.index.names <- rep(names(fits.orig[[country]]$params),each=4)
    param.names <- as.matrix(read.csv("Data/param_names.csv",as.is=T,header=F))
    param.names <- paste0(rep(param.names,each=4)," [",0:3,"]")

    up.funcs <- makeUpFuncs() # get functions that help update parameters
    true.params <-1 - sapply(up.funcs,function(x) all.equal(c(do.call(rbind,x(fits[[country]]$params,-10))),c(do.call(rbind,fits[[country]]$params))) == TRUE)
    true.param.index <- which(true.params == 1)
    original.values <- c(c(t(do.call(rbind,fits[[country]]$params)))[true.param.index],tx.cost.pc[country],tx.cost.mdr.pc[country])

    pdf(sprintf("Figures/oneway_sens_%s_%.fyr_%.fusd.pdf",country,analytic.horizon,cost.per.case),width=5,height=4)
    out <- sens.mat
    run <- runNYearACF(country,pct.incidence = 0.5,
                       case.dt.dif=case.dt.dif,int.dur = 2,total.dur = 10,fits=fits.orig)
    icer.orig <- calcICERFixedCosts(out=run,
                                    eval.times = 1:(10*analytic.horizon + 1),
                                    dtx.cost=case.dt.dif[country]*cost.per.case,
                                    tx.suc=c(1),
                                    tx.cost = tx.cost.pc[country],
                                    tx.cost.partial = tx.cost.partial.pc[country],
                                    tx.cost.mdr = tx.cost.mdr.pc[country],
                                    pct.mdr= pct.mdr.pc[country],
                                    tx.cost.partial.mdr = tx.cost.partial.mdr[country],
                                    params=fits.orig[[country]]$params)[2]
    cat(print(icer.orig))

    layout(matrix(c(1,1,1,2,2,2,2,2,1,1,1,2,2,2,2,2),nrow=2,byrow=T))
    par(mar=c(4.5,1,0,0))
    xlims <-  c(min(out),max(out))

    plot(-100,-100,xlim=xlims,ylim=c(0,1),bty="n",yaxt="n",ylab="",xlab="Cost per DALY Averted (USD)")#ncol(param.array)))
    abline(v=icer.orig,col="grey",lty=2)

    ## sort by extremes
    param.order <- order(apply(out,2,function(x) range(x)[2] - range(x)[1]))
    sorted.out <- out[,param.order]
    y.increment <- 1/min(ncol(out),top.n.params)
    start.iter <- ifelse(top.n.params > ncol(out),1,ncol(out) - top.n.params) # do we start the below iterations from the lowest params?
    for (param in start.iter:ncol(out)){
        tmp.out <- sorted.out[,param]
        greater.than.orig <- param.array[,param] > original.values[param]
        extremes <- range(tmp.out)
        print(range(tmp.out))
        max.col <- ifelse(greater.than.orig[which.max(tmp.out)],"red","blue")
        min.col <- ifelse(max.col == "red","blue","red")
        lines(x=c(extremes[1],icer.orig),y=c((param-start.iter)*y.increment,(param-start.iter)*y.increment),lwd=lwd,lend="butt",col=min.col)
        lines(x=c(icer.orig,extremes[2]),y=c((param-start.iter)*y.increment,(param-start.iter)*y.increment),lwd=lwd,lend="butt",col=max.col)
        }

    text(par("usr")[2]-par("usr")[2]*.1,.1,"High Value",col="red",cex=1)
    text(par("usr")[2]-par("usr")[2]*.1,.14,"Low Value",col="blue",cex=1)

    ## plot ranges for each
    ## plot(-100,-100,axes=F,bty="n",xlim=c(-1,1),ylim=c(0,1),xlab="",ylab="")
    ranges <- apply(param.array,2,range)
    ranges <- apply(ranges,2,function(x) sprintf("(%.2f,%.2f)",x[1],x[2]))
    ## for  (param in 1:ncol(out)) text(.5,(param-start.iter)*y.increment,ranges[param.order[param]],cex=1.1)


    ## plot names of each
    par(mar=c(4.5,0,0,0))
    plot(-100,-100,axes=F,bty="n",xlim=c(-1,1),ylim=c(0,1),xlab="",ylab="")
    for  (param in 1:ncol(out)) text(1,(param-start.iter)*y.increment,
                                     sprintf("%s %s",param.names[true.param.index[param.order[param]]],ranges[param.order[param]]),cex=.9,pos=4,offset=-22)
   dev.off()
}


##' @param nsims
##' @param country
##' @param param_range_file
##' @param output_file
##' @return saves (1) list of run outputs and (2) list of parameters lists
##' @author Andrew Azman
runLHS <- function(nsims=10,
                   country="sa",
                   param_range_prefix="uncer_ranges_",
                   output_file_prefix="uncer_out",
                   case.dt.dif=case.dt.dif,
                   orig.fits=fits,
                   per.person.dx.cost=seq(1000,35000,length=300)
                   ){
    require(tgp)
    ## load in transformation functiosn that deal with dependent params
    up.funcs <- makeUpFuncs()
    params.minmax <- as.matrix(read.csv(paste0("Data/",param_range_prefix,country,".csv"),row.names=1),ncol=4)
    true.params <-1 -sapply(up.funcs,function(x) all.equal(
        unlist(x(orig.fits[[country]]$params,-10)),
        unlist(orig.fits[[country]]$params)) == TRUE)

    true.param.index <- which(true.params == 1)
    param.names <- paste0(rep(names(orig.fits[[country]]$params),each=4),
                          rep(c("_n","_h","_hArt","_hNoART"),
                              length(orig.fits[[country]]$params)))

    ## make the lhs draws
    lhs.draws <- lhs(n=nsims,
                      params.minmax[,2:3],
                      shape=rep(3,nrow(params.minmax)),
                      mode=params.minmax[,1])

    runs <- list("vector",nsims)
    new.params <- list("vector",nsims)
    ## Run a two year ACF and store the results only if
    ## I don't think we are doing the following anymore but left the comment in:
    ## incidence in baseline scenario at year 10 is orig.I <= I_10 <= orig.I*.5
    for (i in 1:nrow(lhs.draws)){
        if (i %% 100 == 0) cat(".")
        ## make the parameter list
        new.params[[i]] <- updateParams(new.values=lhs.draws[i,],
                                         param.indices=true.param.index,
                                         countr=country,
                                         fits=orig.fits)
        tmp.fits <- orig.fits
        (tmp.fits[[country]]$params <- new.params[[i]])
        runs[[i]] <- runNYearACF(country,
                                 pct.incidence=.15,
                                 case.dt.dif=case.dt.dif,
                                 int.dur = 2,
                                 total.dur = 10,
                                 fits=tmp.fits)
    }
    ## going to store as a list of runs

    unix.time.stamp <- sprintf("%.0f",as.numeric(Sys.time()))

    save(runs,file=paste0(output_file_prefix,"_",country,"_runs_",unix.time.stamp,".rda"))
    save(new.params,file=paste0(output_file_prefix,"_",country,"_params_",unix.time.stamp,".rda"))
    save(lhs.draws,file=paste0(output_file_prefix,"_",country,"_lhsdraws_",unix.time.stamp,".rda")) #this is a matrix of the LHS samples and includes the cost
    ## save(runs,file=paste0(output_file_prefix,"_",country,"_runs_",Sys.Date(),".rda"))
    ## save(new.params,file=paste0(output_file_prefix,"_",country,"_params_",Sys.Date(),".rda"))
    ## save(lhs.draws,file=paste0(output_file_prefix,"_",country,"_lhsdraws_",Sys.Date(),".rda")) #this is a matrix of the LHS samples and includes the cost

    horizons <- c(2,5,10)
    out <- array(dim=c(300,3,nsims))

    print(" \n post-processing \n")
    for (i in 1:nsims){
        cat("*")
        for (h in seq_along(horizons)){
            for (t in seq_along(per.person.dx.cost)){
                out[t,h,i] <-
                calcICERFixedCosts(out=runs[[i]],
                                   eval.times = 1:(horizons[h]*10+1),
                                   dtx.cost=case.dt.df[country]*per.person.dx.cost[t],
                                   tx.suc=c(1),
                                   tx.cost = tx.cost.pc[country],
                                   tx.cost.partial = tx.cost.partial.pc[country],
                                   tx.cost.mdr = tx.cost.mdr.pc[country],
                                   pct.mdr= pct.mdr.pc[country],
                                   tx.cost.partial.mdr = tx.cost.partial.mdr[country],
                                   params=new.params[[i]])[2]
            }
        }
    }

    save(out,file=paste0(output_file_prefix,"_",country,"_icers_",unix.time.stamp,".rda"))

}

##' Updates the parameter list for us with a set of new values from LHS
##' @param param.indices
##' @param new.values
##' @param country
##' @param fits
##' @return list of params suitable for model runs
##' @author Andrew Azman
updateParams <- function(new.values,param.indices,country,fits){
    up.funcs <- makeUpFuncs() # get functions that help update parameters

    param.tmp <- fits[[country]]$params
    ## for each parameter we will sequentiually update the parameter list
    ## ineffecient but a function of previous code I wrote.
    for (i in seq_along(param.indices)){
       param.tmp <- up.funcs[[param.indices[i]]](param.tmp,new.values[i])
    }
    param.tmp
}

##' gets the number of of cases that need to be detected for a number of cases equal to pct.first.yr% of either the
##' projected cases detected in the first year (case.det.based == TRUE), or incidence (case.det.base == FALSE).
##' @param case.det.based
##' @param pct.first.yr
##' @return named vector with number of cases for each country
##' @author Andrew Azman
getIncreasedCasesDetected <- function(case.det.based=TRUE,pct.first.yr=0.25){

    if (case.det.based){
	## let's try increasing the number of cases detected by x% of the modeled steady state / first year
	sa.trial <- runTBHIVMod(fit.sa.2011$params,fit.sa.2011$state,1,var.beta=F)
	india.trial <- runTBHIVMod(fit.india.2011$params,fit.india.2011$state,1,var.beta=F)
	china.trial <- runTBHIVMod(fit.china.2011$params,fit.china.2011$state,1,var.beta=F)

	case.dt.dif <- c("china"=round(sum(tail(china.trial[,grep("N.", colnames(india.trial))],1))*pct.first.yr,0),
                         "india"=round(sum(tail(india.trial[,grep("N.", colnames(india.trial))],1))*pct.first.yr,0),
                         "sa"=round(sum(tail(sa.trial[,grep("N.", colnames(india.trial))],1))*pct.first.yr,0))
    } else {
        ## incidence based
	case.dt.dif <- c("china"=round(getWHOStats("China",2011)[,"e_inc_100k"]*pct.first.yr,0),
                         "india"=round(getWHOStats("India",2011)[,"e_inc_100k"]*pct.first.yr,0),
                         "sa"=round(getWHOStats("South Africa",2011)[,"e_inc_100k"]*pct.first.yr,0))

    }

    return(case.dt.dif)
}
