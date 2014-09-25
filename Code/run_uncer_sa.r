## Uncertainty Analysis for Casefinding Paper
library(tgp)
source("Code/ACF-base.R")
set.seed(21062013)
######################################
## import current fits/params
######################################

fit.sa.2011 <- dget("Data/Fits/sa2011fit_run10.rda")
fit.sa.2011$state[37:length(fit.sa.2011$state)] <- 0 ## reset cumulative states at the end
fit.india.2011 <- dget("Data/Fits/india2011fit.rda")
fit.china.2011 <- dget("Data/Fits/china2011fit.rda")
fits <- list("sa"=fit.sa.2011,"india"=fit.india.2011,"china"=fit.china.2011)

## costs per country (pc)
tx.cost.pc <- c("sa"=1029,"india"=81,"china"=232)
tx.cost.partial.pc <- tx.cost.pc*0.75
tx.cost.mdr.pc <-  c("sa"=4672,"india"=4396,"china"=4396)
tx.cost.partial.mdr <- 0.75*tx.cost.mdr.pc

## pct MDR per country
pct.mdr.pc <- c("sa"=0.018,"india"=0.021,"china"=0.057) # WHO 2012 REPORT TABLE 2.3

##2011 gdp/capita for each country --> http://unstats.un.org/unsd/snaama/dnllist.asp
gdp.pc <- c("sa"=8090,"india"=1528,"china"=5439)

## countours to asdd per country
## data from Vassella et al PloS MEd but china is made up for now
## TO DO: are these really fair comparisons?
contours.pc <- list("sa"=c(100,250,500,1000,3000),
                    "india"=c(100,250,500,1000),
                    "china"=c(100,250,500,1000,2000))
## common vars
icer.min <- 1e-6
icer.max <- 10000

sa.trial <- runTBHIVMod(fit.sa.2011$params,fit.sa.2011$state,1,var.beta=F)
india.trial <- runTBHIVMod(fit.india.2011$params,fit.india.2011$state,1,var.beta=F)
china.trial <- runTBHIVMod(fit.china.2011$params,fit.china.2011$state,1,var.beta=F)

pct.first.yr <- 0.25 # pct increase in cases detected in the first year
case.dt.df <- c("china"=round(
                    sum(tail(china.trial[,grep("N.", colnames(india.trial))],1))*pct.first.yr,0),
                "india"=round(sum(tail(india.trial[,grep("N.", colnames(india.trial))],1))*pct.first.yr,0),
                "sa"=round(sum(tail(sa.trial[,grep("N.", colnames(india.trial))],1))*pct.first.yr,0))

nsims <- 5000
runLHS(nsims=nsims,
       country="sa",
       case.dt.dif=case.dt.df,
       orig.fits=fits,
       per.person.dx.cost=seq(1000,35000,length=300))

## horizons <- c(2,5,10)
## per.person.dx.cost.sa <- seq(1000,35000,length=300) # only for sa
## out <- array(dim=c(300,3,nsims))

## print("post-processing")

## for (i in 1:nsims){
##     cat("*")
##     for (h in seq_along(horizons)){
##         for (t in seq_along(per.person.dx.cost.sa)){
##             out[t,h,i] <-
##             calcICERFixedCosts(out=runs[[i]],
##                                eval.times = 1:(horizons[h]*10+1),
##                                dtx.cost=case.dt.df["sa"]*per.person.dx.cost.sa[t],
##                                tx.suc=c(1),
##                                tx.cost = tx.cost.pc["sa"],
##                                tx.cost.partial = tx.cost.partial.pc["sa"],
##                                tx.cost.mdr = tx.cost.mdr.pc["sa"],
##                                pct.mdr= pct.mdr.pc["sa"],
##                                tx.cost.partial.mdr = tx.cost.partial.mdr["sa"],
##                                params=new.params[[i]])[2]
##         }
##     }
## }

## saveRDS(out,"GeneratedData/uncer_out_sa_icers_2014-09-22.rda") #  runs

#source("Code/post_proc_uncer_sa.r")
