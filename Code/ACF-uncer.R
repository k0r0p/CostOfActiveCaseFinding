## Uncertainty Analysis for Casefinding Paper
set.seed(21062013)
library(sensitivity)
######################################
## import current fits/params
######################################
source("Code/common-params.r")

######################################
## import multivariate uncertainty runs
######################################
## India
load("Data/uncer_out_india_2013-07-19.rda")
out.india <- out

## Bring in data from LHS (files run_uncer_country.r)
load("Data/uncer_out_india_runs_2013-07-07.rda") #  runs
runs.india <- runs
load("Data/uncer_out_india_params_2013-07-07.rda") # new.params
new.params.india <- new.params
load("Data/uncer_out_india_lhsdraws_2013-07-07.rda") # new.params
lhs.india <- lhs.draws

# CHINA UNCER
load("Data/uncer_out_china_2013-07-19.rda")
out.china <- out

load("Data/uncer_out_china_runs_2013-07-07.rda") #  runs
runs.china <- runs
load("Data/uncer_out_china_params_2013-07-07.rda") # new.params
new.params.china <- new.params
load("Data/uncer_out_china_lhsdraws_2013-07-07.rda") # new.params
lhs.china <- lhs.draws

# South Africa UNCER
load("Data/uncer_out_sa_2013-07-19.rda")
out.sa <- out

load("Data/uncer_out_sa_runs_2013-07-11.rda") #  runs
runs.sa <- runs
load("Data/uncer_out_sa_params_2013-07-11.rda") # new.params
new.params.sa <- new.params
load("Data/uncer_out_sa_lhsdraws_2013-07-11.rda") # new.params
lhs.sa <- lhs.draws

#########################################################
## Cost Effectiveness Thresholds by Country and Horizon #
#########################################################

## India
per.person.dx.cost.india <- seq(50,20000,length=350) # per person detection cost
hce.india10 <- apply(out.india,3,function(x) approx(x[,3],per.person.dx.cost.india,xout=gdp.pc["india"])$y)
hce.india5 <- apply(out.india,3,function(x) approx(x[,2],per.person.dx.cost.india,xout=gdp.pc["india"])$y)
hce.india2 <- apply(out.india,3,function(x) approx(x[,1],per.person.dx.cost.india,xout=gdp.pc["india"])$y)

(i10 <- quantile(hce.india10,c(0.025,0.975))) # 983.0887 2420.1954
(i5 <- quantile(hce.india5,c(0.025,0.975))) # 794.3113 1336.0736
(i2 <- quantile(hce.india2,c(0.025,0.975))) # 294.6613 401.5656

## now for china
per.person.dx.cost.china <- seq(50,20000,length=350) # per person detection cost
hce.china10 <- apply(out.china,3,function(x) approx(x[,3],per.person.dx.cost.china,xout=gdp.pc["china"])$y)
hce.china5 <- apply(out.china,3,function(x) approx(x[,2],per.person.dx.cost.china,xout=gdp.pc["china"])$y)
hce.china2 <- apply(out.china,3,function(x) approx(x[,1],per.person.dx.cost.china,xout=gdp.pc["china"])$y)

(c10 <- quantile(hce.china10,c(0.025,0.975))) #3181.682 7511.125
(c5 <- quantile(hce.china5,c(0.025,0.975))) # 2661.295 4391.536
(c2 <- quantile(hce.china2,c(0.025,0.975))) # 1181.318 1520.379

## South Africa
per.person.dx.cost.sa <- seq(2000,35000,length=500) # only for sa
hce.sa10 <- apply(out.sa,3,function(x) approx(x[,3],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y)
hce.sa5 <- apply(out.sa,3,function(x) approx(x[,2],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y)
hce.sa2 <- apply(out.sa,3,function(x) approx(x[,1],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y)

(s10 <- quantile(hce.sa10,c(0.025,0.975),na.rm=F)) #11422.91 24174.83
(s5 <- quantile(hce.sa5,c(0.025,0.975),na.rm=F)) # 9046.703 16420.433
(s2 <- quantile(hce.sa2,c(0.025,0.975),na.rm=F)) # 4130.252 6007.599

# export table for supplement
xtable(cbind(rbind(i10,i5,i2),rbind(c10,c5,c2),rbind(s10,s5,s2)))

#######################################################
## calculation of partial rank correlation coeffecients
#######################################################

## load the parameter names
pnames <- read.csv("Data/param_names_old.csv",as.is=T,header=F)
pnames <- sapply(1:nrow(pnames),function(x) paste(pnames[x,2],pnames[x,3],sep="."))
pnames <-gsub(",|-","",pnames)

orig.fits <- fits
## get the list of functions that update correlated parameters (hiv classes etc)
up.funcs <- MakeUpFuncs()
## inidicator for which params are true params and which are simply linked to others by ART multioplier or some assumption
true.params <-1 - sapply(up.funcs,function(x) all.equal(
    unlist(x(orig.fits[["india"]]$params,-10)),
    unlist(orig.fits[["india"]]$params)) == TRUE)

## PCC - India
country = "india"

param.draws<- sapply(new.params.india,function(x) unlist(x))
param.draws <- data.frame(cbind(t(param.draws[true.params == 1,]),lhs.india[,(ncol(lhs.india)-1):ncol(lhs.india)]))
colnames(param.draws) <- pnames
## run pcc
india.pcc <- pcc(X=param.draws,y=hce.india10,rank=T,nboot=0)
## sort results
pcc.ind.sorted <- india.pcc$PRCC[order(abs(india.pcc$PRCC),decreasing=T),]
names(pcc.ind.sorted) <- rownames(india.pcc$PRCC)[order(abs(india.pcc$PRCC),decreasing=T)]
xtable(data.frame(pcc.ind.sorted))

## PCC - China
country = "china"
param.draws <- sapply(new.params.china,function(x) unlist(x))
param.draws <- data.frame(cbind(t(param.draws[true.params == 1,]),lhs.china[,(ncol(lhs.china)-1):ncol(lhs.china)]))
colnames(param.draws) <- pnames
china.pcc <- pcc(X=param.draws,y=hce.china10,rank=T)
## order same as india
pcc.china.sorted <- china.pcc$PRCC[order(abs(india.pcc$PRCC),decreasing=T),]
names(pcc.china.sorted) <- rownames(china.pcc$PRCC)[order(abs(india.pcc$PRCC),decreasing=T)]
xtable(cbind(data.frame(pcc.ind.sorted),data.frame(pcc.china.sorted)))

## PCC - South Africa
country = "sa"
param.draws <- sapply(new.params.sa,function(x) unlist(x))
param.draws <- data.frame(cbind(t(param.draws[true.params == 1,]),lhs.sa[,(ncol(lhs.sa)-1):ncol(lhs.sa)]))
colnames(param.draws) <- pnames
sa.pcc <- pcc(X=param.draws,y=hce.sa10,rank=T)
## order same as india
pcc.sa.sorted <- sa.pcc$PRCC[order(abs(india.pcc$PRCC),decreasing=T),]
names(pcc.sa.sorted) <- rownames(sa.pcc$PRCC)[order(abs(india.pcc$PRCC),decreasing=T)]
xtable(cbind(data.frame(pcc.ind.sorted),data.frame(pcc.china.sorted),data.frame(pcc.sa.sorted)))





## library(epiR)  # just wanted to make sure we get the same results
## test2 <- epi.prcc(cbind(param.draws,hce.india10), sided.test = 2)

## Running one way sens alanyses
india.one.way <- RunOneWaySens(country="india",
                               fits=fits,
                               max.pct.change = 0.5,
                               num.points=5,
                               cost.per.case = 1200,
                               analytic.horizon = 10,
                               min.tx.costs=tx.cost.pc["india"]*0.5,
                               max.tx.costs=tx.cost.pc["india"]*1.5,
                               min.mdr.tx.costs=tx.cost.mdr.pc["india"]*0.5,
                               max.mdr.tx.costs=tx.cost.mdr.pc["india"]*1.5
                               )

param.array <- GenParamSeqs(fits=fits,country="india",p=0.5,seq.lengths=5,true.param.index=which(true.params == 1))
MakeTornadoPlot(india.one.way[[1]],
                param.array=india.one.way[[2]],
                country="india",
                fits.orig=fits,
                analytic.horizon = 10,
                cost.per.case = 1200,
                top.n.params = 10,
                lwd=20)

## CHINA
china.one.way <- RunOneWaySens(country="china",
                               fits=fits,
                               max.pct.change = 0.5,
                               num.points=5,cost.per.case = 3800,
                               analytic.horizon = 10,
                               min.tx.costs=tx.cost.pc["china"]*0.5,
                               max.tx.costs=tx.cost.pc["china"]*1.5,
                               min.mdr.tx.costs=tx.cost.mdr.pc["china"]*0.5,
                               max.mdr.tx.costs=tx.cost.mdr.pc["china"]*1.5)

MakeTornadoPlot(china.one.way[[1]],param.array=china.one.way[[2]],country="china",
                fits.orig=fits,analytic.horizon = 10,
                cost.per.case = 3800,
                top.n.params = 10,
                lwd=20)

## South Africa
sa.one.way <- RunOneWaySens(country="sa",fits=fits,max.pct.change = 0.5,num.points=5,cost.per.case = 11200,analytic.horizon = 10,
                            min.tx.costs=tx.cost.pc["china"]*0.5,
                            max.tx.costs=tx.cost.pc["china"]*1.5,
                            min.mdr.tx.costs=tx.cost.mdr.pc["china"]*0.5,
                            max.mdr.tx.costs=tx.cost.mdr.pc["china"]*1.5)

MakeTornadoPlot(sa.one.way[[1]],param.array=sa.one.way[[2]],country="sa",fits.orig=fits,analytic.horizon = 10,cost.per.case = 11200,top.n.params = 10,lwd=20)

