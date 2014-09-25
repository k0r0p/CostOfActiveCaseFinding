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
load("GeneratedData/uncer_out_india_icers_1411591204.rda")
out.india.1 <- out
load("GeneratedData/uncer_out_india_icers_1411592012.rda")
out.india.2 <- out
load("GeneratedData/uncer_out_india_icers_1411592241.rda")
out.india.3 <- out
load("GeneratedData/uncer_out_india_icers_1411592964.rda")
out.india.4 <- out

load("GeneratedData/uncer_out_india_runs_1411591204.rda")
runs.india.1 <- runs
load("GeneratedData/uncer_out_india_runs_1411592012.rda")
runs.india.2 <- runs
load("GeneratedData/uncer_out_india_runs_1411592241.rda")
runs.india.3 <- runs
load("GeneratedData/uncer_out_india_runs_1411592964.rda")
runs.india.4 <- runs

load("GeneratedData/uncer_out_india_params_1411591204.rda")
new.params.india.1 <- new.params
load("GeneratedData/uncer_out_india_params_1411592012.rda")
new.params.india.2 <- new.params
load("GeneratedData/uncer_out_india_params_1411592241.rda")
new.params.india.3 <- new.params
load("GeneratedData/uncer_out_india_params_1411592964.rda")
new.params.india.4 <- new.params

load("GeneratedData/uncer_out_india_lhsdraws_1411591204.rda")
lhs.india.1 <- lhs.draws
load("GeneratedData/uncer_out_india_lhsdraws_1411592012.rda")
lhs.india.2 <- lhs.draws
load("GeneratedData/uncer_out_india_lhsdraws_1411592241.rda")
lhs.india.3 <- lhs.draws
load("GeneratedData/uncer_out_india_lhsdraws_1411592964.rda")
lhs.india.4 <- lhs.draws


# CHINA UNCER
# need to combine 4 runs of 5k each
load("GeneratedData/uncer_out_china_icers_1411586898.rda")
out.china.1 <- out
load("GeneratedData/uncer_out_china_icers_1411586878.rda")
out.china.2 <- out
load("GeneratedData/uncer_out_china_icers_1411586597.rda")
out.china.3 <- out
load("GeneratedData/uncer_out_china_icers_1411586225.rda")
out.china.4 <- out

load("GeneratedData/uncer_out_china_runs_1411586898.rda")
runs.china.1 <- runs
load("GeneratedData/uncer_out_china_runs_1411586878.rda")
runs.china.2 <- runs
load("GeneratedData/uncer_out_china_runs_1411586597.rda")
runs.china.3 <- runs
load("GeneratedData/uncer_out_china_runs_1411586225.rda")
runs.china.4 <- runs

load("GeneratedData/uncer_out_china_params_1411586898.rda")
new.params.china.1 <- new.params
load("GeneratedData/uncer_out_china_params_1411586878.rda")
new.params.china.2 <- new.params
load("GeneratedData/uncer_out_china_params_1411586597.rda")
new.params.china.3 <- new.params
load("GeneratedData/uncer_out_china_params_1411586225.rda")
new.params.china.4 <- new.params

load("GeneratedData/uncer_out_china_lhsdraws_1411586898.rda")
lhs.china.1 <- lhs.draws
load("GeneratedData/uncer_out_china_lhsdraws_1411586878.rda")
lhs.china.2 <- lhs.draws
load("GeneratedData/uncer_out_china_lhsdraws_1411586597.rda")
lhs.china.3 <- lhs.draws
load("GeneratedData/uncer_out_china_lhsdraws_1411586225.rda")
lhs.china.4 <- lhs.draws

# South Africa UNCER
#out.sa <- out
load("GeneratedData/uncer_out_sa_icers_1411589299.rda")
out.sa.1 <- out
load("GeneratedData/uncer_out_sa_icers_1411590009.rda")
out.sa.2 <- out
load("GeneratedData/uncer_out_sa_icers_1411590038.rda")
out.sa.3 <- out
load("GeneratedData/uncer_out_sa_icers_1411589211.rda")
out.sa.4 <- out

load("GeneratedData/uncer_out_sa_runs_1411589299.rda")
runs.sa.1 <- runs
load("GeneratedData/uncer_out_sa_runs_1411590009.rda")
runs.sa.2 <- runs
load("GeneratedData/uncer_out_sa_runs_1411590038.rda")
runs.sa.3 <- runs
load("GeneratedData/uncer_out_sa_runs_1411589211.rda")
runs.sa.4 <- runs

load("GeneratedData/uncer_out_sa_params_1411589299.rda")
new.params.sa.1 <- new.params
load("GeneratedData/uncer_out_sa_params_1411590009.rda")
new.params.sa.2 <- new.params
load("GeneratedData/uncer_out_sa_params_1411590038.rda")
new.params.sa.3 <- new.params
load("GeneratedData/uncer_out_sa_params_1411589211.rda")
new.params.sa.4 <- new.params

load("GeneratedData/uncer_out_sa_lhsdraws_1411589299.rda")
lhs.sa.1 <- lhs.draws
load("GeneratedData/uncer_out_sa_lhsdraws_1411590009.rda")
lhs.sa.2 <- lhs.draws
load("GeneratedData/uncer_out_sa_lhsdraws_1411590038.rda")
lhs.sa.3 <- lhs.draws
load("GeneratedData/uncer_out_sa_lhsdraws_1411589211.rda")
lhs.sa.4 <- lhs.draws


#########################################################
## Cost Effectiveness Thresholds by Country and Horizon #
#########################################################

## India
per.person.dx.cost.india <- seq(50,20000,length=300) # per person detection cost
hce.india10 <- c(apply(out.india.1,3,function(x) approx(x[,3],per.person.dx.cost.india,xout=gdp.pc["india"])$y),
                 apply(out.india.2,3,function(x) approx(x[,3],per.person.dx.cost.india,xout=gdp.pc["india"])$y),
                 apply(out.india.3,3,function(x) approx(x[,3],per.person.dx.cost.india,xout=gdp.pc["india"])$y),
                 apply(out.india.4,3,function(x) approx(x[,3],per.person.dx.cost.india,xout=gdp.pc["india"])$y))

hce.india5 <- c(apply(out.india.1,3,function(x) approx(x[,2],per.person.dx.cost.india,xout=gdp.pc["india"])$y),
                apply(out.india.2,3,function(x) approx(x[,2],per.person.dx.cost.india,xout=gdp.pc["india"])$y),
                apply(out.india.3,3,function(x) approx(x[,2],per.person.dx.cost.india,xout=gdp.pc["india"])$y),
                apply(out.india.4,3,function(x) approx(x[,2],per.person.dx.cost.india,xout=gdp.pc["india"])$y))

hce.india2 <- c(apply(out.india.1,3,function(x) approx(x[,1],per.person.dx.cost.india,xout=gdp.pc["india"])$y),
                apply(out.india.2,3,function(x) approx(x[,1],per.person.dx.cost.india,xout=gdp.pc["india"])$y),
                apply(out.india.3,3,function(x) approx(x[,1],per.person.dx.cost.india,xout=gdp.pc["india"])$y),
                apply(out.india.4,3,function(x) approx(x[,1],per.person.dx.cost.india,xout=gdp.pc["india"])$y))

(i10 <- quantile(hce.india10,c(0.025,0.975))) # 850.3 2042.7
(i5 <- quantile(hce.india5,c(0.025,0.975))) # 690.5 1151.5
(i2 <- quantile(hce.india2,c(0.025,0.975))) #274.8 342.5

## now for china
per.person.dx.cost.china <- seq(50,20000,length=300) # per person detection cost
hce.china10 <- c(apply(out.china.1,3,function(x) approx(x[,3],per.person.dx.cost.china,xout=gdp.pc["china"])$y),
                 apply(out.china.2,3,function(x) approx(x[,3],per.person.dx.cost.china,xout=gdp.pc["china"])$y),
                 apply(out.china.3,3,function(x) approx(x[,3],per.person.dx.cost.china,xout=gdp.pc["china"])$y),
                 apply(out.china.4,3,function(x) approx(x[,3],per.person.dx.cost.china,xout=gdp.pc["china"])$y))

hce.china5 <- c(apply(out.china.1,3,function(x) approx(x[,2],per.person.dx.cost.china,xout=gdp.pc["china"])$y),
                apply(out.china.2,3,function(x) approx(x[,2],per.person.dx.cost.china,xout=gdp.pc["china"])$y),
                apply(out.china.3,3,function(x) approx(x[,2],per.person.dx.cost.china,xout=gdp.pc["china"])$y),
                apply(out.china.4,3,function(x) approx(x[,2],per.person.dx.cost.china,xout=gdp.pc["china"])$y))

hce.china2 <- c(apply(out.china.1,3,function(x) approx(x[,1],per.person.dx.cost.china,xout=gdp.pc["china"])$y),
                apply(out.china.2,3,function(x) approx(x[,1],per.person.dx.cost.china,xout=gdp.pc["china"])$y),
                apply(out.china.3,3,function(x) approx(x[,1],per.person.dx.cost.china,xout=gdp.pc["china"])$y),
                apply(out.china.4,3,function(x) approx(x[,1],per.person.dx.cost.china,xout=gdp.pc["china"])$y))


(c10 <- quantile(hce.china10,c(0.025,0.975))) #2706  6392
(c5 <- quantile(hce.china5,c(0.025,0.975))) # 2257  3731
(c2 <- quantile(hce.china2,c(0.025,0.975))) # 1030  1258

## South Africa
per.person.dx.cost.sa <- seq(1000,35000,length=300) # only for sa
hce.sa10 <- c(apply(out.sa.1,3,function(x) approx(x[,3],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y),
              apply(out.sa.2,3,function(x) approx(x[,3],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y),
              apply(out.sa.3,3,function(x) approx(x[,3],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y),
              apply(out.sa.4,3,function(x) approx(x[,3],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y))


hce.sa5 <- c(apply(out.sa.1,3,function(x) approx(x[,2],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y),
             apply(out.sa.2,3,function(x) approx(x[,2],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y),
             apply(out.sa.3,3,function(x) approx(x[,2],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y),
             apply(out.sa.4,3,function(x) approx(x[,2],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y))

hce.sa2 <- c(apply(out.sa.1,3,function(x) approx(x[,1],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y),
             apply(out.sa.2,3,function(x) approx(x[,1],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y),
             apply(out.sa.3,3,function(x) approx(x[,1],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y),
             apply(out.sa.4,3,function(x) approx(x[,1],per.person.dx.cost.sa,xout=gdp.pc["sa"])$y))

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
up.funcs <- makeUpFuncs()
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
param.draws <- sapply(new.params.china.1,function(x) unlist(x))
param.draws <- data.frame(cbind(t(param.draws[true.params == 1,]),
                                lhs.china.1[,(ncol(lhs.china.1)-1):ncol(lhs.china.1)]))
colnames(param.draws) <- pnames
china.pcc <- pcc(X=param.draws,y=hce.china10[1:5000],rank=T)
## order same as india
`pcc.china.sorted <- china.pcc$PRCC[order(abs(india.pcc$PRCC),decreasing=T),]
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
india.one.way <- runOneWaySens(country="india",
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

param.array <- genParamSeqs(fits=fits,country="india",p=0.5,seq.lengths=5,true.param.index=which(true.params == 1))
makeTornadoPlot(india.one.way[[1]],
                param.array=india.one.way[[2]],
                country="india",
                fits.orig=fits,
                analytic.horizon = 10,
                cost.per.case = 1200,
                top.n.params = 10,
                lwd=20)

## CHINA
china.one.way <- runOneWaySens(country="china",
                               fits=fits,
                               max.pct.change = 0.5,
                               num.points=5,cost.per.case = 3800,
                               analytic.horizon = 10,
                               min.tx.costs=tx.cost.pc["china"]*0.5,
                               max.tx.costs=tx.cost.pc["china"]*1.5,
                               min.mdr.tx.costs=tx.cost.mdr.pc["china"]*0.5,
                               max.mdr.tx.costs=tx.cost.mdr.pc["china"]*1.5)

makeTornadoPlot(china.one.way[[1]],
                param.array=china.one.way[[2]],
                country="china",
                fits.orig=fits,
                analytic.horizon = 10,
                cost.per.case = 3800,
                top.n.params = 10,
                lwd=20)

## South Africa
sa.one.way <- runOneWaySens(country="sa",
                            fits=fits,max.pct.change = 0.5,
                            num.points=5,
                            cost.per.case = 9408,
                            analytic.horizon = 10,
                            min.tx.costs=tx.cost.pc["china"]*0.5,
                            max.tx.costs=tx.cost.pc["china"]*1.5,
                            min.mdr.tx.costs=tx.cost.mdr.pc["china"]*0.5,
                            max.mdr.tx.costs=tx.cost.mdr.pc["china"]*1.5)

makeTornadoPlot(sa.one.way[[1]],
                param.array=sa.one.way[[2]],
                country="sa",
                fits.orig=fits,
                analytic.horizon = 10,
                cost.per.case = 9400,
                top.n.params = 10,
                lwd=20)

