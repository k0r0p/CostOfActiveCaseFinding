source("Code/ACF-base.R")

######################
## Fitting of India###
######################

## lets look at indias data
tb.burd <- read.csv("Data/TB_burden_countries_2012-12-10.csv")
head(tb.burd)
india.dat <- subset(tb.burd,country == "India")
plotWHOData("India")

##lets look at incidence only
slope.incidence <- with(subset(india.dat,year>2003),{
    plot(year,e_inc_100k)
    fit <- lm(e_inc_100k ~ year)
    lines(year,predict(fit),col=1)
    coef(fit)[2]
    })


## to make sure we are using the most updated version of the parameters
## this is a little app from an apple script that converts sheets to csv
system("open Code/excelsheets_to_csv.app/")

## going to fit to 2003 as steady state - for now assuming no population growth
params.india <- make.params("india")
ss <- c(1e5-184-300,100,100,100,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        rep(0,25))

## setting hiv bits to zero to get an initial fit
params.india$foi.hiv <- c(0,0,0,0)

##init.fit.india <- FitIncPrev(ss,params.india,212,target.prev.tb=383,lowers=c(4.7,.1),uppers=c(20,1.2))
init.fit.india <- fitIncCDR(ss,params.india,212,target.cdr=0.65)

## Fitting to 2004 data
fit.india.2004 <- iterativeHIVTBFit(init.fit.india$ss,
                                    params.india,
                                    target.ci=212,
                                    target.cdr=0.65,
                                        # Taking the upper end of the CI for 2011 to represent the large private sector
                                        # this should be a reasonable approximation for the situation then...
                                        # other options is to try and fit to prev and incidence
                                    target.prev.hiv=0.004,
                                    target.art=0.4,
                                    epsilon.target=1e-2)

## 2004 prevalence from India HIV estiamtes 2006
dput(fit.india.2004,"Data/Fits/india2004fit_run010.rda")
fit.india.2004 <- dget("Data/Fits/india2004fit_run010.rda")
test.run.india.2004 <- runTBHIVMod(fit.india.2004$params,fit.india.2004$state,10,FALSE)
tb.hiv.stats.2004 <- c(getTBStats(test.run.india.2004),hIVStats(addColNames(test.run.india.2004,ext=T)))
xtable(t(as.matrix(tb.hiv.stats.2004)))
GetWHOStats("India",2004)

## now what change in annual pct change in beta can give us the incidence we want in 2011
fit.india.2011 <- fitAnnualBetaDelta(fit.india.2004$params,fit.india.2004$state,target.ci=181,years=7)
params.india.2004.to.2011 <- fit.india.2004$params
params.india.2004.to.2011$beta.delta <- rep(fit.india.2011$par,4)
test.run.india.2004.to.2011 <- runTBHIVMod(params.india.2004.to.2011,fit.india.2004$state,7,var.beta=T)
(tb.hiv.stats.2011 <- c(getTBStats(test.run.india.2004.to.2011),
                        hIVStats(addColNames(test.run.india.2004.to.2011,ext=T))))
xtable(t(as.matrix(tb.hiv.stats.2011)))
getWHOStats("India",2011)

test.run.india.2004.to.2011 <- addColNames(test.run.india.2004.to.2011,ext=T)
start.state.2011 <- tail(test.run.india.2004.to.2011,10)[1,-1]
start.state.2011[37:length(start.state.2011)] <- 0

## need to update beta for 2011
params.india.2004.to.2011$beta.sp <- params.india.2004.to.2011$beta.sp*exp(params.india.2004.to.2011$beta.delta*7)
india2011_params <- params.india.2004.to.2011
dput(params.india.2004.to.2011,"Data/Fits/india2011_params.rda")
india2011_params <- dget("Data/Fits/india2011_params.rda")

## save fit
dput(list(params=params.india.2004.to.2011,state=start.state.2011),"Data/Fits/india2011fit.rda")

######################
## Fitting of China###
######################

## lets look at chinas data
tb.burd <- read.csv("Data/TB_burden_countries_2012-12-10.csv")
china.dat <- subset(tb.burd,country == "China")
# PlotWHOData("China")

##lets look at incidence only
slope.incidence <- with(subset(china.dat,year>2003),{
    plot(year,e_inc_100k)
    fit <- lm(e_inc_100k ~ year)
    lines(year,predict(fit),col=1)
    coef(fit)[2]
    })

## to make sure we are using the most updated version of the parameters
## this is a little app from an apple script that converts sheets to csv
system("open Code/excelsheets_to_csv.app/")

## going to fit to 2004 as steady state - for now assuming no population growth
params.china <- make.params("china")
ss <- c(1e5-184-300,100,100,100,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        rep(0,25))

## setting hiv bits to zero to get an initial fit
params.china$foi.hiv <- c(0,0,0,0)
getWHOStats("China",2004)
init.fit.china <- fitIncCDR(ss,params.china,95,target.cdr=.64)
## 2009 estimate for ART coverage UNAID 31 to 67%
fit.china.2004 <- iterativeHIVTBFit(init.fit.china$ss,
                                    params.china,
                                    target.ci=95,
                                    target.cdr=0.67,
                                    target.prev.hiv=0.00058,
                                    target.art=0.5,
                                    epsilon.target=1e-2)
## assuming that art coverage was less in 2004.


## 2004 prevalence from China HIV estiamtes 2006
dput(fit.china.2004,"Data/Fits/china2004fit_run010.rda")
fit.china.2004 <- dget("Data/Fits/china2004fit_run010.rda")
test.run.china.2004 <- runTBHIVMod(fit.china.2004$params,fit.china.2004$state,10,FALSE)
tb.hiv.stats.2004 <- c(getTBStats(test.run.china.2004),hIVStats(AddColNames(test.run.china.2004,ext=T)))
xtable(t(as.matrix(tb.hiv.stats.2004)))
getWHOStats("China",2004)

## now what change in annual pct change in beta can give us the incidence we want in 2011
fit.china.2011 <- fitAnnualBetaDelta(fit.china.2004$params,fit.china.2004$state,target.ci=75,years=7)
params.china.2004.to.2011 <- fit.china.2004$params
params.china.2004.to.2011$beta.delta <- rep(fit.china.2011$par,4)
test.run.china.2004.to.2011 <- runTBHIVMod(params.china.2004.to.2011,fit.china.2004$state,7,var.beta=T)
(tb.hiv.stats.2011 <- c(getTBStats(test.run.china.2004.to.2011),hIVStats(addColNames(test.run.china.2004.to.2011,ext=T))))
xtable(t(as.matrix(tb.hiv.stats.2011)))
getWHOStats("China",2011)

test.run.china.2004.to.2011 <- addColNames(test.run.china.2004.to.2011,ext=T)
start.state.2011 <- tail(test.run.china.2004.to.2011,10)[1,-1]
start.state.2011[37:length(start.state.2011)] <- 0

## need to update beta for 2011
params.china.2004.to.2011$beta.sp <- params.china.2004.to.2011$beta.sp*exp(params.china.2004.to.2011$beta.delta*7)
china2011_params <- params.china.2004.to.2011
dput(params.china.2004.to.2011,"Data/Fits/china2011_params.rda")
china2011_params <- dget("Data/Fits/china2011_params.rda")

dput(list(params=params.china.2004.to.2011,state=start.state.2011),"Data/Fits/china2011fit.rda")

#############################
## Fitting of South Africa###
#############################

## South Africa fit for case detection work
## PlotWHOData("South Africa")
(stats <- getWHOStats("South Africa",2011))

## since we have increasing incidence we will assume that the counteractual scenario is that they are able to maintain incdidence and prevalance at current levels. In senstivity analsyis we may want to introduce changes to HIV and TB forces of infection.

## to make sure we are using the most updated version of the parameters
## this is a little app from an apple script that converts sheets to csv
system("open Code/excelsheets_to_csv.app/")

params.sa <- make.params("sa")
# getting a starting point for the population
ss <- c(1e5-184-300,100,100,100,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        20,1,1,1,
        rep(0,25))

prelim.fit.sa <- fitIncCDR(initial.state = ss,params= params.sa,target.ci = 993,target.cdr = 0.69)
fit.sa.2011 <- iterativeHIVTBFit(prelim.fit.sa$ss,
                                 prelim.fit.sa$final.pars,
                                 target.ci=993,
                                 target.cdr=0.69,
                                 target.prev.hiv=0.173, #unaids website
                                 target.art=0.752, # country progress report 2012
                                 epsilon.target=1e-2)

dput(fit.sa.2011,file="Data/Fits/sa2011fit_run10.rda")
fit.sa.2011 <- dget("Data/Fits/sa2011fit_run10.rda")
test.run.sa <- runTBHIVMod(fit.sa.2011$params,fit.sa.2011$state,1,var.beta=F)
(tb.hiv.stats.2011 <- c(getTBStats(test.run.sa),
                        hIVStats(addColNames(test.run.sa,ext=T))))
