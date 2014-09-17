# these are common params shared by a lot of the anlayses for paper
## this file should be loaded by most ofther ACF paper files
source("Code/ACF-base.R")
## import current fits/params
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

## countours to add  per country
contours.pc <- list("sa"=c(100,250,500,1000,3000),
                    "india"=c(100,250,500,1000),
                    "china"=c(100,250,500,1000,2000))

## common vars
icer.min <- 1e-6
icer.max <- 10000

## get number of cases detected
case.dt.dif <- getIncreasedCasesDetected(case.det.base=TRUE,pct.first.yr=0.25)
