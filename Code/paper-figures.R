## this script makes the figures and generates numbers used in ACF paper
source("Code/common-params.r")  # bring in a bunch of common params and load base functions


##' Makes short ACF plots
##' @title
##' @param country
##' @return
##' @author Andrew Azman
make.fig.2 <- function(country="india"){
    require(Hmisc)

    case.dt.dif <- getIncreasedCasesDetected(TRUE,0.25)[country]

    ## labels to be used in plot
    country.label.name <- switch(country,
                                 "india" = "India",
                                 "china" = "China",
                                 "sa" = "South Africa")

    fit.tmp <- fitIncreasedDetectionRate(target.detection.increase = case.dt.dif,
                                         duration = 1,
                                         params = fits[[country]]$params,
                                         starting.state = fits[[country]]$state,
                                         ep.sn.multiplier = 1,
                                         var.beta=FALSE)

    theta.reduction <- fit.tmp$par

    year2_int_only<- runIntCont(ss=fits[[country]]$state,
                                params=fits[[country]]$params,
                                time=10,
                                int.theta.sp=theta.reduction,
                                int.theta.sn=theta.reduction*1,
                                int.theta.ep=theta.reduction*1,
                                intervention.duration = 2)

    out <- array(dim=c(99,5))
    colnames(out) <- c(strsplit("txcost.int,ICER,DALYs.averted,DALYs.int,DALYs.cont",",")[[1]])
    for (t in 3:101){
        out[t-2,] <-
            calcICERFixedCosts(out=year2_int_only,
                               eval.times = 1:t,
                               dtx.cost=50*100,
                               tx.suc=(1),
                               tx.cost = tx.cost.pc[country],
                               tx.cost.partial = tx.cost.partial.pc[country],
                               tx.cost.mdr = tx.cost.mdr.pc[country],
                               pct.mdr= pct.mdr.pc[country],
                               tx.cost.partial.mdr = tx.cost.partial.mdr[country],
                               params=fits[[country]]$params)
    }

    col1 <- rgb(0,90,0,maxColorValue=255) #for incdience
    col2 <- rgb(153,0,13,maxColorValue=255) # for mortality
    col3 <- rgb(140,81,10,maxColorValue=255)
                                        #    quartz("",width=5,height=6)
    pdf(sprintf("Figures/Figure2_%scases_%s.pdf",case.dt.dif,country),width=5,height=6)
    par(mfrow=c(3,1),mar=c(.5,2,.5,1.5),mgp=c(1,.2,0),tck=-.01,cex.axis=0.8,oma=c(2,0,1.5,0))

    plotTBIncMort(year2_int_only,xaxt="n",legend=F,cd=T,ylab="count per 100,000 / yr",col1=col1,col2=col2)
    polygon(x=c(0,2,2,0),y=c(-1e6,-1e6,20000,20000),bty="n",border=FALSE,col=AddAlpha("grey",.2))
    abline(h=seq(0,1500,by=50),col=AddAlpha("grey",.15))
    abline(v=seq(0,10,by=.5),col=AddAlpha("grey",.15))

    if (country != "sa") legend("topright",c("Baseline","Intervention"),lty=c(2,1),bty="n")
                                        #text(9.5,185,expression(underline(Intervention)),cex=.85)
    #text(6.5,185,expression(underline(Baseline)),cex=.85)
    #legend(9,185,c("","",""),lty=1,col=c(col1,col2,5),bty="n",cex=0.85)
    if (country == "india"){
        text(5,170,"Incidence")
        text(5,120,"Cases Detected")
        text(5,35,"Mortality")

    }

    if (country == "china"){
        text(5,65,"Incidence")
        text(5,47,"Cases Detected")
        text(5,15,"Mortality")

    }

    if (country == "sa"){
        legend(8,950,c("Baseline","Intervention"),lty=c(1,2),bty="n")

        text(5,1050,"Incidence")
        text(5,750,"Cases Detected")
        text(5,260,"Mortality")

    }

    text(0,par("usr")[4]*0.9,"A",cex=1.3)


    #text(8,165,c("Incidence \n Mortality \n Cases Detected"),cex=.85)
    #legend(6,185,c("","",""),lty=2,col=c(col1,col2,5),bty="n",cex=0.85)

    plotCumTBIncMort(year2_int_only,xaxt="n",legend=F,ylab="count per 100,000",col1=col1,col2=col2,poly=T)
    new.out <- out[seq((11-3),(101-3),by=10),]

    #lines(0:10,c(0,new.out[,"DALYs.cont"]),col=col3,lty=2)
    #lines(0:10,c(0,new.out[,"DALYs.int"]),col=col3,lty=1)
    #polygon(x=c(0:10,10:0),y=c(0,new.out[,"DALYs.cont"],rev(new.out[,"DALYs.int"]),0),col=AddAlpha(3,.2),border=FALSE)

    ## plot(3:101,out[,"DALYs.cont"],xlim=c(0,101),col=3,lty=2,type="l",xaxt="n",ylab="")
    ## lines(3:101,out[,"DALYs.int"],col=3,lty=1)

    polygon(x=c(0,2,2,0),y=c(-1e6,-1e6,2e6,2e6),bty="n",border=FALSE,col=AddAlpha("grey",.2))
    abline(v=seq(0,10,by=.5),col=AddAlpha("grey",.15))

    if (country == "india"){
        abline(h=seq(0,15000,by=500),col=AddAlpha("grey",.15))

        text(x=1,y=8000,"Intervention \n Period",cex=1)
        arrows(x0=.5,x1=0,y0=8000,y1=8000,length=.05)
        arrows(x0=1.5,x1=2,y0=8000,y1=8000,length=.05)

        text(5,11000,"Cumulative Incidence \n (Baseline)",cex=1)
        text(8,9000,"Cumulative Incidence \n (Intervention)",cex=1)
    }

    if (country == "china"){
        abline(h=seq(0,15000,by=500),col=AddAlpha("grey",.15))

        text(x=1,y=3000,"Intervention \n Period",cex=1)
        arrows(x0=.5,x1=0,y0=3000,y1=3000,length=.05)
        arrows(x0=1.5,x1=2,y0=3000,y1=3000,length=.05)

        text(5,4300,"Cumulative Incidence \n (Baseline)",cex=1)
        text(6,2400,"Cumulative Incidence \n (Intervention)",cex=1)
    }

    if (country == "sa"){
        abline(h=seq(0,12e4,by=5000),col=AddAlpha("grey",.15))

        text(x=1,y=4e4,"Intervention \n Period",cex=1)
        arrows(x0=.5,x1=0,y0=4e4,y1=4e4,length=.05)
        arrows(x0=1.5,x1=2,y0=4e4,y1=4e4,length=.05)

        text(5,6.5e4,"Cumulative Incidence \n (Baseline)",cex=1)
        text(6,3e4,"Cumulative Incidence \n (Intervention)",cex=1)
    }

    text(0,par("usr")[4]*0.9,"B",cex=1.3)


    ## text(6,15000,expression(underline(Intervention)),cex=.85)
    ## text(3,15000,expression(underline(Baseline)),cex=.85)
    ## legend(5.5,15000,c("",""),lty=1,col=c(col1,col3),bty="n",cex=0.85)
    ## text(4.5,13400,c("Cum. Incidence \n Cum. DALYs"),cex=.85)
    ## legend(2.5,15000,c("",""),lty=2,col=c(col1,col3),bty="n",cex=0.85)

    plotCumTBIncMort(year2_int_only,xaxt="n",legend=F,diffs=T,ylab="count per 100,000",type="l",col1=col1)
    lines(3:101,out[,"DALYs.cont"] - out[,"DALYs.int"],type="l",col=col3,lty=6)
    axis(1,at=seq(1,101,by=10),label=0:10)
    abline(h=seq(0,2000,by=250),col=AddAlpha("grey",.15))
    abline(v=seq(1,101,by=5),col=AddAlpha("grey",.15))
    polygon(x=c(1,21,21,1),y=c(-1e6,-1e6,2e6,2e6),bty="n",border=FALSE,col=AddAlpha("grey",.2))

    if (country=="india"){
        text(45,650,"Cases Averted \n by Intervention",cex=1)
        text(45,200,"DALYs Averted \n by Intervention",cex=1)
    }

    if (country == "china"){
        text(45,240,"Cases Averted \n by Intervention",cex=1)
        text(55,75,"DALYs Averted \n by Intervention",cex=1)

    }

    if (country == "sa"){
        text(45,6400,"Cases Averted \n by Intervention",cex=1)
        text(55,2000,"DALYs Averted \n by Intervention",cex=1)
    }

    text(0,par("usr")[4]*0.9,"C",cex=1.3)

                                        #    legend(21,2000,paste0(c("DALYs Averted"),c(" - Intervention"," - Baseline")),lty=1:2,col=3,bty="n",cex=0.85)
    mtext("Year",side=1,outer=T,adj=0.5,cex=.8,line=.6)
    mtext(sprintf("2-Year Active Case Finding Campaign in %s",country.label.name),side=3,outer=T,adj=0.5,cex=.85,line=-.4)

    dev.off()
}

## Short Term Campaigns
##' @title
##' @return
##' @author Andrew Azman
make.fig.3 <- function(){

    gdp.pc <- c("sa"=8090,"india"=1528,"china"=5439)
    pct.first.yr <- 0.25

    sa.trial <- runTBHIVMod(fit.sa.2011$params,fit.sa.2011$state,1,var.beta=F)
    india.trial <- runTBHIVMod(fit.india.2011$params,fit.india.2011$state,1,var.beta=F)
    china.trial <- runTBHIVMod(fit.china.2011$params,fit.china.2011$state,1,var.beta=F)

    case.dt.df <- c("china"=round(sum(tail(china.trial[,grep("N.", colnames(india.trial))],1))*pct.first.yr,0),
                    "india"=round(sum(tail(india.trial[,grep("N.", colnames(india.trial))],1))*pct.first.yr,0),
                    "sa"=round(sum(tail(sa.trial[,grep("N.", colnames(india.trial))],1))*pct.first.yr,0))

    country <- c("china","india","sa")


    per.person.dx.cost <- seq(50,15000,length=300) # per person detection cost
    horizons <- c(2,5,10) # time horizons to consider
    out <- array(dim=c(length(per.person.dx.cost),length(horizons),3))

    for (c in seq_along(country)){
        fit.tmp <- fitIncreasedDetectionRate(target.detection.increase = case.dt.df[c],
                                             duration = 1,
                                             params = fits[[country[c]]]$params,
                                             starting.state = fits[[country[c]]]$state,
                                             ep.sn.multiplier = 1,
                                             var.beta=FALSE)

        theta.reduction <- fit.tmp$par

        ## run intervention for 2 years but run model for 10
        year2_int_only<- runIntCont(ss=fits[[country[c]]]$state,
                                    params=fits[[country[c]]]$params,
                                    time=10,
                                    int.theta.sp=theta.reduction,
                                    int.theta.sn=theta.reduction,
                                    int.theta.ep=theta.reduction,
                                    intervention.duration = 2)

        for (h in seq_along(horizons)){
            for (t in seq_along(per.person.dx.cost)){
                out[t,h,c] <- calcICERFixedCosts(out=year2_int_only,
                                                 eval.times = 1:(horizons[h]*10+1),
                                                 dtx.cost=case.dt.df[c]*per.person.dx.cost[t],
                                                 tx.suc=c(1),
                                                 tx.cost = tx.cost.pc[country[c]],
                                                 tx.cost.partial = tx.cost.partial.pc[country[c]],
                                                 tx.cost.mdr = tx.cost.mdr.pc[country[c]],
                                                 pct.mdr= pct.mdr.pc[country[c]],
                                                 tx.cost.partial.mdr = tx.cost.partial.mdr[country[c]],
                                                 params=fits[[country[c]]]$params)[2]
            }
        }
    }

    pdf("Figures/Figure3.pdf",width=6.8,height=3)
    #quartz("",width=6.8,height=3)

    ## make everything per 1000 USD
    out <- out/1000 # make it per 1000 USD
    per.person.dx.cost <- per.person.dx.cost/1000
    gdp.pc <- gdp.pc/1000

    par(mfrow=c(1,3),mar=c(0.1,0.1,0.1,0.5),oma=c(4,3.1,2,0.5),mgp=c(1,.1,0),tck=0.01)
    icer.range <- range(out) #apply(out,3,function(x) range(x[,"ICER"])))
    plot(-1000,-1000,xlim=c(0,10),#range(per.person.dx.cost),
         ylim=c(0,10000/1000),xlab="Cost per Case Detected (USD)",ylab="ICER (10-yr Horizon)")
    abline(v=seq(0,15000/1000,by=1000/1000),col=AddAlpha("grey",.2))
    abline(h=seq(0,15000/1000,by=1000/1000),col=AddAlpha("grey",.2))

    ## first for 2 year horizon
    lines(per.person.dx.cost,out[,1,1],col=1)
    lines(c(-1000,approx(out[,1,1],per.person.dx.cost,xout=gdp.pc["china"])$y),rep(gdp.pc["china"],2),col=1,lty=2)
    ## text(600,gdp.pc["china"]+200,"GDP China",cex=.5,col=1)

    lines(per.person.dx.cost,out[,1,2],col=2)
    lines(c(-1000,approx(out[,1,2],
                         per.person.dx.cost,
                         xout=gdp.pc["india"])$y),
          rep(gdp.pc["india"],2),col=2,lty=2)
    ## text(500,gdp.pc["india"]+200,"GDP India",cex=.5,col=2)

    lines(per.person.dx.cost,out[,1,3],col=3)
    lines(c(-1000,approx(out[,1,3],per.person.dx.cost,xout=gdp.pc["sa"])$y),rep(gdp.pc["sa"],2),col=3,lty=2)
    ## text(2500,gdp.pc["sa"]+200,"GDP South Africa",cex=.5,col=3)
    ## legend("topright",paste0("Per Capita GDP: ",c("China","India","S. Africa")),col=1:3,lty=2,bty="n",cex=0.8)

    # draw lines outside of the plot box to indiate the cost per case detected
    par(xpd=NA)
    lines(rep(approx(out[,1,1],per.person.dx.cost,xout=gdp.pc["china"])$y,2),c(gdp.pc["china"],-1),col=1,lty=2)
    text(x=approx(out[,1,1],per.person.dx.cost,xout=gdp.pc["china"])$y+.5,y=-1.3,
         paste0("$",round(1000*approx(out[,1,1],per.person.dx.cost,xout=gdp.pc["china"])$y,0)),col=1,cex=1.2)

    ##india
    lines(rep(approx(out[,1,2],per.person.dx.cost,xout=gdp.pc["india"])$y,2),c(gdp.pc["india"],-1.5),col=2,lty=2)
        text(x=approx(out[,1,2],per.person.dx.cost,xout=gdp.pc["india"])$y+.5,y=-2,
         paste0("$",round(1000*approx(out[,1,2],per.person.dx.cost,xout=gdp.pc["india"])$y,0)),col=2,cex=1.2)

    lines(rep(approx(out[,1,3],per.person.dx.cost,xout=gdp.pc["sa"])$y,2),c(gdp.pc["sa"],-1.5),col=3,lty=2)
    text(x=approx(out[,1,3],per.person.dx.cost,xout=gdp.pc["sa"])$y+.5,y=-2,paste0("$",round(1000*approx(out[,1,3],per.person.dx.cost,xout=gdp.pc["sa"])$y,0)),col=3,cex=1.2)
    par(xpd=FALSE)

    text(200/1000,par("usr")[4]*0.95,"A",cex=1.5)
    legend("topright",c("India","China","South Africa"),text.col=c(2,1,3),lty=-1,pch=-1,bty="n")
    ## 5 year horzon
    plot(-1000,-1000,xlim=c(0,10),#range(per.person.dx.cost),
         ylim=c(0,10000/1000),xlab="",ylab="",yaxt="n")
    abline(v=seq(0,15000/1000,by=1000/1000),col=AddAlpha("grey",.2))
    abline(h=seq(0,15000/1000,by=1000/1000),col=AddAlpha("grey",.2))

    lines(per.person.dx.cost,out[,2,1],col=1)
    lines(c(-1000,approx(out[,2,1],per.person.dx.cost,xout=gdp.pc["china"])$y),rep(gdp.pc["china"],2),col=1,lty=2)
    #text(2000,gdp.pc["china"]+200,"GDP China",cex=.5,col=1)

    lines(per.person.dx.cost,out[,2,2],col=2)
    lines(c(-1000,approx(out[,2,2],per.person.dx.cost,xout=gdp.pc["india"])$y),rep(gdp.pc["india"],2),col=2,lty=2)
    #text(500,gdp.pc["india"]+200,"GDP India",cex=.5,col=2)

    lines(per.person.dx.cost,out[,2,3],col=3)
    lines(c(-1000,approx(out[,2,3],per.person.dx.cost,xout=gdp.pc["sa"])$y),rep(gdp.pc["sa"],2),col=3,lty=2)
    #text(8000,gdp.pc["sa"]+200,"GDP South Africa",cex=.5,col=3)

    par(xpd=NA)
    ## china
    lines(rep(approx(out[,2,1],per.person.dx.cost,xout=gdp.pc["china"])$y,2),c(gdp.pc["china"],-1),col=1,lty=2)
    text(x=approx(out[,2,1],per.person.dx.cost,xout=gdp.pc["china"])$y+0.5,y=-1.3,
         paste0("$",round(1000*approx(out[,2,1],per.person.dx.cost,xout=gdp.pc["china"])$y,0)),col=1,cex=1.2)

    ##india
    lines(rep(approx(out[,2,2],per.person.dx.cost,xout=gdp.pc["india"])$y,2),c(gdp.pc["india"],-1.5),col=2,lty=2)
    text(x=approx(out[,2,2],per.person.dx.cost,xout=gdp.pc["india"])$y+0.5,y=-2,
         paste0("$",round(1000*approx(out[,2,2],per.person.dx.cost,xout=gdp.pc["india"])$y,0)),col=2,cex=1.2)

    lines(rep(approx(out[,2,3],per.person.dx.cost,xout=gdp.pc["sa"])$y,2),c(gdp.pc["sa"],-1.5),col=3,lty=2)
    text(x=approx(out[,2,3],per.person.dx.cost,xout=gdp.pc["sa"])$y+0.5,y=-2,paste0("$",round(1000*approx(out[,2,3],per.person.dx.cost,xout=gdp.pc["sa"])$y,0))
         ,col=3,cex=1.2)

    par(xpd=FALSE)

    text(200/1000,par("usr")[4]*0.95,"B",cex=1.5)

    ## first for 10 year horizon
    plot(-1000,-1000,xlim=c(0,10),#range(per.person.dx.cost),
         ylim=c(0,10000/1000),xlab="",ylab="",yaxt="n")
    abline(v=seq(0,15000/1000,by=1000/1000),col=AddAlpha("grey",.2))
    abline(h=seq(0,15000/1000,by=1000/1000),col=AddAlpha("grey",.2))

    lines(per.person.dx.cost,out[,3,1],col=1)
    lines(c(-1000,approx(out[,3,1],per.person.dx.cost,xout=gdp.pc["china"])$y),rep(gdp.pc["china"],2),col=1,lty=2)
    #text(3000,gdp.pc["china"]+200,"GDP China",cex=.5,col=1)

    lines(per.person.dx.cost,out[,3,2],col=2)
    lines(c(-1000,approx(out[,3,2],per.person.dx.cost,xout=gdp.pc["india"])$y),rep(gdp.pc["india"],2),col=2,lty=2)
    #text(500,gdp.pc["india"]+200,"GDP India",cex=.5,col=2)

    lines(per.person.dx.cost,out[,3,3],col=3)
    lines(c(-1000,approx(out[,3,3],per.person.dx.cost,xout=gdp.pc["sa"])$y),rep(gdp.pc["sa"],2),col=3,lty=2)
    #text(12000,gdp.pc["sa"]+200,"GDP South Africa",cex=.5,col=3)

    par(xpd=NA)
    lines(rep(approx(out[,3,1],per.person.dx.cost,xout=gdp.pc["china"])$y,2),c(gdp.pc["china"],-1.0),col=1,lty=2)
    text(x=approx(out[,3,1],per.person.dx.cost,xout=gdp.pc["china"])$y+0.5,y=-1.3,
         paste0("$",round(1000*approx(out[,3,1],per.person.dx.cost,xout=gdp.pc["china"])$y,0)),col=1,cex=1.2)

    lines(rep(approx(out[,3,2],per.person.dx.cost,xout=gdp.pc["india"])$y,2),c(gdp.pc["india"],-1.5),col=2,lty=2)
    text(x=approx(out[,3,2],per.person.dx.cost,xout=gdp.pc["india"])$y+0.5,y=-2,
     paste0("$",round(1000*approx(out[,3,2],per.person.dx.cost,xout=gdp.pc["india"])$y,0)),col=2,cex=1.2)

    lines(rep(approx(out[,3,3],per.person.dx.cost,xout=gdp.pc["sa"])$y,2),c(gdp.pc["sa"],-1.5),col=3,lty=2)
    text(x=approx(out[,3,3],per.person.dx.cost,xout=gdp.pc["sa"])$y+0.5,y=-2,paste0("$",round(1000*approx(out[,3,3],per.person.dx.cost,xout=gdp.pc["sa"])$y,0))
     ,col=3,cex=1.2)
    par(xpd=FALSE)

    text(200/1000,par("usr")[4]*0.95,"C",cex=1.5)

    mtext("2-year horizon",side=3,out=T,adj=0.10)
    mtext("5-year horizon",side=3,out=T,adj=0.50)
    mtext("10-year horizon",side=3,out=T,adj=0.90)
    mtext("Cost per DALY Averted (1000 USD)",side=2,out=T,adj=0.5,line=1.5,cex=.8)
    mtext("Cost per Case Detected (1000 USD)",side=1,out=T,adj=0.5,line=2.8,cex=0.8)

    dev.off()
}

##' Generates ICER image plots
##' @param pdf
##' @param fn
##' @return
##' @author Andrew Azman
make.fig.4 <- function(pdf=T,fn="Figures/Figure4.pdf",
                       sa.file="Data/sa_icers_4July.rda",
                       china.file="Data/china_icers_4July.rda",
                       india.file="Data/india_icers_4July.rda"){
    if (is.null(sa.file)){
        ## get ICER data for sustained plots
        sa <- makeICERPlotFixed(icer.min = icer.min,
                                icer.max = icer.max,
                                case.dt.dif = case.dt.dif[3],
                                plot.params = fits[["sa"]]$params,
                                start.state = fits[["sa"]]$state,
                                tx.cost = tx.cost.pc["sa"],
                                tx.cost.partial = tx.cost.partial.pc["sa"],
                                tx.cost.mdr = tx.cost.mdr.pc["sa"],
                                pct.mdr= pct.mdr.pc["sa"],
                                tx.cost.partial.mdr = tx.cost.partial.mdr["sa"],
                                my.title = "",
                                gdp=gdp.pc["sa"],
                                contours=contours.pc["sa"][[1]])
        ## save(sa,file="Data/sa_icers_4July.rda")
    } else {
        load(file=sa.file)
    }

    if (is.null(china.file)){

        china <- makeICERPlotFixed(icer.min = icer.min,
                                   icer.max = icer.max,
                                   case.dt.dif = case.dt.dif[1],
                                   plot.params = fits[["china"]]$params,
                                   start.state = fits[["china"]]$state,
                                   tx.cost = tx.cost.pc["china"],
                                   tx.cost.partial = tx.cost.partial.pc["china"],
                                   tx.cost.mdr = tx.cost.mdr.pc["china"],
                                   pct.mdr= pct.mdr.pc["china"],
                                   tx.cost.partial.mdr = tx.cost.partial.mdr["china"],
                                   my.title = "",
                                   gdp=gdp.pc["china"],
                                   contours=contours.pc["china"][[1]])
                                        #save(china,file="Data/china_icers_4July.rda")
    } else {
        load(file=china.file)
    }

    if (is.null(india.file)){
        india <- makeICERPlotFixed(icer.min = icer.min,
                                   icer.max = icer.max,
                                   case.dt.dif = case.dt.dif[2],
                                   plot.params = fits[["india"]]$params,
                                   start.state = fits[["india"]]$state,
                                   tx.cost = tx.cost.pc["india"],
                                   tx.cost.partial = tx.cost.partial.pc["india"],
                                   tx.cost.mdr = tx.cost.mdr.pc["india"],
                                   pct.mdr= pct.mdr.pc["india"],
                                   tx.cost.partial.mdr = tx.cost.partial.mdr["india"],
                                   my.title = "",
                                   gdp=gdp.pc["india"],
                                   contours=contours.pc["india"][[1]])
                                        #save(india,file="Data/india_icers_4July.rda")
    } else {
        load(file=india.file)
    }


    col1 <- rgb(0,90,0,maxColorValue=255)
    col2 <- rgb(153,0,13,maxColorValue=255)
    if (pdf) {
        pdf(fn,width=6.8,height=6)
    } else {
        quartz("",width=6.8,height=6)
    }

    par(mfrow=c(2,3),mar=c(1,1.5,1,1.5),mgp=c(2,.2,0),tck=-.01,cex.axis=0.8,oma=c(2,2,2,4))
    plotTBIncMort(sa[[2]],legend=F,col1=col1,col2=col2)
    abline(v=1:10,col=AddAlpha("grey",.2))
    abline(h=seq(100,1000,by=100),col=AddAlpha("grey",.2))

    text(par("usr")[2]*.9,par("usr")[4]*0.95,"A",cex=1.5)


    plotTBIncMort(china[[2]],legend=F,col1=col1,col2=col2)
    abline(v=1:10,col=AddAlpha("grey",.2))
    abline(h=seq(0,1000,by=100),col=AddAlpha("grey",.2))

    text(par("usr")[2]*.9,par("usr")[4]*0.95,"B",cex=1.5)


    plotTBIncMort(india[[2]],legend=F,col1=col1,col2=col2)
    abline(v=1:10,col=AddAlpha("grey",.2))
    abline(h=seq(0,1000,by=100),col=AddAlpha("grey",.2))

    text(par("usr")[2]*.9,par("usr")[4]*0.95,"C",cex=1.5)


    null <- makeICERPlotFixed(icer.min = icer.min,
                              icer.max = icer.max,
                              case.dt.dif = case.dt.dif["sa"],
                              plot.params = fits[["sa"]]$params,
                              start.state = fits[["sa"]]$state,
                              tx.cost = tx.cost.pc["sa"],
                              tx.cost.partial = tx.cost.partial.pc["sa"],
                              tx.cost.mdr = tx.cost.mdr.pc["sa"],
                              pct.mdr= pct.mdr.pc["sa"],
                              tx.cost.partial.mdr = tx.cost.partial.mdr["sa"],
                              my.title = "",
                              gdp=gdp.pc["sa"],
                              contours=contours.pc["sa"][[1]],
                              ICERS = sa[[1]],
                              intcont.run = sa[[2]])

    text(par("usr")[2]*.9,par("usr")[4]*0.95,"D",cex=1.5)

    box()

    null <- makeICERPlotFixed(icer.min = icer.min,
                              icer.max = icer.max,
                              case.dt.dif = case.dt.dif["china"],
                              plot.params = fits[["china"]]$params,
                              start.state = fits[["china"]]$state,
                              tx.cost = tx.cost.pc["china"],
                              tx.cost.partial = tx.cost.partial.pc["china"],
                              tx.cost.mdr = tx.cost.mdr.pc["china"],
                              pct.mdr= pct.mdr.pc["china"],
                              tx.cost.partial.mdr = tx.cost.partial.mdr["china"],
                              my.title = "",
                              gdp=gdp.pc["china"],
                              contours=contours.pc["china"][[1]],
                              ICERS = china[[1]],
                              intcont.run = china[[2]])

    text(par("usr")[2]*.9,par("usr")[4]*0.95,"E",cex=1.5)

    box()

    null <- makeICERPlotFixed(icer.min = icer.min,
                              icer.max = icer.max,
                              case.dt.dif = case.dt.dif["india"],
                              plot.params = fits[["india"]]$params,
                              start.state = fits[["india"]]$state,
                              tx.cost = tx.cost.pc["india"],
                              tx.cost.partial = tx.cost.partial.pc["india"],
                              tx.cost.mdr = tx.cost.mdr.pc["india"],
                              pct.mdr= pct.mdr.pc["india"],
                              tx.cost.partial.mdr = tx.cost.partial.mdr["india"],
                              my.title = "",
                              gdp=gdp.pc["india"],
                              contours=contours.pc["india"][[1]],
                              ICERS = india[[1]],
                              intcont.run = india[[2]],leg=TRUE)

    text(par("usr")[2]*.9,par("usr")[4]*0.95,"F",cex=1.5)

    box()
    mtext("ICER ($)",outer=T,side=4,at=.47,las=2,cex=.7,adj=.1)

    mtext("South Africa",side=3,outer=T,at=.15)
    mtext("China",side=3,outer=T,at=.5)
    mtext("India",side=3,outer=T,at=.85)
    mtext("Year",side=1,outer=T,at=0.5,cex=.75,line=.7)
    mtext("count per 100,000",side=2,outer=T,at=0.75,cex=.75)
    mtext("cost per case detected, year 1 (USD)",side=2,outer=T,at=0.25,cex=.75)
    if (pdf) dev.off()
}


#################################
## Numbers mentioned in Paper####
#################################
red.sa <- fitIncreasedDetectionRate(target.detection.increase = case.dt.dif["sa"],
                                    duration = 1,
                                    params = fits[["sa"]]$params,
                                    starting.state = fits[["sa"]]$state,
                                    ep.sn.multiplier = 1,
                                    var.beta=FALSE)

red.india <- fitIncreasedDetectionRate(target.detection.increase = case.dt.dif["india"],
                                       duration = 1,
                                       params = fits[["india"]]$params,
                                       starting.state = fits[["india"]]$state,
                                       ep.sn.multiplier = 1,
                                       var.beta=FALSE)

red.china <- fitIncreasedDetectionRate(target.detection.increase = case.dt.dif["china"],
                                       duration = 1,
                                       params = fits[["china"]]$params,
                                       starting.state = fits[["china"]]$state,
                                       ep.sn.multiplier = 1,
                                       var.beta=FALSE)

## pct increase in theta
(pct.increase.sa <- (red.sa$par+fits[["sa"]]$params$theta.sp[1])/fits[["sa"]]$params$theta.sp[1] - 1)
(pct.increase.india <- (red.india$par+fits[["india"]]$params$theta.sp[1])/fits[["india"]]$params$theta.sp[1] - 1)
(pct.increase.china <- (red.china$par+fits[["china"]]$params$theta.sp[1])/fits[["china"]]$params$theta.sp[1] - 1)

## duration of active TB in each country
## pre duration of aTB in SA
1/(fits[["sa"]]$params$theta.sp[1])*12 + (1/fits[["sa"]]$params$rho.ps[1])*12 # number of months
## post
1/(red.sa$par+fits[["sa"]]$params$theta.sp[1])*12 + (1/fits[["sa"]]$params$rho.ps[1])*12 # number of months
## pre -India
1/(fits[["india"]]$params$theta.sp[1])*12+ (1/fits[["india"]]$params$rho.ps[1])*12 # number of months
## post - India
1/(red.india$par+fits[["india"]]$params$theta.sp[1])*12+ (1/fits[["india"]]$params$rho.ps[1])*12 # number of months
## Pre - China
1/(fits[["china"]]$params$theta.sp[1])*12+ (1/fits[["china"]]$params$rho.ps[1])*12 # number of months
## Post - China
1/(red.china$par+fits[["china"]]$params$theta.sp[1])*12+ (1/fits[["china"]]$params$rho.ps[1])*12 # number of months

## Now running two year ACF campaigns in each country
two.year.india <- runNYearACF("india",
                              case.dt.dif=case.dt.dif,
                              int.dur=2,
                              total.dur = 10,
                              fits=fits)

two.year.china <- runNYearACF("china",
                              case.dt.dif=case.dt.dif,
                              int.dur=2,
                              total.dur = 10,
                              fits=fits)

two.year.sa <- runNYearACF("sa",
                           case.dt.dif=case.dt.dif,
                           int.dur=2,
                           total.dur = 10,
                           fits=fits)

## ## what is the maximum reduction in incidence achieved and timeing
## GetPctReductionAndTime <- function(sim){
##     ci1 <- diff(sim[[1]][,"CIall"])*10 #intervention
##     ci2 <- diff(sim[[2]][,"CIall"])*10 #baseline
##     c("year"=which.max(1- ci1/ci2)/10,
##       "max.pct.diff"=max(1- ci1/ci2))
## }

## GetPctReductionAndTime(two.year.india)
## GetPctReductionAndTime(two.year.china)
## GetPctReductionAndTime(two.year.sa)

## what are the prevalence ratios over time?
prev.ratio.india <- rowSums(getPrevCols(two.year.india[[1]]))/rowSums(getPrevCols(two.year.india[[2]]))
prev.ratio.china <- rowSums(getPrevCols(two.year.china[[1]]))/rowSums(getPrevCols(two.year.china[[2]]))
prev.ratio.sa <- rowSums(getPrevCols(two.year.sa[[1]]))/rowSums(getPrevCols(two.year.sa[[2]]))

## Figure S? prevalance ratios after 2 year campaign
pdf("Figures/prev_ratio_2yrint.pdf")
par(mar=c(3,3,1,1),mgp=c(1.7,0.5,0))
plot(seq(0,10,length=101),prev.ratio.sa,type="l",col=3,ylab="Prevalance Ratio",xlab="Year")
lines(seq(0,10,length=101),prev.ratio.china,col=1)
lines(seq(0,10,length=101),prev.ratio.india,col=2)
abline(v=seq(0,10,by=.25),col=AddAlpha("grey",.2))
abline(h=seq(0.85,1,by=.025),col=AddAlpha("grey",.2))
legend("topright",c("China","India","South Africa"),col=1:3,lty=1,bty="n")
dev.off()

## cumulative deaths averted in population of 1,000,000
10*(sum(tail(two.year.india[[2]][,grep("Mtb",colnames(two.year.india[[2]]))],1)) -
    sum(tail(two.year.india[[1]][,grep("Mtb",colnames(two.year.india[[1]]))],1)))
10*(sum(tail(two.year.china[[2]][,grep("Mtb",colnames(two.year.china[[2]]))],1)) -
    sum(tail(two.year.china[[1]][,grep("Mtb",colnames(two.year.china[[1]]))],1)))
10*(sum(tail(two.year.sa[[2]][,grep("Mtb",colnames(two.year.sa[[2]]))],1)) -
    sum(tail(two.year.sa[[1]][,grep("Mtb",colnames(two.year.sa[[1]]))],1)))

## Uncertainty Ranges for cum deaths
## WARNING: need to load data from uncer file first (see ACF-uncer script)
deaths.averted.india <- c(sapply(runs.india.1,function(x)
                                 sum(tail(x[[2]][,grep("Mtb",colnames(x[[2]]))],1)) -
                                 sum(tail(x[[1]][,grep("Mtb",colnames(x[[1]]))],1)),simplify=T),
                          sapply(runs.india.2,function(x)
                                 sum(tail(x[[2]][,grep("Mtb",colnames(x[[2]]))],1)) -
                                 sum(tail(x[[1]][,grep("Mtb",colnames(x[[1]]))],1)),simplify=T),
                          sapply(runs.india.3,function(x)
                                 sum(tail(x[[2]][,grep("Mtb",colnames(x[[2]]))],1)) -
                                 sum(tail(x[[1]][,grep("Mtb",colnames(x[[1]]))],1)),simplify=T),
                          sapply(runs.india.4,function(x)
                                 sum(tail(x[[2]][,grep("Mtb",colnames(x[[2]]))],1)) -
                                 sum(tail(x[[1]][,grep("Mtb",colnames(x[[1]]))],1)),simplify=T))
10*quantile(deaths.averted.india,c(0.025,.975)) # in a population of a million

deaths.averted.china <- c(sapply(runs.china.1,function(x)
                                 sum(tail(x[[2]][,grep("Mtb",colnames(x[[2]]))],1)) -
                                 sum(tail(x[[1]][,grep("Mtb",colnames(x[[1]]))],1)),simplify=T),
                          sapply(runs.china.2,function(x)
                                 sum(tail(x[[2]][,grep("Mtb",colnames(x[[2]]))],1)) -
                                 sum(tail(x[[1]][,grep("Mtb",colnames(x[[1]]))],1)),simplify=T),
                          sapply(runs.china.3,function(x)
                                 sum(tail(x[[2]][,grep("Mtb",colnames(x[[2]]))],1)) -
                                 sum(tail(x[[1]][,grep("Mtb",colnames(x[[1]]))],1)),simplify=T),
                          sapply(runs.china.4,function(x)
                                 sum(tail(x[[2]][,grep("Mtb",colnames(x[[2]]))],1)) -
                                 sum(tail(x[[1]][,grep("Mtb",colnames(x[[1]]))],1)),simplify=T))

10*quantile(deaths.averted.china,c(0.025,.975)) # in a population of a million

deaths.averted.sa <- c(sapply(runs.sa.1,function(x)
                              sum(tail(x[[2]][,grep("Mtb",colnames(x[[2]]))],1)) -
                              sum(tail(x[[1]][,grep("Mtb",colnames(x[[1]]))],1)),
                              simplify=T),
                       sapply(runs.sa.2,function(x)
                              sum(tail(x[[2]][,grep("Mtb",colnames(x[[2]]))],1)) -
                              sum(tail(x[[1]][,grep("Mtb",colnames(x[[1]]))],1)),
                              simplify=T),
                       sapply(runs.sa.3,function(x)
                              sum(tail(x[[2]][,grep("Mtb",colnames(x[[2]]))],1)) -
                              sum(tail(x[[1]][,grep("Mtb",colnames(x[[1]]))],1)),
                              simplify=T),
                       sapply(runs.sa.4,function(x)
                              sum(tail(x[[2]][,grep("Mtb",colnames(x[[2]]))],1)) -
                              sum(tail(x[[1]][,grep("Mtb",colnames(x[[1]]))],1)),
                              simplify=T))

10*quantile(deaths.averted.sa,c(0.025,.975))

## pecentage of deaths occuring over the time period of the intervention
deaths.averted.india <- rowSums(two.year.india[[2]][,grep("Mtb",colnames(two.year.india[[2]]))]) - rowSums(two.year.india[[1]][,grep("Mtb",colnames(two.year.india[[1]]))])
deaths.averted.india[21]/deaths.averted.india[101]
deaths.averted.china <- rowSums(two.year.china[[2]][,grep("Mtb",colnames(two.year.china[[2]]))]) - rowSums(two.year.china[[1]][,grep("Mtb",colnames(two.year.china[[1]]))])
deaths.averted.china[21]/deaths.averted.china[101]
deaths.averted.sa <- rowSums(two.year.sa[[2]][,grep("Mtb",colnames(two.year.sa[[2]]))]) - rowSums(two.year.sa[[1]][,grep("Mtb",colnames(two.year.sa[[1]]))])
deaths.averted.sa[21]/deaths.averted.sa[101]


## What is the pct reduction in all-cause mortality during the intervention it self?
diff(rowSums(two.year.sa[[2]][,grep("Mtb",colnames(two.year.sa[[2]]))]))
min(diff(rowSums(two.year.sa[[1]][,grep("Mtb",colnames(two.year.sa[[1]]))]))[1:21])/max(diff(rowSums(two.year.sa[[2]][,grep("Mtb",colnames(two.year.sa[[2]]))])))

## how many cases averted and when
## what percentage are durign the intervention itself?
cases.averted <- two.year.india[[2]][,grep("CIall",colnames(two.year.india[[2]]))] - two.year.india[[1]][,grep("CIall",colnames(two.year.india[[1]]))]
cases.averted[21]/max(cases.averted)

## CE thresholds after 10-yrs
(cum.cases.averted <- tail(two.year.india[[2]][,"CIall"],1) - tail(two.year.india[[1]][,"CIall"],1))*10
10*(deaths.averted <- rowSums(two.year.india[[2]][,grep("Mtb",colnames(two.year.india[[2]]))]) - rowSums(two.year.india[[1]][,grep("Mtb",colnames(two.year.india[[1]]))]))

#deaths.averted <- diff(rowSums(two.year.india[[2]][,grep("Mtb",colnames(two.year.india[[2]]))])) - diff(rowSums(two.year.india[[1]][,grep("Mtb",colnames(two.year.india[[1]]))]))

## ho much per daly after 2 years
calcICERFixedCosts(out=two.year.india,
                   eval.times = 1:(10*2+1),
                   dtx.cost=case.dt.dif["india"]*500,
                   tx.suc=c(1),
                   tx.cost = tx.cost.pc["india"],
                   tx.cost.partial = tx.cost.partial.pc["india"],
                   tx.cost.mdr = tx.cost.mdr.pc["india"],
                   pct.mdr= pct.mdr.pc["india"],
                   tx.cost.partial.mdr = tx.cost.partial.mdr["india"],
                   params=fits[["india"]]$params)

## and after 10 years
calcICERFixedCosts(out=two.year.india,
                   eval.times = 1:(10*10+1),
                   dtx.cost=case.dt.dif["india"]*500,
                   tx.suc=c(1),
                   tx.cost = tx.cost.pc["india"],
                   tx.cost.partial = tx.cost.partial.pc["india"],
                   tx.cost.mdr = tx.cost.mdr.pc["india"],
                   pct.mdr= pct.mdr.pc["india"],
                   tx.cost.partial.mdr = tx.cost.partial.mdr["india"],
                   params=fits[["india"]]$params)



par(mfrow=c(3,1))
plot(seq(0,10,length=101),deaths.averted,type="l",main="Deaths Averted in India",xlab="year")

# pct of deaths averted by time
plot(seq(0,10,length=101),(rowSums(two.year.india[[2]][,grep("Mtb",colnames(two.year.india[[2]]))]) - rowSums(two.year.india[[1]][,grep("Mtb",colnames(two.year.india[[1]]))]))/tail(deaths.averted,1),ylab="pct deaths averted by")

plot(seq(0,10,length=101),rowSums(two.year.india[[2]][,grep("Mtb",colnames(two.year.india[[2]]))]),type="l",col=2)
lines(seq(0,10,length=101),rowSums(two.year.india[[1]][,grep("Mtb",colnames(two.year.india[[1]]))]),col=3)


## what percent would occur during the intervention itself?
deaths.during.intervention.baseline <- rowSums(two.year.india[[2]][,grep("Mtb",colnames(two.year.india[[2]]))])[21]
deaths.during.intervention.int <- rowSums(two.year.india[[1]][,grep("Mtb",colnames(two.year.india[[1]]))])[21]
(rowSums(two.year.india[[2]][,grep("Mtb",colnames(two.year.india[[2]]))])[21] - rowSums(two.year.india[[1]][,grep("Mtb",colnames(two.year.india[[1]]))])[21]) / (rowSums(two.year.india[[2]][,grep("Mtb",colnames(two.year.india[[2]]))])[101] - rowSums(two.year.india[[1]][,grep("Mtb",colnames(two.year.india[[1]]))])[101])

deaths.during.10yrs.baseline <- rowSums(two.year.india[[2]][,grep("Mtb",colnames(two.year.india[[2]]))])[101]
deaths.during.10yrs.int <- rowSums(two.year.india[[1]][,grep("Mtb",colnames(two.year.india[[1]]))])[101]

(deaths.during.intervention.baseline - deaths.during.intervention.int)/(deaths.during.10yrs.baseline- deaths.during.10yrs.int)

sum(tail(two.year.india[[2]][,grep("Mtb",colnames(two.year.india[[2]]))],1)) -
    sum(tail(two.year.india[[1]][,grep("Mtb",colnames(two.year.india[[1]]))],1))

sum(tail(two.year.china[[2]][,grep("Mtb",colnames(two.year.china[[2]]))],1)) -
    sum(tail(two.year.china[[1]][,grep("Mtb",colnames(two.year.china[[1]]))],1))

(rowSums(two.year.china[[2]][,grep("Mtb",colnames(two.year.china[[2]]))])[21] - rowSums(two.year.china[[1]][,grep("Mtb",colnames(two.year.china[[1]]))])[21]) / (rowSums(two.year.china[[2]][,grep("Mtb",colnames(two.year.china[[2]]))])[101] - rowSums(two.year.china[[1]][,grep("Mtb",colnames(two.year.china[[1]]))])[101])



## DALYS averted
calcICERFixedCosts(out=two.year.india,
                   eval.times = 1:(10*10+1),
                   dtx.cost=case.dt.dif["india"]*1200,
                   tx.suc=c(1),
                   tx.cost = tx.cost.pc["india"],
                   tx.cost.partial = tx.cost.partial.pc["india"],
                   tx.cost.mdr = tx.cost.mdr.pc["india"],
                   pct.mdr= pct.mdr.pc["india"],
                   tx.cost.partial.mdr = tx.cost.partial.mdr["india"],
                   params=fits[["india"]]$params)


#####################################
#### 10 year sustained ##############
#####################################

## India
theta.reduction.india <- red.india$par
run.india  <- runIntCont(fits[["india"]]$state,
                         fits[["india"]]$params,10,
                         int.theta.sp= theta.reduction.india,
                         int.theta.sn = theta.reduction.india*1,
                         int.theta.ep = theta.reduction.india*1)
out.india <- run.india
## pct recution in cum incidence
cum.inc.diff <- tail(out.india[[2]][,"CIall"],1) - tail(out.india[[1]][,"CIall"],1)
(pct.red.cum.inc <- 1 -
 (tail(out.india[[1]][,"CIall"],1) - out.india[[1]][1,"CIall"]) /
 (tail(out.india[[2]][,"CIall"],1) -  out.india[[2]][1,"CIall"])) #15.28%

## in final year
inc.diff.final.yr <- (tail(out.india[[2]][,"CIall"],1) -  tail(out.india[[2]][,"CIall"],11)[1]) -
    (tail(out.india[[1]][,"CIall"],1) -  tail(out.india[[1]][,"CIall"],11)[1])
(inc.pct.diff.final.yr <- 1 - (tail(out.india[[1]][,"CIall"],1) -  tail(out.india[[1]][,"CIall"],11)[1])/
 (tail(out.india[[2]][,"CIall"],1) -  tail(out.india[[2]][,"CIall"],11)[1]))


##pct reduction in cumulative 10-yr mortality
mtb.cols <- grep("Mtb",colnames(out.india[[1]]))
(pct.red.cum.mort <- 1 - (sum(tail(out.india[[1]][,mtb.cols],1)) /sum(tail(out.india[[2]][,mtb.cols],1)))) #41.68% reduction
total <- sum(tail(out.india[[2]][,mtb.cols],1)) - sum(tail(out.india[[1]][,mtb.cols],1))
(rowSums(out.india[[2]][,mtb.cols]) - rowSums(out.india[[1]][,mtb.cols]))[21]/total
#13.38%

## incidence and mortality rates in final year
1 - diff(tail(out.india[[1]][,"CIall"],2) ) / diff(tail(out.india[[2]][,"CIall"],2) )
1 - (diff(rowSums(tail(out.india[[1]][,mtb.cols],2)))/diff(rowSums(tail(out.india[[2]][,mtb.cols],2))))


## China
theta.reduction.china <- red.china$par
run.china  <- runIntCont(fits[["china"]]$state,
                         fits[["china"]]$params,10,
                         int.theta.sp= theta.reduction.china,
                         int.theta.sn = theta.reduction.china*1,
                         int.theta.ep = theta.reduction.china*1)
out.china <- run.china
## pct reduction in cum incidence
cum.inc.diff <- tail(out.china[[2]][,"CIall"],1) - tail(out.china[[1]][,"CIall"],1)
(pct.red.cum.inc <- 1 - tail(out.china[[1]][,"CIall"],1)/tail(out.china[[2]][,"CIall"],1) ) #19%

## in final year
inc.diff.final.yr <- (tail(out.china[[2]][,"CIall"],1) -  tail(out.china[[2]][,"CIall"],11)[1]) -
    (tail(out.china[[1]][,"CIall"],1) -  tail(out.china[[1]][,"CIall"],11)[1])
(inc.pct.diff.final.yr <-1 - (tail(out.china[[1]][,"CIall"],1) -  tail(out.china[[1]][,"CIall"],11)[1]) /
 (tail(out.china[[2]][,"CIall"],1) -  tail(out.china[[2]][,"CIall"],11)[1])) #32% reduction

##pct reduction in cumulative 10-yr mortality
mtb.cols <- grep("Mtb",colnames(out.china[[1]]))
(pct.red.cum.mort <- 1 - (sum(tail(out.china[[1]][,mtb.cols],1)) /sum(tail(out.china[[2]][,mtb.cols],1)))) #41.68% reduction

## cum.inc.dif <- sum(inc.diffs)
## final.year.inc.dif <- sum(tail(inc.diffs,10))

## incidence and mortality rates in final year
1 - diff(tail(out.china[[1]][,"CIall"],2) ) / diff(tail(out.china[[2]][,"CIall"],2) )
1 - (diff(rowSums(tail(out.china[[1]][,mtb.cols],2)))/diff(rowSums(tail(out.china[[2]][,mtb.cols],2))))


## South Africa
theta.reduction.sa <- red.sa$par
run.sa  <- runIntCont(fits[["sa"]]$state,
                      fits[["sa"]]$params,10,
                      int.theta.sp= theta.reduction.sa,
                      int.theta.sn = theta.reduction.sa*1,
                      int.theta.ep = theta.reduction.sa*1)
out.sa <- run.sa
    ## pct recution in cum incidence
cum.inc.diff <- tail(out.sa[[2]][,"CIall"],1) - tail(out.sa[[1]][,"CIall"],1)
(pct.red.cum.inc <- 1 -
 (tail(out.sa[[1]][,"CIall"],1) - out.sa[[1]][1,"CIall"]) /
 (tail(out.sa[[2]][,"CIall"],1) -  out.sa[[2]][1,"CIall"])) #3.5%

    ## in final year
inc.diff.final.yr <- (tail(out.sa[[2]][,"CIall"],1) -  tail(out.sa[[2]][,"CIall"],11)[1]) -
(tail(out.sa[[1]][,"CIall"],1) -  tail(out.sa[[1]][,"CIall"],11)[1])
(inc.pct.diff.final.yr <- 1 - (tail(out.sa[[1]][,"CIall"],1) -  tail(out.sa[[1]][,"CIall"],11)[1]) /
 (tail(out.sa[[2]][,"CIall"],1) -  tail(out.sa[[2]][,"CIall"],11)[1])) #5% reduction

                                        #pct reduction in cumulative 10-yr mortality
mtb.cols <- grep("Mtb",colnames(out.sa[[1]]))
(pct.red.cum.mort <- 1 - ((sum(tail(out.sa[[1]][,mtb.cols],1)) - sum(out.sa[[1]][1,mtb.cols]))  /
                          (sum(tail(out.sa[[2]][,mtb.cols],1)) - sum(out.sa[[2]][1,mtb.cols])))) #41.68% reduction

1 - sum(tail(out.sa[[1]][,mtb.cols],1))/sum(tail(out.sa[[2]][,mtb.cols],1))

    ## incidence and mortality rates in final year
1 - diff(tail(out.sa[[1]][,"CIall"],2) ) / diff(tail(out.sa[[2]][,"CIall"],2) )
1 - (diff(rowSums(tail(out.sa[[1]][,mtb.cols],2)))/diff(rowSums(tail(out.sa[[2]][,mtb.cols],2))))


# prevalance ratio
prev.ratio.india10 <- rowSums(getPrevCols(out.india[[1]]))/rowSums(getPrevCols(out.india[[2]]))
prev.ratio.china10 <- rowSums(getPrevCols(out.china[[1]]))/rowSums(getPrevCols(out.china[[2]]))
prev.ratio.sa10 <- rowSums(getPrevCols(out.sa[[1]]))/rowSums(getPrevCols(out.sa[[2]]))

pdf("Figures/prevratio_10yrintervention.pdf")
par(mgp = c(1.5,.5,0),mar=c(3,3,2,.5))
plot(seq(0,10,length=101),
     prev.ratio.sa10,
     type="l",col=3,
     ylim=c(0.65,1),
     ylab="prevalence ratio",xlab="year")
lines(seq(0,10,length=101),prev.ratio.india10,col=2)
lines(seq(0,10,length=101),prev.ratio.china10,col=1)
                                        #abline(v=2,col="red",lty=2)
abline(h=seq(0,1,by=0.1),col="grey")
abline(v=seq(0,10,by=1),col="grey")
legend("topright",c("China","India","South Africa"),bty="n",col=1:3,lty=1)
dev.off()
## ICERS

getICER(horiz=2.5,
        cost=5000*case.dt.dif["sa"],
        params= fits[["sa"]]$params,
        out=two.year.sa,
        tx.cost=tx.cost.pc["sa"],
        tx.cost.partial=tx.cost.partial.pc["sa"],
        tx.cost.mdr=tx.cost.mdr.pc["sa"],
        tx.cost.partial.mdr=tx.cost.partial.mdr["sa"],
        pct.mdr=pct.mdr.pc["sa"],
        fixed=TRUE)

getICER(horiz=6.2,
        cost=5000*case.dt.dif["china"],
        params= fits[["china"]]$params,
        out=china[[2]],
        tx.cost=tx.cost.pc["china"],
        tx.cost.partial=tx.cost.partial.pc["china"],
        tx.cost.mdr=tx.cost.mdr.pc["china"],
        tx.cost.partial.mdr=tx.cost.partial.mdr["china"],
        pct.mdr=pct.mdr.pc["china"],
        fixed=TRUE)

getICER(horiz=9.9,
        cost=5000*case.dt.dif["india"],
        params= fits[["india"]]$params,
        out=india[[2]],
        tx.cost=tx.cost.pc["india"],
        tx.cost.partial=tx.cost.partial.pc["india"],
        tx.cost.mdr=tx.cost.mdr.pc["india"],
        tx.cost.partial.mdr=tx.cost.partial.mdr["india"],
        pct.mdr=pct.mdr.pc["india"],
        fixed=TRUE)


getICER(horiz=2,
        cost=1500*case.dt.dif["sa"],
        params= fits[["sa"]]$params,
        out=sa[[2]],
        tx.cost=tx.cost.pc["sa"],
        tx.cost.partial=tx.cost.partial.pc["sa"],
        tx.cost.mdr=tx.cost.mdr.pc["sa"],
        tx.cost.partial.mdr=tx.cost.partial.mdr["sa"],
        pct.mdr=pct.mdr.pc["sa"],
        fixed=TRUE)

getICER(horiz=10,
        cost=1500*case.dt.dif["sa"],
        params= fits[["sa"]]$params,
        out=sa[[2]],
        tx.cost=tx.cost.pc["sa"],
        tx.cost.partial=tx.cost.partial.pc["sa"],
        tx.cost.mdr=tx.cost.mdr.pc["sa"],
        tx.cost.partial.mdr=tx.cost.partial.mdr["sa"],
        pct.mdr=pct.mdr.pc["sa"],
        fixed=TRUE)

getICER(horiz=10,
        cost=1500*case.dt.dif["china"],
        params= fits[["china"]]$params,
        out=china[[2]],
        tx.cost=tx.cost.pc["china"],
        tx.cost.partial=tx.cost.partial.pc["china"],
        tx.cost.mdr=tx.cost.mdr.pc["china"],
        tx.cost.partial.mdr=tx.cost.partial.mdr["china"],
        pct.mdr=pct.mdr.pc["china"],
        fixed=TRUE)

getICER(horiz=10,
        cost=1500*case.dt.dif["india"],
        params= fits[["india"]]$params,
        out=india[[2]],
        tx.cost=tx.cost.pc["india"],
        tx.cost.partial=tx.cost.partial.pc["india"],
        tx.cost.mdr=tx.cost.mdr.pc["india"],
        tx.cost.partial.mdr=tx.cost.partial.mdr["india"],
        pct.mdr=pct.mdr.pc["india"],
        fixed=TRUE)


## Let's look at how the ICERs change with different intervention sizes
##' @param pct.dif pct difference from baseline
##' @param country
##' @param cost.per.case
##' @return returns matrix with column of pct differnces and ICERS for each one (with a single analystic horizon)
##' @author Andrew Azman
int.size.sens <- function(pct.dif=seq(.05,.5,length=10),
                          country="india",
                          cost.per.case=1500,
                          case.dt.dif,
                          fits){


    sa.trial <- runTBHIVMod(fit.sa.2011$params,fit.sa.2011$state,1,var.beta=F)
    india.trial <- runTBHIVMod(fit.india.2011$params,fit.india.2011$state,1,var.beta=F)
    china.trial <- runTBHIVMod(fit.china.2011$params,fit.china.2011$state,1,var.beta=F)

    runs <- sapply(seq_along(pct.dif),function(p) {

        case.difs <- c("china"=
                       round(sum(tail(
                           china.trial[,grep("N.", colnames(india.trial))],1))
                             *pct.dif[p],0),
                       "india"=round(sum(tail(
                           india.trial[,grep("N.", colnames(india.trial))],1))
                           *pct.dif[p],0),
                       "sa"=round(sum(tail(
                           sa.trial[,grep("N.", colnames(india.trial))],1))
                           *pct.dif[p],0))
        print(case.difs)
        run <- runNYearACF(country,
                           pct.incidence = 0,
                           case.dt.dif=case.difs,
                           int.dur = 2,
                           total.dur = 10,
                           fits=fits)

        calcICERFixedCosts(out=run,
                           eval.times = 1:(10*10+1),
                           dtx.cost=round(case.difs[country])*cost.per.case,
                           tx.suc=c(1),
                           tx.cost = tx.cost.pc[country],
                           tx.cost.partial = tx.cost.partial.pc[country],
                           tx.cost.mdr = tx.cost.mdr.pc[country],
                           pct.mdr= pct.mdr.pc[country],
                           tx.cost.partial.mdr = tx.cost.partial.mdr[country],
                           params=fits[[country]]$params)[2]
    })

    cbind(pct.dif,runs)
}


## Here we make the figure for the supplement

india.mag.sense <- int.size.sens(country="india",fits=fits)
sa.mag.sense <- int.size.sens(country="sa",fits=fits)
china.mag.sense <- int.size.sens(country="china",fits=fits)

## get cost per DALY at 1500/case
countries <- c("india","china","sa")
base.case <- sapply(countries,function(country){
    run <- runNYearACF(country,pct.incidence = 0.15,
                       case.dt.dif=case.dt.dif,
                       int.dur = 2,total.dur = 10,fits=fits)
    calcICERFixedCosts(out=run,
                       eval.times = 1:(10*10+1),
                       dtx.cost=case.dt.dif[country]*1500,
                       tx.suc=c(1),
                       tx.cost = tx.cost.pc[country],
                       tx.cost.partial = tx.cost.partial.pc[country],
                       tx.cost.mdr = tx.cost.mdr.pc[country],
                       pct.mdr= pct.mdr.pc[country],
                       tx.cost.partial.mdr = tx.cost.partial.mdr[country],
                       params=fits[[country]]$params)[2]
})

pdf("Figures/change_with_intervention_mag.pdf")
plot(india.mag.sense[,1],1-india.mag.sense[,2]/base.case[1],col=2,xlab="Percent Increase in Number of Cases Detected in First Year",ylab="Percent Change in Cost per DALY Averted",ylim=c(-.2,.2),pch=16)
points(jitter(china.mag.sense[,1]),1-china.mag.sense[,2]/base.case[2],col=1,pch=16)
points(jitter(sa.mag.sense[,1]),1-sa.mag.sense[,2]/base.case[3],col=3,pch=16)
abline(v=.25,col="grey",lty=2)
legend("topright",c("India","China","South Africa"),col=c(2,1,3),pch=16,bty="n")
dev.off()



## explore impact of time horizon on Cost Effectiveness of Intervention
make.time.horiz.ce <- function(pct.first.yr=.25){
    gdp.pc <- c("sa"=8090,"india"=1528,"china"=5439)
    case.dt.df <- getIncreasedCasesDetected(case.det.based=TRUE,pct.first.yr=0.25)

    country <- c("china","india","sa")

    per.person.dx.cost <- c(50,500,1000,5000,10000)
    horizons <- seq(1,10,by=0.1)

    out <- array(dim=c(length(per.person.dx.cost),length(horizons),3))

    for (c in seq_along(country)){
        ## fit the new case detection rate
        fit.tmp <- fitIncreasedDetectionRate(target.detection.increase = case.dt.df[c],
                                             duration = 1,
                                             params = fits[[country[c]]]$params,
                                             starting.state = fits[[country[c]]]$state,
                                             ep.sn.multiplier = 1,
                                             var.beta=FALSE)

        theta.reduction <- fit.tmp$par

        ## run intervention for 2 years but run model for 10
        year2_int_only<- runIntCont(ss=fits[[country[c]]]$state,
                                    params=fits[[country[c]]]$params,
                                    time=10,
                                    int.theta.sp=theta.reduction,
                                    int.theta.sn=theta.reduction,
                                    int.theta.ep=theta.reduction,
                                    intervention.duration = 2)


        for (t in seq_along(per.person.dx.cost)){
            for (h in seq_along(horizons)){
                out[t,h,c] <- calcICERFixedCosts(out=year2_int_only,
                                                 eval.times = 1:((horizons[h]*10)+1),
                                                 dtx.cost=case.dt.df[c]*per.person.dx.cost[t],
                                                 tx.suc=c(1),
                                                 tx.cost = tx.cost.pc[country[c]],
                                                 tx.cost.partial = tx.cost.partial.pc[country[c]],
                                                 tx.cost.mdr = tx.cost.mdr.pc[country[c]],
                                                 pct.mdr= pct.mdr.pc[country[c]],
                                                 tx.cost.partial.mdr = tx.cost.partial.mdr[country[c]],
                                                 params=fits[[country[c]]]$params)[2]
            }
        }
    }

    pdf("Figures/short_icers_horizon_3panel_new.pdf",6,5)
    par(mfrow=c(1,3),mar=c(.1,1,.1,1),oma=c(4,4,4,2),mgp=c(2.5,1,0))

    plot(-1000,-1000,xlim=range(horizons),ylim=c(-.1,10),xlab="",ylab="",axes=FALSE)
    for (i in 1:5) lines(horizons,out[i,,1]/gdp.pc["china"],col=1,lty=i,lwd=1)
    axis(1,at=1:10)
    axis(2,at=0:10)
    abline(h=1,col="grey")
    abline(h=3,col="grey")
    text(8.4,3.1,"Cost Effective",cex=.8)
    text(7.7,1.1,"Highly Cost Effective",cex=.8)

    abline(v=1:10,col=AddAlpha("grey",0.2))
    abline(h=0:10,col=AddAlpha("grey",0.2))


    plot(-1000,-1000,xlim=range(horizons),ylim=c(-.1,10),xlab="",ylab="",axes=FALSE)
    for (i in 1:5) lines(horizons,out[i,,2]/gdp.pc["india"],col=2,lty=i,lwd=1)
    axis(1,at=1:10)
    axis(2,at=0:10,labels=NA)
    abline(h=1,col="grey")
    abline(h=3,col="grey")
    text(8.4,3.1,"Cost Effective",cex=.8)
    text(7.7,1.1,"Highly Cost Effective",cex=.8)

    abline(v=1:10,col=AddAlpha("grey",0.2))
    abline(h=0:10,col=AddAlpha("grey",0.2))


    plot(-1000,-1000,xlim=range(horizons),ylim=c(-.1,10),xlab="Time Horizon (years)",ylab="ICER (USD)",axes=FALSE)
    for (i in 1:5) lines(horizons,out[i,,3]/gdp.pc["sa"],col=3,lty=i,lwd=1)
    axis(1,at=1:10)
    axis(2,at=0:10,labels=NA)

    abline(h=1,col="grey")
    abline(h=3,col="grey")
    text(8.4,3.1,"Cost Effective",cex=.8)
    text(7.7,1.1,"Highly Cost Effective",cex=.8)

    abline(v=1:10,col=AddAlpha("grey",0.2))
    abline(h=0:10,col=AddAlpha("grey",0.2))


    legend("topright",legend=paste0("$",rev(per.person.dx.cost)," per Case"),col="black",lty=5:1,bty="n")

    mtext("China",side=3,outer=T,adj=.15)
    mtext("India",side=3,outer=T,adj=.5)
    mtext("South Africa",side=3,outer=T,adj=.9)
    mtext("Time Horizon (years)",side=1,outer=T,line=2.5)
    mtext("ICER / GDP per capita",side=2,outer=T,line=2)

    dev.off()
}


## % mortality
deaths.sa.2010 <- 543856
deaths.sa.2009 <- 579711
deaths.in.sa.model.09.10 <- 10*(sum(two.year.sa[[2]][,grep("Mtb",colnames(two.year.sa[[2]]))][21,]) - sum(two.year.sa[[1]][,grep("Mtb",colnames(two.year.sa[[1]]))][21,]))
(change.in.all.cause.mortality <- deaths.in.sa.model.09.10/(deaths.sa.2010+deaths.sa.2009)*100)
