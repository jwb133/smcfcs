pkgname <- "Publish"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('Publish')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CiTable")
### * CiTable

flush(stderr()); flush(stdout())

### Name: CiTable
### Title: CiTable data
### Aliases: CiTable
### Keywords: datasets

### ** Examples


data(CiTable)
labellist <- split(CiTable[,c("Dose","Mean","SD","n")],CiTable[,"Drug"])
labellist
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")], labels=labellist)





cleanEx()
nameEx("Diabetes")
### * Diabetes

flush(stderr()); flush(stdout())

### Name: Diabetes
### Title: Diabetes data of Dr John Schorling
### Aliases: Diabetes
### Keywords: datasets

### ** Examples


data(Diabetes)




cleanEx()
nameEx("SpaceT")
### * SpaceT

flush(stderr()); flush(stdout())

### Name: SpaceT
### Title: A study was made of all 26 astronauts on the first eight space
###   shuttle flights (Bungo et.al., 1985). On a voluntary basis 17
###   astronauts consumed large quantities of salt and fluid prior to
###   landing as a countermeasure to space deconditioning, while nine did
###   not.
### Aliases: SpaceT

### ** Examples

data(SpaceT)



cleanEx()
nameEx("Units")
### * Units

flush(stderr()); flush(stdout())

### Name: Units
### Title: Add units to data set
### Aliases: Units

### ** Examples

data(Diabetes)
Diabetes <- Units(Diabetes,list(BMI="kg/m^2"))
Units(Diabetes)
Diabetes <- Units(Diabetes,list(bp.1s="mm Hg",bp.2s="mm Hg"))
Units(Diabetes)



cleanEx()
nameEx("coxphSeries")
### * coxphSeries

flush(stderr()); flush(stdout())

### Name: coxphSeries
### Title: Run a series of Cox regression models
### Aliases: coxphSeries

### ** Examples

library(survival)
data(pbc)
## collect hazard ratios from three univariate Cox regression analyses
pbc$edema <- factor(pbc$edema,levels=c("0","0.5","1"),labels=c("0","0.5","1"))
uni.hr <- coxphSeries(Surv(time,status==2)~1,vars=c("edema","bili","protime"),data=pbc)
uni.hr

## control the logistic regression analyses for age and gender
## but collect only information on the variables in `vars'.
controlled.hr <- coxphSeries(Surv(time,status==2)~age+sex,vars=c("edema","bili","protime"),data=pbc)
controlled.hr




cleanEx()
nameEx("followupTable")
### * followupTable

flush(stderr()); flush(stdout())

### Name: followupTable
### Title: Summary tables for a given followup time point.
### Aliases: followupTable

### ** Examples

library(survival)
data(pbc)
pbc$edema <- factor(pbc$edema,levels=c("0","0.5","1"),labels=c("0","0.5","1"))
pbc$sex <- factor(pbc$sex,levels=c("m","f"),labels=c("m","f"))
followupTable(Hist(time,status)~age+edema+sex,data=pbc,followup.time=1000)




cleanEx()
nameEx("formatCI")
### * formatCI

flush(stderr()); flush(stdout())

### Name: formatCI
### Title: Formatting confidence intervals
### Aliases: formatCI

### ** Examples


x=ci.mean(rnorm(10))
formatCI(lower=x[3],upper=x[4])
formatCI(lower=c(0.001,-2.8413),upper=c(1,3.0008884))
# change format
formatCI(lower=c(0.001,-2.8413),upper=c(1,3.0008884),format="(l, u)")
# show x
formatCI(x=x$mean,lower=x$lower,upper=x$upper,format="(l, u)",show.x=TRUE)

# if the first lower limit is shorter than the second (missing negative sign),
# then, option trim will make a difference:
formatCI(lower=c(0.001,-2.8413),upper=c(1,3.0008884),format="l--u",trim=FALSE)
formatCI(lower=c(0.001,-2.8413),upper=c(1,3.0008884),format="l--u",trim=TRUE)

# change of handler function
l <- c(-0.0890139,0.0084736,144.898333,0.000000001)
u <- c(0.03911392,0.3784706,3338944.8821221,0.00001)
cbind(format=formatCI(lower=l,upper=u,format="[l;u)",digits=2,nsmall=2,handler="format"),
      prettyNum=formatCI(lower=l,upper=u,format="[l;u)",digits=2,nsmall=2,handler="prettyNum"),
      sprintf=formatCI(lower=l,upper=u,format="[l;u)",digits=2,nsmall=2,handler="sprintf"))




cleanEx()
nameEx("glmSeries")
### * glmSeries

flush(stderr()); flush(stdout())

### Name: glmSeries
### Title: Run a series of generalized linear regression analyses
### Aliases: glmSeries

### ** Examples


data(Diabetes)
Diabetes$hyper1 <- factor(1*(Diabetes$bp.1s>140))
## collect odds ratios from three univariate logistic regression analyses
uni.odds <- glmSeries(hyper1~1,vars=c("chol","hdl","location"),data=Diabetes,family=binomial)
uni.odds
## control the logistic regression analyses for age and gender
## but collect only information on the variables in `vars'.
controlled.odds <- glmSeries(hyper1~age+gender,
                             vars=c("chol","hdl","location"),
                             data=Diabetes, family=binomial)
controlled.odds



cleanEx()
nameEx("labelUnits")
### * labelUnits

flush(stderr()); flush(stdout())

### Name: labelUnits
### Title: labelUnits
### Aliases: labelUnits

### ** Examples


data(Diabetes)
tab <- summary(univariateTable(gender~AgeGroups+chol+waist,data=Diabetes))
publish(tab)
ltab <- labelUnits(tab,"chol"="Cholesterol (mg/dL)","<40"="younger than 40")
publish(ltab)

## pass labels immediately to utable
utable(gender~AgeGroups+chol+waist,data=Diabetes,
      "chol"="Cholesterol (mg/dL)","<40"="younger than 40")

## sometimes useful to state explicitly which variables value
## should be re-labelled
utable(gender~AgeGroups+chol+waist,data=Diabetes,
      "chol"="Cholesterol (mg/dL)","AgeGroups.<40"="younger than 40")



cleanEx()
nameEx("lazyFactorCoding")
### * lazyFactorCoding

flush(stderr()); flush(stdout())

### Name: lazyFactorCoding
### Title: Efficient coding of factor levels
### Aliases: lazyFactorCoding

### ** Examples

data(Diabetes)
lazyFactorCoding(Diabetes)




cleanEx()
nameEx("parseInteractionTerms")
### * parseInteractionTerms

flush(stderr()); flush(stdout())

### Name: parseInteractionTerms
### Title: Parse interaction terms
### Aliases: parseInteractionTerms

### ** Examples


tt <- terms(formula(SBP~age+sex*BMI))
xlev <- list(sex=c("male","female"),BMI=c("normal","overweight","obese"))
parseInteractionTerms(terms=tt,xlevels=xlev)
parseInteractionTerms(terms=tt,xlevels=xlev,format.factor="var level")
parseInteractionTerms(terms=tt,xlevels=xlev,format.contrast="var(level:ref)")

tt2 <- terms(formula(SBP~age*factor(sex)+BMI))
xlev2 <- list("factor(sex)"=c("male","female"))
parseInteractionTerms(terms=tt2,xlevels=xlev2)
parseInteractionTerms(terms=tt2,xlevels=xlev2,units=list(age="yrs"))


data(Diabetes)
fit <- glm(bp.2s~age*factor(gender)+BMI,data=Diabetes)
parseInteractionTerms(terms=terms(fit$formula),xlevels=fit$xlevels,
                      format.scale="var -- level:ref",units=list("age"='years'))
parseInteractionTerms(terms=terms(fit$formula),xlevels=fit$xlevels,
                      format.scale.unit="var -- level:ref",units=list("age"='years'))
it <- parseInteractionTerms(terms=terms(fit$formula),xlevels=fit$xlevels)
ivars <- unlist(lapply(it,function(x)attr(x,"variables")))
lava::estimate(fit,function(p)lapply(unlist(it),eval,envir=sys.parent(-1)))





cleanEx()
nameEx("plot.ci")
### * plot.ci

flush(stderr()); flush(stdout())

### Name: plot.ci
### Title: Plot confidence intervals
### Aliases: plot.ci

### ** Examples


data(Diabetes)
x=ci.mean(bp.2s~AgeGroups,data=Diabetes)
plot(x,title.labels="Age groups",xratio=c(0.4,0.3))
x=ci.mean(bp.2s/500~AgeGroups+gender,data=Diabetes)
plot(x,xratio=c(0.4,0.2))
plot(x,xratio=c(0.4,0.2),
     labels=split(x$labels[,"AgeGroups"],x$labels[,"gender"]),
     title.labels="Age groups")
## Not run: 
##D plot(x, leftmargin=0, rightmargin=0)
##D plotConfidence(x, leftmargin=0, rightmargin=0)
##D 
##D data(CiTable)
##D with(CiTable,plotConfidence(x=list(HazardRatio),
##D                                lower=lower,
##D                                upper=upper,
##D                                labels=CiTable[,2:6],
##D                                factor.reference.pos=c(1,10,19),
##D                                format="(u-l)",
##D                                points.col="blue",
##D                                digits=2))
##D 
##D with(CiTable,Publish::plot.ci(x=list(HazardRatio),
##D                                lower=lower,
##D                                upper=upper,
##D                                labels=CiTable[,2:6],
##D                                factor.reference.pos=c(1,10,19),
##D                                format="(u-l)",
##D                                points.col="blue",
##D                                digits=2,
##D                                leftmargin=-2,
##D                                title.labels.cex=1.1,
##D                                labels.cex=0.8,values.cex=0.8))
## End(Not run)



cleanEx()
nameEx("plot.regressionTable")
### * plot.regressionTable

flush(stderr()); flush(stdout())

### Name: plot.regressionTable
### Title: Plotting regression coefficients with confidence limits
### Aliases: plot.regressionTable

### ** Examples

## linear regression
data(Diabetes)
f <- glm(bp.1s~AgeGroups+chol+gender+location,data=Diabetes)
rtf <- regressionTable(f,factor.reference = "inline")
plot(rtf,cex=1.3)

## logistic regression
data(Diabetes)
f <- glm(I(BMI>25)~bp.1s+AgeGroups+chol+gender+location,data=Diabetes,family=binomial)
rtf <- regressionTable(f,factor.reference = "inline")
plot(rtf,cex=1.3)

## Poisson regression
data(trace)
fit <- glm(dead ~ smoking+ sex+ age+Time+offset(log(ObsTime)), family = poisson,data=trace)
rtab <- regressionTable(fit,factor.reference = "inline")
plot(rtab,xlim=c(0.85,1.15),cex=1.8,xaxis.cex=1.5)

## Cox regression
library(survival)
data(pbc)
coxfit <- coxph(Surv(time,status!=0)~age+log(bili)+log(albumin)+factor(edema)+sex,data=pbc)
pubcox <- publish(coxfit)
plot(pubcox,cex=1.5,xratio=c(0.4,0.2))




cleanEx()
nameEx("plotConfidence")
### * plotConfidence

flush(stderr()); flush(stdout())

### Name: plotConfidence
### Title: Plot confidence intervals
### Aliases: plotConfidence

### ** Examples


library(Publish)
data(CiTable) 

## A first draft version of the plot is obtained as follows
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper","p")],
          labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")])

## if argument labels is a named list the table is subdivided:
labellist <- split(CiTable[,c("Dose","Mean","SD","n")],CiTable[,"Drug"])
labellist
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")], labels=labellist)


## The graph consist of at most three columns:
##
## column 1: labels
## column 2: printed values of the confidence intervals
## column 3: graphical presentation of the confidence intervals
##
## NOTE: column 3 appears always, the user decides if also
##       column 1, 2 should appear
##
## The columns are arranged with the function layout
## and the default order is 1,3,2 such that the graphical
## display of the confidence intervals appears in the middle
##
## the order of appearance of the three columns can be changed as follows
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               order=c(1,3,2))
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               order=c(2,3,1))
## if there are only two columns the order is 1, 2
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               values=FALSE,
               order=c(2,1))
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               values=FALSE,
               order=c(1,2))



## The relative size of the columns needs to be controlled manually
## by using the argument xratio. If there are only two columns
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               xratio=c(0.4,0.15))

## The amount of space on the left and right margin can be controlled
## as follows:
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               xratio=c(0.4,0.15),
               leftmargin=0.1,rightmargin=0.00)

## The actual size of the current graphics device determines
## the size of the figures and the space between them.
## The sizes and line widths are increased as follows:
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               xlab="Hazard ratio",
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               points.cex=3,
               cex=2,
               lwd=3,
               xaxis.lwd=1.3,
               xaxis.cex=1.3)
## Note that 'cex' of axis ticks is controlled via 'par' but
## cex of the label via argument 'cex' of 'mtext'.
## The sizes and line widths are decreased as follows:
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               cex=0.8,
               lwd=0.8,
               xaxis.lwd=0.8,
               xaxis.cex=0.8)

## Another good news is that all figures can be controlled separately

## The size of the graphic device can be controlled in the usual way, e.g.:
## Not run: 
##D     pdf("~/tmp/testCI.pdf",width=8,height=8)
##D     plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
##D                    labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")])
##D     dev.off()
## End(Not run)

## More control of the x-axis and confidence intervals that
## stretch outside the x-range end in an arrow. 
## the argument xlab.line adjusts the distance of the x-axis
## label from the graph
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               xlab="Hazard ratio",
               xlab.line=1.8,
               xaxis.at=c(0.8,1,1.3),
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               xlim=c(0.8,1.3))

## log-scale
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               xlab="Hazard ratio",
               xlab.line=1.8,
               xaxis.at=c(0.8,1,1.3),
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               xlim=c(0.8,1.3),plot.log="x")
## More pronounced arrows
## Coloured xlab expression
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               xlab=expression(HR[1](s)),
               xlab.line=1.8,
               xlab.col="darkred",
               extremearrows.angle=50,
               extremearrows.length=0.1,
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               xlim=c(0.8,1.3))

## Controlling the labels and their titles
## and the values and their titles
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               xlab="Hazard ratio",
               title.values=expression(bold(HR (CI[95]))),
               title.labels=c("Drug/Time","Dose","Mean","St.dev.","N"),
               factor.reference.pos=c(1,10,19),
               factor.reference.pch=16,
               cex=1.3,
               xaxis.at=c(0.75,1,1.25,1.5,2))

## For factor reference groups, one may want to replace the
## confidence intervals by the word Reference, as in the previous example.
## To change the word 'Reference' we use the argument factor.reference.label:
## To change the plot symbol for the reference lines factor.reference.pch
## To remove the plot symbol in the reference lines use 'NA' as follows:
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               xlab="Hazard ratio",
               factor.reference.label="Ref",
               title.values=expression(bold(HR (CI[95]))),
               title.labels=c("Drug/Time","Dose","Mean","St.dev.","N"),
               factor.reference.pos=c(1,10,19),
               factor.reference.pch=NA,
               cex=1.3,
               xaxis.at=c(0.75,1,1.25,1.5,2))


## changing the style of the graphical confidence intervals
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               xlab="Hazard ratio",
               factor.reference.pos=c(1,10,19),
               points.pch=15,
               points.col=rainbow(27),
               points.cex=2,
               arrows.col="darkblue",
               cex=1.3,
               order=c(1,3,2),
               xaxis.at=c(0.75,1,1.25,1.5))

## the values column of the graph can have multiple columns as well
## to illustrate this we create the confidence intervals
## before calling the function and then cbind them
## to the pvalues
HR <- pubformat(CiTable[,6])
CI95 <- formatCI(lower=CiTable[,7],upper=CiTable[,8],format="(l-u)")
pval <- format.pval(CiTable[,9],digits=3,eps=10^{-3})
pval[pval=="NA"] <- ""
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               values=list("HR"=HR,"CI-95"=CI95,"P-value"=pval),
               cex=1.2,
               xratio=c(0.5,0.3))

## Finally, vertical columns can be delimited with background color
## NOTE: this may slow things down and potentially create
##       large figures (many bytes)
col1 <- rep(c(prodlim::dimColor("green",density=22),
              prodlim::dimColor("green")),length.out=9)
col2 <- rep(c(prodlim::dimColor("orange",density=22),
              prodlim::dimColor("orange")),length.out=9)
col3 <- rep(c(prodlim::dimColor("blue",density=22),
              prodlim::dimColor("blue")),length.out=9)
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               stripes=c(1,0,1),
               stripes.col=c(col1,col2,col3))
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               stripes=c(1,1,1),
               stripes.col=c(col1,col2,col3))

threegreens <- c(prodlim::dimColor("green",density=55),
                 prodlim::dimColor("green",density=33),
                 prodlim::dimColor("green",density=22))
plotConfidence(x=CiTable[,c("HazardRatio","lower","upper")],
               labels=CiTable[,c("Drug.Time","Dose","Mean","SD","n")],
               values=FALSE,
               xlim=c(0.75,1.5),
               stripes=c(1,1,1),
               xratio=c(0.5,0.15),
               stripes.horizontal=c(0,9,18,27)+0.5,
               stripes.col=threegreens)






cleanEx()
nameEx("print.ci")
### * print.ci

flush(stderr()); flush(stdout())

### Name: print.ci
### Title: Print confidence intervals
### Aliases: print.ci

### ** Examples

library(lava)
m <- lvm(Y~X)
m <- categorical(m,Y~X,K=4)
set.seed(4)
d <- sim(m,24)
ci.mean(Y~X,data=d)
x <- ci.mean(Y~X,data=d)
print(x,format="(l,u)")



cleanEx()
nameEx("print.table2x2")
### * print.table2x2

flush(stderr()); flush(stdout())

### Name: print.table2x2
### Title: print results of 2x2 contingency table analysis
### Aliases: print.table2x2

### ** Examples

table2x2(table("marker"=rbinom(100,1,0.4),"response"=rbinom(100,1,0.1)))
table2x2(matrix(c(71,18,38,8),ncol=2),stats="table")
table2x2(matrix(c(71,18,38,8),ncol=2),stats=c("rr","fisher"))



cleanEx()
nameEx("pubformat")
### * pubformat

flush(stderr()); flush(stdout())

### Name: pubformat
### Title: Format numbers for publication
### Aliases: pubformat

### ** Examples


pubformat(c(0.000143,12.8,1))
pubformat(c(0.000143,12.8,1),handler="format")
pubformat(c(0.000143,12.8,1),handler="format",trim=TRUE)
pubformat(c(0.000143,12.8,1),handler="prettyNum")



cleanEx()
nameEx("publish.CauseSpecificCox")
### * publish.CauseSpecificCox

flush(stderr()); flush(stdout())

### Name: publish.CauseSpecificCox
### Title: Tabulizing cause-specific hazard ratio from all causes with
###   confidence limits and Wald test p-values.
### Aliases: publish.CauseSpecificCox

### ** Examples

library(riskRegression)
library(prodlim)
library(pec)
library(survival)
data(Melanoma,package="riskRegression")
fit1 <- CSC(list(Hist(time,status)~sex,Hist(time,status)~invasion+epicel+age),
            data=Melanoma)
publish(fit1)
publish(fit1,pvalue.stars=TRUE)
publish(fit1,factor.reference="inline",units=list("age"="years"))



cleanEx()
nameEx("publish.MIresult")
### * publish.MIresult

flush(stderr()); flush(stdout())

### Name: publish.MIresult
### Title: Present logistic regression and Cox regression obtained with
###   mitools::MIcombine based on smcfcs::smcfcs multiple imputation
###   analysis
### Aliases: publish.MIresult

### ** Examples


## Not run: 
##D ## continuous outcome: linear regression
##D # lava some data with missing values
##D library(riskRegression)
##D set.seed(7)
##D d=sampleData(78)
##D ## generate missing values
##D d[X1==1,X6:=NA] 
##D d[X2==1,X3:=NA]
##D d=d[,.(X8,X4,X3,X6,X7)]
##D sapply(d,function(x)sum(is.na(x)))
##D 
##D # multiple imputation (should set m to a large value)
##D 
##D library(smcfcs)
##D library(mitools)
##D set.seed(17)
##D f= smcfcs(d,smtype="lm",
##D            smformula=X8~X4+X3+X6+X7,
##D            method=c("","","logreg","norm",""),m=3)
##D ccfit=lm(X8~X4+X3+X6+X7,data=d)
##D mifit=MIcombine(with(imputationList(f$impDatasets),
##D                 lm(X8~X4+X3+X6+X7)))
##D publish(mifit,fit=ccfit,data=d)
##D publish(ccfit)
##D 
##D ## binary outcome
##D # lava some data with missing values
##D library(riskRegression)
##D set.seed(7)
##D db=sampleData(78,outcome="binary")
##D ## generate missing values
##D db[X1==1,X6:=NA] 
##D db[X2==1,X3:=NA]
##D db=db[,.(Y,X4,X3,X6,X7)]
##D sapply(db,function(x)sum(is.na(x)))
##D 
##D # multiple imputation (should set m to a large value)
##D library(smcfcs)
##D library(mitools)
##D set.seed(17)
##D fb= smcfcs(db,smtype="logistic",
##D            smformula=Y~X4+X3+X6+X7,
##D            method=c("","","logreg","norm",""),m=2)
##D ccfit=glm(Y~X4+X3+X6+X7,family="binomial",data=db)
##D mifit=MIcombine(with(imputationList(fb$impDatasets),
##D                 glm(Y~X4+X3+X6+X7,family="binomial")))
##D publish(mifit,fit=ccfit)
##D publish(ccfit)
##D 
##D ## survival: Cox regression
##D library(smcfcs)
##D library(mitools)
##D library(survival)
##D # lava some data with missing values
##D library(riskRegression)
##D set.seed(7)
##D ds=sampleData(78,outcome="survival")
##D ## generate missing values
##D ds[X5==1,X6:=NA] 
##D ds[X2==1,X3:=NA]
##D ds=ds[,.(time,event,X4,X3,X6,X7)]
##D sapply(ds,function(x)sum(is.na(x)))
##D 
##D set.seed(17)
##D fs= smcfcs(ds,smtype="coxph",
##D            smformula="Surv(time,event)~X4+X3+X6+X7",
##D            method=c("","","","logreg","norm",""),m=2)
##D ccfit=coxph(Surv(time,event)~X4+X3+X6+X7,data=ds)
##D mifit=MIcombine(with(imputationList(fs$impDatasets),
##D                 coxph(Surv(time,event)~X4+X3+X6+X7)))
##D publish(mifit,fit=ccfit,data=ds)
##D publish(ccfit)
##D 
##D ## competing risks: Cause-specific Cox regression 
##D library(survival)
##D library(smcfcs)
##D library(mitools)
##D # lava some data with missing values
##D library(riskRegression)
##D set.seed(7)
##D dcr=sampleData(78,outcome="competing.risks")
##D ## generate missing values
##D dcr[X5==1,X6:=NA] 
##D dcr[X2==1,X3:=NA]
##D dcr=dcr[,.(time,event,X4,X3,X6,X7)]
##D sapply(dcr,function(x)sum(is.na(x)))
##D 
##D set.seed(17)
##D fcr= smcfcs(dcr,smtype="compet",
##D            smformula=c("Surv(time,event==1)~X4+X3+X6+X7",
##D                        "Surv(time,event==2)~X4+X3+X6+X7"),
##D            method=c("","","","logreg","norm",""),m=2)
##D ## cause 2 
##D ccfit2=coxph(Surv(time,event==2)~X4+X3+X6+X7,data=dcr)
##D mifit2=MIcombine(with(imputationList(fcr$impDatasets),
##D                 coxph(Surv(time,event==2)~X4+X3+X6+X7)))
##D publish(mifit2,fit=ccfit2,data=dcr)
##D publish(ccfit2)
## End(Not run) 




cleanEx()
nameEx("publish.Score")
### * publish.Score

flush(stderr()); flush(stdout())

### Name: publish.Score
### Title: Publish predictive accuracy results
### Aliases: publish.Score

### ** Examples

library(riskRegression)
library(survival)
learn = sampleData(100)
val= sampleData(100)
f1=CSC(Hist(time,event)~X1+X8,data=learn)
f2=CSC(Hist(time,event)~X1+X5+X6+X8,learn)
xs=Score(list(f1,f2),data=val,formula=Hist(time,event)~1)
publish(xs)




cleanEx()
nameEx("publish.ci")
### * publish.ci

flush(stderr()); flush(stdout())

### Name: publish.ci
### Title: Publish tables with confidence intervals
### Aliases: publish.ci

### ** Examples


data(Diabetes)
publish(ci.mean(chol~location+gender,data=Diabetes),org=TRUE)




cleanEx()
nameEx("publish.coxph")
### * publish.coxph

flush(stderr()); flush(stdout())

### Name: publish.coxph
### Title: Tabulize hazard ratios with confidence intervals and p-values.
### Aliases: publish.coxph

### ** Examples

library(survival)
data(pbc)
pbc$edema <- factor(pbc$edema,
             levels=c("0","0.5","1"), labels=c("0","0.5","1"))
fit = coxph(Surv(time,status!=0)~age+sex+edema+log(bili)+log(albumin),
            data=na.omit(pbc))
publish(fit)
## forest plot
plot(publish(fit),cex=1.3)

publish(fit,ci.digits=2,pvalue.eps=0.01,pvalue.digits=2,pvalue.stars=TRUE)
publish(fit,ci.digits=2,ci.handler="prettyNum",pvalue.eps=0.01,
        pvalue.digits=2,pvalue.stars=TRUE)
publish(fit, ci.digits=2, ci.handler="sprintf", pvalue.eps=0.01,
        pvalue.digits=2,pvalue.stars=TRUE, ci.trim=FALSE)

fit2 = coxph(Surv(time,status!=0)~age+sex+edema+log(bili,base=2)+log(albumin)+log(protime),
    data=na.omit(pbc))
publish(fit2)

# with cluster variable
fit3 = coxph(Surv(time,status!=0)~age+cluster(sex)+edema+log(bili,base=2)
                                    +log(albumin)+log(protime),
    data=na.omit(pbc))
publish(fit3)

# with strata and cluster variable
fit4 = coxph(Surv(time,status!=0)~age+cluster(sex)+strata(edema)+log(bili,base=2)
                 +log(albumin)+log(protime),
    data=pbc)
publish(fit4)




cleanEx()
nameEx("publish.glm")
### * publish.glm

flush(stderr()); flush(stdout())

### Name: publish.glm
### Title: Tabulize regression coefficients with confidence intervals and
###   p-values.
### Aliases: publish.glm

### ** Examples

data(Diabetes)
## Linear regression
f = glm(bp.2s~frame+gender+age,data=Diabetes)
publish(f)
publish(f,factor.reference="inline")
publish(f,pvalue.stars=TRUE)
publish(f,ci.format="(l,u)")

### interaction
fit = glm(bp.2s~frame+gender*age,data=Diabetes)
summary(fit)
publish(fit)

Fit = glm(bp.2s~frame*gender+age,data=Diabetes)
publish(Fit)

## Logistic regression
Diabetes$hyper1 <- factor(1*(Diabetes$bp.1s>140))
lrfit <- glm(hyper1~frame+gender+age,data=Diabetes,family=binomial)
publish(lrfit)

### interaction
lrfit1 <- glm(hyper1~frame+gender*age,data=Diabetes,family=binomial)
publish(lrfit1)

lrfit2 <- glm(hyper1~frame*gender+age,data=Diabetes,family=binomial)
publish(lrfit2)

## Poisson regression
data(trace)
trace <- Units(trace,list("age"="years"))
fit <- glm(dead ~ smoking+sex+age+Time+offset(log(ObsTime)), family="poisson",data=trace)
rtf <- regressionTable(fit,factor.reference = "inline")
summary(rtf)
publish(fit)

## gls regression
library(nlme)
library(lava)
m <- lvm(Y ~ X1 + gender + group + Interaction)
distribution(m, ~gender) <- binomial.lvm()
distribution(m, ~group) <- binomial.lvm(size = 2)
constrain(m, Interaction ~ gender + group) <- function(x){x[,1]*x[,2]}
d <- sim(m, 1e2)
d$gender <- factor(d$gender, labels = letters[1:2])
d$group <- factor(d$group)

e.gls <- gls(Y ~ X1 + gender*group, data = d,
             weights = varIdent(form = ~1|group))
publish(e.gls)

## lme
library(nlme)
fm1 <- lme(distance ~ age*Sex, 
            random = ~1|Subject,
            data = Orthodont) 
res <- publish(fm1)



cleanEx()
nameEx("publish.htest")
### * publish.htest

flush(stderr()); flush(stdout())

### Name: publish.htest
### Title: Pretty printing of test results.
### Aliases: publish.htest

### ** Examples

data(Diabetes)
publish(t.test(bp.2s~gender,data=Diabetes))
publish(wilcox.test(bp.2s~gender,data=Diabetes))
publish(with(Diabetes,t.test(bp.2s,bp.1s,paired=TRUE)))
publish(with(Diabetes,wilcox.test(bp.2s,bp.1s,paired=TRUE)))




cleanEx()
nameEx("publish.matrix")
### * publish.matrix

flush(stderr()); flush(stdout())

### Name: publish.matrix
### Title: Publishing a matrix in raw, org, latex, or muse format
### Aliases: publish.matrix

### ** Examples


x <- matrix(1:12,ncol=3)
publish(x)

# rounding the numeric part of data mixtures 
y <- cbind(matrix(letters[1:12],ncol=3),x,matrix(rnorm(12),ncol=3))
publish(y,digits=1)

publish(x,inter.lines=list("1"="text between line 1 and line 2",
                          "3"="text between line 3 and line 4"))




cleanEx()
nameEx("publish.riskRegression")
### * publish.riskRegression

flush(stderr()); flush(stdout())

### Name: publish.riskRegression
### Title: Publishing results of riskRegression
### Aliases: publish.riskRegression

### ** Examples

library(prodlim)
library(riskRegression)
library(lava)
library(survival)
set.seed(20)
d <- SimCompRisk(20)
f <- ARR(Hist(time,event)~X1+X2,data=d,cause=1)
publish(f)
publish(f,digits=c(1,3))



cleanEx()
nameEx("publish.summary.aov")
### * publish.summary.aov

flush(stderr()); flush(stdout())

### Name: publish.summary.aov
### Title: Format summary table of aov results
### Aliases: publish.summary.aov

### ** Examples

 
data(Diabetes)
f <- glm(bp.1s~age+chol+gender+location,data=Diabetes)
publish(summary(aov(f)),digits=c(1,2))




cleanEx()
nameEx("publish.survdiff")
### * publish.survdiff

flush(stderr()); flush(stdout())

### Name: publish.survdiff
### Title: Alternative summary of survdiff results
### Aliases: publish.survdiff

### ** Examples

library(survival)
data(pbc)
sd <- survdiff(Surv(time,status!=0)~sex,data=pbc)
publish(sd)
publish(sd,digits=c(3,2))




cleanEx()
nameEx("regressionTable")
### * regressionTable

flush(stderr()); flush(stdout())

### Name: regressionTable
### Title: Regression table
### Aliases: regressionTable

### ** Examples

# linear regression
data(Diabetes)
f1 <- glm(bp.1s~age+gender+frame+chol,data=Diabetes)
summary(regressionTable(f1))
summary(regressionTable(f1,units=list("chol"="mmol/L","age"="years")))
## with interaction
f2 <- glm(bp.1s~age*gender+frame+chol,data=Diabetes)
summary(regressionTable(f2))
#Add reference values
summary(regressionTable(f2))
f3 <- glm(bp.1s~age+gender*frame+chol,data=Diabetes)
publish(f3)
regressionTable(f3)

# logistic regression
Diabetes$hyp1 <- factor(1*(Diabetes$bp.1s>140))
l1 <- glm(hyp1~age+gender+frame+chol,data=Diabetes,family="binomial")
regressionTable(l1)
publish(l1)
plot(regressionTable(l1))

## with interaction
l2 <- glm(hyp1~age+gender+frame*chol,data=Diabetes,family="binomial")
regressionTable(l2)
l3 <- glm(hyp1~age*gender+frame*chol,data=Diabetes,family="binomial")
regressionTable(l3)

# Cox regression
library(survival)
data(pbc)
pbc$edema <- factor(pbc$edema,levels=c("0","0.5","1"),labels=c("0","0.5","1"))
c1 <- coxph(Surv(time,status!=0)~log(bili)+age+protime+sex+edema,data=pbc)
regressionTable(c1)
# with interaction
c2 <- coxph(Surv(time,status!=0)~log(bili)+age+protime*sex+edema,data=pbc)
regressionTable(c2)
c3 <- coxph(Surv(time,status!=0)~edema*log(bili)+age+protime+sex+edema+edema:sex,data=pbc)
regressionTable(c3)


## gls regression
library(nlme)
library(lava)
m <- lvm(Y ~ X1 + gender + group + Interaction)
distribution(m, ~gender) <- binomial.lvm()
distribution(m, ~group) <- binomial.lvm(size = 2)
constrain(m, Interaction ~ gender + group) <- function(x){x[,1]*x[,2]}
d <- sim(m, 1e2)
d$gender <- factor(d$gender, labels = letters[1:2])
d$group <- factor(d$group)

e.gls <- gls(Y ~ X1 + gender*group, data = d,
             weights = varIdent(form = ~1|group))
regressionTable(e.gls)
summary(regressionTable(e.gls))



cleanEx()
nameEx("spaghettiogram")
### * spaghettiogram

flush(stderr()); flush(stdout())

### Name: spaghettiogram
### Title: Spaghettiogram
### Aliases: spaghettiogram Spaghettiogram

### ** Examples


data(SpaceT)
Spaghettiogram(HR~Status+id(ID),
               data=SpaceT)



cleanEx()
nameEx("specialFrame")
### * specialFrame

flush(stderr()); flush(stdout())

### Name: specialFrame
### Title: Special frame
### Aliases: specialFrame

### ** Examples


## Here are some data with an event time and no competing risks
## and two covariates X1 and X2.
## Suppose we want to declare that variable X1 is treated differently
## than variable X2. For example, X1 could be a cluster variable, or
## X1 should have a proportional effect on the outcome.
d <- data.frame(y=1:7,
                X2=c(2.24,3.22,9.59,4.4,3.54,6.81,5.05),
                X3=c(1,1,1,1,0,0,1),
                X4=c(44.69,37.41,68.54,38.85,35.9,27.02,41.84),
                X1=factor(c("a","b","a","c","c","a","b"),
                    levels=c("c","a","b")))
## define special functions prop and cluster
prop <- function(x)x
cluster <- function(x)x
## We pass a formula and the data
e <- specialFrame(y~prop(X1)+X2+cluster(X3)+X4,
                  data=d,
                  specials=c("prop","cluster"))
## The first element is the response
e$response
## The other elements are the design, i.e., model.matrix for the non-special covariates
e$design
## and a data.frame for the special covariates
e$prop
## The special covariates can be returned as a model.matrix 
e2 <- specialFrame(y~prop(X1)+X2+cluster(X3)+X4,
                   data=d,
                   specials=c("prop","cluster"),
                   specials.design=TRUE)
e2$prop
## and the non-special covariates can be returned as a data.frame
e3 <- specialFrame(y~prop(X1)+X2+cluster(X3)+X4,
                   data=d,
                   specials=c("prop","cluster"),
                   specials.design=TRUE,
                   unspecials.design=FALSE)
e3$design



cleanEx()
nameEx("splinePlot.lrm")
### * splinePlot.lrm

flush(stderr()); flush(stdout())

### Name: splinePlot.lrm
### Title: Plot predictions of logistic regression
### Aliases: splinePlot.lrm

### ** Examples

data(Diabetes)
Diabetes$hypertension=  1*(Diabetes$bp.1s>140)
library(rms)
uu <- datadist(Diabetes)
options(datadist="uu")
fit=lrm(hypertension~rcs(age)+gender+hdl,data=Diabetes)
splinePlot.lrm(fit,xvar="age",xvalues=seq(30,50,1))



cleanEx()
nameEx("stripes")
### * stripes

flush(stderr()); flush(stdout())

### Name: stripes
### Title: Background and grid color control.
### Aliases: stripes
### Keywords: survival

### ** Examples



plot(0,0)
backGround(bg="beige",fg="red",vertical=0,horizontal=0)

plot(0,0)
stripes(col=c("yellow","green"),gridcol="red",xlim=c(-1,1),horizontal=seq(0,1,.1))
stripes(col=c("yellow","green"),gridcol="red",horizontal=seq(0,1,.1))




cleanEx()
nameEx("summary.ci")
### * summary.ci

flush(stderr()); flush(stdout())

### Name: summary.ci
### Title: Summarize confidence intervals
### Aliases: summary.ci

### ** Examples

library(lava)
m <- lvm(Y~X)
m <- categorical(m,Y~X,K=4)
set.seed(4)
d <- sim(m,24)
ci.mean(Y~X,data=d)
x <- summary(ci.mean(Y~X,data=d),digits=2)
x
x <- summary(ci.mean(Y~X,data=d),format="(u,l)",digits=2)
x <- summary(ci.mean(Y~X,data=d),format="(u,l)",digits=1,se=TRUE)
x <- summary(ci.mean(Y~X,data=d),format="(u,l)",digits=1,handler="format")
x <- summary(ci.mean(Y~X,data=d),format="(u,l)",digits=1,handler="prettyNum")



cleanEx()
nameEx("summary.regressionTable")
### * summary.regressionTable

flush(stderr()); flush(stdout())

### Name: summary.regressionTable
### Title: Formatting regression tables
### Aliases: summary.regressionTable print.summary.regressionTable

### ** Examples

library(survival)
data(pbc)
pbc$edema <- factor(pbc$edema,levels=c("0","0.5","1"),labels=c("0","0.5","1"))
fit = coxph(Surv(time,status!=0)~age+sex+edema+log(bili)+log(albumin)+log(protime),
            data=pbc)
u=summary(regressionTable(fit))
u$regressionTable
u$rawTable
summary(regressionTable(fit),handler="prettyNum")
summary(regressionTable(fit),handler="format")
summary(regressionTable(fit),handler="sprintf",digits=c(2,2),pValue.stars=TRUE)



cleanEx()
nameEx("summary.univariateTable")
### * summary.univariateTable

flush(stderr()); flush(stdout())

### Name: summary.univariateTable
### Title: Preparing univariate tables for publication
### Aliases: summary.univariateTable

### ** Examples

data(Diabetes)
u <- univariateTable(gender~age+location+Q(BMI)+height+weight,
                data=Diabetes)
summary(u)
summary(u,n=NULL)
summary(u,pvalue.digits=2,"age"="Age (years)","height"="Body height (cm)")

u2 <- univariateTable(location~age+AgeGroups+gender+height+weight,
                data=Diabetes)
summary(u2)
summary(u2,drop.reference=TRUE)
## same but more flexible
summary(u2,drop.reference=c("binary"))
## same but even more flexible
summary(u2,drop.reference=c("gender"))





cleanEx()
nameEx("sutable")
### * sutable

flush(stderr()); flush(stdout())

### Name: sutable
### Title: Fast summary of a univariate table
### Aliases: sutable

### ** Examples

data(Diabetes)
sutable(gender~age+location+Q(BMI)+height+weight,data=Diabetes,BMI="Body mass index (kg/m^2)")



cleanEx()
nameEx("table2x2")
### * table2x2

flush(stderr()); flush(stdout())

### Name: table2x2
### Title: 2x2 table calculus for teaching
### Aliases: table2x2

### ** Examples

table2x2(table("marker"=rbinom(100,1,0.4),"response"=rbinom(100,1,0.1)))
table2x2(matrix(c(71,18,38,8),ncol=2),stats="table")
table2x2(matrix(c(71,18,38,8),ncol=2),stats=c("rr","fisher"))



cleanEx()
nameEx("trace")
### * trace

flush(stderr()); flush(stdout())

### Name: trace
### Title: trace data
### Aliases: trace
### Keywords: datasets

### ** Examples


data(trace)
Units(trace,list("age"="years"))
fit <- glm(dead ~ smoking+sex+age+Time+offset(log(ObsTime)), family="poisson",data=trace)
rtf <- regressionTable(fit,factor.reference = "inline")
summary(rtf)
publish(fit)




cleanEx()
nameEx("univariateTable")
### * univariateTable

flush(stderr()); flush(stdout())

### Name: univariateTable
### Title: Univariate table
### Aliases: univariateTable utable

### ** Examples

data(Diabetes)
univariateTable(~age,data=Diabetes)
univariateTable(~gender,data=Diabetes)
univariateTable(~age+gender+ height+weight,data=Diabetes)
## same thing but less typing
utable(~age+gender+ height+weight,data=Diabetes)

## summary by location: 
univariateTable(location~Q(age)+gender+height+weight,data=Diabetes)
## continuous variables marked with Q() are (by default) summarized
## with median (IQR) and kruskal.test (with two groups equivalent to wilcox.test)
## variables not marked with Q() are (by default) summarized
## with mean (sd) and anova.glm(...,test="Chisq")
## the p-value of anova.glm with only two groups is similar
## but not exactly equal to that of a t.test
## categorical variables are (by default) summarized by count
## (percent) and anova.glm(...,family=binomial,test="Chisq")

## export result to csv
table1 = summary(univariateTable(location~age+gender+height+weight,data=Diabetes),
show.pvalues=FALSE)
# write.csv(table1,file="~/table1.csv",rownames=FALSE)

## change labels and values
utable(location~age+gender+height+weight,data=Diabetes,
       age="Age (years)",gender="Sex",
       gender.female="Female",
       gender.male="Male",
       height="Body height (inches)",
       weight="Body weight (pounds)")

## Use quantiles and rank tests for some variables and mean and standard deviation for others
univariateTable(gender~Q(age)+location+Q(BMI)+height+weight,
                data=Diabetes)

## Factor with more than 2 levels
Diabetes$AgeGroups <- cut(Diabetes$age,
                          c(19,29,39,49,59,69,92),
                          include.lowest=TRUE)
univariateTable(location~AgeGroups+gender+height+weight,
                data=Diabetes)

## Row percent
univariateTable(location~gender+age+AgeGroups,
                data=Diabetes,
                column.percent=FALSE)

## change of frequency format
univariateTable(location~gender+age+AgeGroups,
                data=Diabetes,
                column.percent=FALSE,
                freq.format="percent(x) (n=count(x))")

## changing Labels
u <- univariateTable(location~gender+AgeGroups+ height + weight,
                     data=Diabetes,
                     column.percent=TRUE,
                     freq.format="count(x) (percent(x))")
summary(u,"AgeGroups"="Age (years)","height"="Height (inches)")

## more than two groups
Diabetes$frame=factor(Diabetes$frame,levels=c("small","medium","large"))
univariateTable(frame~gender+BMI+age,data=Diabetes)

Diabetes$sex=as.numeric(Diabetes$gender)
univariateTable(frame~sex+gender+BMI+age,
                data=Diabetes,freq.format="count(x) (percent(x))")

## multiple summary formats
## suppose we want for some reason mean (range) for age
## and median (range) for BMI.
## method 1:
univariateTable(frame~Q(age)+BMI,
                data=Diabetes,
                Q.format="mean(x) (range(x))",
                summary.format="median(x) (range(x))")
## method 2:
u1 <- summary(univariateTable(frame~age,
                              data=na.omit(Diabetes),
                              summary.format="mean(x) (range(x))"))
u2 <- summary(univariateTable(frame~BMI,
                              data=na.omit(Diabetes),
                              summary.format="median(x) (range(x))"))
publish(rbind(u1,u2),digits=2)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
