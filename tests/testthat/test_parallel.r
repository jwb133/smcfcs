library(smcfcs)
library(survival)
context("Testing parallel computation")

test_that("Cox imputation runs in parallel", {
  expect_error({
    set.seed(1234)
    n <- 10000
    z <- rnorm(n)
    x <- z+rnorm(n)
    t <- -log(runif(n))/(1*exp(x+z))
    d <- 1*(t<10)
    t[d==0] <- 10
    x[(runif(n)<0.5)] <- NA

    simData <- data.frame(t,d,x,z)

    imps <- smcfcs.parallel(originaldata=simData, smtype="coxph", smformula="Surv(t, d)~x+z",
                   method=c("", "", "norm", ""), m=5, n_cores=2)
  }, NA)
})

test_that("smcfcs.parallel works with 1 imputation per core", {
  expect_error({
    set.seed(1234)
    n <- 100
    z <- rnorm(n)
    x <- z+rnorm(n)
    t <- -log(runif(n))/(1*exp(x+z))
    d <- 1*(t<10)
    t[d==0] <- 10
    x[(runif(n)<0.5)] <- NA

    simData <- data.frame(t,d,x,z)

    imps <- smcfcs.parallel(originaldata=simData, smtype="coxph", smformula="Surv(t, d)~x+z",
                            method=c("", "", "norm", ""), m=5, n_cores=5)
  }, NA)
})

test_that("smcfcs.parallel throws error if you specify more cores than necessary for 1 imp per core", {
  expect_error({
    set.seed(1234)
    n <- 100
    z <- rnorm(n)
    x <- z+rnorm(n)
    t <- -log(runif(n))/(1*exp(x+z))
    d <- 1*(t<10)
    t[d==0] <- 10
    x[(runif(n)<0.5)] <- NA

    simData <- data.frame(t,d,x,z)

    imps <- smcfcs.parallel(originaldata=simData, smtype="coxph", smformula="Surv(t, d)~x+z",
                            method=c("", "", "norm", ""), m=5, n_cores=10)
  }, NA)
})

test_that("Test case cohort imputation runs", {
  expect_error({
    set.seed(1234)

    #run simulation study
    n <- 10000

    z <- rnorm(n)
    x <- 1*(runif(n)<(exp(z)/(1+exp(z))))
    t <- -log(runif(n))/(0.01*exp(x+z))
    d <- 1*(t<2)
    t[d==0] <- 2
    x[(runif(n)<0.5)] <- NA

    fullcohortdata <- data.frame(t,d,x,z)
    fullcohortdata$in.subco <- 0
    #we sample a 10% subcohort
    fullcohortdata$in.subco[sample(n, size=n*0.1)] <- 1
    fullcohortdata$id <- 1:n

    ccdata <- fullcohortdata[(fullcohortdata$in.subco==1) | (fullcohortdata$d==1),]
    ccdata$entertime <- 0
    ccdata$entertime[ccdata$in.subco==0] <- ccdata$t[ccdata$in.subco==0] - 0.000001
    imps <- smcfcs.parallel(originaldata=ccdata, smformula="Surv(entertime, t, d)~x+z", sampfrac=0.1,
                              in.subco="in.subco", method=c("", "", "logreg", "", "", "", ""),
                            smcfcs_func = "smcfcs.casecohort", m=5, n_cores=2)

  }, NA)
})

test_that("Nested case control imputation runs", {
  expect_error({
    set.seed(1234)

    n <- 10000

    z <- rnorm(n)
    x <- z+rnorm(n)
    t <- -log(runif(n))/(0.01*exp(x+z))
    d <- 1*(t<1)
    t[d==0] <- 1
    x[(runif(n)<0.5)] <- NA

    fullcohortdata <- data.frame(t,d,x,z)
    fullcohortdata$id <- 1:n

    # Compute number at risk at each event time using the full cohort data
    nrisk.fit <- survfit(Surv(t,d)~1,data=fullcohortdata)
    ord.t.d1 <- order(fullcohortdata$t[fullcohortdata$d==1])

    m=1 #1 control per case
    ncc=NULL

    no.sample=0
    for (i in which(fullcohortdata$d==1))
    {
      #select control(s) for nested case-control
      possible.controls <- which(fullcohortdata$t>=fullcohortdata$t[i])
      #remove the case from this vector
      possible.controls <- possible.controls[which(possible.controls!=i)]

      if (length(possible.controls)>=m){
        controls <- sample(possible.controls,m)
        numAtRisk <- 1+length(possible.controls)
        ncc <- rbind(ncc,cbind(fullcohortdata[i,],numrisk=numAtRisk))
        ncc <- rbind(ncc,cbind(fullcohortdata[controls,], numrisk=numAtRisk))
        no.sample <- no.sample+1}
    }

    ncc$setno <- rep(1:no.sample,each=m+1)
    ncc$case <- rep(c(1,rep(0,m)),no.sample)

    predictorMatrix <- matrix(0,nrow=dim(ncc)[2],ncol=dim(ncc)[2])
    predictorMatrix[which(colnames(ncc)=="x"),c(which(colnames(ncc)=="z"))] <- 1

    imps <- smcfcs.parallel(originaldata=ncc,set="setno",nrisk="numrisk",event="d",smformula="Surv(t,case)~x+z+strata(setno)",
                            method=c("", "", "norm", "", "", "", "", ""),predictorMatrix=predictorMatrix,
                            smcfcs_func = "smcfcs.nestedcc", m=5, n_cores=2)

    }, NA)
})

test_that("DTSAM imputation runs", {
  expect_error({
    set.seed(1234)
    n <- 1000
    x1 <- 1*(runif(n)<0.5)
    x2 <- x1+rnorm(n)
    T <- 10
    #define vector of intercepts
    alpha <- -3+0.2*(1:T)
    beta <- c(1,-1)
    yMat <- array(0, dim=c(n,T))
    for (i in 1:T) {
      yMat[,i] <- 1*(runif(n)<expit(alpha[i]+beta[1]*x1+beta[2]*x2))
    }
    failtime <- apply(yMat,1,function(x) which(x==1)[1])
    #event indicator
    d <- rep(1,n)
    d[is.na(failtime)] <- 0
    failtime[is.na(failtime)] <- 10
    mean(d)

    simData <- data.frame(x1=x1, x2=x2, failtime=failtime, d=d)
    table(simData$failtime,simData$d)

    simData$x1[runif(n)<expit((x2-mean(x2)/sd(x2)))] <- NA

    #impute using smcfcs.dtsam
    M <- 5
    imps <- smcfcs.parallel(originaldata=simData, smformula="Surv(failtime,d)~x1+x2",
                         method=c("logreg","", "", ""),
                         smcfcs_func = "smcfcs.dtsam", m=M, n_cores=2)
  }, NA)
})
