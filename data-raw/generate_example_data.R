#linear substantive model with quadratic covariate effect
set.seed(1234)
n <- 1000
x <- rnorm(n)
z <- rnorm(n)
#v is an auxiliary variable
v <- x+rnorm(n)
y <- 1+z+x+x^2+rnorm(n)

#make some x values missing
xobsxb <- (y-mean(y))/sd(y)
xobspr <- exp(xobsxb)/(1+exp(xobsxb))
x[runif(n)>xobspr] <- NA

ex_linquad <- data.frame(y,z,x,xsq=x^2,v)

devtools::use_data(ex_linquad, overwrite=TRUE)

#linear substantive model with interaction
n <- 1000
x1 <- rnorm(n)
x2 <- 1*(runif(n)<0.5)
y <- 1+x1+x2+x1*x2+rnorm(n)

#make some x1 and x2 values missing
x1[runif(n)>0.5] <- NA
x2[runif(n)>0.5] <- NA

ex_lininter <- data.frame(y,x1,x2)
ex_lininter$x2 <- factor(ex_lininter$x2)

devtools::use_data(ex_lininter, overwrite=TRUE)

#logistic substantive model with quadratic covariate effect
n <- 1000
x <- rnorm(n)
z <- rnorm(n)
#v is an auxiliary variable
v <- x+rnorm(n)
xb <- z+x+x^2
pr <- exp(xb)/(1+exp(xb))
y <- 1*(runif(n)<pr)

#make some x values missing
xobsxb <- (z-mean(z))/sd(z)
xobspr <- exp(xobsxb)/(1+exp(xobsxb))
x[runif(n)>xobspr] <- NA

ex_logisticquad <- data.frame(y,z,x,xsq=x^2,v)

devtools::use_data(ex_logisticquad, overwrite=TRUE)

#Cox substantive model with quadratic covariate effect
n <- 1000
x <- rnorm(n)
z <- rnorm(n)
#v is an auxiliary variable
v <- x+rnorm(n)
xb <- z+x+x^2
c <- -log(runif(n))/1
t <- -log(runif(n))/exp(xb)
d <- 1*(t<c)
t[d==0] <- c[d==0]

#make some x values missing
xobsxb <- (t-mean(t))/sd(t)
xobspr <- exp(xobsxb)/(1+exp(xobsxb))
x[runif(n)>xobspr] <- NA

ex_coxquad <- data.frame(t,d,z,x,xsq=x^2,v)

devtools::use_data(ex_coxquad, overwrite=TRUE)

#competing risks data
n <- 1000
x1 <- rbinom(n,1,0.5)
x2 <- x1+rnorm(n)
xb1 <- x1+x2
xb2 <- -x1
t1 <- -log(runif(n))/(0.002*exp(xb1))
t2 <- -log(runif(n))/(0.002*exp(xb2))
c <- -log(runif(n))/0.002
d <- rep(0,n)
d[(t1<c) & (t1<t2)] <- 1
d[(t2<c) & (t2<t1)] <- 2
t <- c
t[d==1] <- t1[d==1]
t[d==2] <- t2[d==2]
x1[runif(n)>0.5] <- NA
x2[runif(n)>0.5] <- NA

ex_compet <- data.frame(t,d,x1,x2)

devtools::use_data(ex_compet, overwrite=TRUE)

#Poisson regression
n <- 1000
x <- rnorm(n)
z <- 1*(runif(n)<0.5)
y <- rpois(n, lambda=exp(x-z))

#make some x1 and x2 values missing
x[runif(n)>(exp(y)/(1+exp(y)))] <- NA

ex_poisson <- data.frame(y,x,z)

devtools::use_data(ex_poisson, overwrite=TRUE)

#case cohort
n <- 10000
z <- rnorm(n)
x <- z+rnorm(n)
t <- -log(runif(n))/(0.01*exp(x+z))
d <- 1*(t<1)
t[d==0] <- 1
x[(runif(n)<0.5)] <- NA

fullcohortdata <- data.frame(t,d,x,z)
fullcohortdata$in.subco <- 0
#we sample a 10% subcohort
fullcohortdata$in.subco[sample(n, size=n*0.1)] <- 1
fullcohortdata$id <- 1:n

ex_cc <- fullcohortdata[(fullcohortdata$in.subco==1) | (fullcohortdata$d==1),]
ex_cc$entertime <- 0
ex_cc$entertime[ex_cc$in.subco==0] <- ex_cc$t[ex_cc$in.subco==0] - 0.000001

devtools::use_data(ex_cc, overwrite=TRUE)

#nested case control
n <- 10000
z <- rnorm(n)
x <- rbinom(n,1,exp(z)/(1+exp(z)))
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
ex_ncc=NULL

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
    ex_ncc <- rbind(ex_ncc,cbind(fullcohortdata[i,],numrisk=numAtRisk))
    ex_ncc <- rbind(ex_ncc,cbind(fullcohortdata[controls,], numrisk=numAtRisk))
    no.sample <- no.sample+1}
}

ex_ncc$setno <- rep(1:no.sample,each=m+1)
ex_ncc$case <- rep(c(1,rep(0,m)),no.sample)

devtools::use_data(ex_ncc, overwrite=TRUE)


# #covariate measurement error
# n <- 10000
# x <- rnorm(n)
# xb <- x
# pr <- exp(xb)/(1+exp(xb))
# y <- 1*(runif(n)<pr)
#
# w1 <- x+rnorm(n)
# w2 <- x+rnorm(n)
# w2[runif(1000)<0.9] <- NA
#
# ex_coverr <- data.frame(y,x=NA,w1,w2)
#
# devtools::use_data(ex_coverr, overwrite=TRUE)
