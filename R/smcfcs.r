updatePassiveVars <- function(data, method, passivecols) {
  for (i in passivecols) {
    data[,i] <- with(data, eval(parse(text=method[i])))
  }
  data
}

expit <- function(x) {
  exp(x)/(1+exp(x))
}

catdraw <- function(prob) {
  (1:length(prob))[rmultinom(1,size=1,prob=prob)==1]
}

modPostDraw <- function(modobj) {
  beta <- modobj$coef
  varcov <- vcov(modobj)
  beta + mvrnorm(1, mu=rep(0,ncol(varcov)), Sigma=varcov)
}

#' Substantive model compatible fully conditional specification imputation of covariates.
#'
#' Multiple imputes missing covariate values using substantive model compatible
#' fully conditional specification.
#'
#' smcfcs imputes missing values of covariates using the Substantive Model Compatible
#' Fully Conditional Specification multiple imputation approach proposed by
#' Bartlett et al 2014.
#'
#' Currently imputation is supported for linear regression ("lm"), logistic
#' regression ("logistic"), Cox regression for time to event
#' data ("coxph"), and Cox models for competing risks data. For the latter, a Cox
#' model is assumed for each cause of failure, and the event indicator should be integer
#' coded with 0 corresponding to censoring, 1 corresponding to failure from the first
#' cause etc.
#'
#' At present, in the case of linear or logistic substantive models the outcome  must
#' be fully observed, but this will be relaxed in a future version so that missing outcomes
#' are also imputed.
#'
#' The function returns a list of imputated datasets. Models (e.g. the substantive model)
#' can be fitted to each and results combined using Rubin's rules using the mitools
#' package, as illustrated in the examples.
#'
#' @param originaldata The original data frame with missing values.
#' @param smtype A string specifying the type of substantive model. Possible
#' values are "lm", "logistic", "coxph" and "compet".
#' @param smformula The formula of the substantive model. For "coxph" substantive
#' models the left hand side should be of the form "Surv(t,delta)". For "compet"
#' substantive models, a list should be passed consisting of the Cox models
#' for each cause of failure (see example).
#' @param method A required vector of strings specifying for each variable either
#' that it does not need to be imputed (""), the type of regression model to be
#' be used to impute. Possible values are "norm" (normal linear regression),
#' "logreg" (logistic regression), "poisson" (Poisson regression),
#' "podds" (proportional odds regression for ordered categorical variables),
#' "mlogit" (multinomial logistic regression for unordered categorical variables),
#' or a custom expression which defines a passively imputed variable, e.g.
#' "x^2" or "x1*x2".
#' @param predictorMatrix An optional predictor matrix. If specified, the matrix defines which
#' covariates will be used as predictors in the imputation models
#' (the outcome must not be included). The i'th row of the matrix should consist of
#' 0s and 1s, with a 1 in the j'th column indicating the j'th variable be used
#' as a covariate when imputing the i'th variable. If not specified, when
#' imputing a given variable, the imputation model covariates are the other
#' covariates of the substantive model which are partially observed
#' (but which are not passively imputed) and any fully observed variables (if present).
#' Note that the outcome variable is implicitly conditioned on by the rejection
#' sampling scheme used by smcfcs, and should not be specified as a predictor
#' in the predictor matrix.
#' @param m The number of imputed datasets to generate. The default is 5.
#' @param numit The number of iterations to run when generating each imputation.
#' In a (limited) range of simulations good performance was obtained with the
#' default of 10 iterations. However, particularly when the proportion of missingness
#' is large, more iterations may be required for convergence to stationarity.
#' @param rjlimit Specifies the maximum number of attempts which should be made
#' when using rejection sampling to draw from imputation models. If the limit is reached
#' when running a warning will be issued. In this case it is usually advisable to
#' increase the rjlimit.
#' @param noisy logical value (default FALSE) indicating whether output should be noisy, which can
#' be useful for debugging or checking that models being used are as desired.
#'
#' @return a list of data frames containing the multiply imputed datasets.
#'
#' @example data-raw/examples.r
#'
#' @references Bartlett JW, Seaman SR, White IR, Carpenter JR. Multiple imputation of covariates
#' by fully conditional specification: accommodating the substantive model. Statistical Methods
#' in Medical Research 2014. \url{http://doi.org/10.1177/0962280214521348}

#' @export
smcfcs <- function(originaldata,smtype,smformula,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE) {
  library("MASS")
  stopifnot(is.data.frame(originaldata))
  if (ncol(originaldata)!=length(method)) stop("Method argument must have the same length as the number of columns in the data frame.")

  n <- dim(originaldata)[1]
  #find column numbers of partially observed, fully observed variables, and outcome
  if (smtype=="coxph") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[,dCol]

    nullMod <- coxph(Surv(originaldata[,timeCol],originaldata[,dCol])~1)
    basehaz <- basehaz(nullMod)
    H0indices <- match(originaldata[,timeCol], basehaz[,2])
    rm(nullMod)
  }
  else if (smtype=="compet") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula[[1]])[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula[[1]])[[2]][[3]][[2]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[,dCol]
    numCauses <- length(smformula)
    H0 <- vector("list", numCauses)
    H0indices <- vector("list", numCauses)
    outcomeModBeta <- vector("list", numCauses)
    linpred <- vector("list", numCauses)
    for (cause in 1:numCauses) {
      nullMod <- coxph(as.formula(paste(strsplit(smformula[[cause]],"~")[[1]][1],"~1")), originaldata)
      basehaz <- basehaz(nullMod)
      H0[[cause]] <- basehaz[,1]
      H0indices[[cause]] <- match(originaldata[,timeCol], basehaz[,2])
      linpred[[cause]] <- as.formula(smformula[[cause]])
    }
    rm(nullMod)
  }
  else {
    outcomeCol <- which(colnames(originaldata)==as.formula(smformula)[[2]])
  }

  #create matrix of response indicators
  r <- 1*(is.na(originaldata)==0)

  if (length(c(which(method=="podds"),which(method=="mlogit")))>0) {
    library("VGAM")
  }

  if (smtype=="compet") {
    smcovnames <- attr(terms(as.formula(smformula[[1]])), "term.labels")
    for (cause in 2:numCauses) {
      smcovnames <- c(smcovnames, attr(terms(as.formula(smformula[[cause]])), "term.labels"))
    }
    smcovnames <- unique(smcovnames)
  }
  else {
    smcovnames <- attr(terms(as.formula(smformula)), "term.labels")
  }
  smcovcols <- (1:ncol(originaldata))[colnames(originaldata) %in% smcovnames]

  #partial vars are those variables for which an imputation method has been specified among the available regression types
  partialVars <- which((method=="norm") | (method=="logreg") | (method=="poisson") | (method=="podds") | (method=="mlogit"))

  #fully observed vars are those that are fully observed and are covariates in the substantive model
  fullObsVars <- which((colSums(r)==n) & (colnames(originaldata) %in% smcovnames))

  #passive variables
  passiveVars <- which((method!="") & (method!="norm") & (method!="logreg") & (method!="poisson") & (method!="podds") & (method!="mlogit"))

  print(paste("Outcome variable(s):", paste(colnames(originaldata)[outcomeCol],collapse=',')))
  print(paste("Passive variables:", paste(colnames(originaldata)[passiveVars],collapse=',')))
  print(paste("Partially obs. variables:", paste(colnames(originaldata)[partialVars],collapse=',')))
  print(paste("Fully obs. variables:", paste(colnames(originaldata)[fullObsVars],collapse=',')))

  imputations <- list()
  for (imp in 1:m) {
    imputations[[imp]] <- originaldata
  }

  for (imp in 1:m) {

    print(c("Imputation ",imp))

    #initial imputation of each partially observed variable based on observed values
    for (var in 1:length(partialVars)) {
      imputations[[imp]][r[,partialVars[var]]==0,partialVars[var]] <- sample(imputations[[imp]][r[,partialVars[var]]==1,partialVars[var]], size=sum(r[,partialVars[var]]==0), replace=TRUE)
    }

    for (cyclenum in 1:numit) {

      if (noisy==TRUE) {
        print(paste("Iteration ",cyclenum))
      }
      #update passive variable(s)
      imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)

      for (var in 1:length(partialVars)) {
        targetCol <- partialVars[var]
        if (is.null(predictorMatrix)) {
          predictorCols <- c(partialVars[! partialVars %in% targetCol], fullObsVars)
        }
        else {
          predictorCols <- which(predictorMatrix[targetCol,]==1)
          #ensure that user has not included outcome variable(s) here
          predictorCols <- predictorCols[! predictorCols %in% outcomeCol]
        }
        if ((imp==1) & (cyclenum==1)) {
          print(paste("Imputing: ",colnames(imputations[[imp]])[targetCol]," using ",colnames(imputations[[imp]])[predictorCols],collapse=','))
        }

        xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], "~", paste(colnames(imputations[[imp]])[predictorCols], collapse="+"),sep=""))
        #print(xmodformula)
        if (method[targetCol]=="norm") {
          #estimate parameters of covariate model
          xmod <- lm(xmodformula, data=imputations[[imp]])
          #take draw from posterior of covariate model parameters
          beta <- xmod$coef
          sigmasq <- summary(xmod)$sigma^2
          newsigmasq <- (sigmasq*xmod$df) / rchisq(1,xmod$df)
          covariance <- (newsigmasq/sigmasq)*vcov(xmod)
          #newbeta = beta + chol(covariance) %*% rnorm(length(beta))
          newbeta = beta + mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
          #calculate fitted values
          xfitted <- model.matrix(xmod) %*% newbeta
        }
        else if (method[targetCol]=="logreg") {
          xmod <- glm(xmodformula, family="binomial",data=imputations[[imp]])
          newbeta = modPostDraw(xmod)
          xfitted <- expit(model.matrix(xmod) %*% newbeta)
        }
        else if (method[targetCol]=="poisson") {
          xmod <- glm(xmodformula, family="poisson", data=imputations[[imp]])
          newbeta = modPostDraw(xmod)
          xfitted <- exp(model.matrix(xmod) %*% newbeta)
        }
        else if (method[targetCol]=="podds") {
          xmod <- vglm(xmodformula, propodds, data=imputations[[imp]])
          newbeta <- coefficients(xmod) + mvrnorm(1, mu=rep(0,ncol(vcov(xmod))), Sigma=vcov(xmod))
          linpreds <- matrix(model.matrix(xmod) %*% newbeta, byrow=TRUE, ncol=(nlevels(imputations[[imp]][,targetCol])-1))
          cumprobs <- cbind(1/(1+exp(linpreds)), rep(1,nrow(linpreds)))
          xfitted <- cbind(cumprobs[,1] ,cumprobs[,2:ncol(cumprobs)] - cumprobs[,1:(ncol(cumprobs)-1)])
        }
        else if (method[targetCol]=="mlogit") {
          xmod <- vglm(xmodformula, multinomial(refLevel=1), data=imputations[[imp]])
          newbeta <- coef(xmod) + mvrnorm(1, mu=rep(0,ncol(vcov(xmod))), Sigma=vcov(xmod))
          linpreds <- matrix(model.matrix(xmod) %*% newbeta, byrow=TRUE, ncol=(nlevels(imputations[[imp]][,targetCol])-1))
          denom <- 1+rowSums(exp(linpreds))
          xfitted <-cbind(1/denom, exp(linpreds) / denom)
        }
        if (noisy==TRUE) {
          print(summary(xmod))
        }

        #estimate parameters of substantive model
        if (smtype=="lm") {
          ymod <- lm(as.formula(smformula),imputations[[imp]])
          beta <- ymod$coef
          sigmasq <- summary(ymod)$sigma^2
          varcov <- vcov(ymod)
          outcomeModResVar <- (sigmasq*ymod$df) / rchisq(1,ymod$df)
          outcomeModBeta = beta + chol((outcomeModResVar/sigmasq)*varcov) %*% rnorm(length(beta))
        }
        else if (smtype=="logistic") {
          ymod <- glm(as.formula(smformula),family="binomial",imputations[[imp]])
          outcomeModBeta = modPostDraw(ymod)
        }
        else if (smtype=="coxph") {
          ymod <- coxph(as.formula(smformula), imputations[[imp]])
          outcomeModBeta <- modPostDraw(ymod)
          ymod$coefficients <- outcomeModBeta
          basehaz <- basehaz(ymod, centered=FALSE)[,1]
          H0 <- basehaz[H0indices]
        }
        else if (smtype=="compet") {
          for (cause in 1:numCauses) {
            ymod <- coxph(as.formula(smformula[[cause]]), imputations[[imp]])
            outcomeModBeta[[cause]] <- modPostDraw(ymod)
            ymod$coefficients <- outcomeModBeta[[cause]]
            basehaz <- basehaz(ymod, centered=FALSE)[,1]
            H0[[cause]] <- basehaz[H0indices[[cause]]]
          }
        }
        if (noisy==TRUE) {
          print(summary(ymod))
        }

        #impute x, either directly where possibly, or using rejection sampling otherwise
        imputationNeeded <- (1:n)[r[,targetCol]==0]

        if ((method[targetCol]=="logreg") | (method[targetCol]=="podds") | (method[targetCol]=="mlogit")) {
          #directly sample
          if (method[targetCol]=="logreg") {
            numberOutcomes <- 2
            fittedMean <- cbind(1-xfitted, xfitted)
          }
          else {
            numberOutcomes <- nlevels(imputations[[imp]][,targetCol])
            fittedMean <- xfitted
          }

          outcomeDensCovDens = array(dim=c(length(imputationNeeded),numberOutcomes),0)

          for (xMisVal in 1:numberOutcomes) {
            if (method[targetCol]=="logreg") {
              valToImpute <- xMisVal-1
            }
            else {
              valToImpute <- levels(imputations[[imp]][,targetCol])[xMisVal]
            }
            imputations[[imp]][imputationNeeded,targetCol] <- valToImpute

            #update passive variables
            imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)

            if (smtype=="lm") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              deviation <- imputations[[imp]][imputationNeeded,outcomeCol] - outmodxb[imputationNeeded]
              outcomeDens <- dnorm(deviation, mean=0, sd=1)
              outcomeDensCovDens[,xMisVal] <- outcomeDens * fittedMean[imputationNeeded,xMisVal]
            }
            else if (smtype=="logistic") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              prob <- expit(outmodxb[imputationNeeded])
              outcomeDens <- prob*imputations[[imp]][imputationNeeded,outcomeCol] + (1-prob)*(1-imputations[[imp]][imputationNeeded,outcomeCol])
            }
            else if (smtype=="coxph") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]])
              outmodxb <- outmodxb[,2:dim(outmodxb)[2]] %*% outcomeModBeta
              outcomeDens <- exp(-H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))* (exp(outmodxb[imputationNeeded])^d[imputationNeeded])
            }
            else if (smtype=="compet") {
              outcomeDens <- rep(1,length(imputationNeeded))
              for (cause in 1:numCauses) {
                outmodxb <-  model.matrix(linpred[[cause]],imputations[[imp]])
                outmodxb <- outmodxb[,2:dim(outmodxb)[2]] %*% outcomeModBeta[[cause]]
                outcomeDens <- outcomeDens * exp(-H0[[cause]][imputationNeeded] * exp(outmodxb[imputationNeeded]))* (exp(outmodxb[imputationNeeded])^(d[imputationNeeded]==cause))
              }
            }
            outcomeDensCovDens[,xMisVal] <- outcomeDens * fittedMean[imputationNeeded,xMisVal]
          }
          directImpProbs = outcomeDensCovDens / rowSums(outcomeDensCovDens)

          if (method[targetCol]=="logreg") {
            directImpProbs = directImpProbs[,2]
            imputations[[imp]][imputationNeeded,targetCol] <- rbinom(length(imputationNeeded),1,directImpProbs)
          }
          else {
            imputations[[imp]][imputationNeeded,targetCol] <- levels(imputations[[imp]][,targetCol])[apply(directImpProbs, 1, catdraw)]
          }

          #update passive variables
          imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
        }
        else {
          #use rejection sampling
          #first draw for all subjects who need imputing, using a small number of attempts
          firstTryLimit <- 25
          j <- 1

          while ((length(imputationNeeded)>0) & (j<firstTryLimit)) {
            #sample from covariate model
            if (method[targetCol]=="norm") {
              imputations[[imp]][imputationNeeded,targetCol] <- rnorm(length(imputationNeeded),xfitted[imputationNeeded],newsigmasq^0.5)
            }
            else if (method[targetCol]=="logreg") {
              imputations[[imp]][imputationNeeded,targetCol] <- rbinom(length(imputationNeeded),size=1,prob=xfitted[imputationNeeded])
            }
            else if (method[targetCol]=="poisson") {
              imputations[[imp]][imputationNeeded,targetCol] <- rpois(length(imputationNeeded),xfitted[imputationNeeded])
            }

            #update passive variables
            imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)

            #accept/reject
            uDraw <- runif(length(imputationNeeded))
            if (smtype=="lm") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              deviation <- imputations[[imp]][imputationNeeded,outcomeCol] - outmodxb[imputationNeeded]
              reject = 1*(log(uDraw) > -(deviation^2) / (2*array(outcomeModResVar,dim=c(length(imputationNeeded),1))))
            }
            else if (smtype=="logistic") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              prob = expit(outmodxb[imputationNeeded])
              prob = prob*imputations[[imp]][imputationNeeded,outcomeCol] + (1-prob)*(1-imputations[[imp]][imputationNeeded,outcomeCol])
              reject = 1*(uDraw>prob)
            }
            else if (smtype=="coxph") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]])
              outmodxb <- outmodxb[,2:dim(outmodxb)[2]] %*% outcomeModBeta
              s_t = exp(-H0[imputationNeeded]* exp(outmodxb[imputationNeeded]))
              prob = exp(1 + outmodxb[imputationNeeded] - (H0[imputationNeeded]* exp(outmodxb[imputationNeeded])) ) * H0[imputationNeeded]
              prob = d[imputationNeeded]*prob + (1-d[imputationNeeded])*s_t
              reject = 1*(uDraw > prob )
            }
            else if (smtype=="compet") {
              prob <- rep(1,length(imputationNeeded))
              for (cause in 1:numCauses) {
                outmodxb <-  model.matrix(linpred[[cause]],imputations[[imp]])
                outmodxb <- outmodxb[,2:dim(outmodxb)[2]] %*% outcomeModBeta[[cause]]
                prob = prob * exp(-H0[[cause]][imputationNeeded] * exp(outmodxb[imputationNeeded]))* (H0[[cause]][imputationNeeded]*exp(1+outmodxb[imputationNeeded]))^(d[imputationNeeded]==cause)
              }
              reject = 1*(uDraw > prob )
            }
            imputationNeeded <- imputationNeeded[reject==1]

            j <- j+1
          }

          #now, for those remaining, who must have low acceptance probabilities, sample by subject
          for (i in imputationNeeded) {
            impFound <- 0
            j <- 0
            while ((impFound<1) & (j<rjlimit)) {
              tempData <- imputations[[imp]][i,]
              tempData <- tempData[rep(1,rjlimit),]
              if (method[targetCol]=="norm") {
                tempData[,targetCol] <- rnorm(rjlimit,xfitted[i],newsigmasq^0.5)
              }
              else if (method[targetCol]=="logreg") {
                tempData[,targetCol] <- rbinom(rjlimit,size=1,xfitted[i])
              }
              else if (method[targetCol]=="poisson") {
                tempData[,targetCol] <- rpois(rjlimit,xfitted[i])
              }

              #passively impute
              tempData <- updatePassiveVars(tempData, method, passiveVars)

              #accept reject
              uDraw <- runif(rjlimit)
              if (smtype=="lm") {
                outmodxb <-  model.matrix(as.formula(smformula),tempData) %*% outcomeModBeta
                deviation <- tempData[,outcomeCol] - outmodxb
                reject = 1*(log(uDraw) > -(deviation^2) / (2*array(outcomeModResVar,dim=c(rjlimit,1))))
              }
              else if (smtype=="logistic") {
                outmodxb <-  model.matrix(as.formula(smformula),tempData) %*% outcomeModBeta
                prob = expit(outmodxb)
                prob = prob*tempData[,outcomeCol] + (1-prob)*(1-tempData[,outcomeCol])
                reject = 1*(uDraw>prob)
              }
              else if (smtype=="coxph") {
                outmodxb <-  model.matrix(as.formula(smformula),tempData)
                outmodxb <- outmodxb[,2:dim(outmodxb)[2]] %*% outcomeModBeta
                s_t = exp(-H0[i]* exp(outmodxb))
                prob = exp(1 + outmodxb - (H0[i]* exp(outmodxb)) ) * H0[i]
                prob = d[i]*prob + (1-d[i])*s_t
                reject = 1*(uDraw > prob )
              }
              else if (smtype=="compet") {
                prob <- rep(1,rjlimit)
                for (cause in 1:numCauses) {
                  outmodxb <-  model.matrix(linpred[[cause]],tempData)
                  outmodxb <- outmodxb[,2:dim(outmodxb)[2]] %*% outcomeModBeta[[cause]]
                  prob = prob * exp(-H0[[cause]][i] * exp(outmodxb))* (H0[[cause]][i]*exp(1+outmodxb))^(d[i]==cause)
                }
                reject = 1*(uDraw > prob )
              }

              if (sum(reject)<rjlimit) {
                impFound <- 1
                imputations[[imp]][i,targetCol] <- tempData[reject==0,targetCol][1]
              }

              j <- j+1
            }
            if (j==rjlimit) {
              print("Rejection sampling has failed for one record. You may want to increase the rejecton sampling limit.")
            }
          }
          #update passive variables
          imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
        }
      }


    }

  }

  imputations

}

