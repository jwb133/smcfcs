#' Substantive model compatible fully conditional specification imputation of covariates.
#'
#' Multiply imputes missing covariate values using substantive model compatible
#' fully conditional specification.
#'
#' smcfcs imputes missing values of covariates using the Substantive Model Compatible
#' Fully Conditional Specification multiple imputation approach proposed by
#' Bartlett \emph{et al} 2015 (see references).
#'
#' Currently imputation is supported for linear regression (\code{"lm"}),
#' logistic regression (\code{"logistic"}), Poisson regression
#' (\code{"poisson"}), Cox regression for time to event data (\code{"coxph"}),
#' and Cox models for competing risks data (\code{"compet"}). For the latter, a
#' Cox model is assumed for each cause of failure, and the event indicator
#' should be integer coded with 0 corresponding to censoring, 1 corresponding to
#' failure from the first cause etc.
#'
#' The function returns a list. The first element \code{impDataset} of the list is a list of the imputed
#' datasets. Models (e.g. the substantive model) can be fitted to each and results
#' combined using Rubin's rules using the mitools package, as illustrated in the examples.
#'
#' The second element \code{smCoefIter} is a three dimensional array containing the values
#' of the substantive model parameters obtained at the end of each iteration of the algorithm.
#' The array is indexed by: imputation number, parameter number, iteration.
#'
#' If the substantive model is linear, logistic or Poisson regression,
#' \code{smcfcs} will automatically impute missing outcomes, if present, using
#' the specified substantive model. However, even in this case, the user should
#' specify "" in the element of method corresponding to the outcome variable.
#'
#'
#' The development of this package was supported by a UK Medical Research Council
#' Fellowship (MR/K02180X/1). Part of its development took place while the author was
#' kindly hosted by the University of Michigan's Department of Biostatistics & Institute for
#' Social Research.
#'
#' The structure of many of the arguments to \code{smcfcs} are based on those of
#' the excellent \code{mice} package.
#'
#' @param originaldata The original data frame with missing values.
#' @param smtype A string specifying the type of substantive model. Possible
#' values are \code{"lm"}, \code{"logistic"}, \code{"poisson"}, \code{"coxph"}
#'  and \code{"compet"}.
#' @param smformula The formula of the substantive model. For \code{"coxph"} substantive
#' models the left hand side should be of the form \code{"Surv(t,d)"}. For \code{"compet"}
#' substantive models, a list should be passed consisting of the Cox models
#' for each cause of failure (see example).
#' @param method A required vector of strings specifying for each variable either
#' that it does not need to be imputed (""), the type of regression model to be
#' be used to impute. Possible values are \code{"norm"} (normal linear regression),
#' \code{"logreg"} (logistic regression), \code{"poisson"} (Poisson regression),
#' \code{"podds"} (proportional odds regression for ordered categorical variables),
#' \code{"mlogit"} (multinomial logistic regression for unordered categorical variables),
#' or a custom expression which defines a passively imputed variable, e.g.
#' \code{"x^2"} or \code{"x1*x2"}. \code{"latnorm"} indicates the variable is a latent
#' normal variable which is measured with error. If this is specified for a variable,
#' the \code{"errorProneMatrix"} argument should also be used.
#' @param predictorMatrix An optional predictor matrix. If specified, the matrix defines which
#' covariates will be used as predictors in the imputation models
#' (the outcome must not be included). The i'th row of the matrix should consist of
#' 0s and 1s, with a 1 in the j'th column indicating the j'th variable be used
#' as a covariate when imputing the i'th variable. If not specified, when
#' imputing a given variable, the imputation model covariates are the other
#' covariates of the substantive model which are partially observed
#' (but which are not passively imputed) and any fully observed covariates (if present)
#' in the substantive model. Note that the outcome variable is implicitly conditioned
#' on by the rejection sampling scheme used by smcfcs, and should not be specified as a predictor
#' in the predictor matrix.
#' @param m The number of imputed datasets to generate. The default is 5.
#' @param numit The number of iterations to run when generating each imputation.
#' In a (limited) range of simulations good performance was obtained with the
#' default of 10 iterations. However, particularly when the proportion of missingness
#' is large, more iterations may be required for convergence to stationarity.
#' @param rjlimit Specifies the maximum number of attempts which should be made
#' when using rejection sampling to draw from imputation models. If the limit is reached
#' when running a warning will be issued. In this case it is probably advisable to
#' increase the \code{rjlimit} until the warning does not appear.
#' @param noisy logical value (default FALSE) indicating whether output should be noisy, which can
#' be useful for debugging or checking that models being used are as desired.
#' @param errorProneMatrix An optional matrix which if specified indicates that some variables
#' are measured with classical measurement error. If the i'th variable is measured with error
#' by variables j and k, then the (i,j) and (i,k) entries of this matrix should be 1, with the
#' remainder of entries 0. The i'th element of the method argument should then be specified
#' as \code{"latnorm"}. See the measurement error vignette for more details.
#'
#' @return A list containing:
#'
#' \code{impDatasets} a list containing the imputed datasets
#'
#' \code{smCoefIter} a three dimension matrix containing the substantive model parameter
#' values. The matrix is indexed by [imputation,parameter number,iteration]
#'
#' @author Jonathan Bartlett \email{j.w.bartlett@@bath.ac.uk} \url{http://www.missingdata.org.uk}
#' \url{http://thestatsgeek.com}
#'
#' @example data-raw/examples.r
#'
#' @references Bartlett JW, Seaman SR, White IR, Carpenter JR. Multiple imputation of covariates
#' by fully conditional specification: accommodating the substantive model. Statistical Methods
#' in Medical Research 2015; 24(4): 462-487. \url{http://doi.org/10.1177/0962280214521348}

#' @import stats

#' @export
smcfcs <- function(originaldata,smtype,smformula,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,errorProneMatrix=NULL) {
  #call core  smcfcs function, passing through arguments
  smcfcs.core(originaldata,smtype,smformula,method,predictorMatrix,m,numit,rjlimit,noisy,errorProneMatrix=errorProneMatrix)
}

#' Substantive model compatible fully conditional specification imputation of covariates for case cohort studies
#'
#' Multiply imputes missing covariate values using substantive model compatible
#' fully conditional specification for case cohort studies.
#'
#' This version of \code{smcfcs} is designed for use with case cohort studies but where the analyst does not wish to,
#' or cannot (due to not having the necessary data) impute the full cohort. The function's arguments are the same
#' as for the main smcfcs function, except for \code{smformula}, \code{in.subco}, and \code{sampfrac} - see above
#' for details on how these should be specified.
#'
#' @author Ruth Keogh \email{ruth.keogh@@lshtm.ac.uk}
#' @author Jonathan Bartlett \email{j.w.bartlett@@bath.ac.uk}
#'
#' @param originaldata The case-cohort data set (NOT a full cohort data set with a case-cohort substudy within it)
#' @param smformula A formula of the form "Surv(entertime,t,d)~x", where d is the event (d=1) or censoring (d=0) indicator, t is the event or censoring time and entertime is equal to the time origin (typically 0) for individuals in the subcohort and is equal to (t-0.001) for cases outside the subcohort [this sets cases outside the subcohort to enter follow-up just before their event time. The value 0.001 may need to be modified depending on the time scale.]
#' @param in.subco The name of a column in the dataset with 0/1s that indicates whether the subject is in the subcohort
#' @param sampfrac The proportion of individuals from the underlying full cohort who are in the subcohort
#'
#' @inheritParams smcfcs
#'
#' @example data-raw/cc_example.r
#'
#' @export
smcfcs.casecohort <- function(originaldata,smformula,sampfrac,in.subco,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,
                              errorProneMatrix=NULL) {
  smcfcs.core(originaldata,smtype="casecohort",smformula,method,predictorMatrix,m,numit,rjlimit,noisy,sampfrac=sampfrac,in.subco=in.subco,
              errorProneMatrix=errorProneMatrix)
}

#' Substantive model compatible fully conditional specification imputation of covariates for nested case control
#' studies
#'
#' Multiply imputes missing covariate values using substantive model compatible
#' fully conditional specification for nested case control studies.
#'
#' This version of \code{smcfcs} is designed for use with nested case control studies. The function's arguments are the same
#' as for the main smcfcs function, except for \code{smformula}, \code{set}, \code{event} and \code{nrisk} - see above
#' for details on how these should be specified.
#'
#' @author Ruth Keogh \email{ruth.keogh@@lshtm.ac.uk}
#' @author Jonathan Bartlett \email{j.w.bartlett@@bath.ac.uk}
#'
#' @param originaldata The nested case-control data set (NOT a full cohort data set with a case-cohort substudy within it)
#' @param smformula A formula of the form "Surv(t,case)~x+strata(set)", where case is case-control indicator, t is the event or censoring time. Note that t could be set to the case's event time for the matched controls in a given set. The right hand side should include the case control set as a strata term (see example).
#' @param set variable identifying matched sets in nested case-control study
#' @param event variable which indicates who is a case/control in the nested case-control sample. Note that this is distinct from d.
#' @param nrisk variable which is the number at risk (in the underlying full cohort) at the event time for the case in each matched set (i.e. nrisk is the same for all individuals in a matched set).
#'
#' @inheritParams smcfcs
#' @example data-raw/ncc_example.r
#' @export
smcfcs.nestedcc <- function(originaldata,smformula,set,event,nrisk,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,errorProneMatrix=NULL) {
  smcfcs.core(originaldata,smtype="nestedcc",smformula,method,predictorMatrix,m,numit,rjlimit,noisy,set=set,event=event,nrisk=nrisk,
              errorProneMatrix=errorProneMatrix)
}

#this is the core of the smcfcs function, called by wrapper functions for certain different substantive models
smcfcs.core <- function(originaldata,smtype,smformula,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,errorProneMatrix=NULL,
                        ...) {
  #get extra arguments passed in ...
  extraArgs <- list(...)

  stopifnot(is.data.frame(originaldata))
  if (ncol(originaldata)!=length(method)) stop("Method argument must have the same length as the number of columns in the data frame.")

  n <- dim(originaldata)[1]

  #create matrix of response indicators
  r <- 1*(is.na(originaldata)==0)

  if ((smtype %in% c("lm", "logistic", "poisson", "coxph", "compet", "casecohort","nestedcc"))==FALSE) #RUTH 21/03/2017: ADDED 'nestedcc' OPTION HERE
    stop(paste("Substantive model type ",smtype," not recognised.",sep=""))

  #find column numbers of partially observed, fully observed variables, and outcome
  if (smtype=="coxph") {

    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[,dCol]

    nullMod <- survival::coxph(Surv(originaldata[,timeCol],originaldata[,dCol])~1)
    basehaz <- survival::basehaz(nullMod)
    H0indices <- match(originaldata[,timeCol], basehaz[,2])
    rm(nullMod)
  } else if (smtype=="compet") {

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
      nullMod <- survival::coxph(as.formula(paste(strsplit(smformula[[cause]],"~")[[1]][1],"~1")), originaldata)
      basehaz <- survival::basehaz(nullMod)
      H0[[cause]] <- basehaz[,1]
      H0indices[[cause]] <- match(originaldata[,timeCol], basehaz[,2])
      linpred[[cause]] <- as.formula(smformula[[cause]])
    }
    rm(nullMod)
  } else if (smtype=="casecohort") {

    subcoCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% extraArgs$in.subco]
    #subcoMembers is a vector of row numbers of those in the subcohort
    subcoMembers <- which(originaldata[,subcoCol]==1)

    #generate weights for use in later analysis which we use to obtain baseline cumulative hazard
    #assign a weight of /samp.frac to individuals in the subcohort and 0 to those outside the subcohort
    subco.weight<-ifelse(originaldata[,subcoCol]==1,1/extraArgs$sampfrac,0)


    entertimeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[4]])]
    outcomeCol <- c(entertimeCol,timeCol, dCol)
    d <- originaldata[,dCol]

    #list of unique event times - used in calculation of baseline cumulative hazard
    list.times=sort(unique(originaldata[,timeCol][originaldata[,dCol]==1]))  #RUTH 21/03/17: ADDED THIS LINE

    smcfcsid <- 1:n
    smformula2<-paste(smformula,"+cluster(smcfcsid)",sep="")
  }
  else if (smtype=="nestedcc") {

    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    setCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(extraArgs$set)]
    nriskCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(extraArgs$nrisk)]
    eventCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(extraArgs$event)] #this is distinct from the dCol

    #the below command creates "Surv(t,case)~x" from "Surv(t,case)~x+strata(setno)" (for example) (i.e. it removes the strata part of the formula)
    #this is used when obtaining outmodxb
    #note this is done slightly oddly, but this is because as.formula does not work well for long formulas as it splits across lines
    exp1=as.formula(paste(smformula))[[2]]
    exp2=as.formula(smformula)[[3]][[2]]
    smformula2<-paste(deparse(exp1),"~",deparse(exp2,width.cutoff = 500L))

    #This is the indicator of whether an individual ever has the event (regardless of whether they are sometimes used as a control and sometimes (one) as a case)
    d <- originaldata[,eventCol]

    #noncases is a vector of row numbers of those who never have the event (which is a subset of the controls)
    noncases <- which(originaldata[,eventCol]==0)

    #number of individuals in each sampled risk set (matched set)
   num.sampriskset<-ave(rep(1,dim(originaldata)[1]), originaldata[,setCol], FUN = function(x) sum(x))

  }
  else {
    outcomeCol <- which(colnames(originaldata)==as.formula(smformula)[[2]])
  }

  if (smtype=="logistic") {
    if (is.numeric(originaldata[,outcomeCol])==FALSE) {
      stop("For logistic substantive models the outcome variable must be numeric 0/1.")
    } else {
      if (all.equal(unique(originaldata[,outcomeCol]),c(0,1))==FALSE) {
        stop("For logistic substantive models the outcome variable must be coded 0/1.")
      }
    }
  }

  if (smtype=="compet") {
    smcovnames <- attr(terms(as.formula(smformula[[1]])), "term.labels")
    for (cause in 2:numCauses) {
      smcovnames <- c(smcovnames, attr(terms(as.formula(smformula[[cause]])), "term.labels"))
    }
    smcovnames <- unique(smcovnames)
  } else  if (smtype=="nestedcc") {
    smcovnames <- attr(terms(as.formula(smformula)), "term.labels")[-length(attr(terms(as.formula(smformula)), "term.labels"))]
  }
  else {
    smcovnames <- attr(terms(as.formula(smformula)), "term.labels")
  }
  smcovcols <- (1:ncol(originaldata))[colnames(originaldata) %in% smcovnames]

  #partial vars are those variables for which an imputation method has been specified among the available regression types
  partialVars <- which((method=="norm") | (method=="latnorm") | (method=="logreg") | (method=="poisson") | (method=="podds") | (method=="mlogit"))

  if (length(partialVars)==0) stop("You have not specified any valid imputation methods in the method argument.")

  #check that methods are given for each partially observed column, and not given for fully observed columns
  for (colnum in 1:ncol(originaldata)) {
    if (method[colnum]!="") {
      #an imputation method has been specified
      if (colnum %in% outcomeCol) {
        stop(paste("An imputation method has been specified for ",colnames(originaldata)[colnum],
                   ". Elements of the method argument corresponding to the outcome variable(s) should be empty.",sep=""))
      }
      else {
        if (sum(r[,colnum])==n) {
          stop(paste("An imputation method has been specified for ",colnames(originaldata)[colnum],
                     ", but it appears to be fully observed.",sep=""))
        }
      }
    }
    else {
      #no imputation method has been specified
      if (sum(r[,colnum])<n) {
        #some values are missing
        if (((colnum %in% outcomeCol)==FALSE) & (sum(errorProneMatrix[,colnum])==0)) {
          stop(paste("Variable ",colnames(originaldata)[colnum], " does not have an imputation method specified, yet appears to have missing values.",sep=""))
        }
        }
      }

    }

  #fully observed vars are those that are fully observed and are covariates in the substantive model
  fullObsVars <- which((colSums(r)==n) & (colnames(originaldata) %in% smcovnames))

  #passive variables
  passiveVars <- which((method!="") & (method!="norm") & (method!="logreg") & (method!="poisson") & (method!="podds") & (method!="mlogit") & (method!="latnorm"))

  print(paste("Outcome variable(s):", paste(colnames(originaldata)[outcomeCol],collapse=',')))
  print(paste("Passive variables:", paste(colnames(originaldata)[passiveVars],collapse=',')))
  print(paste("Partially obs. variables:", paste(colnames(originaldata)[partialVars],collapse=',')))
  print(paste("Fully obs. substantive model variables:", paste(colnames(originaldata)[fullObsVars],collapse=',')))

  imputations <- list()
  for (imp in 1:m) {
    imputations[[imp]] <- originaldata
  }

  rjFailCount <- 0

  for (imp in 1:m) {

    print(paste("Imputation ",imp))

    #initial imputation of each partially observed variable based on observed values
    for (var in 1:length(partialVars)) {
      targetCol <- partialVars[var]
      if (method[targetCol]=="latnorm") {
        #first impute any missing replicate error-prone measurements of this variable by a randomly chosen observed value
        errorProneCols <- which(errorProneMatrix[targetCol,]==1)
        for (measure in 1:length(errorProneCols)) {
          if (sum(r[,errorProneCols[measure]])<n) {
            imputations[[imp]][r[,errorProneCols[measure]]==0,errorProneCols[measure]] <- sample(imputations[[imp]][r[,errorProneCols[measure]]==1,
                errorProneCols[measure]], size=sum(r[,errorProneCols[measure]]==0), replace=TRUE)
          }
        }

        #initialize latent predictors with mean of their error-prone measurements
        imputations[[imp]][,targetCol] <- apply(imputations[[imp]][,errorProneCols], 1, mean)
      }
      else {
        imputations[[imp]][r[,targetCol]==0,targetCol] <- sample(imputations[[imp]][r[,targetCol]==1,targetCol], size=sum(r[,targetCol]==0), replace=TRUE)
      }
    }

    #initial imputations of missing outcomes, if present (using improper imputation)
    if ((smtype=="lm") | (smtype=="logistic") | (smtype=="poisson")) {
      if (sum(r[,outcomeCol])<n) {
        if (imp==1) {
          print("Imputing missing outcomes using specified substantive model.")
        }
        #update passive variable(s)
        imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)

        imputationNeeded <- (1:n)[r[,outcomeCol]==0]
        #estimate parameters of substantive model
        if (smtype=="lm") {
          ymod <- stats::lm(as.formula(smformula),imputations[[imp]])
          beta <- ymod$coef
          sigmasq <- summary(ymod)$sigma^2
          #fill out missing values so that model.matrix works for all rows
          imputations[[imp]][imputationNeeded,outcomeCol] <- 0
          outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% beta
          imputations[[imp]][imputationNeeded,outcomeCol] <- rnorm(length(imputationNeeded),outmodxb[imputationNeeded], sigmasq^0.5)
        }
        else if (smtype=="logistic") {
          ymod <- glm(as.formula(smformula),family="binomial",imputations[[imp]])
          beta <- ymod$coef
          imputations[[imp]][imputationNeeded,outcomeCol] <- 0
          outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% beta
          prob <- expit(outmodxb[imputationNeeded])
          imputations[[imp]][imputationNeeded,outcomeCol] <- rbinom(length(imputationNeeded),1,prob)
        }
        else if (smtype=="poisson") {
          ymod <- glm(as.formula(smformula),family="poisson",imputations[[imp]])
          beta <- ymod$coef
          imputations[[imp]][imputationNeeded,outcomeCol] <- 0
          outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% beta
          imputations[[imp]][imputationNeeded,outcomeCol] <- rpois(length(imputationNeeded),exp(outmodxb[imputationNeeded]))
        }
      }
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
          if (method[targetCol]=="latnorm") {
            print(paste("Imputing: ",colnames(imputations[[imp]])[targetCol]," using ",paste(colnames(imputations[[imp]])[c(predictorCols,which(errorProneMatrix[targetCol,]==1))],collapse=',')," plus outcome",collapse=','))
          }
          else {
            print(paste("Imputing: ",colnames(imputations[[imp]])[targetCol]," using ",paste(colnames(imputations[[imp]])[predictorCols],collapse=',')," plus outcome",collapse=','))
          }
        }
        if (length(predictorCols)>0) {
          xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], "~", paste(colnames(imputations[[imp]])[predictorCols], collapse="+"),sep=""))
        }
        else {
          xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], "~1",sep=""))
        }
        if (smtype=="casecohort") {
          xmoddata <- imputations[[imp]][subcoMembers,]
        } else if (smtype=="nestedcc"){
          xmoddata <- imputations[[imp]][noncases,]
        } else {
          xmoddata <- imputations[[imp]]
        }
        if (method[targetCol]=="norm") {
          #estimate parameters of covariate model
          xmod <- lm(xmodformula, data=xmoddata)
          #take draw from posterior of covariate model parameters
          beta <- xmod$coef
          sigmasq <- summary(xmod)$sigma^2
          newsigmasq <- (sigmasq*xmod$df) / rchisq(1,xmod$df)
          covariance <- (newsigmasq/sigmasq)*vcov(xmod)
          newbeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
          #calculate fitted values
          if ((smtype=="casecohort")|(smtype=="nestedcc")) {
            xfitted <- model.matrix(xmodformula, data=imputations[[imp]]) %*% newbeta
          } else {
            xfitted <- model.matrix(xmod) %*% newbeta
          }
        } else if  (method[targetCol]=="latnorm") {
            #estimate parameters of covariate model
            xmod <- lm(xmodformula, data=xmoddata)
            #take draw from posterior of covariate model parameters
            beta <- xmod$coef
            sigmasq <- summary(xmod)$sigma^2
            #newsigmasq <- (sigmasq*xmod$df) / rchisq(1,xmod$df)
            newsigmasq <- 1/rgamma(1,shape=((n+1)/2), rate=((n*sigmasq+1)/2))
            covariance <- (newsigmasq/sigmasq)*vcov(xmod)
            newbeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
            #calculate fitted values
            if ((smtype=="casecohort")|(smtype=="nestedcc")) {
              xfitted <- model.matrix(xmodformula, data=imputations[[imp]]) %*% newbeta
            } else {
              xfitted <- model.matrix(xmod) %*% newbeta
            }

            #estimate error variance and draw new value of error variance
            errorProneCols <- which(errorProneMatrix[targetCol,]==1)
            xmat <- matrix(imputations[[imp]][,targetCol], nrow=nrow(imputations[[imp]]), ncol=length(errorProneCols))
            uVec <- c(as.matrix(imputations[[imp]][,errorProneCols] - xmat))
            sigmausq <- mean(uVec^2)
            #take draw from posterior of error variance
            sum_ni <- n*length(errorProneCols)
            sigmausq <- 1/rgamma(1,shape=((sum_ni+1)/2), rate=((sum_ni*sigmausq+1)/2))

            #re-impute any originally missing error-prone measurements, based on classical error model assumption
            for (measure in 1:length(errorProneCols)) {
              nToImpute <- n-sum(r[,errorProneCols[measure]])
              if (nToImpute>0) {
                #then some values need imputing
                imputations[[imp]][r[,errorProneCols[measure]]==0,errorProneCols[measure]] <- imputations[[imp]][r[,errorProneCols[measure]]==0,targetCol] +
                  rnorm(nToImpute, 0, sd=sqrt(sigmausq))
              }
            }

            #calculate conditional mean and variance of X|everything else except outcome
            wmean <- rowMeans(imputations[[imp]][,errorProneCols])
            lambda <- newsigmasq/(newsigmasq+sigmausq/length(errorProneCols))
            xfitted <- xfitted + lambda * (wmean - xfitted)
            newsigmasq <- rep(newsigmasq*(1-lambda), n)

        } else if (method[targetCol]=="logreg") {
          xmod <- glm(xmodformula, family="binomial",data=xmoddata)
          newbeta = modPostDraw(xmod)
          if ((smtype=="casecohort")|(smtype=="nestedcc")) {
            xfitted <- expit(model.matrix(xmodformula, data=imputations[[imp]]) %*% newbeta)
          } else {
            xfitted <- expit(model.matrix(xmod) %*% newbeta)
          }
        } else if (method[targetCol]=="poisson") {
          xmod <- glm(xmodformula, family="poisson", data=xmoddata)
          newbeta = modPostDraw(xmod)
          if ((smtype=="casecohort")|(smtype=="nestedcc")) {
            xfitted <- exp(model.matrix(xmodformula, data=imputations[[imp]]) %*% newbeta)
          } else {
            xfitted <- exp(model.matrix(xmod) %*% newbeta)
          }
        } else if (method[targetCol]=="podds") {
          if (is.ordered(imputations[[imp]][,targetCol])==FALSE) stop("Variables to be imputed using method podds must be stored as ordered factors.")
          xmod <- VGAM::vglm(xmodformula, VGAM::propodds, data=xmoddata)
          xmod.dummy <- VGAM::vglm(xmodformula, VGAM::propodds, data=imputations[[imp]])
          newbeta <- VGAM::coef(xmod) + MASS::mvrnorm(1, mu=rep(0,ncol(VGAM::vcov(xmod))), Sigma=VGAM::vcov(xmod))
          linpreds <- matrix((VGAM::model.matrix(xmod.dummy)) %*% newbeta, byrow=TRUE, ncol=(nlevels(imputations[[imp]][,targetCol])-1))
          cumprobs <- cbind(1/(1+exp(linpreds)), rep(1,nrow(linpreds)))
          xfitted <- cbind(cumprobs[,1] ,cumprobs[,2:ncol(cumprobs)] - cumprobs[,1:(ncol(cumprobs)-1)])
        } else if (method[targetCol]=="mlogit") {
          if (is.factor(imputations[[imp]][,targetCol])==FALSE) stop("Variables to be imputed using method modds must be stored as factors.")
          xmod <- VGAM::vglm(xmodformula, VGAM::multinomial(refLevel=1), data=xmoddata)
          xmod.dummy <- VGAM::vglm(xmodformula, VGAM::multinomial(refLevel=1), data=imputations[[imp]])
          newbeta <- VGAM::coef(xmod) + MASS::mvrnorm(1, mu=rep(0,ncol(VGAM::vcov(xmod))), Sigma=VGAM::vcov(xmod))
          linpreds <- matrix((VGAM::model.matrix(xmod.dummy)) %*% newbeta, byrow=TRUE, ncol=(nlevels(imputations[[imp]][,targetCol])-1))
          denom <- 1+rowSums(exp(linpreds))
          xfitted <-cbind(1/denom, exp(linpreds) / denom)
        }
        if (noisy==TRUE) {
          print(summary(xmod))
        }

        #if latent normal, estimate error variance and calculate required conditional expectation and variance
        #         if (method[targetCol]=="latnorm") {
        #           errorProneCols <- which(errorProneMatrix[targetCol,]==1)
        #           wmean <- rowMeans(imputations[[imp]][,errorProneCols], na.rm=TRUE)
        #           n_i <- apply(imputations[[imp]][,errorProneCols], 1, sumna)
        #           sum_ni <- sum(n_i)
        #           #estimate error variance
        #           if (cyclenum==1) {
        #             xmat <- matrix(wmean, nrow=nrow(imputations[[imp]]), ncol=length(errorProneCols))
        #             uVec <- c(as.matrix(imputations[[imp]][,errorProneCols] - xmat))
        #             sigmausq <- sum(uVec^2, na.rm=TRUE) / (sum_ni - n)
        #           }
        #           else {
        #             xmat <- matrix(imputations[[imp]][,targetCol], nrow=nrow(imputations[[imp]]), ncol=length(errorProneCols))
        #             uVec <- c(as.matrix(imputations[[imp]][,errorProneCols] - xmat))
        #             sigmausq <- sum(uVec^2, na.rm=TRUE) / sum_ni
        #           }
        #           #take draw from posterior of error variance
        #           sigmausq <- sigmausq*sum_ni/rchisq(1,sum_ni)
        #           #calculate conditional mean and variance
        #           lambda <- newsigmasq/(newsigmasq+sigmausq/n_i)
        #           xfitted <- xfitted + lambda * (wmean - xfitted)
        #           newsigmasq <- newsigmasq*(1-lambda)
        #         }

        #estimate parameters of substantive model
        if (smtype=="lm") {
          ymod <- lm(as.formula(smformula),imputations[[imp]])
          beta <- ymod$coef
          sigmasq <- summary(ymod)$sigma^2
          varcov <- vcov(ymod)
          outcomeModResVar <- (sigmasq*ymod$df) / rchisq(1,ymod$df)
          covariance <- (outcomeModResVar/sigmasq)*vcov(ymod)
          outcomeModBeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        else if (smtype=="logistic") {
          ymod <- glm(as.formula(smformula),family="binomial",imputations[[imp]])
          outcomeModBeta = modPostDraw(ymod)
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        else if (smtype=="poisson") {
          ymod <- glm(as.formula(smformula),family="poisson",imputations[[imp]])
          outcomeModBeta = modPostDraw(ymod)
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        else if (smtype=="coxph") {
          ymod <- survival::coxph(as.formula(smformula), imputations[[imp]])
          outcomeModBeta <- modPostDraw(ymod)
          ymod$coefficients <- outcomeModBeta
          basehaz <- survival::basehaz(ymod, centered=FALSE)[,1]
          H0 <- basehaz[H0indices]
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        else if (smtype=="compet") {
          for (cause in 1:numCauses) {
            ymod <- survival::coxph(as.formula(smformula[[cause]]), imputations[[imp]])
            outcomeModBeta[[cause]] <- modPostDraw(ymod)
            ymod$coefficients <- outcomeModBeta[[cause]]
            basehaz <- survival::basehaz(ymod, centered=FALSE)[,1]
            H0[[cause]] <- basehaz[H0indices[[cause]]]
            if (noisy==TRUE) {
              print(summary(ymod))
            }
          }
        }
        else if (smtype=="casecohort") {
          ymod <- survival::coxph(as.formula(smformula2), imputations[[imp]])
          outcomeModBeta <- modPostDraw(ymod)

          cumhaz.denom.elements=exp(model.matrix(as.formula(smformula),imputations[[imp]])[,-1] %*% outcomeModBeta)
          cumhaz.denom=sapply(list.times,function(x){sum(cumhaz.denom.elements[which(originaldata[,timeCol]>=x)]*subco.weight[which(originaldata[,timeCol]>=x)])})
          exp.func.denom=cumsum(1/cumhaz.denom)
          H0.fun=stepfun(list.times,c(0,exp.func.denom))
          H0=H0.fun(originaldata[,timeCol])

          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        else if (smtype=="nestedcc") {
          ymod <- survival::coxph(as.formula(smformula), imputations[[imp]])
          outcomeModBeta <- modPostDraw(ymod)

          explan.matrix<-model.matrix(ymod)
          cumbasehaz.denom<-exp(matrix(outcomeModBeta,nrow=1)%*%t(explan.matrix))*originaldata[,nriskCol]/num.sampriskset
          cumbasehaz.denom<-ave(cumbasehaz.denom, originaldata[,setCol], FUN = sum)[originaldata[,dCol]==1] #this is the denominator of the contribution to the cumulative baseline hazard at each event time
          cumbasehaz.t<-originaldata[,timeCol][originaldata[,dCol]==1] #times to which the baseline cumulative hazards refer
          H0<-unlist(lapply(originaldata[,timeCol],function(x) {sum((1/cumbasehaz.denom)[cumbasehaz.t<=x])}))

          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }

        if ((imp==1) & (cyclenum==1) & (var==1)) {
          if (smtype=="compet") {
            totalCoefVec <- outcomeModBeta[[1]]
            for (cause in 2:numCauses) {
              totalCoefVec <- c(totalCoefVec, outcomeModBeta[[cause]])
            }
            smCoefIter <- array(0, dim=c(m, length(totalCoefVec), numit))
          }
          else {
            smCoefIter <- array(0, dim=c(m, length(outcomeModBeta), numit))
          }
        }

        if (var==length(partialVars)) {
          #then we have reached end of a cycle
          if (smtype=="compet") {
            totalCoefVec <- outcomeModBeta[[1]]
            for (cause in 2:numCauses) {
              totalCoefVec <- c(totalCoefVec, outcomeModBeta[[cause]])
            }
            smCoefIter[imp,,cyclenum] <- totalCoefVec
          }
          else {
            smCoefIter[imp,,cyclenum] <- outcomeModBeta
          }
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
              if (is.factor(imputations[[imp]][,targetCol])==TRUE) {
                valToImpute <- levels(imputations[[imp]][,targetCol])[xMisVal]
              }
              else {
                valToImpute <- xMisVal-1
              }
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
              outcomeDens <- dnorm(deviation, mean=0, sd=outcomeModResVar^0.5)
            }
            else if (smtype=="logistic") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              prob <- expit(outmodxb[imputationNeeded])
              outcomeDens <- prob*imputations[[imp]][imputationNeeded,outcomeCol] + (1-prob)*(1-imputations[[imp]][imputationNeeded,outcomeCol])
            }
            else if (smtype=="poisson") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              outcomeDens <- dpois(imputations[[imp]][imputationNeeded,outcomeCol], exp(outmodxb[imputationNeeded]))
            }
            else if ((smtype=="coxph") | (smtype=="casecohort")) {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]])
              outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              outcomeDens <- exp(-H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))* (exp(outmodxb[imputationNeeded])^d[imputationNeeded])
            }
            else if (smtype=="nestedcc") {
              outmodxb <-  model.matrix(as.formula(smformula2),imputations[[imp]])
              outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              outcomeDens <- exp(-H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))* (exp(outmodxb[imputationNeeded])^d[imputationNeeded])
            }
            else if (smtype=="compet") {
              outcomeDens <- rep(1,length(imputationNeeded))
              for (cause in 1:numCauses) {
                outmodxb <-  model.matrix(linpred[[cause]],imputations[[imp]])
                outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta[[cause]])
                outcomeDens <- outcomeDens * exp(-H0[[cause]][imputationNeeded] * exp(outmodxb[imputationNeeded]))* (exp(outmodxb[imputationNeeded])^(d[imputationNeeded]==cause))
              }
            }
            outcomeDensCovDens[,xMisVal] <- outcomeDens * fittedMean[imputationNeeded,xMisVal]
          }
          directImpProbs = outcomeDensCovDens / rowSums(outcomeDensCovDens)

          if (method[targetCol]=="logreg") {
            directImpProbs = directImpProbs[,2]
            if (is.factor(imputations[[imp]][,targetCol])==TRUE) {
              imputations[[imp]][imputationNeeded,targetCol] <- levels(imputations[[imp]][,targetCol])[1]
              imputations[[imp]][imputationNeeded,targetCol][rbinom(length(imputationNeeded),1,directImpProbs)==1] <- levels(imputations[[imp]][,targetCol])[2]
            }
            else {
              imputations[[imp]][imputationNeeded,targetCol] <- rbinom(length(imputationNeeded),1,directImpProbs)
            }
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
            if ((method[targetCol]=="norm") | (method[targetCol]=="latnorm")) {
              imputations[[imp]][imputationNeeded,targetCol] <- rnorm(length(imputationNeeded),xfitted[imputationNeeded],newsigmasq^0.5)
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
            else if (smtype=="poisson") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              prob = dpois(imputations[[imp]][imputationNeeded,outcomeCol], exp(outmodxb[imputationNeeded]))
              reject = 1*(uDraw>prob)
            }
            else if ((smtype=="coxph") | (smtype=="casecohort")) {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]])
              outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              s_t = exp(-H0[imputationNeeded]* exp(outmodxb[imputationNeeded]))
              prob = exp(1 + outmodxb[imputationNeeded] - (H0[imputationNeeded]* exp(outmodxb[imputationNeeded])) ) * H0[imputationNeeded]
              prob = d[imputationNeeded]*prob + (1-d[imputationNeeded])*s_t
              reject = 1*(uDraw > prob )
            }
            else if (smtype=="nestedcc") {
              outmodxb <-  model.matrix(as.formula(smformula2),imputations[[imp]])
              outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              s_t = exp(-H0[imputationNeeded]* exp(outmodxb[imputationNeeded]))
              prob = exp(1 + outmodxb[imputationNeeded] - (H0[imputationNeeded]* exp(outmodxb[imputationNeeded])) ) * H0[imputationNeeded]
              prob = d[imputationNeeded]*prob + (1-d[imputationNeeded])*s_t
              reject = 1*(uDraw > prob )
            }
            else if (smtype=="compet") {
              prob <- rep(1,length(imputationNeeded))
              for (cause in 1:numCauses) {
                outmodxb <-  model.matrix(linpred[[cause]],imputations[[imp]])
                outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta[[cause]])
                prob = prob * exp(-H0[[cause]][imputationNeeded] * exp(outmodxb[imputationNeeded]))* (H0[[cause]][imputationNeeded]*exp(1+outmodxb[imputationNeeded]))^(d[imputationNeeded]==cause)
              }
              reject = 1*(uDraw > prob )
            }
            imputationNeeded <- imputationNeeded[reject==1]

            j <- j+1
          }

          #now, for those remaining, who must have low acceptance probabilities, sample by subject
          for (i in imputationNeeded) {

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
            else if (method[targetCol]=="latnorm") {
              tempData[,targetCol] <- rnorm(rjlimit,xfitted[i],newsigmasq[i]^0.5)
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
            else if (smtype=="poisson") {
              outmodxb <-  model.matrix(as.formula(smformula),tempData) %*% outcomeModBeta
              prob = dpois(tempData[,outcomeCol], exp(outmodxb))
              reject = 1*(uDraw>prob)
            }
            else if ((smtype=="coxph") | (smtype=="casecohort")) {
              outmodxb <-  model.matrix(as.formula(smformula),tempData)
              outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              s_t = exp(-H0[i]* exp(outmodxb))
              prob = exp(1 + outmodxb - (H0[i]* exp(outmodxb)) ) * H0[i]
              prob = d[i]*prob + (1-d[i])*s_t
              reject = 1*(uDraw > prob )
            }
            else if (smtype=="nestedcc") {
              outmodxb <-  model.matrix(as.formula(smformula2),tempData)
              outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              s_t = exp(-H0[i]* exp(outmodxb))
              prob = exp(1 + outmodxb - (H0[i]* exp(outmodxb)) ) * H0[i]
              prob = d[i]*prob + (1-d[i])*s_t
              reject = 1*(uDraw > prob )
            }
            else if (smtype=="compet") {
              prob <- rep(1,rjlimit)
              for (cause in 1:numCauses) {
                outmodxb <-  model.matrix(linpred[[cause]],tempData)
                outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta[[cause]])
                prob = prob * exp(-H0[[cause]][i] * exp(outmodxb))* (H0[[cause]][i]*exp(1+outmodxb))^(d[i]==cause)
              }
              reject = 1*(uDraw > prob )
            }

            if (sum(reject)<rjlimit) {
              imputations[[imp]][i,targetCol] <- tempData[reject==0,targetCol][1]
            } else {
              rjFailCount <- rjFailCount + 1
            }
          }
          #update passive variables
          imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
        }
      }

      #imputations of missing outcomes, if present (using proper imputation), for regression and logistic
      #substantive models
      if ((smtype=="lm") | (smtype=="logistic")) {
        if (sum(r[,outcomeCol])<n) {
          imputationNeeded <- (1:n)[r[,outcomeCol]==0]
          #estimate parameters of substantive model using those with outcomes observed
          if (smtype=="lm") {
            ymod <- lm(as.formula(smformula),imputations[[imp]][r[,outcomeCol]==1,])
            beta <- ymod$coef
            sigmasq <- summary(ymod)$sigma^2
            varcov <- vcov(ymod)
            outcomeModResVar <- (sigmasq*ymod$df) / rchisq(1,ymod$df)
            covariance <- (outcomeModResVar/sigmasq)*vcov(ymod)
            outcomeModBeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
            outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
            imputations[[imp]][imputationNeeded,outcomeCol] <- rnorm(length(imputationNeeded),outmodxb[imputationNeeded], sigmasq^0.5)
          }
          else if (smtype=="logistic") {
            ymod <- glm(as.formula(smformula),family="binomial",imputations[[imp]][r[,outcomeCol]==1,])
            outcomeModBeta = modPostDraw(ymod)
            outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
            prob <- expit(outmodxb[imputationNeeded])
            imputations[[imp]][imputationNeeded,outcomeCol] <- rbinom(length(imputationNeeded),1,prob)
          }
        }
      }
    }

  }

  if (rjFailCount>0) {
    warning(paste("Rejection sampling failed ",rjFailCount," times (across all variables, iterations, and imputations). You may want to increase the rejection sampling limit.",sep=""))
  }

  list(impDatasets=imputations, smCoefIter=smCoefIter)

}

updatePassiveVars <- function(data, method, passivecols) {
  for (i in passivecols) {
    data[,i] <- with(data, eval(parse(text=method[i])))
  }
  data
}

expit <- function(x) {
  exp(x)/(1+exp(x))
}

sumna <- function(x) {
  sum(is.na(x)==FALSE)
}

#returns first non missing entry of x
firstnonna <- function(x) {
  x[is.na(x)==FALSE][1]
}

catdraw <- function(prob) {
  (1:length(prob))[rmultinom(1,size=1,prob=prob)==1]
}

modPostDraw <- function(modobj) {
  beta <- modobj$coef
  varcov <- vcov(modobj)
  beta + MASS::mvrnorm(1, mu=rep(0,ncol(varcov)), Sigma=varcov)
}

