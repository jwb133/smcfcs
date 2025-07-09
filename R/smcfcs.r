#' Substantive model compatible fully conditional specification imputation of covariates.
#'
#' Multiply imputes missing covariate values using substantive model compatible
#' fully conditional specification.
#'
#' smcfcs imputes missing values of covariates using the Substantive Model Compatible
#' Fully Conditional Specification multiple imputation approach proposed by
#' Bartlett \emph{et al} 2015 (see references).
#'
#' Imputation is supported for linear regression (\code{"lm"}),
#' logistic regression (\code{"logistic"}), bias reduced logistic regression (\code{"brlogistic"}),
#' Poisson regression (\code{"poisson"}), Weibull (\code{"weibull"}) and Cox regression
#' for time to event data (\code{"coxph"}),
#' and Cox models for competing risks data (\code{"compet"}). For \code{"coxph"},
#' the event indicator should be integer coded with 0 for censoring and 1 for event.
#' For \code{"compet"}, a Cox model is assumed for each cause specific hazard function,
#' and the event indicator
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
#' The bias reduced methods make use of the \code{\link[brglm2]{brglm2}} package to fit the corresponding glms
#' using Firth's bias reduced approach. These may be particularly useful to use in case
#' of perfect prediction, since the resulting model estimates are always guaranteed to be
#' finite, even in the case of perfect prediction.
#'
#' The development of this package was supported by the UK Medical Research Council
#' (Fellowship MR/K02180X/1 and grant MR/T023953/1). Part of its development took place while Bartlett was
#' kindly hosted by the University of Michigan's Department of Biostatistics & Institute for
#' Social Research.
#'
#' The structure of many of the arguments to \code{smcfcs} are based on those of
#' the excellent \code{mice} package.
#'
#' @param originaldata The original data frame with missing values.
#' @param smtype A string specifying the type of substantive model. Possible
#' values are \code{"lm"}, \code{"logistic"}, \code{"brlogistic"}, \code{"poisson"}, \code{"weibull"},
#' \code{"coxph"}, \code{"compet"}.
#' @param smformula The formula of the substantive model. For \code{"weibull"} and \code{"coxph"}
#' substantive models the left hand side should be of the form \code{"Surv(t,d)"}. For \code{"compet"}
#' substantive models, a list should be passed consisting of the Cox models
#' for each cause of failure (see example).
#' @param method A required vector of strings specifying for each variable either
#' that it does not need to be imputed (""), the type of regression model to be
#' be used to impute. Possible values are \code{"norm"} (normal linear regression),
#' \code{"logreg"} (logistic regression), \code{"brlogreg"} (bias reduced logistic regression),
#' \code{"poisson"} (Poisson regression),
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
#' on by smcfcs, and should not be specified as a predictor
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
#' as \code{"latnorm"}. See \code{vignette("coverror", package = "smcfcs")} for more details.
#' @param restrictions Optional string which specifies restrictions for handling
#' coarsened factor level covariates. This is where for a factor variable for
#' some individuals we do not their value of the variable but we do know
#' it belongs to some subset of the sample space. For further details
#' on how to specify this argument, see \code{vignette("coarsening", package = "smcfcs")}.
#'
#' @return A list containing:
#'
#' \code{impDatasets} a list containing the imputed datasets
#'
#' \code{smCoefIter} a three dimension matrix containing the substantive model parameter
#' values. The matrix is indexed by [imputation,parameter number,iteration]
#'
#' @author Jonathan Bartlett \email{jonathan.bartlett1@@lshtm.ac.uk}
#'
#' @example data-raw/examples.r
#'
#' @references Bartlett JW, Seaman SR, White IR, Carpenter JR. Multiple imputation of covariates
#' by fully conditional specification: accommodating the substantive model. Statistical Methods
#' in Medical Research 2015; 24(4): 462-487. \doi{10.1177/0962280214521348}
#' @import stats
#' @importFrom survival Surv
#' @export
smcfcs <- function(originaldata, smtype, smformula, method, predictorMatrix = NULL, m = 5, numit = 10, rjlimit = 1000, noisy = FALSE, errorProneMatrix = NULL, restrictions = NULL) {
  # call core  smcfcs function, passing through arguments
  smcfcs.core(originaldata, smtype, smformula, method, predictorMatrix, m, numit, rjlimit, noisy, errorProneMatrix = errorProneMatrix, restrictions = restrictions)
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
#' @author Jonathan Bartlett \email{jonathan.bartlett1@@lshtm.ac.uk}
#'
#' @param originaldata The case-cohort data set (NOT a full cohort data set with a case-cohort substudy within it)
#' @param smformula A formula of the form "Surv(entertime,t,d)~x", where d is the event (d=1) or censoring (d=0) indicator, t is the event or censoring time and entertime is equal to the time origin (typically 0) for individuals in the subcohort and is equal to (t-0.001) for cases outside the subcohort [this sets cases outside the subcohort to enter follow-up just before their event time. The value 0.001 may need to be modified depending on the time scale.]
#' @param in.subco The name of a column in the dataset with 0/1s that indicates whether the subject is in the subcohort
#' @param sampfrac The proportion of individuals from the underlying full cohort who are in the subcohort
#' @param ... Additional arguments to pass on to \link[smcfcs]{smcfcs}
#'
#' @inheritParams smcfcs
#'
#' @example data-raw/cc_example.r
#'
#' @export
smcfcs.casecohort <- function(originaldata, smformula, method, sampfrac, in.subco, ...) {
  smcfcs.core(originaldata=originaldata,
              smformula=smformula,
              method=method,
              smtype = "casecohort",
              sampfrac = sampfrac,
              in.subco = in.subco,
              ...)
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
#' @author Jonathan Bartlett \email{jonathan.bartlett1@@lshtm.ac.uk}
#'
#' @param originaldata The nested case-control data set (NOT a full cohort data set with a case-cohort substudy within it)
#' @param smformula A formula of the form "Surv(t,case)~x+strata(set)", where case is case-control indicator, t is the event or censoring time. Note that t could be set to the case's event time for the matched controls in a given set. The right hand side should include the case control set as a strata term (see example).
#' @param set variable identifying matched sets in nested case-control study
#' @param event variable which indicates who is a case/control in the nested case-control sample. Note that this is distinct from d.
#' @param nrisk variable which is the number at risk (in the underlying full cohort) at the event time for the case in each matched set (i.e. nrisk is the same for all individuals in a matched set).
#' @param ... Additional arguments to pass on to \link[smcfcs]{smcfcs}
#'
#' @inheritParams smcfcs
#' @example data-raw/ncc_example.r
#' @export
smcfcs.nestedcc <- function(originaldata, smformula, method, set, event, nrisk, ...) {
  smcfcs.core(originaldata=originaldata,
              smformula=smformula,
              method=method,
              smtype = "nestedcc",
              set = set,
              event = event,
              nrisk = nrisk,
              ...)
}

#' Substantive model compatible fully conditional specification imputation of covariates for
#' discrete time survival analysis
#'
#' Multiply imputes missing covariate values using substantive model compatible
#' fully conditional specification for discrete time survival analysis.
#'
#' For this substantive model type, like for the other substantive model types, \code{smcfcs} expects the \code{originaldata} to have
#' one row per subject. Variables indicating the discrete time of failure/censoring
#' and the event indicator should be passed in \code{smformula}, as described.
#'
#' The default is to model the effect of time as a factor. This will not work in datasets where
#' there is not at least one observed event in each time period. In such cases you must specify
#' a simpler parametric model for the effect of time. At the moment you can specify either a linear or quadratic
#' effect of time (on the log odds scale).
#'
#' @author Jonathan Bartlett \email{jonathan.bartlett1@@lshtm.ac.uk}
#'
#' @param originaldata The data in wide form (i.e. one row per subject)
#' @param smformula A formula of the form "Surv(t,d)~x1+x2+x3", where t is the discrete time variable, d is the binary event
#' indicator, and the covariates should not include time. The time variable should be
#' an integer coded numeric variable taking values from 1 up to the final time period.
#' @param timeEffects Specifies how the effect of time is modelled. \code{timeEffects="factor"} (the default) models time as a
#' factor variable. \code{timeEffects="linear"} and \code{timeEffects="quad"} specify that time be modelled as a continuous
#' linear or quadratic effect on the log odds scale respectively.
#' @param ... Additional arguments to pass on to \link[smcfcs]{smcfcs}
#'
#' @inheritParams smcfcs
#' @example data-raw/dtsam_example.r
#' @export
smcfcs.dtsam <- function(originaldata, smformula, method, timeEffects = "factor",...) {
  smcfcs.core(originaldata=originaldata,
              smformula=smformula,
              method=method,
              smtype = "dtsam",
              timeEffects = timeEffects,
              ...)
}

# this is the core of the smcfcs function, called by wrapper functions for certain different substantive models
smcfcs.core <- function(originaldata, smtype, smformula, method, predictorMatrix = NULL, m = 5, numit = 10, rjlimit = 1000, noisy = FALSE, errorProneMatrix = NULL,restrictions = NULL,
                        ...) {
  # get extra arguments passed in ...
  extraArgs <- list(...)

  stopifnot(is.data.frame(originaldata))
  if (ncol(originaldata) != length(method)) stop("Method argument must have the same length as the number of columns in the data frame.")

  n <- dim(originaldata)[1]

  # create matrix of response indicators
  r <- 1 * (is.na(originaldata) == 0)

  if ((smtype %in% c(
    "lm", "logistic", "brlogistic", "poisson", "coxph", "compet", "casecohort", "nestedcc",
    "weibull", "dtsam", "flexsurv"
  )) == FALSE) {
    stop(paste("Substantive model type ", smtype, " not recognised.", sep = ""))
  }

  # find column numbers of partially observed, fully observed variables, and outcome
  if (smtype == "flexsurv") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[, dCol]
    if (!(all(sort(unique(d)) == c(0, 1))) & !(all(unique(d) == 1))) {
      stop("Event indicator for flexsurv must be coded 0/1 for censoring/event.")
    }
    if ((sum(is.na(d))+sum(is.na(originaldata[,timeCol])))>0) {
      stop("Event indicator and time variables should not have NAs.")
    }
    # stop if time-varying effects of a variable for which rejection
    # sampling would be used has been specified, since bound is not correct
    # in this case
    if (grepl("gamma",smformula)==TRUE) {
      timeVaryingVars <- extract_gamma_values(smformula)
      for (i in 1:length(timeVaryingVars)) {
        if ((method[which(colnames(originaldata)==timeVaryingVars[i])] %in% c("norm", "latnorm", "poisson"))==TRUE) {
          stop("You cannot use method norm/latnorm/poisson for covariates with time-varying effects (yet).")
        }
      }
    }
  } else if (smtype == "coxph") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[, dCol]
    if (!(all(sort(unique(d)) == c(0, 1))) & !(all(unique(d) == 1))) {
      stop("Event indicator for coxph must be coded 0/1 for censoring/event.")
    }
    if ((sum(is.na(d))+sum(is.na(originaldata[,timeCol])))>0) {
      stop("Event indicator and time variables should not have NAs.")
    }

    nullMod <- survival::coxph(survival::Surv(originaldata[, timeCol], originaldata[, dCol]) ~ 1,
      control = survival::coxph.control(timefix = FALSE)
    )
    basehaz <- survival::basehaz(nullMod)
    H0indices <- match(originaldata[, timeCol], basehaz[, 2])
    rm(nullMod)
  } else if (smtype == "weibull") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[, dCol]
  } else if (smtype == "dtsam") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[, dCol]
    # determine cut points
    cutPoints <- 1:max(originaldata[, timeCol])
    nTimePoints <- length(cutPoints)

    # check all time points are integers
    if (!all(unique(originaldata[, timeCol]) == floor(unique(originaldata[, timeCol])))) {
      stop("Your time variable must only take positive integer values.")
    }
    # check all time points are positive
    if (any(unique(originaldata[, timeCol]) <= 0)) {
      stop("Your time variable must only take positive integer values.")
    }

    # if factor time effects, check there are some events at each integer value
    # so that model can fit
    if (extraArgs$timeEffects == "factor") {
      if (!identical(sort(unique(originaldata[, timeCol][originaldata[, dCol] == 1])), as.numeric(cutPoints))) {
        stop("You cannot fit a dtsam model with factor time effects since there are some periods with no events. See documentation for timeEffects argument for parametric alternatives")
      }
    }
  } else if (smtype == "compet") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula[[1]])[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula[[1]])[[2]][[3]][[2]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[, dCol]
    numCauses <- length(smformula)
    H0 <- vector("list", numCauses)
    H0indices <- vector("list", numCauses)
    outcomeModBeta <- vector("list", numCauses)
    linpred <- vector("list", numCauses)
    for (cause in 1:numCauses) {
      nullMod <- survival::coxph(as.formula(paste(strsplit(smformula[[cause]], "~")[[1]][1], "~1")), originaldata,
        control = survival::coxph.control(timefix = FALSE)
      )
      basehaz <- survival::basehaz(nullMod)
      H0[[cause]] <- basehaz[, 1]
      H0indices[[cause]] <- match(originaldata[, timeCol], basehaz[, 2])
      linpred[[cause]] <- as.formula(smformula[[cause]])
    }
    rm(nullMod)
  } else if (smtype == "casecohort") {
    subcoCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% extraArgs$in.subco]
    # subcoMembers is a vector of row numbers of those in the subcohort
    subcoMembers <- which(originaldata[, subcoCol] == 1)

    # generate weights for use in later analysis which we use to obtain baseline cumulative hazard
    # assign a weight of /samp.frac to individuals in the subcohort and 0 to those outside the subcohort
    subco.weight <- ifelse(originaldata[, subcoCol] == 1, 1 / extraArgs$sampfrac, 0)


    entertimeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[4]])]
    outcomeCol <- c(entertimeCol, timeCol, dCol)
    d <- originaldata[, dCol]

    # list of unique event times - used in calculation of baseline cumulative hazard
    list.times <- sort(unique(originaldata[, timeCol][originaldata[, dCol] == 1])) # RUTH 21/03/17: ADDED THIS LINE

    smcfcsid <- 1:n
    smformula2 <- paste(smformula, "+cluster(smcfcsid)", sep = "")
  } else if (smtype == "nestedcc") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    setCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(extraArgs$set)]
    nriskCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(extraArgs$nrisk)]
    eventCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(extraArgs$event)] # this is distinct from the dCol

    # the below command creates "Surv(t,case)~x" from "Surv(t,case)~x+strata(setno)" (for example) (i.e. it removes the strata part of the formula)
    # this is used when obtaining outmodxb
    # note this is done slightly oddly, but this is because as.formula does not work well for long formulas as it splits across lines
    exp1 <- as.formula(paste(smformula))[[2]]
    exp2 <- as.formula(smformula)[[3]][[2]]
    smformula2 <- paste(deparse(exp1), "~", deparse(exp2, width.cutoff = 500L))

    # This is the indicator of whether an individual ever has the event (regardless of whether they are sometimes used as a control and sometimes (one) as a case)
    d <- originaldata[, eventCol]

    # noncases is a vector of row numbers of those who never have the event (which is a subset of the controls)
    noncases <- which(originaldata[, eventCol] == 0)

    # number of individuals in each sampled risk set (matched set)
    num.sampriskset <- ave(rep(1, dim(originaldata)[1]), originaldata[, setCol], FUN = function(x) sum(x))
  } else {
    outcomeCol <- which(colnames(originaldata) == as.formula(smformula)[[2]])
  }

  if ((smtype == "logistic") | (smtype == "brlogistic")) {
    if (is.numeric(originaldata[, outcomeCol]) == FALSE) {
      stop("For logistic substantive models the outcome variable must be numeric 0/1.")
    } else {
      if (all.equal(sort(unique(originaldata[, outcomeCol])), c(0, 1)) == FALSE) {
        stop("For logistic substantive models the outcome variable must be coded 0/1.")
      }
    }
  }

  if (smtype == "compet") {
    smcovnames <- attr(terms(as.formula(smformula[[1]])), "term.labels")
    for (cause in 2:numCauses) {
      smcovnames <- c(smcovnames, attr(terms(as.formula(smformula[[cause]])), "term.labels"))
    }
    smcovnames <- unique(smcovnames)
  } else if (smtype == "nestedcc") {
    smcovnames <- attr(terms(as.formula(smformula)), "term.labels")[-length(attr(terms(as.formula(smformula)), "term.labels"))]
  } else {
    smcovnames <- attr(terms(as.formula(smformula)), "term.labels")
  }
  smcovcols <- (1:ncol(originaldata))[colnames(originaldata) %in% smcovnames]

  # partial vars are those variables for which an imputation method has been specified among the available regression types
  partialVars <- which((method == "norm") | (method == "latnorm") | (method == "logreg") | (method == "poisson") |
    (method == "podds") | (method == "mlogit") | (method == "brlogreg"))

  if (length(partialVars) == 0) {
    if (smtype!="flexsurv") {
      stop("You have not specified any valid imputation methods in the method argument.")
    } else {
      if (extraArgs$imputeTimes==FALSE) {
        stop("You have not specified any valid imputation methods in the method argument.")
      }
    }
  }

  # check that methods are given for each partially observed column, and not given for fully observed columns
  for (colnum in 1:ncol(originaldata)) {
    if (method[colnum] != "") {
      # an imputation method has been specified
      if (colnum %in% outcomeCol) {
        stop(paste("An imputation method has been specified for ", colnames(originaldata)[colnum],
          ". Elements of the method argument corresponding to the outcome variable(s) should be empty.",
          sep = ""
        ))
      } else {
        if (sum(r[, colnum]) == n) {
          stop(paste("An imputation method has been specified for ", colnames(originaldata)[colnum],
            ", but it appears to be fully observed.",
            sep = ""
          ))
        }
      }
    } else {
      # no imputation method has been specified
      if (sum(r[, colnum]) < n) {
        # some values are missing
        if (((colnum %in% outcomeCol) == FALSE) & (sum(errorProneMatrix[, colnum]) == 0)) {
          stop(paste("Variable ", colnames(originaldata)[colnum], " does not have an imputation method specified, yet appears to have missing values.", sep = ""))
        }
      }
    }
  }

  # check that if any variables have latnorm specified as method, that at least two error-prone measures
  # are specified for it in errorProneMatrix
  if ("latnorm" %in% method) {
    if (is.null(errorProneMatrix) == TRUE) {
      stop("If you specify method latnorm you must specify the errorProneMatrix argument.")
    } else {
      # check errorProneMatrix is of correct dimensional
      if (identical(dim(errorProneMatrix), c(length(originaldata), length(originaldata))) == FALSE) {
        stop("The errorProneMatrix should be a square matrix with number of rows equal to the number of variables in the dataset.")
      }
      # check entries of errorProneMatrix only consists of 0s and 1s
      if (identical(sort(unique(as.vector(errorProneMatrix))), c(0, 1)) == FALSE) {
        stop("The errorProneMatrix should only consist of 0s and 1s.")
      }
      # check each latnorm variable has at least 2 error-prones
      for (varNum in 1:length(method)) {
        if (method[varNum] == "latnorm") {
          if (sum(errorProneMatrix[varNum, ]) < 2) {
            stop("Each latnorm variable must have two or more error prone measurements specified in the errorProneMatrix argument.")
          }
        }
      }
      # check no error-prone measurement is allocated to more than one latnorm
      if (sum(colSums(errorProneMatrix) > 1) > 0) {
        stop("Each error-prone measurement should be allocated to exactly one latnorm variable.")
      }
    }
  }
  if (is.null(errorProneMatrix) == FALSE) {
    if (("latnorm" %in% method) == FALSE) {
      stop("If you specify errorProneMatrix then at least one variable must be imputed using latnorm.")
    }
  }

  # fully observed vars are those that are fully observed and are covariates in the substantive model
  fullObsVars <- which((colSums(r) == n) & (colnames(originaldata) %in% smcovnames))

  # passive variables
  passiveVars <- which((method != "") & (method != "norm") & (method != "logreg") & (method != "poisson") & (method != "podds") &
    (method != "mlogit") & (method != "latnorm") & (method != "brlogreg"))

  print(paste("Outcome variable(s):", paste(colnames(originaldata)[outcomeCol], collapse = ",")))
  print(paste("Passive variables:", paste(colnames(originaldata)[passiveVars], collapse = ",")))
  print(paste("Partially obs. variables:", paste(colnames(originaldata)[partialVars], collapse = ",")))
  print(paste("Fully obs. substantive model variables:", paste(colnames(originaldata)[fullObsVars], collapse = ",")))

  imputations <- list()
  for (imp in 1:m) {
    imputations[[imp]] <- originaldata
  }

  rjFailCount <- 0
  flexsurvFailCount <- 0

  if(!is.null(restrictions)){
    if (!requireNamespace("stringr", quietly = TRUE)) {
      stop("The package 'stringr' is required if using restrictions. Please install it.", call. = FALSE)
    }

    index_restrictions = restrictions_index = which(unlist(lapply(restrictions, function(sub){TRUE %in% (sub != "")}))); len_restrictions = restrictions_len = length(restrictions_index)

    if(restrictions_len == 0){
      stop("No additional information is supplied")
    }
    if(length(restrictions) != ncol(originaldata)) {
      stop("Restrictions should be a list with length equal to the number of columns/variables.")
    }

    restrictions_work = restrictions_work_conds = rep(list(NULL), length(restrictions))
    for(i in restrictions_index){
      for(j in 1:length(restrictions[[i]])){
        restrictions_work[[i]][[j]] = as.list(stringr::str_trim(unlist(strsplit(restrictions[[i]][[j]], "[=~]+")), side = "both"))

        if(length(restrictions_work[[i]][[j]]) == 2){
          restrictions_work[[i]][[j]] = append(restrictions_work[[i]][[j]], stringr::str_extract_all(restrictions_work[[i]][[j]][[2]], "[0-9A-Za-z/]+"))
        } else if(length(restrictions_work[[i]][[j]]) == 3){
          restrictions_work[[i]][[j]][[3]] = unlist(stringr::str_extract_all(restrictions_work[[i]][[j]][[3]], "[0-9A-Za-z/]+"))
        }
      }
    }

    for (var in 1:length(partialVars)){
      restrictions_work_conds[[var]] = unlist(lapply(restrictions_work[[var]], function(sub_rest){sub_rest[[2]]}))
    }
  }

  for (imp in 1:m) {
    print(paste("Imputation ", imp))

    # initial imputation of each partially observed variable based on observed values
    if (length(partialVars)>0) {
      for (var in 1:length(partialVars)) {
        targetCol <- partialVars[var]
        if (method[targetCol] == "latnorm") {
          # first impute any missing replicate error-prone measurements of this variable by a randomly chosen observed value
          errorProneCols <- which(errorProneMatrix[targetCol, ] == 1)
          for (measure in 1:length(errorProneCols)) {
            if (sum(r[, errorProneCols[measure]]) < n) {
              imputations[[imp]][r[, errorProneCols[measure]] == 0, errorProneCols[measure]] <- sample(imputations[[imp]][
                r[, errorProneCols[measure]] == 1,
                errorProneCols[measure]
              ], size = sum(r[, errorProneCols[measure]] == 0), replace = TRUE)
            }
          }

          # initialize latent predictors with mean of their error-prone measurements
          imputations[[imp]][, targetCol] <- apply(imputations[[imp]][, errorProneCols], 1, mean)
        } else {
          if(length(imputations[[imp]][r[,targetCol]==1,targetCol]) != 0){
            imputations[[imp]][r[, targetCol] == 0, targetCol] <- sample(imputations[[imp]][r[, targetCol] == 1, targetCol], size = sum(r[, targetCol] == 0), replace = TRUE)
          }

          if(!is.null(restrictions) && var %in% index_restrictions){
            len_restrictions = length(restrictions_work[[var]])

            for(i in 1:len_restrictions){
              restrictions_work_vari = restrictions_work[[var]][[i]]
              res_var = restrictions_work_vari[[1]]; res_cond = restrictions_work_vari[[2]]; res_opts = restrictions_work_vari[[3]]

              ## To determine which NAs must be re-sampled in this run
              replace_index = which(is.na(originaldata[, targetCol]) & originaldata[, res_var] == res_cond); len_replace_index = length(replace_index)

              if(method[var] == "mlogit"){

                ## Index on the individuals that match the individuals we want to re-sample
                observed_index = originaldata[r[, targetCol] == 1, targetCol] %in% res_opts

                if(sum(observed_index) != 0){
                  sub_rest_i = originaldata[r[, targetCol] == 1, targetCol][observed_index]
                } else {
                  sub_rest_i = res_opts
                }

                ## Re-sample for the individuals that match the supplied restrictions (e.g., those that are present for C)
                imputations[[imp]][replace_index, targetCol] = sample(x = sub_rest_i, size = len_replace_index, replace = TRUE)


              } else if(method[var] == "norm"){

                res_opts = as.numeric(res_opts)

                observed_index = (originaldata[r[, targetCol] == 1, targetCol] > res_opts[1] & originaldata[r[, targetCol] == 1, targetCol] <= res_opts[2])

                if(sum(observed_index) != 0){
                  sub_rest_i = originaldata[r[, targetCol] == 1, targetCol][observed_index]

                  ## Re-sample for the individuals that match the supplied restrictions (e.g., those that are present for C)
                  imputations[[imp]][replace_index, targetCol] = sample(x = sub_rest_i, size = len_replace_index, replace = TRUE)

                } else {

                  ## Sample from uniform distribution for coarsened observations
                  imputations[[imp]][replace_index, targetCol] = runif(len_replace_index, min = min(res_opts), max = max(res_opts))

                }
              }
            }

            if(length(imputations[[imp]][r[, targetCol] == 1, targetCol]) == 0 && method[var] == "norm"){
              ## Sample missing observations from previously imputed coarsened observations if there are no fully observed observations
              index_fully_missing = is.na(imputations[[imp]][, targetCol])
              imputations[[imp]][index_fully_missing, targetCol] = sample(imputations[[imp]][!index_fully_missing, targetCol], size = sum(index_fully_missing), replace = TRUE)
            }
          }
        }
      }
    }

    # initial imputations of missing outcomes, if present (using improper imputation)
    if ((smtype == "lm") | (smtype == "logistic") | (smtype == "brlogistic") | (smtype == "poisson")) {
      if (sum(r[, outcomeCol]) < n) {
        if (imp == 1) {
          print("Imputing missing outcomes using specified substantive model.")
        }
        # update passive variable(s)
        imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)

        imputationNeeded <- (1:n)[r[, outcomeCol] == 0]
        # estimate parameters of substantive model
        if (smtype == "lm") {
          ymod <- stats::lm(as.formula(smformula), imputations[[imp]])
          beta <- ymod$coef
          sigmasq <- summary(ymod)$sigma^2
          # fill out missing values so that model.matrix works for all rows
          imputations[[imp]][imputationNeeded, outcomeCol] <- 0
          outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]]) %*% beta
          imputations[[imp]][imputationNeeded, outcomeCol] <- rnorm(length(imputationNeeded), outmodxb[imputationNeeded], sigmasq^0.5)
        } else if ((smtype == "logistic") | (smtype == "brlogistic")) {
          if (smtype == "logistic") {
            ymod <- glm(as.formula(smformula), family = "binomial", imputations[[imp]])
          } else {
            ymod <- glm(as.formula(smformula),
              family = "binomial", imputations[[imp]],
              method = brglm2::brglmFit
            )
          }
          beta <- ymod$coef
          imputations[[imp]][imputationNeeded, outcomeCol] <- 0
          outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]]) %*% beta
          prob <- plogis(outmodxb[imputationNeeded])
          imputations[[imp]][imputationNeeded, outcomeCol] <- rbinom(length(imputationNeeded), 1, prob)
        } else if (smtype == "poisson") {
          ymod <- glm(as.formula(smformula), family = "poisson", imputations[[imp]])
          beta <- ymod$coef
          imputations[[imp]][imputationNeeded, outcomeCol] <- 0
          outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]]) %*% beta
          imputations[[imp]][imputationNeeded, outcomeCol] <- rpois(length(imputationNeeded), exp(outmodxb[imputationNeeded]))
        }
      }
    }

    for (cyclenum in 1:numit) {
      if (noisy == TRUE) {
        print(paste("Iteration ", cyclenum))
      }
      # update passive variable(s)
      imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)

      if (length(partialVars)>0) {
        for (var in 1:length(partialVars)) {
        targetCol <- partialVars[var]
        if (is.null(predictorMatrix)) {
          predictorCols <- c(partialVars[!partialVars %in% targetCol], fullObsVars)
        } else {
          predictorCols <- which(predictorMatrix[targetCol, ] == 1)
          # ensure that user has not included outcome variable(s) here
          predictorCols <- predictorCols[!predictorCols %in% outcomeCol]
        }
        if ((imp == 1) & (cyclenum == 1)) {
          if (method[targetCol] == "latnorm") {
            print(paste("Imputing: ", colnames(imputations[[imp]])[targetCol], " using ", paste(colnames(imputations[[imp]])[c(predictorCols, which(errorProneMatrix[targetCol, ] == 1))], collapse = ","), " plus outcome", collapse = ","))
          } else {
            print(paste("Imputing: ", colnames(imputations[[imp]])[targetCol], " using ", paste(colnames(imputations[[imp]])[predictorCols], collapse = ","), " plus outcome", collapse = ","))
          }
        }
        if (length(predictorCols) > 0) {
          xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], "~", paste(colnames(imputations[[imp]])[predictorCols], collapse = "+"), sep = ""))
        } else {
          xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], "~1", sep = ""))
        }
        if (smtype == "casecohort") {
          xmoddata <- imputations[[imp]][subcoMembers, ]
        } else if (smtype == "nestedcc") {
          xmoddata <- imputations[[imp]][noncases, ]
        } else {
          xmoddata <- imputations[[imp]]
        }
        if (method[targetCol] == "norm") {
          # estimate parameters of covariate model
          xmod <- lm(xmodformula, data = xmoddata)
          # take draw from posterior of covariate model parameters
          beta <- xmod$coef
          sigmasq <- summary(xmod)$sigma^2
          newsigmasq <- (sigmasq * xmod$df) / rchisq(1, xmod$df)
          covariance <- (newsigmasq / sigmasq) * vcov(xmod)
          newbeta <- beta + MASS::mvrnorm(1, mu = rep(0, ncol(covariance)), Sigma = covariance)
          # calculate fitted values
          if ((smtype == "casecohort") | (smtype == "nestedcc")) {
            xfitted <- model.matrix(xmodformula, data = imputations[[imp]]) %*% newbeta
          } else {
            xfitted <- model.matrix(xmod) %*% newbeta
          }
        } else if (method[targetCol] == "latnorm") {
          # estimate parameters of covariate model
          xmod <- lm(xmodformula, data = xmoddata)
          # take draw from posterior of covariate model parameters
          beta <- xmod$coef
          sigmasq <- summary(xmod)$sigma^2
          # draw from sigmasq posterior based on proper inverse gamma prior for sigmasq
          # prior equivalent to 1 observation and guess of sigmasq=1
          newsigmasq <- 1 / rgamma(1, shape = ((n + 1) / 2), rate = ((n * sigmasq + 1) / 2))
          covariance <- (newsigmasq / sigmasq) * vcov(xmod)
          newbeta <- beta + MASS::mvrnorm(1, mu = rep(0, ncol(covariance)), Sigma = covariance)
          # calculate fitted values
          if ((smtype == "casecohort") | (smtype == "nestedcc")) {
            xfitted <- model.matrix(xmodformula, data = imputations[[imp]]) %*% newbeta
          } else {
            xfitted <- model.matrix(xmod) %*% newbeta
          }

          # estimate error variance and draw new value of error variance
          errorProneCols <- which(errorProneMatrix[targetCol, ] == 1)
          xmat <- matrix(imputations[[imp]][, targetCol], nrow = nrow(imputations[[imp]]), ncol = length(errorProneCols))
          uVec <- c(as.matrix(imputations[[imp]][, errorProneCols] - xmat))
          sigmausq <- mean(uVec^2)
          # take draw from posterior of error variance, using proper inverse gamma prior
          sum_ni <- n * length(errorProneCols)
          sigmausq <- 1 / rgamma(1, shape = ((sum_ni + 1) / 2), rate = ((sum_ni * sigmausq + 1) / 2))

          # re-impute any originally missing error-prone measurements, based on classical error model assumption
          for (measure in 1:length(errorProneCols)) {
            nToImpute <- n - sum(r[, errorProneCols[measure]])
            if (nToImpute > 0) {
              # then some values need imputing
              imputations[[imp]][r[, errorProneCols[measure]] == 0, errorProneCols[measure]] <- imputations[[imp]][r[, errorProneCols[measure]] == 0, targetCol] +
                rnorm(nToImpute, 0, sd = sqrt(sigmausq))
            }
          }

          # calculate conditional mean and variance of X|everything else except outcome
          wmean <- rowMeans(imputations[[imp]][, errorProneCols])
          lambda <- newsigmasq / (newsigmasq + sigmausq / length(errorProneCols))
          xfitted <- xfitted + lambda * (wmean - xfitted)
          newsigmasq <- rep(newsigmasq * (1 - lambda), n)
        } else if (method[targetCol] == "logreg") {
          xmod <- glm(xmodformula, family = "binomial", data = xmoddata)
          newbeta <- modPostDraw(xmod)
          if ((smtype == "casecohort") | (smtype == "nestedcc")) {
            xfitted <- plogis(model.matrix(xmodformula, data = imputations[[imp]]) %*% newbeta)
          } else {
            xfitted <- plogis(model.matrix(xmod) %*% newbeta)
          }
        } else if (method[targetCol] == "brlogreg") {
          xmod <- glm(xmodformula, family = "binomial", data = xmoddata, method = brglm2::brglmFit)
          newbeta <- modPostDraw(xmod)
          if ((smtype == "casecohort") | (smtype == "nestedcc")) {
            xfitted <- plogis(model.matrix(xmodformula, data = imputations[[imp]]) %*% newbeta)
          } else {
            xfitted <- plogis(model.matrix(xmod) %*% newbeta)
          }
        } else if (method[targetCol] == "poisson") {
          xmod <- glm(xmodformula, family = "poisson", data = xmoddata)
          newbeta <- modPostDraw(xmod)
          if ((smtype == "casecohort") | (smtype == "nestedcc")) {
            xfitted <- exp(model.matrix(xmodformula, data = imputations[[imp]]) %*% newbeta)
          } else {
            xfitted <- exp(model.matrix(xmod) %*% newbeta)
          }
        } else if (method[targetCol] == "podds") {
          if (is.ordered(imputations[[imp]][, targetCol]) == FALSE) stop("Variables to be imputed using method podds must be stored as ordered factors.")
          xmod <- MASS::polr(xmodformula, data = xmoddata, Hess=TRUE)
          xmod.dummy <- MASS::polr(xmodformula, data = imputations[[imp]], Hess=TRUE)
          polr_cutpts <- xmod$zeta
          # Thanks to Ed Bonneville: polr does maximization using first cutpoint and then log differences
          # in the cut-points. So we extract estimates use this transformation, and then apply a MVN draw
          polr_theta <- c(xmod$coefficients, polr_cutpts[1], log(diff(polr_cutpts)))
          newtheta <- polr_theta + MASS::mvrnorm(1, mu = rep(0, length(polr_theta)), Sigma = MASS::ginv(xmod$Hessian))
          # extract new 'coefficients'
          newbeta <- newtheta[seq_along(xmod$coefficients)]
          # transform back to cutpoints
          newcutpoints <- cumsum(c(newtheta[length(xmod$coefficients)+1], exp(newtheta[-(1:(length(xmod$coefficients)+1))])))
          # calculate fitted probabilities
          polr_cumlogitprob <- array(0, dim=c(dim(imputations[[imp]])[1],length(newcutpoints)))
          for (j in 1:length(newcutpoints)) {
            polr_cumlogitprob[,j] <- model.matrix(xmod) %*% c(-newcutpoints[j],newbeta)
          }
          polr_cumprobs <- cbind(rep(1,dim(imputations[[imp]])[1]),plogis(polr_cumlogitprob))
          xfitted <- cbind(-t(apply(polr_cumprobs, MARGIN=1, FUN=diff)), polr_cumprobs[,dim(polr_cumprobs)[2]])
        } else if (method[targetCol] == "mlogit") {
          if (is.factor(imputations[[imp]][, targetCol]) == FALSE) stop("Variables to be imputed using method mlogit must be stored as factors.")
          xmod <- VGAM::vglm(xmodformula, VGAM::multinomial(refLevel = 1), data = xmoddata)
          xmod.dummy <- VGAM::vglm(xmodformula, VGAM::multinomial(refLevel = 1), data = imputations[[imp]])
          newbeta <- VGAM::coef(xmod) + MASS::mvrnorm(1, mu = rep(0, ncol(VGAM::vcov(xmod))), Sigma = VGAM::vcov(xmod))
          linpreds <- matrix((VGAM::model.matrix(xmod.dummy)) %*% newbeta, byrow = TRUE, ncol = (nlevels(imputations[[imp]][, targetCol]) - 1))
          denom <- 1 + rowSums(exp(linpreds))
          xfitted <- cbind(1 / denom, exp(linpreds) / denom)
        }
        if (noisy == TRUE) {
          print(summary(xmod))
        }

        # estimate parameters of substantive model
        if (smtype == "lm") {
          ymod <- lm(as.formula(smformula), imputations[[imp]])
          beta <- ymod$coef
          sigmasq <- summary(ymod)$sigma^2
          varcov <- vcov(ymod)
          outcomeModResVar <- (sigmasq * ymod$df) / rchisq(1, ymod$df)
          covariance <- (outcomeModResVar / sigmasq) * vcov(ymod)
          outcomeModBeta <- beta + MASS::mvrnorm(1, mu = rep(0, ncol(covariance)), Sigma = covariance)
          if (noisy == TRUE) {
            print(summary(ymod))
          }
        } else if ((smtype == "logistic") | (smtype == "brlogistic")) {
          if (smtype == "logistic") {
            ymod <- glm(as.formula(smformula), family = "binomial", imputations[[imp]])
          } else {
            ymod <- glm(as.formula(smformula),
              family = "binomial", imputations[[imp]],
              method = brglm2::brglmFit
            )
          }
          outcomeModBeta <- modPostDraw(ymod)
          if (noisy == TRUE) {
            print(summary(ymod))
          }
        } else if (smtype == "dtsam") {
          # split data to long form
          longData <- survival::survSplit(as.formula(smformula), data = imputations[[imp]], cut = cutPoints)
          # fit logistic model
          if (extraArgs$timeEffects == "factor") {
            dtsamFormula <- paste(colnames(imputations[[imp]])[dCol], "~-1+factor(tstart)+",
              strsplit(smformula, "~")[[1]][2],
              sep = ""
            )
          } else if (extraArgs$timeEffects == "linear") {
            # linear time effect
            dtsamFormula <- paste(colnames(imputations[[imp]])[dCol], "~tstart+",
              strsplit(smformula, "~")[[1]][2],
              sep = ""
            )
          } else {
            # quadratic effect of time
            dtsamFormula <- paste(colnames(imputations[[imp]])[dCol], "~tstart+I(tstart^2)+",
              strsplit(smformula, "~")[[1]][2],
              sep = ""
            )
          }
          ymod <- glm(as.formula(dtsamFormula), family = "binomial", data = longData)
          outcomeModBeta <- modPostDraw(ymod)
          if (noisy == TRUE) {
            print(summary(ymod))
          }
        } else if (smtype == "poisson") {
          ymod <- glm(as.formula(smformula), family = "poisson", imputations[[imp]])
          outcomeModBeta <- modPostDraw(ymod)
          if (noisy == TRUE) {
            print(summary(ymod))
          }
        } else if (smtype == "coxph") {
          ymod <- survival::coxph(as.formula(smformula), imputations[[imp]],
            control = survival::coxph.control(timefix = FALSE), model = TRUE
          )
          outcomeModBeta <- modPostDraw(ymod)
          ymod$coefficients <- outcomeModBeta
          basehaz <- survival::basehaz(ymod, centered = FALSE)[, 1]
          H0 <- basehaz[H0indices]
          if (noisy == TRUE) {
            print(summary(ymod))
          }
        } else if (smtype == "flexsurv") {
          if (cyclenum==1) {
            tryCatch({
              ymod <- flexsurv::flexsurvspline(as.formula(smformula), imputations[[imp]],
                                    k=extraArgs$k, scale="hazard")
            },error = function(e) {
              # Silently increment the failure counter
              flexsurvFailCount <<- flexsurvFailCount + 1
            })
          } else {
            # use previous estimates as initial values, to try and prevent non-convergence
            tryCatch(
              {
                # first attempt to fit model
                if (extraArgs$originalKnots==FALSE) {
                  ymod <- flexsurv::flexsurvspline(as.formula(smformula), imputations[[imp]],
                                                   k=extraArgs$k, scale="hazard",
                                                   inits=flexsurvEsts)
                } else {
                  ymod <- flexsurv::flexsurvspline(as.formula(smformula), imputations[[imp]],
                                                   scale="hazard",
                                                   knots=ymod$knots[2:(length(ymod$knots)-1)],
                                                   inits=flexsurvEsts)
                }
              },
              error = function(e) {
                # Silently increment the failure counter
                flexsurvFailCount <<- flexsurvFailCount + 1
              }
              )
          }
          # save parameter estimates to use as initial values in subsequent iterations
          # to help prevent convergence problems
          flexsurvEsts <- ymod$res.t[,1]
          outcomeModBeta <- as.numeric(flexsurv::normboot.flexsurvreg(ymod, B=1, raw=TRUE))
          # overwrite model estimates in fitted model object
          ymod$res.t[,1] <- outcomeModBeta
          if (noisy == TRUE) {
            print(ymod)
          }
        } else if (smtype == "weibull") {
          ymod <- survival::survreg(as.formula(smformula), data = imputations[[imp]], dist = "weibull")
          outcomeModBeta <- c(coef(ymod), log(ymod$scale)) +
            MASS::mvrnorm(1, mu = rep(0, ncol(vcov(ymod))), Sigma = vcov(ymod))
          weibullScale <- exp(utils::tail(outcomeModBeta, 1))
          outcomeModBeta <- utils::head(outcomeModBeta, -1)
          if (noisy == TRUE) {
            print(summary(ymod))
          }
        } else if (smtype == "compet") {
          for (cause in 1:numCauses) {
            ymod <- survival::coxph(as.formula(smformula[[cause]]), imputations[[imp]],
              control = survival::coxph.control(timefix = FALSE), model = TRUE
            )
            outcomeModBeta[[cause]] <- modPostDraw(ymod)
            ymod$coefficients <- outcomeModBeta[[cause]]
            basehaz <- survival::basehaz(ymod, centered = FALSE)[, 1]
            H0[[cause]] <- basehaz[H0indices[[cause]]]
            if (noisy == TRUE) {
              print(summary(ymod))
            }
          }
        } else if (smtype == "casecohort") {
          ymod <- survival::coxph(as.formula(smformula2), imputations[[imp]],
            control = survival::coxph.control(timefix = FALSE)
          )
          outcomeModBeta <- modPostDraw(ymod)

          cumhaz.denom.elements <- exp(model.matrix(as.formula(smformula), imputations[[imp]])[, -1] %*% outcomeModBeta)
          cumhaz.denom <- sapply(list.times, function(x) {
            sum(cumhaz.denom.elements[which(originaldata[, timeCol] >= x)] * subco.weight[which(originaldata[, timeCol] >= x)])
          })
          exp.func.denom <- cumsum(1 / cumhaz.denom)
          H0.fun <- stepfun(list.times, c(0, exp.func.denom))
          H0 <- H0.fun(originaldata[, timeCol])

          if (noisy == TRUE) {
            print(summary(ymod))
          }
        } else if (smtype == "nestedcc") {
          ymod <- survival::coxph(as.formula(smformula), imputations[[imp]],
            control = survival::coxph.control(timefix = FALSE)
          )
          outcomeModBeta <- modPostDraw(ymod)

          explan.matrix <- model.matrix(ymod)
          cumbasehaz.denom <- exp(matrix(outcomeModBeta, nrow = 1) %*% t(explan.matrix)) * originaldata[, nriskCol] / num.sampriskset
          cumbasehaz.denom <- ave(cumbasehaz.denom, originaldata[, setCol], FUN = sum)[originaldata[, dCol] == 1] # this is the denominator of the contribution to the cumulative baseline hazard at each event time
          cumbasehaz.t <- originaldata[, timeCol][originaldata[, dCol] == 1] # times to which the baseline cumulative hazards refer
          H0 <- unlist(lapply(originaldata[, timeCol], function(x) {
            sum((1 / cumbasehaz.denom)[cumbasehaz.t <= x])
          }))

          if (noisy == TRUE) {
            print(summary(ymod))
          }
        }

        if ((imp == 1) & (cyclenum == 1) & (var == 1)) {
          if (smtype == "compet") {
            totalCoefVec <- outcomeModBeta[[1]]
            for (cause in 2:numCauses) {
              totalCoefVec <- c(totalCoefVec, outcomeModBeta[[cause]])
            }
            smCoefIter <- array(0, dim = c(m, length(totalCoefVec), numit))
          } else {
            smCoefIter <- array(0, dim = c(m, length(outcomeModBeta), numit))
          }
        }

        if (var == length(partialVars)) {
          # then we have reached end of a cycle
          if (smtype == "compet") {
            totalCoefVec <- outcomeModBeta[[1]]
            for (cause in 2:numCauses) {
              totalCoefVec <- c(totalCoefVec, outcomeModBeta[[cause]])
            }
            smCoefIter[imp, , cyclenum] <- totalCoefVec
          } else {
            smCoefIter[imp, , cyclenum] <- outcomeModBeta
          }
        }

        # impute x, either directly where possibly, or using rejection sampling otherwise
        imputationNeeded <- (1:n)[r[, targetCol] == 0]

        if ((method[targetCol] == "logreg") | (method[targetCol] == "podds") | (method[targetCol] == "mlogit") |
          (method[targetCol] == "brlogreg")) {
          # directly sample
          if ((method[targetCol] == "logreg") | (method[targetCol] == "brlogreg")) {
            numberOutcomes <- 2
            fittedMean <- cbind(1 - xfitted, xfitted)
          } else {
            numberOutcomes <- nlevels(imputations[[imp]][, targetCol])
            fittedMean <- xfitted
          }

          outcomeDensCovDens <- array(dim = c(length(imputationNeeded), numberOutcomes), 0)

          for (xMisVal in 1:numberOutcomes) {
            if ((method[targetCol] == "logreg") | (method[targetCol] == "brlogreg")) {
              if (is.factor(imputations[[imp]][, targetCol]) == TRUE) {
                valToImpute <- levels(imputations[[imp]][, targetCol])[xMisVal]
              } else {
                valToImpute <- xMisVal - 1
              }
            } else {
              valToImpute <- levels(imputations[[imp]][, targetCol])[xMisVal]
            }
            imputations[[imp]][imputationNeeded, targetCol] <- valToImpute

            # update passive variables
            imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)

            if (smtype == "lm") {
              outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]]) %*% outcomeModBeta
              deviation <- imputations[[imp]][imputationNeeded, outcomeCol] - outmodxb[imputationNeeded]
              outcomeDens <- dnorm(deviation, mean = 0, sd = outcomeModResVar^0.5)
            } else if ((smtype == "logistic") | (smtype == "brlogistic")) {
              outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]]) %*% outcomeModBeta
              prob <- plogis(outmodxb[imputationNeeded])
              outcomeDens <- prob * imputations[[imp]][imputationNeeded, outcomeCol] + (1 - prob) * (1 - imputations[[imp]][imputationNeeded, outcomeCol])
            } else if (smtype == "dtsam") {
              outcomeDens <- dtsamOutcomeDens(
                imputations[[imp]], extraArgs$timeEffects, outcomeModBeta, nTimePoints,
                smformula, timeCol, dCol
              )[imputationNeeded]
            } else if (smtype == "poisson") {
              outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]]) %*% outcomeModBeta
              outcomeDens <- dpois(imputations[[imp]][imputationNeeded, outcomeCol], exp(outmodxb[imputationNeeded]))
            } else if (smtype == "weibull") {
              outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]]) %*% outcomeModBeta
              # weibull survival function
              outcomeDens <- (1 - d[imputationNeeded]) * (1 - survival::psurvreg(imputations[[imp]][imputationNeeded, timeCol], mean = outmodxb[imputationNeeded], scale = weibullScale)) +
                d[imputationNeeded] * (survival::dsurvreg(imputations[[imp]][imputationNeeded, timeCol], mean = outmodxb[imputationNeeded], scale = weibullScale))
            } else if ((smtype == "coxph") | (smtype == "casecohort")) {
              outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]])
              outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              outcomeDens <- exp(-H0[imputationNeeded] * exp(outmodxb[imputationNeeded])) * (exp(outmodxb[imputationNeeded])^d[imputationNeeded])
            } else if (smtype == "flexsurv") {
              survEst <- summary(ymod, newdata=imputations[[imp]][imputationNeeded,], type="survival",
                                 ci=FALSE, t=imputations[[imp]][imputationNeeded,timeCol],
                                 cross=FALSE, tidy=TRUE)
              survEst <- as.matrix(survEst)[,"est"]
              hazEst <- summary(ymod, newdata=imputations[[imp]][imputationNeeded,], type="hazard",
                                ci=FALSE, t=imputations[[imp]][imputationNeeded,timeCol],
                                cross=FALSE, tidy=TRUE)
              hazEst <- as.matrix(hazEst)[,"est"]
              outcomeDens <- survEst*(hazEst^imputations[[imp]][imputationNeeded,dCol])
            } else if (smtype == "nestedcc") {
              outmodxb <- model.matrix(as.formula(smformula2), imputations[[imp]])
              outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              outcomeDens <- exp(-H0[imputationNeeded] * exp(outmodxb[imputationNeeded])) * (exp(outmodxb[imputationNeeded])^d[imputationNeeded])
            } else if (smtype == "compet") {
              outcomeDens <- rep(1, length(imputationNeeded))
              for (cause in 1:numCauses) {
                outmodxb <- model.matrix(linpred[[cause]], imputations[[imp]])
                outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta[[cause]])
                outcomeDens <- outcomeDens * exp(-H0[[cause]][imputationNeeded] * exp(outmodxb[imputationNeeded])) * (exp(outmodxb[imputationNeeded])^(d[imputationNeeded] == cause))
              }
            }
            outcomeDensCovDens[, xMisVal] <- outcomeDens * fittedMean[imputationNeeded, xMisVal]
          }

          if(!is.null(restrictions) && var %in% index_restrictions){
            len_restrictions = length(restrictions_work[[var]])

            ## For each information source, a separate step is run
            for(i in 1:len_restrictions){
              restrictions_work_vari = restrictions_work[[var]][[i]]
              res_var = restrictions_work_vari[[1]]; res_cond = restrictions_work_vari[[2]]; res_opts = restrictions_work_vari[[3]]

              ## If a column is not compatible with the restricting information, the probability for that column is set to 0
              outcomeDensCovDens[which(originaldata[imputationNeeded, res_var] == res_cond), !(levels(imputations[[imp]][, targetCol]) %in% res_opts)] = 0
            }
          }

          directImpProbs <- outcomeDensCovDens / rowSums(outcomeDensCovDens)

          if ((method[targetCol] == "logreg") | (method[targetCol] == "brlogreg")) {
            directImpProbs <- directImpProbs[, 2]
            if (is.factor(imputations[[imp]][, targetCol]) == TRUE) {
              imputations[[imp]][imputationNeeded, targetCol] <- levels(imputations[[imp]][, targetCol])[1]
              imputations[[imp]][imputationNeeded, targetCol][rbinom(length(imputationNeeded), 1, directImpProbs) == 1] <- levels(imputations[[imp]][, targetCol])[2]
            } else {
              imputations[[imp]][imputationNeeded, targetCol] <- rbinom(length(imputationNeeded), 1, directImpProbs)
            }
          } else {
            imputations[[imp]][imputationNeeded, targetCol] <- levels(imputations[[imp]][, targetCol])[apply(directImpProbs, 1, catdraw)]
          }

          # update passive variables
          imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
        } else {
          # use rejection sampling
          # first draw for all subjects who need imputing, using a small number of attempts
          firstTryLimit <- 25
          j <- 1

          while ((length(imputationNeeded) > 0) & (j < firstTryLimit)) {
            # sample from covariate model
            if ((method[targetCol] == "norm") | (method[targetCol] == "latnorm")) {
              imputations[[imp]][imputationNeeded, targetCol] <- rnorm(length(imputationNeeded), xfitted[imputationNeeded], newsigmasq^0.5)
            } else if (method[targetCol] == "poisson") {
              imputations[[imp]][imputationNeeded, targetCol] <- rpois(length(imputationNeeded), xfitted[imputationNeeded])
            }

            # update passive variables
            imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)

            # accept/reject
            uDraw <- runif(length(imputationNeeded))
            if (smtype == "lm") {
              outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]]) %*% outcomeModBeta
              deviation <- imputations[[imp]][imputationNeeded, outcomeCol] - outmodxb[imputationNeeded]
              reject <- 1 * (log(uDraw) > -(deviation^2) / (2 * array(outcomeModResVar, dim = c(length(imputationNeeded), 1))))
            } else if ((smtype == "logistic") | (smtype == "brlogistic")) {
              outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]]) %*% outcomeModBeta
              prob <- plogis(outmodxb[imputationNeeded])
              prob <- prob * imputations[[imp]][imputationNeeded, outcomeCol] + (1 - prob) * (1 - imputations[[imp]][imputationNeeded, outcomeCol])
              reject <- 1 * (uDraw > prob)
            } else if (smtype == "dtsam") {
              prob <- dtsamOutcomeDens(
                imputations[[imp]], extraArgs$timeEffects, outcomeModBeta, nTimePoints,
                smformula, timeCol, dCol
              )[imputationNeeded]
              reject <- 1 * (uDraw > prob)
            } else if (smtype == "poisson") {
              outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]]) %*% outcomeModBeta
              prob <- dpois(imputations[[imp]][imputationNeeded, outcomeCol], exp(outmodxb[imputationNeeded]))
              reject <- 1 * (uDraw > prob)
            } else if (smtype == "weibull") {
              outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]]) %*% outcomeModBeta
              s_t <- 1 - survival::psurvreg(imputations[[imp]][imputationNeeded, timeCol], mean = outmodxb[imputationNeeded], scale = weibullScale)
              prob <- -exp(1) * log(s_t) * s_t
              prob <- d[imputationNeeded] * prob + (1 - d[imputationNeeded]) * s_t
              reject <- 1 * (uDraw > prob)
            } else if ((smtype == "coxph") | (smtype == "casecohort")) {
              outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]])
              outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              s_t <- exp(-H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))
              prob <- exp(1 + outmodxb[imputationNeeded] - (H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))) * H0[imputationNeeded]
              prob <- d[imputationNeeded] * prob + (1 - d[imputationNeeded]) * s_t
              reject <- 1 * (uDraw > prob)
            } else if (smtype == "flexsurv") {
              survEst <- summary(ymod, newdata=imputations[[imp]][imputationNeeded,], type="survival",
                                 ci=FALSE, t=imputations[[imp]][imputationNeeded,timeCol],
                                 cross=FALSE, tidy=TRUE)
              survEst <- as.matrix(survEst)[,"est"]
              hazEst <- summary(ymod, newdata=imputations[[imp]][imputationNeeded,], type="hazard",
                                ci=FALSE, t=imputations[[imp]][imputationNeeded,timeCol],
                                cross=FALSE, tidy=TRUE)
              hazEst <- as.matrix(hazEst)[,"est"]
              cumhazEst <- summary(ymod, newdata=imputations[[imp]][imputationNeeded,], type="cumhaz",
                                ci=FALSE, t=imputations[[imp]][imputationNeeded,timeCol],
                                cross=FALSE, tidy=TRUE)
              cumhazEst <- as.matrix(cumhazEst)[,"est"]
              prob <- imputations[[imp]][imputationNeeded,dCol] * (survEst*exp(1)*cumhazEst) +
                (1 - imputations[[imp]][imputationNeeded,dCol]) * survEst
              reject <- 1 * (uDraw > prob)
            } else if (smtype == "nestedcc") {
              outmodxb <- model.matrix(as.formula(smformula2), imputations[[imp]])
              outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              s_t <- exp(-H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))
              prob <- exp(1 + outmodxb[imputationNeeded] - (H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))) * H0[imputationNeeded]
              prob <- d[imputationNeeded] * prob + (1 - d[imputationNeeded]) * s_t
              reject <- 1 * (uDraw > prob)
            } else if (smtype == "compet") {
              prob <- rep(1, length(imputationNeeded))
              for (cause in 1:numCauses) {
                outmodxb <- model.matrix(linpred[[cause]], imputations[[imp]])
                outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta[[cause]])
                prob <- prob * exp(-H0[[cause]][imputationNeeded] * exp(outmodxb[imputationNeeded])) * (H0[[cause]][imputationNeeded] * exp(1 + outmodxb[imputationNeeded]))^(d[imputationNeeded] == cause)
              }
              reject <- 1 * (uDraw > prob)
            }
            imputationNeeded <- imputationNeeded[reject == 1]

            j <- j + 1
          }

          # now, for those remaining, who must have low acceptance probabilities, sample by subject
          for (i in imputationNeeded) {
            tempData <- imputations[[imp]][i, ]
            tempData <- tempData[rep(1, rjlimit), ]
            if (method[targetCol] == "norm") {
              tempData[, targetCol] <- rnorm(rjlimit, xfitted[i], newsigmasq^0.5)
            } else if ((method[targetCol] == "logreg") | (method[targetCol] == "brlogreg")) {
              tempData[, targetCol] <- rbinom(rjlimit, size = 1, xfitted[i])
            } else if (method[targetCol] == "poisson") {
              tempData[, targetCol] <- rpois(rjlimit, xfitted[i])
            } else if (method[targetCol] == "latnorm") {
              tempData[, targetCol] <- rnorm(rjlimit, xfitted[i], newsigmasq[i]^0.5)
            }

            # passively impute
            tempData <- updatePassiveVars(tempData, method, passiveVars)

            # accept reject
            uDraw <- runif(rjlimit)
            if (smtype == "lm") {
              outmodxb <- model.matrix(as.formula(smformula), tempData) %*% outcomeModBeta
              deviation <- tempData[, outcomeCol] - outmodxb
              reject <- 1 * (log(uDraw) > -(deviation^2) / (2 * array(outcomeModResVar, dim = c(rjlimit, 1))))
            } else if ((smtype == "logistic") | (smtype == "brlogistic")) {
              outmodxb <- model.matrix(as.formula(smformula), tempData) %*% outcomeModBeta
              prob <- plogis(outmodxb)
              prob <- prob * tempData[, outcomeCol] + (1 - prob) * (1 - tempData[, outcomeCol])
              reject <- 1 * (uDraw > prob)
            } else if (smtype == "dtsam") {
              prob <- dtsamOutcomeDens(
                tempData, extraArgs$timeEffects, outcomeModBeta, nTimePoints,
                smformula, timeCol, dCol
              )
              reject <- 1 * (uDraw > prob)
            } else if (smtype == "poisson") {
              outmodxb <- model.matrix(as.formula(smformula), tempData) %*% outcomeModBeta
              prob <- dpois(tempData[, outcomeCol], exp(outmodxb))
              reject <- 1 * (uDraw > prob)
            } else if (smtype == "weibull") {
              outmodxb <- model.matrix(as.formula(smformula), tempData) %*% outcomeModBeta
              s_t <- 1 - survival::psurvreg(tempData[, timeCol], mean = outmodxb, scale = weibullScale)
              if (d[i] == 1) {
                prob <- -exp(1) * log(s_t) * s_t
                # the following line fixes a numerical error which occurs if s_t=0 for any draws
                prob[is.na(prob)] <- 0
              } else {
                prob <- s_t
              }
              reject <- 1 * (uDraw > prob)
            } else if ((smtype == "coxph") | (smtype == "casecohort")) {
              outmodxb <- model.matrix(as.formula(smformula), tempData)
              outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              s_t <- exp(-H0[i] * exp(outmodxb))
              prob <- exp(1 + outmodxb - (H0[i] * exp(outmodxb))) * H0[i]
              prob <- d[i] * prob + (1 - d[i]) * s_t
              reject <- 1 * (uDraw > prob)
            } else if (smtype == "flexsurv") {
              survEst <- summary(ymod, newdata=tempData, type="survival",
                                 ci=FALSE, t=tempData[,timeCol],
                                 cross=FALSE, tidy=TRUE)
              survEst <- as.matrix(survEst)[,"est"]
              hazEst <- summary(ymod, newdata=tempData, type="hazard",
                                ci=FALSE, t=tempData[,timeCol],
                                cross=FALSE, tidy=TRUE)
              hazEst <- as.matrix(hazEst)[,"est"]
              cumhazEst <- summary(ymod, newdata=tempData, type="cumhaz",
                                   ci=FALSE, t=tempData[,timeCol],
                                   cross=FALSE, tidy=TRUE)
              cumhazEst <- as.matrix(cumhazEst)[,"est"]
              prob <- imputations[[imp]][i,dCol] * (survEst*exp(1)*cumhazEst) +
                (1 - imputations[[imp]][i,dCol]) * survEst
              reject <- 1 * (uDraw > prob)
            } else if (smtype == "nestedcc") {
              outmodxb <- model.matrix(as.formula(smformula2), tempData)
              outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              s_t <- exp(-H0[i] * exp(outmodxb))
              prob <- exp(1 + outmodxb - (H0[i] * exp(outmodxb))) * H0[i]
              prob <- d[i] * prob + (1 - d[i]) * s_t
              reject <- 1 * (uDraw > prob)
            } else if (smtype == "compet") {
              prob <- rep(1, rjlimit)
              for (cause in 1:numCauses) {
                outmodxb <- model.matrix(linpred[[cause]], tempData)
                outmodxb <- as.matrix(outmodxb[, 2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta[[cause]])
                prob <- prob * exp(-H0[[cause]][i] * exp(outmodxb)) * (H0[[cause]][i] * exp(1 + outmodxb))^(d[i] == cause)
              }
              reject <- 1 * (uDraw > prob)
            }
            if (sum(reject) < rjlimit) {
              imputations[[imp]][i, targetCol] <- tempData[reject == 0, targetCol][1]
            } else {
              rjFailCount <- rjFailCount + 1
            }
          }
          # update passive variables
          imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
        }
      }
      }

      # imputations of missing outcomes, if present (using proper imputation), for regression and logistic
      # substantive models
      if ((smtype == "lm") | (smtype == "logistic") | (smtype == "brlogistic")) {
        if (sum(r[, outcomeCol]) < n) {
          imputationNeeded <- (1:n)[r[, outcomeCol] == 0]
          # estimate parameters of substantive model using those with outcomes observed
          if (smtype == "lm") {
            ymod <- lm(as.formula(smformula), imputations[[imp]][r[, outcomeCol] == 1, ])
            beta <- ymod$coef
            sigmasq <- summary(ymod)$sigma^2
            varcov <- vcov(ymod)
            outcomeModResVar <- (sigmasq * ymod$df) / rchisq(1, ymod$df)
            covariance <- (outcomeModResVar / sigmasq) * vcov(ymod)
            outcomeModBeta <- beta + MASS::mvrnorm(1, mu = rep(0, ncol(covariance)), Sigma = covariance)
            outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]]) %*% outcomeModBeta
            imputations[[imp]][imputationNeeded, outcomeCol] <- rnorm(length(imputationNeeded), outmodxb[imputationNeeded], sigmasq^0.5)
          } else if ((smtype == "logistic") | (smtype == "brlogistic")) {
            if (smtype == "logistic") {
              ymod <- glm(as.formula(smformula), family = "binomial", imputations[[imp]][r[, outcomeCol] == 1, ])
            } else {
              ymod <- glm(as.formula(smformula),
                family = "binomial", imputations[[imp]][r[, outcomeCol] == 1, ],
                method = brglm2::brglmFit
              )
            }
            outcomeModBeta <- modPostDraw(ymod)
            outmodxb <- model.matrix(as.formula(smformula), imputations[[imp]]) %*% outcomeModBeta
            prob <- plogis(outmodxb[imputationNeeded])
            imputations[[imp]][imputationNeeded, outcomeCol] <- rbinom(length(imputationNeeded), 1, prob)
          }
        }
      } else if (smtype == "flexsurv") {
        if (extraArgs$imputeTimes==TRUE) {
          # impute censored times
          if (cyclenum==1) {
            tryCatch({
              ymod <- flexsurv::flexsurvspline(as.formula(smformula), imputations[[imp]],
                                               k=extraArgs$k, scale="hazard")
            },error = function(e) {
              # Silently increment the failure counter
              flexsurvFailCount <<- flexsurvFailCount + 1
            })
          } else {
            # use previous estimates as initial values, to try and prevent non-convergence
            tryCatch(
              {
                if (extraArgs$originalKnots==FALSE) {
                  ymod <- flexsurv::flexsurvspline(as.formula(smformula), imputations[[imp]],
                                                 k=extraArgs$k, scale="hazard",
                                                 inits=flexsurvEsts)
                } else {
                  # use knots based on fit to original dataset where censored times had not been imputed
                  ymod <- flexsurv::flexsurvspline(as.formula(smformula), imputations[[imp]],
                                                   scale="hazard",
                                                   knots=ymod$knots[2:(length(ymod$knots)-1)],
                                                   inits=flexsurvEsts)
                }
              },
              error = function(e) {
                flexsurvFailCount <<- flexsurvFailCount + 1
              }
            )
          }
          flexsurvEsts <- ymod$res.t[,1]
          outcomeModBeta <- as.numeric(flexsurv::normboot.flexsurvreg(ymod, B=1, raw=TRUE))

          if (exists("smCoefIter")==FALSE) {
            smCoefIter <- array(0, dim = c(m, length(outcomeModBeta), numit))
          }
          smCoefIter[imp, , cyclenum] <- outcomeModBeta

          # overwrite model estimates in fitted model object
          ymod$res.t[,1] <- outcomeModBeta
          if (noisy == TRUE) {
            print(ymod)
          }

          # impute censored times
          if (is.null(extraArgs$censtime)) {
            # impute all censored times
            timeImp <- simulate(ymod, nsim=1, newdata=imputations[[imp]][d==0,],
                              start=originaldata[d==0,timeCol])
            imputations[[imp]][d==0,timeCol] <- timeImp$time_1
            imputations[[imp]][,dCol] <- 1
          } else {
            # impute but still impose censoring at specified value(s)
            timeImp <- simulate(ymod, nsim=1, newdata=imputations[[imp]][d==0,],
                                start=originaldata[d==0,timeCol],
                                censtime=extraArgs$censtime)
            imputations[[imp]][d==0,timeCol] <- timeImp$time_1
            imputations[[imp]][d==0,dCol] <- timeImp$event_1
          }
        }
      }
    }
  }

  if (rjFailCount > 0) {
    warning(paste("Rejection sampling failed ", rjFailCount, " times (across all variables, iterations, and imputations). You may want to increase the rejection sampling limit.", sep = ""))
  }
  if (flexsurvFailCount > 0) {
    warning(paste("Flexsurv fit failed ", flexsurvFailCount, " times (across all variables, iterations, and imputations). See documentation for more details.", sep = ""))
  }

  # Added smformula and smtype to metadata, and make "smcfcs class"
  res <- list(
    impDatasets = imputations,
    smCoefIter = smCoefIter,
    smInfo = list("smtype" = smtype, "smformula" = smformula),
    extraArgs = extraArgs # For plot of dtsam
  )
  class(res) <- "smcfcs"

  return(res)
}

updatePassiveVars <- function(data, method, passivecols) {
  for (i in passivecols) {
    data[, i] <- with(data, eval(parse(text = method[i])))
  }
  data
}

sumna <- function(x) {
  sum(is.na(x) == FALSE)
}

# returns first non missing entry of x
firstnonna <- function(x) {
  x[is.na(x) == FALSE][1]
}

catdraw <- function(prob) {
  (1:length(prob))[rmultinom(1, size = 1, prob = prob) == 1]
}

modPostDraw <- function(modobj) {
  beta <- modobj$coef
  varcov <- vcov(modobj)
  beta + MASS::mvrnorm(1, mu = rep(0, ncol(varcov)), Sigma = varcov)
}

# a helper function which calculates the substantive model probability value
# for each individual for the discrete time survival (logistic) substantive model
dtsamOutcomeDens <- function(inputData, timeEffects, outcomeModBeta, nTimePoints, smformula, timeCol, dCol) {
  inputDataN <- dim(inputData)[1]

  # calculate covariate effects
  modelMatrix <- model.matrix(
    as.formula(paste("~", strsplit(smformula, "~")[[1]][2], sep = "")),
    inputData
  )
  # remove intercept column
  modelMatrix <- modelMatrix[,2:ncol(modelMatrix)]

  if (timeEffects == "factor") {
    # first add in time effects on log odds scale
    outmodxb <- matrix(outcomeModBeta[1:nTimePoints], nrow = inputDataN, ncol = nTimePoints, byrow = TRUE)
    covXbEffects <- modelMatrix %*% utils::tail(outcomeModBeta, length(outcomeModBeta) - nTimePoints)
  } else if (timeEffects == "linear") {
    # linear time
    # first add in time effects on log odds scale
    outmodxb <- outcomeModBeta[1] +
      matrix(outcomeModBeta[2] * ((1:nTimePoints) - 1), nrow = inputDataN, ncol = nTimePoints, byrow = TRUE)
    # calculate covariate effects
    covXbEffects <- modelMatrix %*% utils::tail(outcomeModBeta, length(outcomeModBeta) - 2)
  } else {
    # quadratic time
    # first add in time effects on log odds scale
    outmodxb <- outcomeModBeta[1] +
      matrix(outcomeModBeta[2] * ((1:nTimePoints) - 1) + outcomeModBeta[3] * (((1:nTimePoints) - 1)^2), nrow = inputDataN, ncol = nTimePoints, byrow = TRUE)
    # calculate covariate effects
    covXbEffects <- modelMatrix %*% utils::tail(outcomeModBeta, length(outcomeModBeta) - 3)
  }

  # add in covariate effects
  outmodxb <- outmodxb + matrix(covXbEffects, nrow = inputDataN, ncol = nTimePoints)
  # prob is matrix of conditional probabilities/hazard of event in each period
  prob <- expit(outmodxb)
  logSurvProb <- log(1 - prob)
  logSurvProbCumSum <- cbind(rep(0, inputDataN), t(apply(logSurvProb, 1, cumsum)))

  # create vector of last time point each person survived to, +1
  lastSurvPlusOne <- 1 + inputData[, timeCol] - inputData[, dCol]
  # create vector of log of probability of survival to time when each person actually survived to
  logSurvProbIndividual <- logSurvProbCumSum[cbind(1:inputDataN, lastSurvPlusOne)]
  # return vector of outcome density values
  exp(logSurvProbIndividual + inputData[, dCol] * log(prob[cbind(1:inputDataN, inputData[, timeCol])]))
}

# this helper function is used in a check for smcfcs.flexsurv
# it is used to extract names of variables which are being modelled
# using time-varying effects
extract_gamma_values <- function(s) {
  # Define the regex pattern: "gamma" followed by digits, then capture what's inside the parentheses.
  pattern <- "gamma\\d+\\(([^)]+)\\)"

  # Find all matches with capture groups enabled
  m <- gregexpr(pattern, s, perl = TRUE)

  # Get capture group positions and lengths
  capture_starts <- attr(m[[1]], "capture.start")
  capture_lengths <- attr(m[[1]], "capture.length")

  # If no matches are found, return an empty character vector
  if (capture_starts[1] == -1) {
    return(character(0))
  }

  # Extract the captured text for each match
  extracted <- substring(s, capture_starts, capture_starts + capture_lengths - 1)
  return(unique(extracted))
}

