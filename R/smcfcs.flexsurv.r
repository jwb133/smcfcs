#' Substantive model compatible fully conditional specification imputation of covariates and
#' event times using flexible parametric survival models
#'
#' Multiply imputes missing covariate values and event times using substantive model compatible
#' fully conditional specification with a Royston-Parmar flexible parametric survival model.
#'
#' This version of \code{smcfcs} is for time-to-event outcomes which are modelled
#' using a flexible parametric proportional hazards survival model, as proposed
#' by Royston and Parmar (2002). The model is
#' fitted using the \code{\link[flexsurv]{flexsurvspline}} function in the
#' \pkg{flexsurv} package. Specifically it fits models using the hazard scale. The
#' flexibility of the model can be changed by modifying the k argument, which
#' specifies the number of knots.
#'
#' If desired, \code{smcfcs.flexsurv} can be used to impute event times for individuals
#' who are originally censored, by specifying \code{imputeTimes=TRUE}. In the resulting
#' imputed datasets every individual will have an event time and the event indicator will
#' be one for all. Alternatively, you can impute censored times, but setting a larger
#' potential censoring time, which is either a common value used for all or a vector of times,
#' by using the \code{censtime} argument. If some individuals have their time-to-event
#' outcome completely missing and you want to impute this, they should have a time of zero
#' and the event indicator set to zero.
#'
#' \code{smcfcs.flexsurv} will not let you impute using norm, latnorm or poisson methods
#' for variables that are allowed to have time-varying effects, because the usual
#' rejection sampling bound used by smcfcs is not valid in this setting.
#'
#' \code{\link[flexsurv]{flexsurvspline}} sometimes fails during model fitting.
#' If/when this occurs, \code{smcfcs.flexsurv} takes a posterior draw based
#' on the model fit from the preceding iteration, and a warning is printed at
#' the end of the \code{smcfcs.flexsurv} run detailing how many times it occurred.
#'
#' @author Jonathan Bartlett \email{jonathan.bartlett1@@lshtm.ac.uk}
#'
#' @param smformula A formula of the form "Surv(t,d)~x+z"
#' @param k Number of knots to use in the flexible parametric survival model
#' @param imputeTimes If set to TRUE, \code{smcfcs.flexsurv} will impute
#' censored survival times, as well as any missing covariates
#' @param censtime Value(s) to use for censoring of imputed event times. If
#' a vector, it should be of length equal to the number of original censored
#' individuals
#' @param originalKnots If imputing censored event times, setting
#' \code{originalKnots=TRUE} means the automatically chosen knot locations
#' from the model fitted to the observed times are used throughout. If \code{FALSE},
#' knots are chosen automatically at each iteration by \code{flexsurvspline}
#' based on the current observed+imputed event times, according to the chosen
#' value of \code{k}.
#' @param ... Additional arguments to pass on to \link[smcfcs]{smcfcs}
#'
#' @inheritParams smcfcs
#' @example data-raw/flexsurv_example.r
#'
#' @references Royston P, Parmar MKB. Flexible parametric proportional-hazards
#' and proportional-odds models for censored survival data, with application
#' to prognostic modelling and estimation of treatment effects.
#' Statistics in Medicine 2002; 21(15): 2175-2197. \doi{doi:10.1002/sim.1203}
#'
#' @export
smcfcs.flexsurv <- function(originaldata, smformula, method, k=2, imputeTimes=FALSE,
                    censtime = NULL,
                    originalKnots = TRUE,
                    ...) {

  if (!requireNamespace("flexsurv", quietly = TRUE)) {
    stop("The package 'flexsurv' is required for this function. Please install it.", call. = FALSE)
  }

  smcfcs.core(originaldata=originaldata,
              smtype = "flexsurv",
              smformula = smformula,
              method = method,
              k=k,
              imputeTimes=imputeTimes,
              censtime=censtime,
              originalKnots=originalKnots,
              ...)
}
