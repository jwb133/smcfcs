#' Substantive model compatible fully conditional specification imputation of covariates and
#' event times using flexible parametric survival models
#'
#' Multiply imputes missing covariate values and event times using substantive model compatible
#' fully conditional specification with flexible parametric survival models.
#'
#' This version of \code{smcfcs} is for time-to-event outcomes which are modelled
#' using a flexible parametric proportional hazards survival model. The model is
#' fitted using the \code{\link[flexsurv]{flexsurvspline}} function in the
#' \pkg{flexsurv} package. Specifically it fits models using the hazard scale. The
#' flexibility of the model can be changed by modifying the k argument, which
#' specifies the number of knots.
#'
#' \code{\link[flexsurv]{flexsurvspline}} sometimes fails during model fitting.
#' If/when this occurs, \code{smcfcs.flexsurv} uses takes a posterior draw based
#' on the model fit from the preceding iteration, and a warning is printed at
#' the end of the \code{smcfcs.flexsurv} run detailing how many times it occurred.
#'
#' @author Jonathan Bartlett \email{jonathan.bartlett1@@lshtm.ac.uk}
#'
#' @param smformula A formula of the form "Surv(t,d)~x+z"
#' @param k Number of knots to use in the flexible parametric survival model
#'
#' @inheritParams smcfcs
#' @example data-raw/flexsurv_example.r
#' @export
smcfcs.flexsurv <- function(originaldata, smformula, k=2, method, predictorMatrix = NULL, m = 5, numit = 10, rjlimit = 1000, noisy = FALSE, errorProneMatrix = NULL) {
  smcfcs.core(originaldata=originaldata,
              smtype = "flexsurv",
              smformula = smformula,
              method = method,
              predictorMatrix = predictorMatrix,
              m = m,
              numit = numit,
              rjlimit = rjlimit,
              noisy = noisy,
              errorProneMatrix = errorProneMatrix,
              k=k)
}
