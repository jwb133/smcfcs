#' Multiple imputation for confounders in propensity score analysis.
#'
#' Multiply imputes missing covariate values using an adaptation of the substantive
#' model compatible fully conditional specification (SMC-FCS) approach to multiple imputation.
#'
#' \code{mips} imputes missing values in confounders prior to propensity score analysis.
#'
#' @param omtype A string specifying the type of outcome model. Possible
#' values are \code{"lm"}, \code{"logistic"}, \code{"poisson"}, \code{"coxph"}
#'  and \code{"compet"}.
#' @param omformula The formula of the outcome model. For \code{"coxph"} outcome
#' models the left hand side should be of the form \code{"Surv(t,d)"}. For \code{"compet"}
#' outcome models, a list should be passed consisting of the Cox models
#' for each cause of failure.
#' @param psformula The model formula for the propensity score model.
#'
#' @return A list containing the imputed datasets.
#'
#' @inheritParams smcfcs
#'
#' @example data-raw/mips-examples.r

#' @export
mips <- function(originaldata,omtype,omformula,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,psformula) {
  if (missing(psformula)) {
    stop("You must specify a propensity score model formula")
  }
  smcfcs.core(originaldata,omtype,omformula,method,predictorMatrix,m,numit,rjlimit,noisy,psformula)
}
