
#' Simulated example data with continuous outcome and quadratic covariate effects
#'
#' A dataset containing simulated data where the outcome depends quadratically
#' on a partially observed covariate.
#'
#' @format A data frame with 1000 rows and 5 variables:
#' \describe{
#'   \item{y}{Continuous outcome}
#'   \item{z}{Fully observed covariate, with linear effect on outcome}
#'   \item{x}{Partially observed normally distributed covariate, with quadratic effect on outcome}
#'   \item{xsq}{The square of x, which thus has missing values also}
#'   \item{v}{An auxiliary variable (i.e. not contained in the substantive model)}
#' }
#'
"ex_linquad"

#' Simulated example data with continuous outcome and interaction between two partially observed covariates
#'
#' A dataset containing simulated data where the outcome depends on both main
#' effects and interaction of two partially observed covariates.
#'
#' @format A data frame with 1000 rows and 4 variables:
#' \describe{
#'   \item{y}{Continuous outcome}
#'   \item{x1}{Partially observed normally distributed covariate}
#'   \item{x2}{Partially observed binary covariate}
#' }
#'
"ex_lininter"

#' Simulated example data with binary outcome and quadratic covariate effects
#'
#' A dataset containing simulated data where the binary outcome depends quadratically
#' on a partially observed covariate.
#'
#' @format A data frame with 1000 rows and 5 variables:
#' \describe{
#'   \item{y}{Binary outcome}
#'   \item{z}{Fully observed covariate, with linear effect on outcome (on log odds scale)}
#'   \item{x}{Partially observed normally distributed covariate, with quadratic effect on outcome (on log odds scale)}
#'   \item{xsq}{The square of x, which thus has missing values also}
#'   \item{v}{An auxiliary variable (i.e. not contained in the substantive model)}
#' }
#'
"ex_logisticquad"

#' Simulated example data with count outcome, modelled using Poisson regression
#'
#' A dataset containing simulated data where the count outcome depends on two
#' covariates, x and z, with missing values in x. The substantive model is
#' Poisson regression.
#'
#' @format A data frame with 1000 rows and 3 variables:
#' \describe{
#'   \item{y}{Count outcome}
#'   \item{z}{Fully observed covariate, with linear effect on outcome}
#'   \item{x}{Partially observed normally distributed covariate, with linear effect on outcome}
#' }
#'
"ex_poisson"

#' Simulated example data with time to event outcome and quadratic covariate effects
#'
#' A dataset containing simulated data where a time to event outcome depends quadratically
#' on a partially observed covariate.
#'
#' @format A data frame with 1000 rows and 6 variables:
#' \describe{
#'   \item{t}{Time to event or censoring}
#'   \item{d}{Binary indicator of whether event occurred or individual was censored}
#'   \item{z}{Fully observed covariate, with linear effect on outcome (on log hazard scale)}
#'   \item{x}{Partially observed normally distributed covariate, with quadratic effect on outcome (on log hazard scale)}
#'   \item{xsq}{The square of x, which thus has missing values also}
#'   \item{v}{An auxiliary variable (i.e. not contained in the substantive model)}
#' }
#'
"ex_coxquad"

#' Simulated example data with competing risks outcome and partially observed covariates
#'
#' A dataset containing simulated competing risks data. There are two competing risks, and
#' some times are also censored.
#'
#' @format A data frame with 1000 rows and 4 variables:
#' \describe{
#'   \item{t}{Time to event or censoring}
#'   \item{d}{Indicator of whether event 1 occurred (d=1), event 2 occurred (d=2) or individual was censored (d=0)}
#'   \item{x1}{Partially observed binary covariate, with linear effects on log competing risk hazards}
#'   \item{x2}{Partially observed normally distributed (conditional on x1) covariate, with linear effects
#'   on log competing risk hazards}
#' }
#'
"ex_compet"

#' Simulated case cohort data
#'
#' A dataset containing simulated case cohort data, where the sub-cohort was a 10\% random sample of the full cohort.
#'
#' @format A data frame with 1571 rows and 7 variables:
#' \describe{
#'   \item{t}{Time to event or censoring}
#'   \item{d}{Indicator of whether event 1 occurred (d=1), or not (d=0)}
#'   \item{x}{Partially observed continuous covariate}
#'   \item{z}{Fully observed covariate}
#'   \item{in.subco}{A binary indicator of whether the subject is in the sub-cohort}
#'   \item{id}{An id variable}
#'   \item{entertime}{The entry time variable to be used in the analysis}
#' }
#'
"ex_cc"

#' Simulated nested case-control data
#'
#' A dataset containing simulated nested case-control data.
#'
#' @format A data frame with 728 rows and 8 variables:
#' \describe{
#'   \item{t}{Time to event or censoring}
#'   \item{d}{Indicator of whether event 1 occurred (d=1), or not (d=0)}
#'   \item{x}{Partially observed binary covariate}
#'   \item{z}{Fully observed covariate}
#'   \item{id}{An id variable}
#'   \item{numrisk}{Number of patients at risk at time of case's event}
#'   \item{setno}{The case-control set number}
#'   \item{case}{Binary indicator of case (=1) or control (=0)}
#' }
#'
"ex_ncc"

#' Simulated discrete time survival data set
#'
#' A dataset containing simulated discrete time survival data.
#'
#' @format A data frame with 1000 rows and 8 variables:
#' \describe{
#'   \item{x1}{A binary variable with missing values}
#'   \item{x2}{A fully observed continuous variable}
#'   \item{failtime}{The discrete failure/censoring time}
#'   \item{d}{Indicator of failure (=1) or censoring (=0)}
#' }
#'
"ex_ncc"
