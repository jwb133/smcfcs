#' Substantive model compatible fully conditional specification imputation of
#' covariates for a Fine-Gray model
#'
#' Multiply imputes missing covariate values using substantive model compatible
#' fully conditional specification for competing risks outcomes, when the
#' substantive model is a Fine-Gray model for the subdistribution hazard of one event.
#'
#' In the presence of random right censoring, the function first multiply imputes
#' the potential censoring times for those failing from competing events using
#' \link[kmi]{kmi}, and thereafter uses \link[smcfcs]{smcfcs} to impute the missing
#' covariates. See Bonneville \emph{et al.} 2024 for further details on the methodology.
#'
#' The function does not (yet) support parallel computation.
#'
#' @author Edouard F. Bonneville \email{e.f.bonneville@@lumc.nl}
#'
#' @param smformula The formula of the substantive model, given as a string. Needs to
#' be of the form "Surv(t, d) ~ x1 + x2", where t is a vector of competing event
#' times, and d is a (numeric) competing event indicator, where 0 must designate
#' a censored observation.
#' @param cause Numeric, designating the competing event of interest (default is
#' `cause = 1`).
#' @param kmi_args List, containing arguments to be passed on to \link[kmi]{kmi}.
#' The "formula" element is a formula where the right-hand side specifies the
#' covariates used for multiply imputing the potential censoring times for
#' individual's failing from competing events. The default is `formula = ~ 1`,
#' which uses marginal Kaplan-Meier estimator of the censoring distribution.
#' @param ... Additional arguments to pass on to \link[smcfcs]{smcfcs}
#' @inheritParams smcfcs
#'
#' @return An object of type "smcfcs", as would usually be returned from
#' \link[smcfcs]{smcfcs}.
#' @export
#'
#' @references Bonneville EF, Beyersmann J, Keogh RH, Bartlett JW, Morris TP,
#' Polverelli N, de Wreede LC, Putter H. Multiple imputation of missing covariates
#' when using the Fine--Gray model. 2024. Submitted.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' library(kmi)
#'
#' imps <- smcfcs.finegray(
#'   originaldata = ex_finegray,
#'   smformula = "Surv(times, d) ~ x1 + x2",
#'   method = c("", "", "logreg", "norm"),
#'   cause = 1,
#'   kmi_args = list("formula" = ~ 1)
#' )
#'
#' if (requireNamespace("mitools", quietly = TRUE)) {
#'   library(mitools)
#'   impobj <- imputationList(imps$impDatasets)
#'   # Important: use Surv(newtimes, newevent) ~ ... when pooling
#'   # (respectively: subdistribution time and indicator for cause of interest)
#'   models <- with(impobj, coxph(Surv(newtimes, newevent) ~ x1 + x2))
#'   summary(MIcombine(models))
#' }
#' }
#'
smcfcs.finegray <- function(originaldata,
                            smformula,
                            method,
                            cause = 1,
                            m = 5,
                            numit = 10,
                            rjlimit = 5000,
                            kmi_args = list(
                              "formula" = ~ 1,
                              "bootstrap" = FALSE,
                              "nboot" = 10
                            ),
                            ...) {

  if (!requireNamespace("kmi", quietly = TRUE)) {
    stop("Package {kmi} needed for this function to work. Please install it.", call. = FALSE)
  }

  # Locate/sort out outcome variables
  outcome_vars <- all.vars(update(as.formula(smformula), . ~ 1))
  time_var_name <- outcome_vars[[1]] # will need to add warning/error if (tstart, tstop) format
  status_var_name <- outcome_vars[[2]]
  time_var <- originaldata[[time_var_name]]
  status_var <- originaldata[[status_var_name]]

  # Some checks
  checkmate::assert_numeric(x = time_var, lower = 0, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(x = status_var, lower = 0, finite = TRUE, any.missing = FALSE)

  # Get functional form RHS of formula
  smformula_rhs <- unlist(strsplit(smformula, split = "~"))[2]

  # Use names consistent with kmi's names
  smformula_processed <- paste0("Surv(newtimes, newevent) ~", smformula_rhs)
  meths_smcfcs <- c(method, "newtimes" = "", newevent = "")

  # Check number censored, since based on this we use {kmi} or not
  num_censored <- sum(status_var == 0)

  # If no censoring: just pre-process data
  if (num_censored == 0) {

    # Set competing events to event time larger that largest event 1 time
    # (in this case: largest event 1 time + 0.1)
    eps <- 0.1
    max_ev1_time <- max(time_var[status_var == cause])
    newtimes <- time_var
    newtimes[status_var != cause] <- max_ev1_time + eps
    newevent <- as.numeric(status_var == cause)

    # Bind to original data
    originaldata_processed <- cbind.data.frame(originaldata, newtimes, newevent)

    # Run imputations
    smcfcs_obj <- smcfcs(
      originaldata = originaldata_processed,
      smtype = "coxph",
      smformula = smformula_processed,
      method = meths_smcfcs,
      rjlimit = rjlimit,
      numit = numit,
      m = m,
      ...
    )

  } else { # Need to multiply impute censoring times!

    # Prepare kmi() formula (for now default Kaplan-Meier imputation)
    lhs_kmi <- paste0("Surv(", paste(outcome_vars, collapse = ", "), " != 0)")
    form_kmi <- as.formula(paste0(lhs_kmi, deparse1(kmi_args$formula)))
    rhs_cens <- all.vars(kmi_args$formula)

    # Some checks for the censoring
    if (length(rhs_cens) > 0) {
      if (any(!(rhs_cens %in% colnames(originaldata))))
        stop("Variables used to model the censoring distribution need to be in originaldata.")
      if (anyNA(originaldata[, rhs_cens]))
        warning(
          "          One or more of the variables used to model the censoring
          distribution has missing data. This MI method works under the assumption
          that only complete covariates can be associated with the censoring. If
          nothing is changed by the user, {kmi} deals with these missing values
          internally by using single imputation."
        )
    }

    args_cens_imps <- c(
      list(
        "formula" = form_kmi,
        "data" = originaldata,
        "etype" = as.symbol(status_var_name),
        "failcode" = cause,
        "nimp" = m
      ),
      kmi_args[-1] # remove formula
    )

    # Impute missing censoring times in first loop
    kmi_imps <- do.call(kmi::kmi, args = args_cens_imps)

    # And now impute the covariates
    imps_loop <- lapply(kmi_imps$imputed.data, function(new_outcomes) {

      df_imp <- cbind(kmi_imps$original.data, new_outcomes)
      df_imp$newevent <- as.numeric(df_imp$newevent) - 1L # Make numeric for smcfcs

      smcfcs_modif <- smcfcs(
        originaldata = df_imp,
        smtype = "coxph",
        smformula = smformula_processed,
        method = meths_smcfcs,
        rjlimit = rjlimit,
        numit = numit,
        m = 1L, # one imputation per kmi dataset
        ...
      )

      return(smcfcs_modif)
    })

    smcfcs_obj <- combine_smcfcs_objects(imps_loop)
  }

  # Do invisible and capture output bits here?
  return(smcfcs_obj)
}
