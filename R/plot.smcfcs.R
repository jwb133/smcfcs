#' Assess convergence of a smcfcs object
#'
#' Visualises the contents of smCoefIter. Specifically, it plots the parameter
#' estimates of the substantive model against the number of iterations from
#' the imputation procedure. This is done for each regression coefficient,
#' and each line corresponds to an imputed dataset.
#'
#' Requires loading of ggplot2 plotting library.
#'
#' @author Edouard F. Bonneville \email{e.f.bonneville@@lumc.nl}
#'
#' @param x An object of class 'smcfcs'
#' @param include Character vector of coefficient names for which to return the
#' convergence plot. Default is "all" and returns plots for all coefficients in
#' a facetted manner.
#'
#' Recommendation is to plot first with include = "all", and then select
#' coefficient names to zoom in to.
#'
#' For competing risks, the coefficients are indexed by their cause. E.g. for
#' coefficient of a variable x1 in a model for cause 2, will be labelled
#' "x1-cause2".
#' @param ... Additional parameters to pass on to ggplot2::facet_wrap(),
#' eg. nrow = 2
#'
#' @return A ggplot2 object, containing the convergence plots, facetted per
#' covariate in the substantive model
#' @export
#'
#' @examples
#' \dontrun{
#' # Use simulated competing risks example in package
#' imps <- smcfcs(
#'   originaldata = ex_compet,
#'   smtype = "compet",
#'   smformula = list(
#'     "Surv(t, d == 1) ~ x1 + x2",
#'     "Surv(t, d == 2) ~ x1 + x2"
#'   ),
#'   method = c("", "", "norm", "norm")
#' )
#'
#' plot(imps)
#' plot(imps, include = c("x1-cause1", "x2-cause2"))
#' }
#'
#' @importFrom rlang .data
plot.smcfcs <- function(x,
                        include = "all",
                        ...) {
  if (!inherits(x, "smcfcs")) {
    stop("'x' must be a 'smcfcs' object")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  # Prepare data
  df_plot <- prep_iters(x)

  # Choose plots to include
  if (length(include) >= 1 & include[1] != "all") {
    coef_names <- unique(df_plot$covar)
    if (any(!(include %in% coef_names))) {
      mssg <- paste0(
        "include should be character vector containing any combination of: '",
        paste0(coef_names, collapse = "','"),
        "'; or simply 'all'"
      )
      stop(mssg)
    } else {
      df_plot <- df_plot[df_plot$covar %in% include, ]
    }
  }

  # Make plot
  p <- ggplot2::ggplot(
    data = df_plot,
    ggplot2::aes(
      x = .data$iters,
      y = .data$value,
      col = factor(.data$imp)
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = "Iterations", y = "Coefficient") +
    ggplot2::facet_wrap(~covar, ...)

  return(p)
}


# Prepare data for plotting
prep_iters <- function(x) {
  # Extract meta data
  M <- dim(x$smCoefIter)[1]
  smtype <- x$smInfo$smtype
  smformula <- if (inherits(x$smInfo$smformula, "formula")) deparse1(x$smInfo$smformula) else x$smInfo$smformula
  dat <- x$impDatasets[[1]] # for names in model matrix
  numit <- dim(x$smCoefIter)[3]

  if (numit < 2) {
    stop("Re-run smcfcs() with numit >= 2 in order to assess convergence")
  }

  # Check if competing risks
  coef_names <- if (smtype == "compet") {
    K <- length(smformula)
    cause_coef_names <- lapply(X = seq_len(K), FUN = function(k) {
      names_mod <- get_coef_names(smformula[k], dat, intercept = FALSE)
      paste0(names_mod, "-cause", as.character(k))
    })

    unlist(cause_coef_names)
  } else if (smtype %in% c("weibull", "coxph", "casecohort", "nestedcc")) {
    get_coef_names(smformula, dat, intercept = FALSE)
  } else if (smtype == "dtsam") {
    get_dtsam_names(smformula, dat, x$extraArgs$timeEffects)
  } else {
    get_coef_names(smformula, dat, intercept = TRUE)
  }

  # Prepare df for plotting
  ests_list <- lapply(X = seq_len(M), function(m) {
    coef_dat <- as.data.frame(t(x$smCoefIter[m, , ]))
    coef_dat$iters <- seq_len(numit)
    coef_dat$imp <- m

    return(coef_dat)
  })

  ests_df <- do.call(rbind.data.frame, ests_list)
  colnames(ests_df) <- c(coef_names, "iters", "imp")

  # Make long format
  ests_long <- stats::reshape(
    data = ests_df,
    varying = coef_names,
    timevar = "covar",
    v.names = "value",
    idvar = c("imp", "iters"),
    direction = "long",
    times = coef_names
  )

  return(ests_long)
}

get_dtsam_names <- function(smformula,
                            dat,
                            timeEffects) {
  # Get sides of formula
  rhs <- gsub(x = smformula, pattern = ".*~", replacement = "")
  lhs <- gsub(x = smformula, pattern = "~.*", replacement = "")
  surv_obj <- with(dat, eval(parse(text = lhs)))
  time_var <- as.matrix(surv_obj)[, "time"]

  # Make longdata
  cutPoints <- seq_len(max(time_var)) # Probably change to unique timepoints later
  longData <- survival::survSplit(as.formula(smformula), data = dat, cut = cutPoints)

  # Get model matrix formula
  smformula_matrix <- if (timeEffects == "factor") {
    paste0("~ -1 + factor(tstart) + ", rhs)
  } else if (timeEffects == "linear") {
    paste0("~ tstart + ", rhs)
  } else {
    paste0("~ tstart + I(tstart^2) + ", rhs)
  }

  # Make matrix
  model_mat <- stats::model.matrix(
    object = as.formula(smformula_matrix),
    data = longData
  )

  coef_names <- colnames(model_mat)

  return(coef_names)
}

get_coef_names <- function(smformula,
                           dat,
                           intercept) {
  rhs <- gsub(x = smformula, pattern = ".*~", replacement = "")
  smformula_matrix <- as.formula(paste0("~ +", rhs))

  # Check if there is stratification - if so remove from model matrix
  # (stratification means different baseline hazards, coefficients still same)
  if (grepl(x = rhs, pattern = "strata")) {
    strata_var <- gsub(x = rhs, pattern = ".*\\(|\\).*", replacement = "")
    rm_strata <- as.formula(paste0("~ . - strata(", strata_var, ")"))
    smformula_matrix <- update(smformula_matrix, rm_strata)
  }

  model_mat <- stats::model.matrix(
    object = smformula_matrix,
    data = dat
  )

  # For survival models
  if (intercept == FALSE) {
    model_mat <- model_mat[, !(colnames(model_mat) %in% "(Intercept)")]
  }

  coef_names <- colnames(model_mat)

  return(coef_names)
}
