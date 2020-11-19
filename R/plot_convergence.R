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
#' @param contrast Contrast to choose for any ordered categorical covariates
#' in the substantive model, see ?stats::contrasts .
#' Default is "contr.treatment".
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
#' originaldata = ex_compet,
#' smtype = "compet",
#' smformula = list(
#' "Surv(t, d == 1) ~ x1 + x2",
#' "Surv(t, d == 2) ~ x1 + x2"
#' ),
#' method = c("", "", "norm", "norm")
#' )
#'
#' library(ggplot2)
#' plot(imps)
#' plot(imps, include = c("x1-cause1", "x2-cause2"))
#' }
#'
#' @importFrom rlang .data
plot.smcfcs <- function(x,
                        include = "all",
                        contrast = "contr.treatment",
                        #use_ggplot = F,
                        ...) {

  if (!inherits(x, "smcfcs"))
    stop("'x' must be a 'smcfcs' object")

  # Prepare data
  df_plot <- prep_iters(x, contrast = contrast)

  # Choose plots to include
  if (length(include) >= 1 & include[1] != "all") {

    coef_names = unique(df_plot$covar)
    if (any(!(include %in% coef_names))) {

      mssg <- paste0(
        "include should be character vector containing any combination of: '",
        paste0(coef_names, collapse = "','"),
        "'; or simply 'all'"
      )
      stop(mssg)
    } else df_plot <- df_plot[df_plot$covar %in% include, ]
  }

  # Make plot
  p <- ggplot2::ggplot(
    data = df_plot,
    ggplot2::aes(x = .data$iters,
                 y = .data$value,
                 col = factor(.data$imp))
  ) +
    ggplot2::geom_line() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = "Iterations", y = "Coefficient") +
    ggplot2::facet_wrap(~ covar, ...)

  return(p)
}


# Prepare data for plotting
prep_iters <- function(x, contrast) {

  # Extract meta data
  M <- dim(x$smCoefIter)[1]
  smtype <- x$smInfo$smtype
  smformula <- x$smInfo$smformula
  dat <- x$impDatasets[[1]] # for names in model matrix
  numit <- dim(x$smCoefIter)[3]

  if (numit < 2)
    stop("Re-run smcfcs() with numit >= 2 in order to assess convergence")

  # Check if competing risks
  if (smtype == "compet") {

    K <- length(smformula)
    cause_coef_names <- lapply(X = 1:K, FUN = function(k) {
      names_mod <- get_coef_names(smformula[k], dat, intercept = F, contrast)
      paste0(names_mod, "-cause", as.character(k))
    })

    coef_names <- unlist(cause_coef_names)

  } else {

    # No intercept for other survival models
    if (smtype %in% c("weibull", "coxph")) {
      coef_names <- get_coef_names(smformula, dat, intercept = F, contrast)
    } else {
      coef_names <- get_coef_names(smformula, dat, intercept = T, contrast)
    }
  }

  # Prepare df for plotting
  ests_list <- lapply(X = 1:M, function(m) {

    coef_dat <- as.data.frame(t(x$smCoefIter[m, ,]))
    coef_dat$iters <- 1:numit
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

get_coef_names <- function(smformula,
                           dat,
                           intercept,
                           contrast) {

  rhs <- gsub(x = smformula, pattern = ".*~", replacement = "")

  # Check if any ordered categorical
  check_ordered <- lapply(dat, is.ordered)

  if (any(unlist(check_ordered))) {
    contr_list <- check_ordered[check_ordered == 1]
    contr_list <- replace(contr_list, values = contrast)
  } else contr_list <- NULL

  model_mat <- stats::model.matrix(
    object = as.formula(paste0("~ +", rhs)),
    data = dat,
    contrasts.arg = contr_list
  )

  # For survival models
  if (intercept == F) {
    model_mat <- model_mat[, !(colnames(model_mat) %in% "(Intercept)")]
  }

  coef_names <- colnames(model_mat)

  return(coef_names)
}
