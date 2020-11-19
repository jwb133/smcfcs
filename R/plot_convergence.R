#' Plot convergence of a smcfcs object
#'
#' Visualises the contents of smCoefIter.
#'
#' @param x An object of class 'smcfcs'
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
#' plot(imps)
#' }
#'
#' @importFrom rlang .data
plot.smcfcs <- function(x,
                        #use_ggplot = F,
                        #include = "all",
                        ...) {

  if (!inherits(x, "smcfcs"))
    stop("'x' must be a 'smcfcs' object")

  # Prepare data
  df_plot <- prep_iters(x)

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


# Helpers for the plot function
prep_iters <- function(x) {

  # Extract meta data
  M <- dim(x$smCoefIter)[1]
  numit <- dim(x$smCoefIter)[3]
  smtype <- x$smInfo$smtype
  smformula <- x$smInfo$smformula
  dat <- x$impDatasets[[1]] # for names in model matrix

  # Check if competing risks
  if (smtype == "compet") {

    K <- length(smformula)
    cause_coef_names <- lapply(X = 1:K, FUN = function(k) {
      names_mod <- get_coef_names(smformula[k], dat, intercept = F)
      paste0(names_mod, ".cause", as.character(k))
    })

    coef_names <- unlist(cause_coef_names)

  } else {

    # No intercept for other survival models
    if (smtype %in% c("weibull", "coxph")) {
      coef_names <- get_coef_names(smformula, dat, intercept = F)
    } else {
      coef_names <- get_coef_names(smformula, dat, intercept = T)
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
                           intercept) {

  rhs <- gsub(x = smformula, pattern = ".*~", replacement = "")

  model_mat <- stats::model.matrix(
    object = as.formula(paste0("~ +", rhs)),
    data = dat
  )

  # For survival models
  if (intercept == F) {
    model_mat <- model_mat[, !(colnames(model_mat) %in% "(Intercept)")]
  }

  coef_names <- colnames(model_mat)

  return(coef_names)
}
