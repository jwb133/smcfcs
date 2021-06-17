#' Parallel substantive model compatible imputation
#'
#' Runs substantive model compatible imputation using parallel cores
#'
#' This function can be used to call one of the substantive model compatible imputation
#' methods using parallel cores, to reduce computation time. You must specify
#' the arguments required for the standard smcfcs call, and then specify your
#' the arguments for how to use parallel cores.
#'
#' @author Edouard F. Bonneville \email{e.f.bonneville@@lumc.nl}
#' @author Jonathan Bartlett \email{j.w.bartlett@@bath.ac.uk}
#'
#' @param smcfcs_func Specifies which base smcfcs function to call. Possible values
#' are `smcfcs`, `smcfcs.casecohort`, `smcfcs.dtasam`, `smcfcs.nestedcc`. Defaults
#' to `smcfcs`.
#' @param seed Optional seed, set as `set.seed` when `n_cores = 1`,
#' or as `parallel::clusterSetRNGStream` when `n_cores > 1`.
#' @param m Number of imputed datasets to generate.
#' @param n_cores Number of cores over which to split the `m` imputations. If
#' `n_cores` is not divisible exactly by `m`, one of the cores will perform
#' more/less imputations that the rest such that the final result still contains
#' `m` imputed datasets.
#' @param cl_type Either "PSOCK" or "FORK". If running on a Windows system
#' "PSOCK" is recommended, otherwise for Linux/Mac machines "FORK" tends to
#' offer faster computation - see \link[mice]{parlmice}.
#' @param outfile Optional character path to location for
#' output from the workers. Useful to diagnose rejection sampling warnings.
#' File path must be formulated as "path/to/filename.txt".
#' @param ... Additional arguments to pass on to \link[smcfcs]{smcfcs},
#' \link[smcfcs]{smcfcs.casecohort},
#' \link[smcfcs]{smcfcs.dtsam}, or
#' \link[smcfcs]{smcfcs.nestedcc}.
#'
#' @return An object of type "smcfcs", as would usually be returned from
#' \link[smcfcs]{smcfcs}.
#' @export
#'
#' @examples
#' \dontrun{
#' # Detect number of cores
#' parallel::detectCores()
#'
#' imps <- smcfcs.parallel(
#' smcfcs_func="smcfcs",
#' seed = 2021,
#' n_cores = 2,
#' originaldata = smcfcs::ex_compet,
#' m = 10,
#' smtype = "compet",
#' smformula = list(
#' "Surv(t, d == 1) ~ x1 + x2",
#' "Surv(t, d == 2) ~ x1 + x2"
#' ),
#' method = c("", "", "norm", "norm")
#' )
#' }
smcfcs.parallel <- function(smcfcs_func = "smcfcs",
                            seed = NULL,
                            m = 5,
                            n_cores = parallel::detectCores() - 1,
                            cl_type = "PSOCK",
                            outfile = "",
                            ...) {

  checkmate::matchArg(
    x = smcfcs_func,
    choices = c("smcfcs", "smcfcs.casecohort", "smcfcs.dtsam", "smcfcs.nestedcc")
  )

  # Save chosen function as string and and expression for later use
  smcfcs_func_string <- paste0("smcfcs::", smcfcs_func)
  smcfcs_func_expr <- parse(text = smcfcs_func_string)

  # Check smcfcs arguments
  args <- list(...)
  args_smcfcs <- names(formals(eval(smcfcs_func_expr)))
  check_args <- !(names(args) %in% args_smcfcs)

  if (any(check_args)) {
    wrong_args <- paste(names(args)[check_args], collapse = ", ")
    mssg <- paste0("The following are not valid arguments of",
                   smcfcs_func_string, " : ", wrong_args)
    stop(mssg)
  }

  # Check parallel arguments
  checkmate::assert_numeric(x = seed, null.ok = TRUE, any.missing = FALSE, len = 1)
  checkmate::assert_int(x = m, lower = 1)
  #checkmate::assert_int(x = m_per_core, lower = 1, upper = max(1,floor(m / n_cores)), null.ok = TRUE)
  checkmate::matchArg(x = cl_type, choices = c("PSOCK", "FORK"))
  checkmate::assert_int(x = n_cores, lower = 1, upper = min(parallel::detectCores(), m))
  if (outfile != "") checkmate::assert_path_for_output(x = outfile, overwrite = TRUE)

  # Standard smcfcs if n_cores = 1
  if (n_cores == 1) {
    if (!is.null(seed)) set.seed(seed)
    args$m <- m
    res <- do.call(eval(smcfcs_func_expr), args)

  } else {

    # Determine number of imputations per core
    imp_specs <- determine_imp_specs(n_cores, m)

    # Set up the cluster
    cl <- parallel::makeCluster(n_cores, type = cl_type, outfile = outfile)
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)

    # Export necessary objects/functions
    parallel::clusterCall(cl, assign, "Surv", survival::Surv, envir = .GlobalEnv)
    parallel::clusterCall(cl, assign, "strata", survival::strata, envir = .GlobalEnv)
    parallel::clusterExport(
      cl = cl,
      varlist = c(
        "args", "imp_specs", "seed", "m", "n_cores", "cl_type", "smcfcs",
        "smcfcs.casecohort", "smcfcs.dtsam", "smcfcs.nestedcc", "smcfcs_func_expr"
      ),
      envir = environment()
    )

    # Run the imputations
    imps <- parallel::parLapply(cl = cl, X = seq_along(imp_specs), function(x) {
      args$m <- imp_specs[x]
      do.call(eval(smcfcs_func_expr), args)
    })

    parallel::stopCluster(cl)

    # Combine imputations
    res <- combine_smcfcs_objects(imps)
  }

  return(res)
}


# Prepare imputations per core
determine_imp_specs <- function(n_cores,
                                m) {

  imp_specs <- rep(floor(m / n_cores), times = n_cores)
  modul <- m %% n_cores

  # Add remaining imps to add to m
  if (modul != 0) imp_specs[length(imp_specs)] <- imp_specs[length(imp_specs)] + modul
  return(imp_specs)
}


# Helper to combine smcfcs objects
combine_smcfcs_objects <- function(smcfcs_list) {

  # Combine imputed datasets
  ls_impdats <- do.call("c",  lapply(smcfcs_list, "[[", "impDatasets"))

  # Combine monitoring of imputations
  coef_array <- abind::abind(lapply(smcfcs_list, "[[", "smCoefIter"), along = 1)

  # Polish and return
  res <- list(
    impDatasets = ls_impdats,
    smCoefIter = coef_array,
    smInfo = list(
      smtype = smcfcs_list[[1]]$smInfo$smtype,
      smformula = smcfcs_list[[1]]$smInfo$smformula
    )
  )

  class(res) <- "smcfcs"
  return(res)
}
