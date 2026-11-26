# Shared code for the frozen-reference regression test (design note,
# section 4.2). Used by data-raw/make-references.R to generate the
# references and by test_reference.R to check them.

# Directory holding the reference files: one subdirectory per BLAS
# library, with one .rds per case. Results differ between BLAS libraries
# (a rounding difference can flip a rejection-sampling decision, after
# which the random-number stream diverges), so each BLAS needs its own set.
reference_dir <- function(root = ".") {
  file.path(root, "tests", "testthat", "reference")
}

# A short label for the BLAS library in use, used as the subdirectory
# name: "Rblas" for R's reference BLAS, "vecLib" for Apple's Accelerate,
# otherwise derived from the library's file name.
blas_tag <- function(path = extSoftVersion()[["BLAS"]]) {
  p <- normalizePath(path, mustWork = FALSE)
  if (grepl("vecLib|Accelerate", p)) return("vecLib")
  if (grepl("libRblas", p)) return("Rblas")
  b <- sub("[.](dylib|so|dll).*$", "", basename(p))
  if (!nzchar(b)) "unknown" else gsub("[^A-Za-z0-9]+", "_", b)
}

# Run one case: set its seed, call the function, and collect the result
# and the warnings. Printed output and messages are discarded, because
# the refactor deliberately turns print() into message() (section 2.5).
run_reference_case <- function(case) {
  fun <- get(case$fun, envir = asNamespace("smcfcs"))
  warnings <- character(0)
  set.seed(case$seed)
  utils::capture.output(
    result <- suppressMessages(withCallingHandlers(
      do.call(fun, case$args),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )),
    type = "output"
  )
  list(
    impDatasets = result$impDatasets,
    smCoefIter = result$smCoefIter,
    smInfo = result$smInfo,
    warnings = warnings
  )
}

# The environment a reference was generated in. A difference in any of
# these can change results without any change to smcfcs, so the test
# skips a case rather than fail when they differ.
reference_environment <- function() {
  pkgs <- c("MASS", "survival", "VGAM", "brglm2", "flexsurv", "kmi")
  list(
    R = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    RNGkind = RNGkind(),
    BLAS = unname(extSoftVersion()["BLAS"]),
    LAPACK = La_version(),
    packages = vapply(pkgs, function(p) {
      if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
    }, character(1))
  )
}
