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

# A fingerprint of the linear algebra actually in use: an md5 checksum of
# a matrix product, a Cholesky factor and an eigendecomposition of fixed
# numbers. Different BLAS/LAPACK libraries round differently, so they give
# different fingerprints. This is more reliable than the library name
# reported by extSoftVersion(): a session can report Accelerate (because
# an IDE loaded it into the process) while R's arithmetic still runs
# through R's own BLAS. The global random-number state is left unchanged.
linalg_fingerprint <- function() {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
          else rm(".Random.seed", envir = globalenv()))
  set.seed(20261126)
  X <- matrix(rnorm(400 * 60), 400)
  A <- crossprod(X)
  x <- list(gemm = A %*% t(X[1:60, ]), chol = chol(A),
            eigen = eigen(A, symmetric = TRUE))
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  saveRDS(x, f, version = 2, compress = FALSE)
  unname(tools::md5sum(f))
}

# A short label for the BLAS library in use, used only to name the
# subdirectory a new set of references is written to: "Rblas" for R's
# reference BLAS, "vecLib" for Apple's Accelerate, otherwise derived from
# the library's file name. The test chooses a set by linalg_fingerprint(),
# not by this label.
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
    BLAS = unname(extSoftVersion()["BLAS"]),   # recorded for information only
    LAPACK = La_version(),
    fingerprint = linalg_fingerprint(),
    packages = vapply(pkgs, function(p) {
      if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
    }, character(1))
  )
}
