# Generate the frozen references for the smcfcs.cycle() refactor
# (design note design/smcfcs-cycle.md, section 4.2).
#
# Run this ONCE, on the baseline commit (ae95ab3), BEFORE changing any
# code in R/:
#
#   Rscript data-raw/make-references.R [pkg_dir]
#
# Run it from the package root. pkg_dir (default ".") is the package
# source whose R/ code is used; it may be a git worktree checked out at
# the baseline, e.g.
#
#   git worktree add ../smcfcs-baseline ae95ab3
#   Rscript data-raw/make-references.R ../smcfcs-baseline
#
# The script loads that source with pkgload::load_all(), so the installed
# smcfcs is not used. It writes one .rds per case to
# tests/testthat/reference/<blas>/ in the current directory, where <blas>
# labels the BLAS library in use (see blas_tag()). Results differ between
# BLAS libraries, so run the script once under each BLAS the test should
# run with (e.g. R's reference BLAS from the terminal, and Apple's
# Accelerate where an IDE uses it).
#
# Each .rds holds the case (function, arguments including the data, seed),
# the result (impDatasets, smCoefIter, smInfo, warnings), and the
# environment it was generated in. tests/testthat/test_reference.R reruns
# each case from the stored inputs and requires identical() results.
#
# All cases are small (n = 200 or so, m = 2, numit = 3), so the whole grid
# runs in well under a minute.

needed <- c("pkgload", "survival", "flexsurv", "kmi", "brglm2", "VGAM")
missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Install these packages first, so that every case is generated: ",
       paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
pkg_dir <- if (length(args) >= 1) args[1] else "."

pkgload::load_all(pkg_dir, quiet = TRUE)
suppressPackageStartupMessages(library(survival))
source(file.path("tests", "testthat", "helper-reference.R"))

baseline <- tryCatch(
  system2("git", c("-C", pkg_dir, "rev-parse", "--short", "HEAD"), stdout = TRUE),
  error = function(e) NA_character_
)

# ---------------------------------------------------------------------------
# Data generators. Each takes a seed, so the data are fixed independently
# of the seed used for imputation.
# ---------------------------------------------------------------------------

mar <- function(v) runif(length(v)) > plogis(scale(v)[, 1])

sim_lm <- function(seed, n = 200) {
  set.seed(seed)
  z <- rnorm(n)
  x <- z + rnorm(n)
  v <- x + rnorm(n)                      # auxiliary variable
  y <- 1 + z + x + 0.5 * x^2 + rnorm(n)
  x[mar(y)] <- NA
  data.frame(y, z, x, v)
}

sim_lm_two <- function(seed, n = 200) {
  # one continuous and one binary (factor) incomplete covariate
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, plogis(x1))
  y <- 1 + x1 + x2 + x1 * x2 + rnorm(n)
  x1[runif(n) < 0.3] <- NA
  x2[runif(n) < 0.3] <- NA
  data.frame(y, x1, x2 = factor(x2))
}

sim_lm_cat <- function(seed, n = 200, ordered = FALSE) {
  set.seed(seed)
  z <- rnorm(n)
  xt <- z + rlogis(n)
  x <- cut(xt, c(-Inf, -1, 0.5, 1.5, Inf), labels = c("a", "b", "c", "d"))
  y <- z + (x == "b") - (x == "c") + 2 * (x == "d") + rnorm(n)
  x[mar(z)] <- NA
  x <- if (ordered) factor(x, ordered = TRUE) else factor(x, ordered = FALSE)
  data.frame(x, z, y)
}

sim_lm_count <- function(seed, n = 200) {
  set.seed(seed)
  z <- rnorm(n)
  x <- rpois(n, exp(0.5 + 0.3 * z))
  y <- z + 0.5 * x + rnorm(n)
  x[mar(y)] <- NA
  data.frame(y, z, x)
}

sim_lm_brlogreg <- function(seed, n = 200) {
  # binary covariate with near-perfect prediction by z
  set.seed(seed)
  z <- rnorm(n)
  x <- 1 * (z + 0.1 * rnorm(n) > 0)
  y <- x + z + rnorm(n)
  x[runif(n) < 0.4] <- NA
  data.frame(x, z, y)
}

sim_lm_misy <- function(seed, n = 200) {
  # missing outcome and missing covariate
  d <- sim_lm(seed, n)
  set.seed(seed + 1)
  d$y[runif(n) < 0.25] <- NA
  d
}

sim_restrictions <- function(seed, n = 200) {
  set.seed(seed)
  z <- rnorm(n)
  p1 <- 1 / (1 + exp(z) + exp(1 + z))
  p2 <- exp(z) / (1 + exp(z) + exp(1 + z))
  u <- runif(n)
  x <- rep("a", n)
  x[u > p1 & u < p1 + p2] <- "b"
  x[u > p1 + p2] <- "c"
  y <- z + (x == "b") + 2 * (x == "c") + rnorm(n)
  cc <- sample(0:3, n, replace = TRUE)
  xobs <- x
  xobs[cc == 1] <- "NA"
  xobs[cc == 2 & x %in% c("b", "c")] <- "b/c"
  xobs[cc == 3 & x %in% c("a", "c")] <- "a/c"
  x[cc == 1 | (cc == 2 & x %in% c("b", "c")) | (cc == 3 & x %in% c("a", "c"))] <- NA
  data.frame(x = factor(x), xobs, z, y)
}

sim_latnorm <- function(seed, n = 200) {
  set.seed(seed)
  x <- rnorm(n)
  w1 <- x + rnorm(n)
  w2 <- x + rnorm(n)
  w2[runif(n) < 0.5] <- NA               # second replicate on half the sample
  y <- x + rnorm(n)
  data.frame(x = NA_real_, w1, w2, y)
}

sim_logistic <- function(seed, n = 200, misy = FALSE) {
  set.seed(seed)
  z <- rnorm(n)
  x <- z + rnorm(n)
  y <- rbinom(n, 1, plogis(-0.5 + z + x))
  x[mar(z)] <- NA
  if (misy) y[runif(n) < 0.25] <- NA
  data.frame(y, z, x)
}

sim_poisson <- function(seed, n = 200) {
  set.seed(seed)
  z <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- rpois(n, exp(0.5 + 0.5 * x - z))
  x[runif(n) > plogis(y - 1)] <- NA
  y[runif(n) < 0.2] <- NA
  data.frame(y, x, z)
}

sim_surv <- function(seed, n = 200) {
  set.seed(seed)
  z <- rnorm(n)
  x <- z + rnorm(n)
  v <- x + rnorm(n)
  t <- -log(runif(n)) / exp(0.5 * z + 0.5 * x)
  cens <- -log(runif(n)) / 0.5
  d <- 1 * (t < cens)
  t <- pmin(t, cens)
  x[mar(t)] <- NA
  data.frame(t, d, z, x, v)
}

sim_compet <- function(seed, n = 200) {
  set.seed(seed)
  x1 <- rbinom(n, 1, 0.5)
  x2 <- x1 + rnorm(n)
  t1 <- -log(runif(n)) / (0.5 * exp(0.5 * x1 + 0.5 * x2))
  t2 <- -log(runif(n)) / (0.5 * exp(-0.5 * x1))
  cens <- -log(runif(n)) / 0.5
  t <- pmin(t1, t2, cens)
  d <- ifelse(t == cens, 0, ifelse(t == t1, 1, 2))
  x1[runif(n) < 0.3] <- NA
  x2[runif(n) < 0.3] <- NA
  data.frame(t, d, x1, x2)
}

sim_dtsam <- function(seed, n = 300) {
  set.seed(seed)
  x1 <- rbinom(n, 1, 0.5)
  x2 <- x1 + rnorm(n)
  nT <- 5
  alpha <- seq(-1.5, -1, length.out = nT)
  ev <- sapply(1:nT, function(i) rbinom(n, 1, plogis(alpha[i] + x1 - 0.5 * x2)))
  failtime <- apply(ev, 1, function(e) which(e == 1)[1])
  d <- 1 * !is.na(failtime)
  failtime[is.na(failtime)] <- nT
  x1[mar(x2)] <- NA
  data.frame(x1, x2, failtime, d)
}

sim_casecohort <- function(seed, n = 2000) {
  set.seed(seed)
  z <- rnorm(n)
  x <- z + rnorm(n)
  t <- -log(runif(n)) / (0.05 * exp(x + z))
  d <- 1 * (t < 1)
  t[d == 0] <- 1
  x[runif(n) < 0.5] <- NA
  full <- data.frame(t, d, x, z, in.subco = 0, id = 1:n)
  full$in.subco[sample(n, n * 0.1)] <- 1
  cc <- full[full$in.subco == 1 | full$d == 1, ]
  cc$entertime <- 0
  cc$entertime[cc$in.subco == 0] <- cc$t[cc$in.subco == 0] - 1e-6
  rownames(cc) <- NULL
  cc
}

sim_nestedcc <- function(seed, n = 2000) {
  set.seed(seed)
  z <- rnorm(n)
  x <- rbinom(n, 1, plogis(z))
  t <- -log(runif(n)) / (0.05 * exp(x + z))
  d <- 1 * (t < 1)
  t[d == 0] <- 1
  x[runif(n) < 0.5] <- NA
  full <- data.frame(t, d, x, z, id = 1:n)
  sets <- list()
  for (i in which(full$d == 1)) {
    ctrl <- setdiff(which(full$t >= full$t[i]), i)
    if (length(ctrl) >= 1) {
      sets[[length(sets) + 1]] <- cbind(
        rbind(full[i, ], full[ctrl[sample.int(length(ctrl), 1)], ]),
        numrisk = 1 + length(ctrl)
      )
    }
  }
  ncc <- do.call(rbind, sets)
  ncc$setno <- rep(seq_along(sets), each = 2)
  ncc$case <- rep(c(1, 0), length(sets))
  rownames(ncc) <- NULL
  ncc
}

sim_finegray <- function(seed, n = 200) {
  set.seed(seed)
  z <- rnorm(n)
  x <- rnorm(n, mean = z)
  p <- 0.15
  xb1 <- 0.75 * x - 0.5 * z
  xb2 <- -x + 0.5 * z
  dt <- 1 + rbinom(n, 1, (1 - p)^exp(xb1))
  c1 <- dt == 1
  u <- runif(n)
  tt <- numeric(n)
  num <- 1 - (1 - u[c1] * (1 - (1 - p)^exp(xb1[c1])))^(1 / exp(xb1[c1]))
  tt[c1] <- -log(1 - num / p)
  tt[!c1] <- -log(u[!c1]) / exp(xb2[!c1])
  cens <- -log(runif(n)) / 0.1
  t <- pmin(cens, tt)
  d <- ifelse(cens < tt, 0, dt)
  x[runif(n) > 0.5] <- NA
  data.frame(t, d, x, z)
}

# ---------------------------------------------------------------------------
# The grid. Each case: name, function, arguments, imputation seed.
# m = 2 and numit = 3 throughout unless stated.
# ---------------------------------------------------------------------------

small <- list(m = 2, numit = 3)
case <- function(name, fun, seed, ...) {
  list(name = name, fun = fun, seed = seed, args = c(list(...), small))
}

d_lm <- sim_lm(1)
pm_aux <- matrix(0, 4, 4, dimnames = list(names(d_lm), names(d_lm)))
pm_aux["x", c("z", "v")] <- 1

d_res <- sim_restrictions(2)
restr <- append(list(c("xobs = a/c ~ a + c", "xobs = b/c ~ b + c")),
                as.list(c("", "", "")))

d_lat <- sim_latnorm(3)
err_mat <- matrix(0, 4, 4)
err_mat[1, 2:3] <- 1

d_surv <- sim_surv(4)
d_ncc <- sim_nestedcc(5)
pm_ncc <- matrix(0, ncol(d_ncc), ncol(d_ncc))
pm_ncc[which(names(d_ncc) == "x"), which(names(d_ncc) == "z")] <- 1

cases <- list(
  # lm: covariate methods, formula terms, passive variables, predictors
  case("lm_norm_quad", "smcfcs", 101, originaldata = d_lm, smtype = "lm",
       smformula = "y ~ z + x + I(x^2)", method = c("", "", "norm", "")),
  case("lm_passive", "smcfcs", 102,
       originaldata = transform(d_lm, xsq = x^2), smtype = "lm",
       smformula = "y ~ z + x + xsq", method = c("", "", "norm", "", "x^2")),
  case("lm_predictormatrix", "smcfcs", 103, originaldata = d_lm, smtype = "lm",
       smformula = "y ~ z + x + I(x^2)", method = c("", "", "norm", ""),
       predictorMatrix = pm_aux),
  case("lm_interaction_logreg", "smcfcs", 104, originaldata = sim_lm_two(6),
       smtype = "lm", smformula = "y ~ x1 * x2", method = c("", "norm", "logreg")),
  case("lm_brlogreg", "smcfcs", 105, originaldata = sim_lm_brlogreg(7),
       smtype = "lm", smformula = "y ~ x + z", method = c("brlogreg", "", "")),
  case("lm_poisson_cov", "smcfcs", 106, originaldata = sim_lm_count(8),
       smtype = "lm", smformula = "y ~ z + x", method = c("", "", "poisson")),
  case("lm_podds", "smcfcs", 107, originaldata = sim_lm_cat(9, ordered = TRUE),
       smtype = "lm", smformula = "y ~ z + x", method = c("podds", "", "")),
  case("lm_mlogit", "smcfcs", 108, originaldata = sim_lm_cat(10),
       smtype = "lm", smformula = "y ~ z + x", method = c("mlogit", "", "")),
  case("lm_restrictions", "smcfcs", 109, originaldata = d_res, smtype = "lm",
       smformula = "y ~ z + x", method = c("mlogit", "", "", ""),
       restrictions = restr),
  case("lm_latnorm", "smcfcs", 110, originaldata = d_lat, smtype = "lm",
       smformula = "y ~ x", method = c("latnorm", "", "", ""),
       errorProneMatrix = err_mat),
  case("lm_missing_outcome", "smcfcs", 111, originaldata = sim_lm_misy(11),
       smtype = "lm", smformula = "y ~ z + x + I(x^2)",
       method = c("", "", "norm", "")),
  case("lm_noisy", "smcfcs", 112, originaldata = d_lm, smtype = "lm",
       smformula = "y ~ z + x", method = c("", "", "norm", ""), noisy = TRUE),

  # other outcome models with an outcome step
  case("logistic", "smcfcs", 201, originaldata = sim_logistic(12),
       smtype = "logistic", smformula = "y ~ z + x", method = c("", "", "norm")),
  case("logistic_missing_outcome", "smcfcs", 202,
       originaldata = sim_logistic(13, misy = TRUE), smtype = "logistic",
       smformula = "y ~ z + x", method = c("", "", "norm")),
  case("brlogistic", "smcfcs", 203, originaldata = sim_logistic(14),
       smtype = "brlogistic", smformula = "y ~ z + x", method = c("", "", "norm")),
  case("poisson_missing_outcome", "smcfcs", 204, originaldata = sim_poisson(15),
       smtype = "poisson", smformula = "y ~ x + z", method = c("", "norm", "")),

  # survival models
  case("coxph", "smcfcs", 301, originaldata = d_surv, smtype = "coxph",
       smformula = "Surv(t, d) ~ z + x + I(x^2)",
       method = c("", "", "", "norm", "")),
  case("coxph_rjfail", "smcfcs", 302, originaldata = d_surv, smtype = "coxph",
       smformula = "Surv(t, d) ~ z + x", method = c("", "", "", "norm", ""),
       rjlimit = 5),  # rjlimit = 1 errors at baseline (drop = FALSE missing)
  case("compet", "smcfcs", 303, originaldata = sim_compet(16), smtype = "compet",
       smformula = c("Surv(t, d == 1) ~ x1 + x2", "Surv(t, d == 2) ~ x1 + x2"),
       method = c("", "", "logreg", "norm")),
  case("weibull", "smcfcs", 304, originaldata = d_surv, smtype = "weibull",
       smformula = "Surv(t, d) ~ z + x", method = c("", "", "", "norm", "")),
  case("dtsam_factor", "smcfcs.dtsam", 305, originaldata = sim_dtsam(17),
       smformula = "Surv(failtime, d) ~ x1 + x2", method = c("logreg", "", "", "")),
  case("dtsam_linear", "smcfcs.dtsam", 306, originaldata = sim_dtsam(18),
       smformula = "Surv(failtime, d) ~ x1 + x2", method = c("logreg", "", "", ""),
       timeEffects = "linear"),
  case("casecohort", "smcfcs.casecohort", 307, originaldata = sim_casecohort(19),
       smformula = "Surv(entertime, t, d) ~ x + z", sampfrac = 0.1,
       in.subco = "in.subco", method = c("", "", "norm", "", "", "", "")),
  case("nestedcc", "smcfcs.nestedcc", 308, originaldata = d_ncc,
       set = "setno", nrisk = "numrisk", event = "d",
       smformula = "Surv(t, case) ~ x + z + strata(setno)",
       method = c("", "", "logreg", "", "", "", "", ""), predictorMatrix = pm_ncc),
  case("finegray", "smcfcs.finegray", 309, originaldata = sim_finegray(20),
       smformula = "Surv(t, d) ~ x + z", method = c("", "", "norm", "")),
  case("flexsurv", "smcfcs.flexsurv", 310, originaldata = d_surv[, 1:4],
       smformula = "Surv(t, d) ~ z + x", method = c("", "", "", "norm"), k = 2),
  case("flexsurv_imputetimes", "smcfcs.flexsurv", 311,
       originaldata = d_surv[, 1:4], smformula = "Surv(t, d) ~ z + x",
       method = c("", "", "", "norm"), k = 2, imputeTimes = TRUE, censtime = 3)
)

# ---------------------------------------------------------------------------
# Run and save
# ---------------------------------------------------------------------------

out_dir <- file.path(reference_dir("."), blas_tag())
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
env <- reference_environment()

for (cs in cases) {
  cat(sprintf("%-26s", cs$name))
  t0 <- proc.time()[["elapsed"]]
  res <- run_reference_case(cs)
  saveRDS(
    list(case = cs, result = res, baseline = baseline, environment = env),
    file.path(out_dir, paste0(cs$name, ".rds")),
    version = 2  # readable by R < 3.5.0; avoids raising the R dependency
  )
  cat(sprintf("%5.1fs  %d warning(s)\n", proc.time()[["elapsed"]] - t0,
              length(res$warnings)))
}
cat("Wrote", length(cases), "references to", out_dir, "from", baseline, "\n")
