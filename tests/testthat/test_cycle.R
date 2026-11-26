# Contract tests for smcfcs.cycle() (design/smcfcs-cycle.md, section 4.3)

library(smcfcs)
library(survival)

# Small dataset with a partially observed normal x and a binary z, and an
# outcome of the requested type. For lm and logistic some outcomes are
# missing as well.
cycleData <- function(smtype, n = 200) {
  set.seed(4321)
  z <- rbinom(n, 1, 0.5)
  x <- z + rnorm(n)
  if (smtype == "lm") {
    y <- x + z + rnorm(n)
  } else if (smtype == "logistic") {
    y <- rbinom(n, 1, plogis(x + z))
  } else if (smtype == "poisson") {
    y <- rpois(n, exp(0.5 * x + 0.5 * z))
  } else if (smtype == "coxph") {
    t <- -log(runif(n)) / exp(x + z)
    d <- 1 * (t < 2)
    t[d == 0] <- 2
  } else if (smtype == "compet") {
    t1 <- -log(runif(n)) / exp(x + z)
    t2 <- -log(runif(n)) / exp(0.5 * x)
    t <- pmin(t1, t2, 2)
    d <- ifelse(t == 2, 0, ifelse(t == t1, 1, 2))
  }
  x[runif(n) < 0.4] <- NA
  if (smtype %in% c("lm", "logistic", "poisson")) {
    if (smtype != "poisson") y[runif(n) < 0.1] <- NA
    dat <- data.frame(y, z, x)
    smformula <- "y ~ z + x"
    method <- c("", "", "norm")
  } else {
    dat <- data.frame(t, d, z, x)
    smformula <- if (smtype == "coxph") "Surv(t, d) ~ z + x" else
      c("Surv(t, d == 1) ~ z + x", "Surv(t, d == 2) ~ z + x")
    method <- c("", "", "", "norm")
  }
  list(data = dat, smformula = smformula, method = method)
}

# Runs smcfcs() and records the input and the random number state of the
# first call to smcfcs.cycle() for each imputation
recordFirstCycles <- function(smtype, cd, m = 2, numit = 3) {
  realCycle <- smcfcs:::smcfcs.cycle
  first <- list()
  local_mocked_bindings(smcfcs.cycle = function(data, ...) {
    args <- list(...)
    if (is.null(args$cycleState$flexsurvEsts) && is.null(args$cycleState$smCoef)) {
      first[[length(first) + 1]] <<- list(data = data, seed = .Random.seed)
    }
    realCycle(data, ...)
  })
  set.seed(99)
  imps <- suppressWarnings(smcfcs(cd$data, smtype = smtype, smformula = cd$smformula,
    method = cd$method, m = m, numit = numit, verbose = FALSE
  ))
  list(imps = imps, first = first)
}

# Runs numit cycles from data with cycleState = NULL on the first call
runCycles <- function(data, r, smtype, cd, numit, ...) {
  state <- NULL
  coefs <- NULL
  for (it in 1:numit) {
    res <- smcfcs.cycle(data, r, smtype = smtype, smformula = cd$smformula,
      method = cd$method, cycleState = state, ...
    )
    data <- res$data
    state <- res$cycleState
    coefs <- cbind(coefs, state$smCoef)
  }
  list(data = data, state = state, coefs = coefs)
}

for (smtype in c("lm", "logistic", "poisson", "coxph", "compet")) {
  test_that(paste0("feeding back cycleState reproduces smcfcs() (", smtype, ")"), {
    cd <- cycleData(smtype)
    rec <- recordFirstCycles(smtype, cd)
    r <- !is.na(cd$data)
    expect_length(rec$first, 2)
    for (imp in 1:2) {
      assign(".Random.seed", rec$first[[imp]]$seed, envir = globalenv())
      res <- runCycles(rec$first[[imp]]$data, r, smtype, cd, numit = 3)
      expect_equal(res$data, rec$imps$impDatasets[[imp]])
      expect_identical(unname(res$coefs), unname(rec$imps$smCoefIter[imp, , ]))
    }
  })
}

test_that("cycleState = NULL works and returns data and state", {
  cd <- cycleData("lm")
  r <- !is.na(cd$data)
  dat <- cd$data
  dat$y[!r[, "y"]] <- mean(dat$y, na.rm = TRUE)
  dat$x[!r[, "x"]] <- mean(dat$x, na.rm = TRUE)
  set.seed(1)
  res <- smcfcs.cycle(dat, r, "lm", cd$smformula, cd$method)
  expect_named(res, c("data", "cycleState"))
  expect_false(anyNA(res$data))
  expect_identical(names(res$data), names(dat))
  expect_length(res$cycleState$smCoef, 3)
  expect_equal(res$cycleState$rjFailCount, 0)
  expect_equal(res$cycleState$flexsurvFailCount, 0)
})

test_that("cells with r == 1 are never changed", {
  for (smtype in c("lm", "logistic", "coxph")) {
    cd <- cycleData(smtype)
    r <- !is.na(cd$data)
    dat <- cd$data
    # fill with an observed value, which keeps the type of the column
    for (j in which(colSums(!r) > 0)) {
      dat[!r[, j], j] <- dat[which(r[, j])[1], j]
    }
    set.seed(2)
    state <- NULL
    for (it in 1:3) {
      res <- smcfcs.cycle(dat, r, smtype, cd$smformula, cd$method, cycleState = state)
      for (j in seq_along(dat)) {
        expect_identical(res$data[r[, j], j], cd$data[r[, j], j])
      }
      expect_false(isTRUE(all.equal(res$data, dat)))
      dat <- res$data
      state <- res$cycleState
    }
  }
})

test_that("a logical r and a 0/1 r give identical results", {
  cd <- cycleData("coxph")
  rl <- !is.na(cd$data)
  dat <- cd$data
  dat$x[!rl[, "x"]] <- 0
  set.seed(3)
  resl <- runCycles(dat, rl, "coxph", cd, numit = 3)
  set.seed(3)
  res01 <- runCycles(dat, 1 * rl, "coxph", cd, numit = 3)
  expect_identical(resl$data, res01$data)
  expect_identical(resl$coefs, res01$coefs)
})

test_that("verbose = FALSE is silent and verbose = TRUE gives messages", {
  cd <- cycleData("lm")
  r <- !is.na(cd$data)
  dat <- cd$data
  dat$y[!r[, "y"]] <- 0
  dat$x[!r[, "x"]] <- 0
  set.seed(4)
  expect_silent(smcfcs.cycle(dat, r, "lm", cd$smformula, cd$method))
  set.seed(4)
  expect_message(smcfcs.cycle(dat, r, "lm", cd$smformula, cd$method, verbose = TRUE))
})

test_that("the setup is reused, and redone when the arguments change", {
  cd <- cycleData("lm")
  r <- !is.na(cd$data)
  dat <- cd$data
  dat$y[!r[, "y"]] <- 0
  dat$x[!r[, "x"]] <- 0
  realSetup <- smcfcs:::smcfcs.setup
  nSetup <- 0
  local_mocked_bindings(smcfcs.setup = function(...) {
    nSetup <<- nSetup + 1
    realSetup(...)
  })
  set.seed(5)
  res <- runCycles(dat, r, "lm", cd, numit = 3)
  res <- smcfcs.cycle(res$data, 1 * r, "lm", cd$smformula, cd$method,
    cycleState = res$state
  )
  expect_equal(nSetup, 1)
  # changed smformula
  res2 <- smcfcs.cycle(res$data, r, "lm", "y ~ x", cd$method, cycleState = res$cycleState)
  expect_equal(nSetup, 2)
  expect_length(res2$cycleState$smCoef, 2)
  # changed r: z now has missing values to impute
  r2 <- r
  r2[1:10, "z"] <- FALSE
  res3 <- smcfcs.cycle(res$data, r2, "lm", cd$smformula, c("", "logreg", "norm"),
    cycleState = res$cycleState
  )
  expect_equal(nSetup, 3)
  # the observed values are kept (as double: the logreg draws are stored
  # as double, as in smcfcs())
  expect_equal(res3$data$z[r2[, "z"]], res$data$z[r2[, "z"]])
  # smcfcs() prepares once for all imputations and cycles
  nSetup <- 0
  suppressWarnings(smcfcs(cd$data, "lm", cd$smformula, cd$method, m = 2, numit = 2,
    verbose = FALSE
  ))
  expect_equal(nSetup, 1)
})

test_that("changing the arguments after imputing flexsurv times is an error", {
  skip_if_not_installed("flexsurv")
  set.seed(6)
  n <- 200
  z <- rnorm(n)
  x <- z + rnorm(n)
  t <- -log(runif(n)) / exp(x + z)
  d <- 1 * (t < 1)
  t[d == 0] <- 1
  x[runif(n) < 0.3] <- NA
  dat <- data.frame(t, d, x, z)
  r <- !is.na(dat)
  dat$x[!r[, "x"]] <- 0
  args <- list(smtype = "flexsurv", smformula = "Surv(t, d) ~ x + z",
    method = c("", "", "norm", ""), k = 1, imputeTimes = TRUE, censtime = 2,
    originalKnots = TRUE
  )
  res <- do.call(smcfcs.cycle, c(list(dat, r), args))
  args$k <- 2
  expect_error(
    do.call(smcfcs.cycle, c(list(res$data, r), args, list(cycleState = res$cycleState))),
    "censored times were imputed"
  )
})
