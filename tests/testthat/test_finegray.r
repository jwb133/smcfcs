library(smcfcs)
library(survival)
library(kmi)
context("Fine-Gray model testing")

# In style of test_coxph.r, probably need to re-write with mini simulations
test_that("Fine-Gray imputation is approximately unbiased", {
  skip_on_cran()
  expect_equal(
    {
      set.seed(4321)
      n <- 10000
      z <- rnorm(n = n)
      x <- rnorm(n = n, mean = z)
      p <- 0.15
      betas <- c(0.75, 0.5)
      xb1 <- betas[1] * x - betas[2] * z # Fine-Gray linear predictor
      xb2 <- -x + 0.5 * z # Linear predictor for hazard conditional on cause 2
      d_tilde <- 1 + rbinom(n = n, size = 1, prob = (1 - p)^exp(xb1))
      cause1_ind <- d_tilde == 1
      U <- runif(n = n)
      t_tilde <- numeric(length = n)
      num <- 1 - (1 - U[cause1_ind] * (1 - (1 - p)^xb1[cause1_ind]))^(1 / xb1[cause1_ind])
      t_tilde[cause1_ind] <- -log(1 - num / p)
      t_tilde[!cause1_ind] <- -log(U[!cause1_ind]) / exp(xb2[!cause1_ind])
      cens <- -log(runif(n = n)) / 0.1
      t <- pmin(cens, t_tilde)
      d <- ifelse(cens < t_tilde, 0, d_tilde)
      x[runif(n) > 0.5] <- NA

      simData <- data.frame(t, d, x, z)

      imps <- smcfcs.finegray(
        originaldata = simData,
        smformula = "Surv(t, d) ~ x + z",
        method = c("", "", "norm", "")
      )
      library(mitools)
      impobj <- imputationList(imps$impDatasets)
      models <- with(impobj, coxph(Surv(newtimes, newevent) ~ x + z))
      abs(summary(MIcombine(models))[1, 1] - betas[1]) < 0.1
    },
    TRUE
  )
})

test_that("Fine-Gray imputation works with no censoring", {
  skip_on_cran()
  expect_equal(
    {
      set.seed(4321)
      n <- 10000
      z <- rnorm(n = n)
      x <- rnorm(n = n, mean = z)
      p <- 0.15
      betas <- c(0.75, 0.5)
      xb1 <- betas[1] * x - betas[2] * z # Fine-Gray linear predictor
      xb2 <- -x + 0.5 * z # Linear predictor for hazard conditional on cause 2
      d_tilde <- 1 + rbinom(n = n, size = 1, prob = (1 - p)^exp(xb1))
      cause1_ind <- d_tilde == 1
      U <- runif(n = n)
      t_tilde <- numeric(length = n)
      num <- 1 - (1 - U[cause1_ind] * (1 - (1 - p)^xb1[cause1_ind]))^(1 / xb1[cause1_ind])
      t_tilde[cause1_ind] <- -log(1 - num / p)
      t_tilde[!cause1_ind] <- -log(U[!cause1_ind]) / exp(xb2[!cause1_ind])
      x[runif(n) > 0.5] <- NA

      simData <- data.frame(t_tilde, d_tilde, x, z)

      imps <- smcfcs.finegray(
        originaldata = simData,
        smformula = "Surv(t_tilde, d_tilde) ~ x + z",
        method = c("", "", "norm", "")
      )
      library(mitools)
      impobj <- imputationList(imps$impDatasets)
      models <- with(impobj, coxph(Surv(newtimes, newevent) ~ x + z))
      abs(summary(MIcombine(models))[1, 1] - betas[1]) < 0.1
    },
    TRUE
  )
})
