library(smcfcs)
library(survival)
library(mitools)
context("Nested case control testing")

test_that("Nested case control imputation runs and is approximately unbiased, continuous covariate missing", {
  skip_on_cran()
  expect_equal(
    {
      set.seed(1234)

      # run simulation study
      nSims <- 100
      xEsts <- array(0, dim = nSims)
      n <- 10000

      for (sim in 1:nSims) {
        print(sim)
        z <- rnorm(n)
        x <- z + rnorm(n)
        t <- -log(runif(n)) / (0.01 * exp(x + z))
        d <- 1 * (t < 1)
        t[d == 0] <- 1
        x[(runif(n) < 0.5)] <- NA

        fullcohortdata <- data.frame(t, d, x, z)
        fullcohortdata$id <- 1:n

        # Compute number at risk at each event time using the full cohort data
        nrisk.fit <- survfit(Surv(t, d) ~ 1, data = fullcohortdata)
        ord.t.d1 <- order(fullcohortdata$t[fullcohortdata$d == 1])

        m <- 1 # 1 control per case
        ncc <- NULL

        no.sample <- 0
        for (i in which(fullcohortdata$d == 1))
        {
          # select control(s) for nested case-control
          possible.controls <- which(fullcohortdata$t >= fullcohortdata$t[i])
          # remove the case from this vector
          possible.controls <- possible.controls[which(possible.controls != i)]

          if (length(possible.controls) >= m) {
            controls <- sample(possible.controls, m)
            numAtRisk <- 1 + length(possible.controls)
            ncc <- rbind(ncc, cbind(fullcohortdata[i, ], numrisk = numAtRisk))
            ncc <- rbind(ncc, cbind(fullcohortdata[controls, ], numrisk = numAtRisk))
            no.sample <- no.sample + 1
          }
        }

        ncc$setno <- rep(1:no.sample, each = m + 1)
        ncc$case <- rep(c(1, rep(0, m)), no.sample)

        predictorMatrix <- matrix(0, nrow = dim(ncc)[2], ncol = dim(ncc)[2])
        predictorMatrix[which(colnames(ncc) == "x"), c(which(colnames(ncc) == "z"))] <- 1

        imps <- smcfcs.nestedcc(
          originaldata = ncc, set = "setno", nrisk = "numrisk", event = "d", smformula = "Surv(t,case)~x+z+strata(setno)",
          method = c("", "", "norm", "", "", "", "", ""), predictorMatrix = predictorMatrix
        )

        impobj <- imputationList(imps$impDatasets)
        models <- with(impobj, clogit(case ~ x + z + strata(setno)))

        xEsts[sim] <- summary(MIcombine(models))[1, 1]
      }
      print(mean(xEsts))
      abs(mean(xEsts) - 1) < 0.1
    },
    TRUE
  )
})

test_that("Nested case control imputation runs and is approximately unbiased, binary covariate missing", {
  skip_on_cran()
  expect_equal(
    {
      set.seed(5678)

      # perform simulation study
      nSims <- 100
      xEsts <- array(0, dim = nSims)
      n <- 10000

      for (sim in 1:nSims) {
        print(sim)
        z <- rnorm(n)
        x <- rbinom(n, 1, exp(z) / (1 + exp(z)))
        t <- -log(runif(n)) / (0.01 * exp(x + z))
        d <- 1 * (t < 1)
        t[d == 0] <- 1
        x[(runif(n) < 0.5)] <- NA

        fullcohortdata <- data.frame(t, d, x, z)
        fullcohortdata$id <- 1:n

        # Compute number at risk at each event time using the full cohort data
        nrisk.fit <- survfit(Surv(t, d) ~ 1, data = fullcohortdata)
        ord.t.d1 <- order(fullcohortdata$t[fullcohortdata$d == 1])

        m <- 1 # 1 control per case
        ncc <- NULL

        no.sample <- 0
        for (i in which(fullcohortdata$d == 1))
        {
          # select control(s) for nested case-control
          possible.controls <- which(fullcohortdata$t >= fullcohortdata$t[i])
          # remove the case from this vector
          possible.controls <- possible.controls[which(possible.controls != i)]

          if (length(possible.controls) >= m) {
            controls <- sample(possible.controls, m)
            numAtRisk <- 1 + length(possible.controls)
            ncc <- rbind(ncc, cbind(fullcohortdata[i, ], numrisk = numAtRisk))
            ncc <- rbind(ncc, cbind(fullcohortdata[controls, ], numrisk = numAtRisk))
            no.sample <- no.sample + 1
          }
        }

        ncc$setno <- rep(1:no.sample, each = m + 1)
        ncc$case <- rep(c(1, rep(0, m)), no.sample)

        predictorMatrix <- matrix(0, nrow = dim(ncc)[2], ncol = dim(ncc)[2])
        predictorMatrix[which(colnames(ncc) == "x"), c(which(colnames(ncc) == "z"))] <- 1

        imps <- smcfcs.nestedcc(
          originaldata = ncc, set = "setno", nrisk = "numrisk", event = "d", smformula = "Surv(t,case)~x+z+strata(setno)",
          method = c("", "", "logreg", "", "", "", "", ""), predictorMatrix = predictorMatrix
        )

        impobj <- imputationList(imps$impDatasets)
        models <- with(impobj, clogit(case ~ x + z + strata(setno)))

        xEsts[sim] <- summary(MIcombine(models))[1, 1]
      }
      print(mean(xEsts))
      abs(mean(xEsts) - 1) < 0.1
    },
    TRUE
  )
})
