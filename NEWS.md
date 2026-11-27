Version 2.1.0

* New function `smcfcs.cycle()`, which runs a single SMC-FCS cycle on a completed
  dataset. It is a low-level building block for packages that run their own
  imputation loop, such as a block imputation method for mice; most users want
  `smcfcs()`. `smcfcs()` now calls it for each iteration, with identical results.
  `smcfcs.cycle()` and the associated restructuring were contributed by
  Stef van Buuren.

* Status messages (variable roles, imputation number, imputation models) are now
  issued with `message()` instead of `print()`, so they can be silenced with
  `suppressMessages()` or the new `verbose` argument (default `TRUE`).
  `noisy = TRUE` still prints the model summaries; its "Iteration" lines are now
  messages too.

* Fixed an error ("incorrect number of subscripts") in `smcfcs.flexsurv()` with
  `imputeTimes = TRUE` and no missing covariates when the workspace contained an
  object named `smCoefIter`.
  
* Fixed outcome residual variance used when imputing missing outcomes for
  linear (`"lm"`) substantive models: it was using point estimate of residual
  variance, rather than posterior draw.

Version 2.0.3 - fix to imputation of missing Poisson outcomes (spotted by Stef van Buuren).

Version 2.0.2 - fix to flexsurv imputation to overcome error sometimes occurring with
factor covariates.

Version 2.0.1 - numerical instability fixes - replacing expit with plogis and replacing
manual calculation of fitted values in mlogit imputation with predictvglm from VGAM
package.

Version 2.0.0 - adding new restrictions argument for imputation with coarsening, and
ordered factors now imputed using MASS::polr.

Version 1.9.2 - fixing dtsam imputation for factor covariates

Version 1.9.1 - adding model=TRUE to coxph call to resolve an error some have experienced

Version 1.9.0 - new functionality for flexible parametric survival models

Version 1.8.0 - new functionality for Fine-Gray models for competing risks outcomes

Version 1.7.1 - Allowing parallel runs without seed being specified, and fixing survival data error trap

Version 1.7.0 - Adding bias reduced logistic options to smtype and method
arguments

Version 1.6.2 - fixing an issue with coxph error check added in previous version

Version 1.6.1 - adding check for event indicator values for coxph models

Version 1.6.0 - adding parallel core functionality

Version 1.5.0 - adding discrete time survival analysis functionality

Version 1.4.2 - adding convergence plot functionality

Version 1.4.1 - adding coxph control option to stop error occurring due to recent change in survival package default behaviour handling of ties

Version 1.4.0 - adding covariate measurement error and Weibull substantive model functionality.

Version 1.3.1 - fixing a bug in a check smcfcs does that every variable that should have an imputation method specified does.

Version 1.3.0 - adding functionality for nested case control and case cohort studies

Version 1.2.1 - fixing bug in code for linear regression substantive models

Version 1.2.0 - functionality added to handle Poisson substantive models

Version 1.1.1 - a number of small changes, plus critical bug fix for linear substantive models with categorical missing covariates. The bug which has been fixed meant that if the substantive model was linear regression, imputation of categorical covariates was done incorrectly.

Version 1.1.0 - second release on CRAN, including functionality for competing risks outcomes.

Version 1.0.0 - first release on CRAN.
