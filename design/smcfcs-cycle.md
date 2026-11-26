# Design note: a single-cycle SMC-FCS function for use from mice

Stef van Buuren, draft of 26 November 2026

Status: draft for discussion with Jonathan Bartlett. Based on smcfcs master at commit `ae95ab3`, which includes the fix for re-imputing missing Poisson outcomes (`cee0737`) and the fix for the outcome imputation variance (open question 1). CRAN has 2.0.3. Line numbers refer to `R/smcfcs.r` at `ae95ab3`. Steps 1–5 of the refactor are done on the branch `smcfcs-cycle`; section 6 lists the progress and section 2 describes the code as it now stands.

## 1. Goal

Users should be able to impute a substantive-model-compatible block inside `mice()`, together with ordinary FCS imputation of auxiliary variables, and keep the `mids` → `with()` → `pool()` workflow. smcfcs stays the reference implementation. mice calls it and does not reimplement it.

mice runs its own Gibbs loop and calls a block method once per iteration with the current completed data. To fit in, mice needs a function that runs exactly one SMC-FCS cycle on a given dataset.

## 2. Proposed changes to smcfcs

### 2.1 `smcfcs.cycle()` (exported)

``` r
smcfcs.cycle(data, r,                                         # what to update
             smtype, smformula, method, predictorMatrix = NULL,
             errorProneMatrix = NULL, restrictions = NULL,    # the model
             rjlimit = 1000, noisy = FALSE, verbose = FALSE,  # control
             cycleState = NULL, ...)                          # state
```

The help page will say plainly that this is a low-level function for package developers, and that most users want `smcfcs()`. The name and argument order are discussed in open question 5.

- `data` is the current completed data frame. The partially observed variables and the outcome contain current imputations, with no `NA`s.
- `r` is the original response-indicator matrix, n × p, with the same column order as `data`. It may be 0/1 or logical, so mice's logical `r` passes straight through. Only cells with `r == 0` are re-imputed.
- `smtype`, `smformula`, `method`, `predictorMatrix`, `rjlimit`, `noisy`, `errorProneMatrix`, `restrictions` and `...` have the same meaning as in `smcfcs.core()`.
- `verbose` controls the descriptive messages (outcome, passive variables, "Imputing: x using ..."). `noisy` still controls the model summaries.
- `cycleState` is `NULL` on the first call, and otherwise the value returned by the previous call (section 2.3).

The function returns `list(data = <updated data>, cycleState = <state>)`.

The body is the `for (cyclenum in 1:numit)` loop, moved essentially unchanged: the update of passive variables, the loop over `partialVars`, and the outcome update for lm, logistic, brlogistic and poisson (and for flexsurv time imputation). Inside the body, `imputations[[imp]]` becomes `data`. Conditions on `imp == 1 & cyclenum == 1` become `verbose`. The test `cyclenum == 1` for flexsurv becomes "no initial values from a previous cycle yet" (section 2.3).

**The original data.** The body also needs the data as originally observed, not only `r`: for restrictions, for the case-cohort and nested case-control quantities, and for flexsurv time imputation, which starts from the censoring times. In the last case `data` and `r` are not enough, because the imputed times overwrite the censoring times in `data` while `r` still marks them as observed. So `smcfcs.setup()` (section 2.2) keeps `originaldata`, and the setup is stored in `cycleState`. On a first call without a setup, `smcfcs.cycle()` rebuilds `originaldata` from `data` by setting the cells with `r == 0` to `NA`. That is exact on a first call, because no times have been imputed yet. We checked it on the reference cases (section 4.2): the rebuilt data equal the original data, apart from 0/1 columns stored as integer, which come back as double with the same values.

Implemented in step 4 as an internal function. It is exported, with a help page, in step 6.

**Scope.** `smcfcs.core()` will call `smcfcs.cycle()` for every `smtype`, so the function will run all model types. The "first version" restriction to lm, logistic, brlogistic, poisson, coxph and compet applies to the mice interface, and to what we document and test as supported for direct use. Case-cohort, nested case-control, dtsam, weibull, flexsurv and latnorm stay available through `smcfcs()`. They can be opened to mice later, once their extra requirements are handled: subcohort and risk-set data, flexsurv warm starts across mice iterations, and error-prone replicate measurements.

### 2.2 Shared internal setup: `smcfcs.setup()`

Before imputing anything, `smcfcs.core()` does some preparation (lines 233–536):

- It works out the role of each column (outcome, incomplete covariate, passive, fully observed).
- It checks the inputs.
- For Cox and competing-risks models, it fits a null model that fixes the event times used for the baseline hazard.

We move this preparation into an internal function, `smcfcs.setup()`, which both `smcfcs.core()` and `smcfcs.cycle()` use. It draws no random numbers, so results are unchanged.

**Doing the preparation only once.** mice calls `smcfcs.cycle()` many times on the same data. Repeating the preparation each time is wasteful: for a Cox model it costs about a fifth of a cycle (2.8 ms against 13 ms at n = 1,000). So the first call does the preparation and keeps the result in `cycleState`. Later calls reuse it.

**Guarding against a mismatch.** A caller could change, say, `smformula` and still pass in the old `cycleState`. To catch that, `cycleState` also records the arguments the preparation was based on. If they differ, the preparation is simply redone. The check costs almost nothing.

**No change for `smcfcs()` users.** `smcfcs()` does the preparation once, as it does now. Results, seeds and arguments are unchanged.

The contents of the stored preparation stay internal. If needed later, an optional `setup` argument can be added without breaking existing code. With it, a caller could prepare once and share the result across mice's m chains.

**As implemented (step 3).** `smcfcs.setup(originaldata, smtype, smformula, method, errorProneMatrix, restrictions, verbose, ...)` contains the preparation code unchanged, including the messages about the variables' roles. It returns a named list of the 35 quantities the imputation needs, among them `n`, `r`, the column numbers, `partialVars`, the null-model quantities, the parsed restrictions and `originaldata`. Quantities that do not apply to the substantive model, such as `numCauses` for anything but `compet`, are `NULL`. `smcfcs.core()` and `smcfcs.cycle()` put the list into their own environment with `list2env()`, so the code that uses the quantities is unchanged. Unpacking is needed in any case for one of them: `smformula2` for case-cohort contains `cluster(smcfcsid)`, so `smcfcsid` has to exist where the formula is evaluated, although the code never refers to it by name. `codetools::findGlobals()` confirms that every variable the two functions read without defining it comes from the setup list. So none of them can fall back to an object of the same name in the user's workspace.

**As implemented (step 5).** The key that is compared is built by an internal function, `smcfcs.setupKey()`. It holds `smtype`, `smformula`, `method`, `errorProneMatrix`, `restrictions`, the extra arguments in `...`, and `r`. `r` is included because the setup is prepared from `data` with the cells `r == 0` set to `NA`. It is stored as a logical matrix without dimnames, so a logical `r` and a 0/1 `r` give the same key. `predictorMatrix`, `rjlimit`, `noisy` and `verbose` are left out, because the setup does not use them (`verbose` only controls its messages). `smcfcs.cycle()` builds the key on every call and compares it with `identical()` to the key stored in `cycleState$setupKey`. If they differ, or there is no setup yet, it prepares the setup afresh, as on a first call. `smcfcs.core()` puts the key in the state together with the setup, so its cycles reuse the setup and the preparation still runs once per `smcfcs()` call. Changing the formula, `method` or `r` between calls redoes the setup; feeding back the state with the same arguments does not.

Two details on a redo. The flexsurv initial values and the last flexsurv fit are dropped, because they belong to the old model; the failure counts are kept. And for flexsurv with `imputeTimes = TRUE`, `originaldata` cannot be rebuilt once a cycle has run, because `data` then holds imputed times instead of the censoring times (section 2.1). In that case `smcfcs.cycle()` stops with an error that says to start again with `cycleState = NULL`. A formula is compared with its environment, so a caller that builds the same formula anew inside a function on each call (rather than reusing one object, as mice does) would redo the setup each time. A character `smformula`, or a formula created once, has no such problem. That is wasteful but gives correct results.

### 2.3 Explicit state

These quantities carry over between cycles:

| Quantity | At `ae95ab3` | Now (steps 2 and 4) |
|-----------------|---------------------|----------------------------------|
| `rjFailCount` | counter in `smcfcs.core()` | `cycleState$rjFailCount`, accumulated |
| `flexsurvFailCount` | counter, updated with `<<-` from a `tryCatch()` error handler | `cycleState$flexsurvFailCount`, accumulated; a failed fit returns `NULL` and is counted by the caller of `tryCatch()` |
| `flexsurvEsts` | initial values for the next fit when `cyclenum > 1` | `cycleState$flexsurvEsts`; `NULL` means first cycle |
| last flexsurv fit (`ymod`) | kept in `smcfcs.core()`; reused when a fit fails, and its knots reused with `originalKnots = TRUE` | `cycleState$flexsurvFit` |
| substantive-model coefficients at end of cycle | written to `smCoefIter[imp, , cyclenum]` | `cycleState$smCoef`; `smcfcs.core()` fills `smCoefIter` |
| preparation (setup) | done once per `smcfcs()` call | `cycleState$setup`, redone if the arguments change (2.2, step 5); the arguments it was based on are in `cycleState$setupKey` |

The coefficients are `outcomeModBeta`, except for `compet`, where they are the cause-specific vectors joined into one (`totalCoefVec`). They are taken after the last partially observed variable (lines 978–1000), for every model. For flexsurv with `imputeTimes = TRUE` they are then overwritten with the draw from the time imputation (line 1375). The outcome update for lm, logistic, brlogistic and poisson draws new coefficients too, but those are not recorded. `smcfcs.cycle()` keeps exactly these values in `cycleState$smCoef`.

**Between imputations.** `smcfcs.core()` starts each imputation with a new `cycleState` that holds the setup, the two failure counts and the last flexsurv fit. The flexsurv initial values are reset, so the first cycle of each imputation fits without them, as before. The last flexsurv fit is carried over on purpose, to keep results identical: if the first fit of an imputation fails, the code at `ae95ab3` silently reuses the model fitted in the previous imputation (open question 7).

**A bug found along the way.** At `ae95ab3`, `smCoefIter` was allocated in the time-imputation step if `exists("smCoefIter")` was false. `exists()` also searches the global environment, so an object named `smCoefIter` in the user's workspace made `smcfcs.flexsurv()` with `imputeTimes = TRUE` and no missing covariates stop with "incorrect number of subscripts". Step 2 replaces the test with `is.null()` on a local variable, which fixes this. NEWS mentions it.

`smcfcs.cycle()` does not issue warnings. The caller does: `smcfcs.core()` issues the two existing warnings at the end, unchanged, and mice writes them to `loggedEvents`.

### 2.4 What `smcfcs.core()` becomes

As implemented in step 4, abridged:

``` r
setup <- smcfcs.setup(originaldata, smtype, smformula, method, ...)
list2env(setup, envir = environment())
rjFailCount <- 0; flexsurvFailCount <- 0; flexsurvFit <- NULL; smCoefIter <- NULL
for (imp in 1:m) {
  # initial imputation: unchanged (lines 541-661)
  state <- list(setup = setup, rjFailCount = rjFailCount,
                flexsurvFailCount = flexsurvFailCount, flexsurvFit = flexsurvFit)
  for (cyclenum in 1:numit) {
    if (noisy) message(paste("Iteration ", cyclenum))
    res <- smcfcs.cycle(imputations[[imp]], r, smtype, smformula, method, ...,
                        verbose = verbose && (imp == 1) && (cyclenum == 1),
                        cycleState = state)
    imputations[[imp]] <- res$data; state <- res$cycleState
    if (is.null(smCoefIter)) smCoefIter <- array(0, dim = c(m, length(state$smCoef), numit))
    smCoefIter[imp, , cyclenum] <- state$smCoef
  }
  rjFailCount <- state$rjFailCount; flexsurvFailCount <- state$flexsurvFailCount
  flexsurvFit <- state$flexsurvFit
}
# the two warnings and the return value: unchanged
```

Random numbers are drawn in the same order as now, because neither setup nor the bookkeeping draws any. A regression test (section 4.2) checks this exactly.

### 2.5 Messages

The descriptive `print()` calls (lines 493–496, 539, 624 and 682–684) become `message()` calls controlled by `verbose`. This applies to every caller, including `smcfcs()`, which keeps `verbose = TRUE` by default. The `noisy` model summaries stay as `print()`, because they are printed results rather than status messages. The one status line under `noisy`, `"Iteration <cyclenum>"` (line 665), becomes a `message()` as well, but stays controlled by `noisy`, so the existing behaviour of `noisy = TRUE` is kept. Within `smcfcs.cycle()` it moves to the caller, because the function does not know the iteration number: `smcfcs.core()` emits it before each call.

What changes for users:

- The lines appear without the `[1]` prefix and the quotes.
- `suppressMessages(smcfcs(...))` silences them. Today that takes `invisible(capture.output(...))`.
- In R Markdown and Quarto, the chunk option `message = FALSE` hides them.
- They go to stderr instead of stdout. So `capture.output()` and `sink()` no longer capture them, and a shell redirect needs `2>`.

Wrappers such as `smcfcs.finegray()` (see the comment at its line 185) and mice can then silence the messages cleanly or log them. The PR updates the test that relies on `expect_output()` ("Checking outcome model check for logistic models", `test_errorchecks.R`) to use `expect_message()`, and adds a NEWS entry. Keeping `print()` for `smcfcs()` only would leave those wrappers with the current problem, so we recommend the switch for all callers.

## 3. Proposed changes to mice

### 3.1 A block method, with a name that doesn't clash

Working name: `mice.impute.smc()`. The name is open, but it must not be `smcfcs`. miceadds already exports a `mice.impute.smcfcs()`, which is its own implementation, with different arguments. mice looks up `mice.impute.<method>` from inside its own namespace, so a mice function of that name would always win. Existing miceadds users with `method = "smcfcs"` would then silently get the new function. smcfcs goes in Suggests.

The substantive model can be given in two ways, discussed in open question 6. Option A, which we recommend, uses mice's own `formulas` argument:

``` r
blocks <- make.blocks(list(smc = c("y", "x1", "x2"), "aux1", "aux2"))
meth   <- make.method(data, blocks); meth["smc"] <- "smc"
imp <- mice(data, blocks = blocks, method = meth,
            formulas = list(smc = y ~ x1 + I(x1^2) + x2),
            blots = list(smc = list(smtype = "lm",
                                    smmethod = c(x1 = "norm", x2 = "logreg"),
                                    auxiliary = c("aux1", "aux2"))))
```

Option B passes the model as `blots = list(smc = list(smformula = "y ~ x1 + I(x1^2) + x2", ...))`, gives no formula for the block, and takes auxiliary variables from mice's `predictorMatrix`.

- The block holds the outcome and the incomplete covariates in the substantive model. Variables outside the block are imputed by ordinary mice methods. Inside a cycle they are treated as fixed at their current values. The block must be given explicitly through `blocks`. The wrapper checks that it contains exactly the incomplete outcome and the incomplete covariates in the substantive model, and stops with an error naming the variables to add or remove.
- Squares and interactions are written in the substantive-model formula, e.g. `I(x1^2)` or `x1:x2`, as in the smcfcs vignette. This is the recommended form (open question 4). A passive column is also possible. It must then be in the block, with a smcfcs-style expression without `~`, e.g. `smmethod = c(x1sq = "x1^2")`. A mice passive method (`"~I(x1^2)"`) is not suitable here, because smcfcs would not update it within a cycle.
- The argument names follow smcfcs (`smtype`, `rjlimit`, `predictorMatrix`, and under option B `smformula`) and are passed straight through. The per-variable smcfcs methods go in `smmethod`, because `method` already means something else in mice. Defaults come from each variable's type: numeric → `norm`, binary → `logreg`, ordered factor → `podds`, other factor → `mlogit`.
- The function has a `format` formal argument, which is how mice recognises a block method.
- It builds the smcfcs inputs as follows:
  - `method`: `smmethod` for the incomplete block variables and `""` elsewhere.
  - `r`: mice's `r` for the block columns, set to observed for all other columns.
  - `smformula`: the block's formula, converted to a string (option A), or `blots$smformula` (option B).
  - `predictorMatrix`: taken from `blots$predictorMatrix` if given. Otherwise, under option A, smcfcs's default plus the variables in `blots$auxiliary`, and under option B, derived from mice's predictor row for the block. Auxiliary variables outside the substantive model enter the covariate models only this way.
- It calls `smcfcs::smcfcs.cycle()` once per iteration, then returns the imputations for the block's `r == 0` cells.
- It refuses `smtype` values outside the first-version set with a clear error. In the first version it also requires `where == !r` within the block.

### 3.2 Sampler: give block methods their current imputations

Today \[mice `sampler.R`\] sets the block's missing cells to `NA` before calling a block method (`data[mis] <- NA`). SMC-FCS needs the current imputations, as do the rejection-sampling fallback and the warm starts. The proposal:

- If a block method has an `r` argument (`handles.arg(fm, "r")`), mice passes the completed `data` and `r`, and does not blank the cells. Other block methods behave as now.
- If the method has a `cycleState` argument, mice keeps one state per block and chain, passes it in, and takes the updated state back from `attr(imputes, "cycleState")`. At the end of `mice()`, the fail counts go to `loggedEvents`. With futuremice, each worker keeps its own chain's state.
- The sampler's local names (`newstate`, `yname`, `data`, `blocks`) stay unchanged, so that `miceadds` keeps working.

### 3.3 Differences from calling `smcfcs()` directly

- mice initialises the block with its own starting values (random draws from the observed values). smcfcs initialises a missing outcome with an improper fit of the substantive model. So the two give different results under the same seed. Seed equality is checked only for `smcfcs()` itself.
- `smCoefIter` has no counterpart in mice. Convergence is checked with mice's chain means and variances, as for any other method.
- When the outcome is missing, the outcome update refits the substantive model on rows with an observed outcome only. The covariate updates fit on all rows. This is deliberate, and the mice documentation will say so.

## 4. Tests

### 4.1 Why the current tests are not enough

The existing suite checks two things: that smcfcs runs without error, and that it is approximately unbiased. The bias tests use n = 10,000, accept a pooled estimate within about 0.1 of the truth, and are all `skip_on_cran()`. No test compares imputed values or `smCoefIter` with a known result. So a refactor could change every imputation and still pass, as long as the results stay unbiased. Errors that would slip through include:

- random numbers drawn in a different order, so that results under a fixed seed change, although section 2.4 promises they won't;
- a warm start (`flexsurvEsts`) or a fail count lost between cycles;
- a stale setup reused after an argument has changed;
- state not reset between imputations;
- a mistake in a path that no test runs at all: passive variables, `I()` terms in `smformula`, missing `lm` or `logistic` outcomes, `compet` called directly, the rejection-sampling and flexsurv fail-count warnings, and `noisy = TRUE`.

### 4.2 A frozen-reference regression test

The refactor is meant to change no results at all. That allows an exact test instead of a statistical one: run the old code once, store what it produces, and require the new code to produce the same objects, bit for bit.

**Step 1: generate the references before the refactor.** A script, `data-raw/make-references.R` (not shipped with the package), loads the package source at `ae95ab3` with `pkgload::load_all()`. It runs `smcfcs()` and its wrappers on a grid of 27 small cases (mostly n = 200, m = 2, numit = 3), each under its own fixed seed. For each case it saves one `.rds` file in `tests/testthat/reference/<blas>/`, where `<blas>` names the BLAS library in use (see Limits). The file holds the inputs (function, arguments including the data, seed), the results (`impDatasets`, `smCoefIter`, `smInfo` and any warnings), and the environment they were generated in. The grid covers:

- every `smtype`: lm, logistic, brlogistic, poisson, coxph, compet, weibull and dtsam, and the wrappers `smcfcs.casecohort()`, `smcfcs.nestedcc()`, `smcfcs.finegray()` and `smcfcs.flexsurv()`, including the imputation of censored times;
- every covariate method: norm, latnorm (with `errorProneMatrix`), logreg, brlogreg, poisson, podds and mlogit;
- a missing outcome for lm, logistic and poisson;
- a passive variable, and an `I()` term in `smformula`;
- restrictions on a factor covariate;
- a user-supplied `predictorMatrix`;
- a very small `rjlimit`, so that rejection sampling fails and the warning count becomes part of the reference;
- `noisy = TRUE`.

The data are simulated inside the script, so every case is small and self-contained. The whole grid runs in a few seconds and the references take about 300 KB. Because the inputs are stored with each reference, the test does not depend on the simulation code.

**Step 2: check every refactoring commit against the references.** A test file, `tests/testthat/test_reference.R`, re-runs each case from its stored inputs with the same seed and checks that `impDatasets`, `smCoefIter`, `smInfo` and the warnings are `identical()` to the stored versions. The script and the test share one function, in `tests/testthat/helper-reference.R`, that runs a case. It discards printed output and messages, because the refactor turns `print()` into `message()` on purpose (section 2.5). Because the test runs on every commit of the refactor, a difference points to the commit that caused it.

We tried this on the baseline. All 27 cases reproduce exactly. Run against the code before the fix for open question 1, exactly one case fails, `lm_missing_outcome`, which is the only case that runs the changed line.

All 27 cases still reproduce exactly after each of steps 1–5, and the full test suite passed after steps 1 and 5. No case in the grid has a failed flexsurv fit, so for steps 2 and 4 we also checked that path separately. We made `flexsurv::flexsurvspline()` fail on chosen calls, including the first fit of the second imputation, with and without `imputeTimes`. The code before and after each step gave identical results and identical warnings.

**How this prevents errors.** Any change in which random numbers are drawn, or in their order, changes the imputations, and the test fails. The same holds for lost or stale state, a wrong setup, or a broken code path anywhere in the grid. The comparison is exact, so there is no tolerance in which a small error can hide, and a failure never depends on chance. The bias tests stay and still check that the method is correct. The reference test checks that the refactor has not changed it.

**Limits.** Exact equality holds only in the environment the references were made in. Results can change without any change to smcfcs when one of these differs:

- the R version or `RNGkind()`;
- the BLAS or LAPACK library, because `MASS::mvrnorm()` uses an eigendecomposition, whose rounding depends on the library;
- the versions of the packages smcfcs calls: MASS, survival, VGAM, brglm2, flexsurv and kmi.

The BLAS matters most. Under Apple's Accelerate instead of R's reference BLAS, 26 of the 27 cases differ: some only by rounding (about 1e-15), most by 0.1 to 1.5 in the coefficients, because one rounding difference can flip a rejection-sampling decision, after which the random-number stream takes another path. Both are valid draws, but not the same ones. We therefore keep one set of references per BLAS (`reference/Rblas/` and `reference/vecLib/`, about 300 KB each), and the test picks the set that matches the linear algebra actually in use. It does not go by the library name that R reports, which can be misleading: an IDE session can report Accelerate, because the IDE loaded it into the process, while R's arithmetic still runs through R's own BLAS. Instead, each reference stores a fingerprint, a checksum of a matrix product, a Cholesky factor and an eigendecomposition of fixed numbers, and the test uses the set whose fingerprint matches the current session. Within that set it compares the rest of the stored environment with the current one and skips a case, with a message saying what differs, rather than fail. For the same reason it is `skip_on_cran()`. It is a guard for the refactor, run on the developer's machine, not a test that holds on every platform. After the refactor, the references can be regenerated whenever the environment changes. When results are changed on purpose later, the references are regenerated in the same commit, and the commit message says why. Open question 1 was settled before step 1 for this reason: its fix (`ae95ab3`) is part of the baseline.

### 4.3 Other tests

In smcfcs:

1.  **Frozen-reference regression test** (section 4.2).
2.  **`smcfcs.cycle()` contract.** These tests check that:
    - `cycleState = NULL` works;
    - feeding back `cycleState` over `numit` calls reproduces the `smcfcs.core()` loop;
    - cells with `r == 1` are never changed;
    - a logical `r` and a 0/1 `r` give identical results;
    - `verbose = FALSE` produces no output.
3.  The existing test suite passes, with one deliberate change: the `expect_output()` in `test_errorchecks.R` becomes `expect_message()` (section 2.5).

In mice:

- The block method recovers the true coefficients in a simulation where a covariate is missing at random depending on the outcome, under the first-version `smtype`s.
- Under the same seed, `mice()` with an SMC block plus `pmm` for auxiliary variables runs end-to-end through `with()` and `pool()`.
- An unsupported `smtype` gives a clear error.
- The `miceadds` example still runs.

## 5. Open questions for Jonathan

1.  **Outcome imputation variance (lm). Settled.** At line 1316 a missing outcome was drawn with `sd = sigmasq^0.5`, the residual variance of the fitted model, rather than the posterior draw `outcomeModResVar^0.5` computed on line 1312. This was not intended, and Jonathan fixed it in `ae95ab3`. The regression references (section 4.2) are generated after this fix.

2.  **`print()` → `message()` for all callers** (section 2.5). **Settled.** The change is cosmetic, and it gives users the standard way to silence smcfcs. Jonathan agreed to it, including the NEWS entry and the updated `expect_output()` test.

3.  **Keeping the preparation in `cycleState`** (section 2.2). **Settled.** The first call prepares, later calls reuse the result, and a change of arguments triggers a fresh preparation. This keeps a single exported function, and the stored preparation stays internal, so it can change freely. Jonathan agreed to it.

4.  **Squares and interactions: write them in the formula.** A model with `x^2` can be specified in two ways in smcfcs. One is directly in the formula, `smformula = "y ~ x + I(x^2)"`, as in the vignette. The other is a separate column `xsq` with the passive method `"x^2"`. Both give the same substantive model. For the mice interface we recommend the formula form, for three reasons:

    - The user writes the same formula for imputation and for the analysis (`with(imp, lm(y ~ x + I(x^2)))`), with no extra column to create.
    - A passive column brings no extra benefit. By default it doesn't enter the imputation models of other incomplete covariates either.
    - In smcfcs's default predictor choice, a fully observed `z` that appears only inside a term like `I(z^2)` is not used as a predictor of the incomplete covariates. A passive column `zsq` doesn't help here either. Under option B of question 6, the wrapper takes the predictors from mice's `predictorMatrix`, where `z` is included by default. Under option A, the user adds `z` to `auxiliary`.

    Passive columns would still work for users who want them. Do you agree with recommending the formula form?

5.  **Name and argument order of `smcfcs.cycle()`.** Our proposals:

    - **Name: `smcfcs.cycle()`.** It matches the existing family (`smcfcs.casecohort()`, `smcfcs.flexsurv()`, ...), so it is easy to find. Because it could be mistaken for another way to impute, its help page will make clear it is a low-level building block. It will not be added to the functions accepted by `smcfcs.parallel()`. We considered `smcfcs_cycle()` and `smcfcs.step()`. The first would be the only exported snake_case name, and the second has no real advantage.
    - **Order: grouped by purpose** (section 2.1). First comes what to update (`data`, `r`), then the model (`smtype`, `smformula`, `method`, `predictorMatrix`, `errorProneMatrix`, `restrictions`), then the controls (`rjlimit`, `noisy`, `verbose`), and last the state (`cycleState`). The alternative is to copy the order of `smcfcs()` and add `r` at the end. That would put a required argument after optional ones, only to match positions in a function that callers use with named arguments. Neither choice changes any existing function.
    - **`data`, not `originaldata`**, because the input is the completed data. The return value is `list(data, cycleState)`.
    - **`verbose` next to `noisy`.** `noisy` keeps its meaning (model summaries), and `verbose` controls the status messages. They will be documented side by side.

    Any preferences, given how you name things in smcfcs?

6.  **Where the mice user writes the substantive model.** mice has its own `formulas` argument, with one formula per block. When a block has a formula, mice passes it to the block method unchanged, so it can serve as the substantive model. We checked this with mice 3.20.0, including `Surv(t, d) ~ x + z`. mice then passes no predictor row. When a block has no formula, mice passes the predictor row instead, and the model has to come from `blots$smformula`. This gives two options.

    **Option A: the block's formula is the substantive model.** `formulas = list(smc = y ~ x1 + I(x1^2) + x2)`. Auxiliary variables are listed in `blots$auxiliary`, or given in full through `blots$predictorMatrix`.

    Upsides:

    - The model is written where mice users already write models, as a formula rather than a string, and in the same form as the later analysis, `with(imp, lm(y ~ x1 + I(x1^2) + x2))`.
    - There is one way to specify the model, so there is nothing to check for consistency.
    - Auxiliary variables are an explicit choice. That matters, because an auxiliary variable in the covariate models but not in the substantive model implies that the outcome is independent of it given the covariates. smcfcs makes the same assumption.

    Downsides:

    - mice reads a block formula differently. Its left-hand side normally lists the variables to impute, but here it holds the outcome, while the incomplete covariates are on the right-hand side. Without explicit `blocks`, mice builds the wrong block: `y` alone for `y ~ x1 + x2`, and the fully observed `t, d` for `Surv(t, d) ~ x + z`. So `blocks` is required, and the wrapper must check it.
    - mice still builds a row of `predictorMatrix` for the block, all ones, but never uses it. A user who edits it gets no effect and no warning. The help page has to say so.
    - The default has no auxiliary variables. That is smcfcs's default, but not what mice users expect, since in mice every other variable is a predictor by default.

    **Option B: `blots$smformula`, with predictors from mice's predictor row.** The block has no formula. The wrapper converts mice's predictor row into smcfcs's `predictorMatrix`.

    Upsides:

    - `formulas` keeps its usual mice meaning, and the block is set up like any other mice block.
    - Auxiliary variables follow mice's defaults and can be edited with the usual tools, such as `quickpred()`.
    - The argument name matches smcfcs, so it is passed straight through.

    Downsides:

    - The model is a string inside `blots`, apart from the arguments where mice users normally describe models.
    - By default every other variable becomes an auxiliary variable. That imposes the independence assumption above on all of them, often without the user noticing.
    - `blocks` still has to be given explicitly, because mice cannot infer the block from a string.

    A third possibility is to support both, with an error if the two disagree. That gives the most flexibility, at the cost of two code paths to document and test.

    We recommend option A: one way to specify the model, written in mice's own place for models, with auxiliary variables chosen explicitly. Do you agree, or do you prefer mice's predictor defaults (option B)?

7.  **A failed flexsurv fit reuses the previous model** (section 2.3). When `flexsurvspline()` fails, the code counts the failure and carries on with the previous fitted model. Within an imputation that is the fit from the previous variable or cycle. For the first fit of the second or a later imputation, it is the last fit of the previous imputation, made on other imputed data. For the first fit of the first imputation there is none, and smcfcs stops with an error. At `ae95ab3` the error is "object 'ymod' not found". Since step 4 the previous fit is `NULL` rather than missing, and the error comes from inside flexsurv ("non-numeric matrix extent"), which says even less. Apart from that message, the refactor keeps all of this, so that results stay identical. Is reusing a model from another imputation intended? If not, a later change could refit, or stop with a clear error, and regenerate the references in the same commit.

## 6. Progress

Each step is a separate commit on the branch `smcfcs-cycle` of stefvanbuuren/smcfcs. After each one the frozen-reference test (section 4.2) passes.

| Step | What | Status |
|------|------|--------|
| 1 | Status lines via `message()`, new `verbose` argument (section 2.5) | done, `7205b69` |
| 2 | State carried explicitly instead of via `<<-` and `exists()` (section 2.3) | done, `981297c` |
| 3 | Preparation moved into `smcfcs.setup()` (section 2.2) | done, `260bdec` |
| 4 | One cycle moved into `smcfcs.cycle()` (sections 2.1, 2.4) | done, `ea089a2` |
| 5 | Redo the stored setup when the arguments change (section 2.2) | done |
| 6 | Contract tests (section 4.3), help page, export | to do |