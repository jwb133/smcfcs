# Design note: a single-cycle SMC-FCS function for use from mice

Stef van Buuren, draft of 18 November 2026

Status: draft for discussion with Jonathan Bartlett. Based on smcfcs master at
commit `cee0737`, which includes the fix for re-imputing missing Poisson
outcomes. CRAN has 2.0.3.

## 1. Goal

Users should be able to impute a substantive-model-compatible block inside
`mice()`, together with ordinary FCS imputation of auxiliary variables, and
keep the `mids` → `with()` → `pool()` workflow. smcfcs stays the reference
implementation. mice calls it and does not reimplement it.

mice runs its own Gibbs loop and calls a block method once per iteration with
the current completed data. To fit in, mice needs a function that runs exactly
one SMC-FCS cycle on a given dataset.

## 2. Proposed changes to smcfcs

### 2.1 `smcfcs.cycle()` (exported)

```r
smcfcs.cycle(data, r,                                         # what to update
             smtype, smformula, method, predictorMatrix = NULL,
             errorProneMatrix = NULL, restrictions = NULL,    # the model
             rjlimit = 1000, noisy = FALSE, verbose = FALSE,  # control
             cycleState = NULL, ...)                          # state
```

The help page will say plainly that this is a low-level function for package
developers, and that most users want `smcfcs()`. The name and argument order
are discussed in open question 5.

- `data` is the current completed data frame. The partially observed
  variables and the outcome contain current imputations, with no `NA`s.
- `r` is the original response-indicator matrix, n × p, with the same column
  order as `data`. It may be 0/1 or logical, so mice's logical `r` passes
  straight through. Only cells with `r == 0` are re-imputed.
- `smtype`, `smformula`, `method`, `predictorMatrix`, `rjlimit`, `noisy`,
  `errorProneMatrix`, `restrictions` and `...` have the same meaning as in
  `smcfcs.core()`.
- `verbose` controls the descriptive messages (outcome, passive variables,
  "Imputing: x using ..."). `noisy` still controls the model summaries.
- `cycleState` is `NULL` on the first call, and otherwise the value returned
  by the previous call (section 2.3).

The function returns `list(data = <updated data>, cycleState = <state>)`.

The body is the current `for (cyclenum in 1:numit)` loop, moved essentially
unchanged: the update of passive variables, the loop over `partialVars`, and
the outcome update for lm, logistic, brlogistic and poisson (and for flexsurv
time imputation). Inside the body, `imputations[[imp]]` becomes `data`.
Conditions on `imp == 1 & cyclenum == 1` become `verbose`.

**Scope.** `smcfcs.core()` will call `smcfcs.cycle()` for every `smtype`, so
the function will run all model types. The "first version" restriction to
lm, logistic, brlogistic, poisson, coxph and compet applies to the mice
interface, and to what we document and test as supported for direct use.
Case-cohort, nested case-control, dtsam, weibull, flexsurv and latnorm stay
available through `smcfcs()`. They can be opened to mice later, once their
extra requirements are handled: subcohort and risk-set data, flexsurv warm
starts across mice iterations, and error-prone replicate measurements.

### 2.2 Shared internal setup: `smcfcs.setup()`

Before imputing anything, `smcfcs.core()` does some preparation (lines
233–536):

- It works out the role of each column (outcome, incomplete covariate,
  passive, fully observed).
- It checks the inputs.
- For Cox and competing-risks models, it fits a null model that fixes the
  event times used for the baseline hazard.

We move this preparation into an internal function, `smcfcs.setup()`, which
both `smcfcs.core()` and `smcfcs.cycle()` use. It draws no random numbers, so
results are unchanged.

**Doing the preparation only once.** mice calls `smcfcs.cycle()` many times on
the same data. Repeating the preparation each time is wasteful: for a Cox
model it costs about a fifth of a cycle (2.8 ms against 13 ms at n = 1,000).
So the first call does the preparation and keeps the result in `cycleState`.
Later calls reuse it.

**Guarding against a mismatch.** A caller could change, say, `smformula` and
still pass in the old `cycleState`. To catch that, `cycleState` also records
the arguments the preparation was based on. If they differ, the preparation
is simply redone. The check costs almost nothing.

**No change for `smcfcs()` users.** `smcfcs()` does the preparation once, as
it does now. Results, seeds and arguments are unchanged.

The contents of the stored preparation stay internal. If needed later, an
optional `setup` argument can be added without breaking existing code. With
it, a caller could prepare once and share the result across mice's m chains.

### 2.3 Explicit state

These quantities carry over between cycles today:

| Quantity               | Today                                 | Proposed                                                       |
|------------------------|---------------------------------------|----------------------------------------------------------------|
| `rjFailCount`          | counter in `smcfcs.core()`            | `cycleState$rjFailCount`, accumulated                          |
| `flexsurvFailCount`    | counter, updated with `<<-`           | `cycleState$flexsurvFailCount`                                 |
| `flexsurvEsts`         | warm start when `cyclenum > 1`        | `cycleState$flexsurvEsts`; `NULL` means first cycle            |
| `outcomeModBeta` at end of cycle | written to `smCoefIter[imp, , cyclenum]` | `cycleState$smCoef`; `smcfcs.core()` fills `smCoefIter` |
| preparation (setup)    | done once per `smcfcs()` call         | `cycleState$setup`, redone if the arguments change (2.2)       |

`smcfcs.cycle()` does not issue warnings. The caller does: `smcfcs.core()`
issues the two existing warnings at the end, unchanged, and mice writes them
to `loggedEvents`.

### 2.4 What `smcfcs.core()` becomes

```r
setup <- smcfcs.setup(...)
for (imp in 1:m) {
  # initial imputation: unchanged (lines 541-661)
  state <- list(setup = setup)
  for (cyclenum in 1:numit) {
    res <- smcfcs.cycle(imputations[[imp]], r, ..., verbose = (imp == 1 && cyclenum == 1),
                        cycleState = state)
    imputations[[imp]] <- res$data; state <- res$cycleState
    smCoefIter[imp, , cyclenum] <- state$smCoef   # allocated after first cycle
  }
  # accumulate fail counts
}
```

Random numbers are drawn in the same order as now, because neither setup nor
the bookkeeping draws any. A regression test (section 4) checks this exactly.

### 2.5 Messages

The descriptive `print()` calls (lines 493–496, 539, 624 and 682–684) become
`message()` calls controlled by `verbose`. This applies to every caller,
including `smcfcs()`, which keeps `verbose = TRUE` by default. The `noisy`
model summaries stay as `print()`, because they are printed results rather
than status messages.

What changes for users:

- The lines appear without the `[1]` prefix and the quotes.
- `suppressMessages(smcfcs(...))` silences them. Today that takes
  `invisible(capture.output(...))`.
- In R Markdown and Quarto, the chunk option `message = FALSE` hides them.
- They go to stderr instead of stdout. So `capture.output()` and `sink()` no
  longer capture them, and a shell redirect needs `2>`.

Wrappers such as `smcfcs.finegray()` (see the comment at its line 185) and
mice can then silence the messages cleanly or log them. The PR updates the
test that relies on `expect_output()` ("Checking outcome model check for
logistic models", `test_errorchecks.R`) to use `expect_message()`, and adds a
NEWS entry. Keeping `print()` for `smcfcs()` only would leave those wrappers
with the current problem, so we recommend the switch for all callers.

## 3. Proposed changes to mice

### 3.1 A block method, with a name that doesn't clash

Working name: `mice.impute.smc()`. The name is open, but it must not be
`smcfcs`. miceadds already exports a `mice.impute.smcfcs()`, which is its own
implementation, with different arguments. mice looks up
`mice.impute.<method>` from inside its own namespace, so a mice function of
that name would always win. Existing miceadds users with
`method = "smcfcs"` would then silently get the new function. smcfcs goes in
Suggests.

```r
blocks <- make.blocks(list(smc = c("y", "x1", "x2"), "aux1", "aux2"))
meth   <- make.method(data, blocks); meth["smc"] <- "smc"
imp <- mice(data, blocks = blocks, method = meth,
            blots = list(smc = list(smtype = "lm",
                                    smformula = "y ~ x1 + I(x1^2) + x2",
                                    smmethod = c(x1 = "norm", x2 = "logreg"))))
```

- The block holds the outcome and the incomplete covariates in the
  substantive model. Variables outside the block are imputed by ordinary
  mice methods. Inside a cycle they are treated as fixed at their current
  values.
- Squares and interactions are written in `smformula`, e.g. `I(x1^2)` or
  `x1:x2`, as in the smcfcs vignette. This is the recommended form (open
  question 4). A passive column is also possible. It must then be in the
  block, with a smcfcs-style expression without `~`, e.g.
  `smmethod = c(x1sq = "x1^2")`. A mice passive method (`"~I(x1^2)"`) is not
  suitable here, because smcfcs would not update it within a cycle.
- The argument names follow smcfcs (`smtype`, `smformula`, `rjlimit`,
  `predictorMatrix`) and are passed straight through. The per-variable
  smcfcs methods go in `smmethod`, because `method` already means something
  else in mice. Defaults come from each variable's type: numeric → `norm`,
  binary → `logreg`, ordered factor → `podds`, other factor → `mlogit`.
- The function has a `format` formal argument, which is how mice recognises a
  block method.
- It builds the smcfcs inputs as follows:
  - `method`: `smmethod` for the incomplete block variables and `""`
    elsewhere.
  - `r`: mice's `r` for the block columns, set to observed for all other
    columns.
  - `predictorMatrix`: derived from mice's predictor row for the block, or
    taken from `blots` if given. Auxiliary variables outside the substantive
    model enter the covariate models only this way.
- It calls `smcfcs::smcfcs.cycle()` once per iteration, then returns the
  imputations for the block's `r == 0` cells.
- It refuses `smtype` values outside the first-version set with a clear
  error. In the first version it also requires `where == !r` within the
  block.

### 3.2 Sampler: give block methods their current imputations

Today [mice `sampler.R`] sets the block's missing cells to `NA` before
calling a block method (`data[mis] <- NA`). SMC-FCS needs the current
imputations, as do the rejection-sampling fallback and the warm starts. The
proposal:

- If a block method has an `r` argument (`handles.arg(fm, "r")`), mice passes
  the completed `data` and `r`, and does not blank the cells. Other block
  methods behave as now.
- If the method has a `cycleState` argument, mice keeps one state per block
  and chain, passes it in, and takes the updated state back from
  `attr(imputes, "cycleState")`. At the end of `mice()`, the fail counts go
  to `loggedEvents`. With futuremice, each worker keeps its own chain's state.
- The sampler's local names (`newstate`, `yname`, `data`, `blocks`) stay
  unchanged, so that `miceadds` keeps working.

### 3.3 Differences from calling `smcfcs()` directly

- mice initialises the block with its own starting values (random draws from
  the observed values). smcfcs initialises a missing outcome with an improper
  fit of the substantive model. So the two give different results under the
  same seed. Seed equality is checked only for `smcfcs()` itself.
- `smCoefIter` has no counterpart in mice. Convergence is checked with
  mice's chain means and variances, as for any other method.
- When the outcome is missing, the outcome update refits the substantive
  model on rows with an observed outcome only. The covariate updates fit on
  all rows. This is deliberate, and the mice documentation will say so.

## 4. Tests

In smcfcs:

1. **Regression against the pre-refactor code.** Before touching
   `smcfcs.core()`, generate reference results on master (`cee0737`) for
   small cases: n = 200, m = 2, numit = 3, under fixed seeds. Cover every
   `smtype` and every covariate `method`, including restrictions and latnorm.
   Store them as `.rds` files in `tests/testthat/`. After the refactor, check
   that `impDatasets` and `smCoefIter` are `identical()`.
2. **`smcfcs.cycle()` contract.** These tests check that:
   - `cycleState = NULL` works;
   - feeding back `cycleState` over `numit` calls reproduces the
     `smcfcs.core()` loop;
   - cells with `r == 1` are never changed;
   - a logical `r` and a 0/1 `r` give identical results;
   - `verbose = FALSE` produces no output.
3. The existing test suite passes, with one deliberate change: the
   `expect_output()` in `test_errorchecks.R` becomes `expect_message()`
   (section 2.5).

In mice:

- The block method recovers the true coefficients in a simulation where a
  covariate is missing at random depending on the outcome, under the
  first-version `smtype`s.
- Under the same seed, `mice()` with an SMC block plus `pmm` for auxiliary
  variables runs end-to-end through `with()` and `pool()`.
- An unsupported `smtype` gives a clear error.
- The `miceadds` example still runs.

## 5. Open questions for Jonathan

1. **Outcome imputation variance (lm).** At line 1316 a missing outcome is
   drawn with `sd = sigmasq^0.5`, the residual variance of the fitted model,
   rather than the posterior draw `outcomeModResVar^0.5` computed on line
   1312. Is that intended? If not, fixing it changes results, so it should
   land before the regression references are generated, or in its own
   commit.
2. **`print()` → `message()` for all callers** (section 2.5). We recommend
   it: the change is cosmetic, and it gives users the standard way to silence
   smcfcs. Are you happy with it, including the NEWS entry and the updated
   `expect_output()` test?
3. **Keeping the preparation in `cycleState`** (section 2.2). We propose that
   the first call prepares, later calls reuse the result, and a change of
   arguments triggers a fresh preparation. This keeps a single exported
   function, and the stored preparation stays internal, so you remain free
   to change it. Does that work for you?
4. **Squares and interactions: write them in the formula.** A model with
   `x^2` can be specified in two ways in smcfcs. One is directly in the
   formula, `smformula = "y ~ x + I(x^2)"`, as in the vignette. The other is
   a separate column `xsq` with the passive method `"x^2"`. Both give the
   same substantive model. For the mice interface we recommend the formula
   form, for three reasons:
   - The user writes the same formula for imputation and for the analysis
     (`with(imp, lm(y ~ x + I(x^2)))`), with no extra column to create.
   - A passive column brings no extra benefit. By default it doesn't enter
     the imputation models of other incomplete covariates either.
   - In smcfcs's default predictor choice, a fully observed `z` that appears
     only inside a term like `I(z^2)` is not used as a predictor of the
     incomplete covariates. In mice this doesn't arise, because the wrapper
     takes the predictors from mice's `predictorMatrix` (section 3.1), and
     there `z` is included by default.

   Passive columns would still work for users who want them. Do you agree
   with recommending the formula form?
5. **Name and argument order of `smcfcs.cycle()`.** Our proposals:
   - **Name: `smcfcs.cycle()`.** It matches the existing family
     (`smcfcs.casecohort()`, `smcfcs.flexsurv()`, ...), so it is easy to
     find. Because it could be mistaken for another way to impute, its help
     page will make clear it is a low-level building block. It will not be
     added to the functions accepted by `smcfcs.parallel()`. We considered
     `smcfcs_cycle()` and `smcfcs.step()`. The first would be the only
     exported snake_case name, and the second has no real advantage.
   - **Order: grouped by purpose** (section 2.1). First comes what to update
     (`data`, `r`), then the model (`smtype`, `smformula`, `method`,
     `predictorMatrix`, `errorProneMatrix`, `restrictions`), then the
     controls (`rjlimit`, `noisy`, `verbose`), and last the state
     (`cycleState`). The alternative is to copy the order of `smcfcs()` and
     add `r` at the end. That would put a required argument after optional
     ones, only to match positions in a function that callers use with named
     arguments. Neither choice changes any existing function.
   - **`data`, not `originaldata`**, because the input is the completed data.
     The return value is `list(data, cycleState)`.
   - **`verbose` next to `noisy`.** `noisy` keeps its meaning (model
     summaries), and `verbose` controls the status messages. They will be
     documented side by side.

   Any preferences, given how you name things in smcfcs?
