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

## 2. Lessons from `miceadds::mice.impute.smcfcs()`

miceadds (3.20.10) already has a `mice.impute.smcfcs()`. It does not call
smcfcs. It is a separate implementation with these features:

- It is univariate. Each incomplete variable gets method `"smcfcs"`, and the
  substantive model is passed through mice's `...` as `sm`.
- Each call runs 20 Metropolis–Hastings steps with a proposal width that
  adapts across iterations. smcfcs uses direct sampling or rejection sampling.
- Parameter uncertainty comes from a bootstrap of the rows where the target
  variable is observed.
- The substantive model can be normal, logistic, Box-Cox, Yeo-Johnson or
  lognormal. There is no Cox, Poisson or competing-risks model.
- It finds its state by searching parent frames for the sampler's local
  variables (`newstate`, `yname`, `data`, `blocks`). It keeps its tuning
  factor across iterations by calling `assign()` into mice's frame.

Three lessons for this design.

**(a) Fit the models on all rows.** miceadds fits both the covariate model and
the substantive model only on rows where the covariate is observed. In the
simulation below, x is missing at random depending on y (47% missing). The
estimates are biased towards complete-case analysis. Changing that one line
so it bootstraps all rows removes most of the bias. Each entry is the mean
over 10 datasets (n = 1000, m = 5, 10 iterations).

Model: `y = 0.5 z + x + x^2 + e`, with `P(x missing) = plogis(-1 + y)`.

| Method                                | z    | x    | x²   |
|---------------------------------------|------|------|------|
| Truth                                 | 0.50 | 1.00 | 1.00 |
| Complete cases                        | 0.42 | 0.85 | 0.82 |
| `miceadds::mice.impute.smcfcs`        | 0.48 | 0.86 | 0.88 |
| miceadds, bootstrap over all rows     | 0.51 | 0.96 | 1.03 |
| `smcfcs::smcfcs` 2.0.3                | 0.52 | 1.06 | 1.06 |

In the smcfcs runs, rejection sampling failed 50–150 times per call with
`rjlimit = 1000`. A failed draw keeps the previous value. Counts like these
must reach the mice user, which is why section 3.3 returns them in
`cycleState`.

**(b) Keep state in a supported place.** State that persists across
iterations should be passed in and returned. It should not be stored by
reaching into the caller's frame.

**(c) Avoid a name clash.** mice looks up `mice.impute.<method>` from inside
its own namespace, so a function of that name in mice always wins over the
one in miceadds. If mice exported `mice.impute.smcfcs()`, every existing
miceadds user with `method = "smcfcs"` would silently get the new function,
which takes different arguments. The mice method therefore needs a different
name (section 4.1).

## 3. Proposed changes to smcfcs

### 3.1 `smcfcs.cycle()` (exported)

```r
smcfcs.cycle(data, r, smtype, smformula, method,
             predictorMatrix = NULL, rjlimit = 1000, noisy = FALSE,
             errorProneMatrix = NULL, restrictions = NULL,
             verbose = FALSE, cycleState = NULL, ...)
```

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
  by the previous call (section 3.3).

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

### 3.2 Shared internal setup: `smcfcs.setup()`

This is lines 233–536 of `smcfcs.core()`, moved into a function that returns
a list. The list holds:

- the column roles: `outcomeCol`, `timeCol`, `dCol`, `partialVars`,
  `passiveVars`, `fullObsVars` and `smcovcols`;
- the input checks;
- the Cox and competing-risks null-model baseline hazard indices (`H0indices`)
  and `linpred`;
- the case-cohort and nested case-control quantities;
- the parsed restrictions.

None of this consumes random numbers. Both `smcfcs.core()` and
`smcfcs.cycle()` call it. `smcfcs.cycle()` caches the result in
`cycleState$setup`, so that mice does not refit the Cox null model in every
iteration.

### 3.3 Explicit state

These quantities carry over between cycles today:

| Quantity               | Today                                 | Proposed                                                       |
|------------------------|---------------------------------------|----------------------------------------------------------------|
| `rjFailCount`          | counter in `smcfcs.core()`            | `cycleState$rjFailCount`, accumulated                          |
| `flexsurvFailCount`    | counter, updated with `<<-`           | `cycleState$flexsurvFailCount`                                 |
| `flexsurvEsts`         | warm start when `cyclenum > 1`        | `cycleState$flexsurvEsts`; `NULL` means first cycle            |
| `outcomeModBeta` at end of cycle | written to `smCoefIter[imp, , cyclenum]` | `cycleState$smCoef`; `smcfcs.core()` fills `smCoefIter` |
| setup                  | recomputed per call                   | `cycleState$setup`                                             |

`smcfcs.cycle()` does not issue warnings. The caller does: `smcfcs.core()`
issues the two existing warnings at the end, unchanged, and mice writes them
to `loggedEvents`.

### 3.4 What `smcfcs.core()` becomes

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
the bookkeeping draws any. A regression test (section 5) checks this exactly.

### 3.5 Messages

The descriptive `print()` calls (lines 493–496, 539, 624 and 682–684) become
`message()` calls controlled by `verbose`. `smcfcs()` keeps
`verbose = TRUE` by default, so existing users still see the messages, but on
stderr instead of stdout. If that change of stream matters, we can keep
`print()` for `smcfcs()` and use `message()` only in `smcfcs.cycle()` (open
question 2).

## 4. Proposed changes to mice

### 4.1 A block method, with a name that doesn't clash

Working name: `mice.impute.smc()`. The name is open; it must not be `smcfcs`
(section 2c). smcfcs goes in Suggests.

```r
blocks <- make.blocks(list(smc = c("y", "x1", "x2"), "aux1", "aux2"))
meth   <- make.method(data, blocks); meth["smc"] <- "smc"
imp <- mice(data, blocks = blocks, method = meth,
            blots = list(smc = list(smtype = "lm",
                                    smformula = "y ~ x1 + I(x1^2) + x2",
                                    smmethod = c(x1 = "norm", x2 = "logreg"))))
```

- The block holds the outcome and the incomplete covariates in the
  substantive model. Transformations and interactions can be written inside
  `smformula` (e.g. `I(x1^2)`, `x1:x2`), as in the smcfcs vignette. smcfcs
  re-evaluates `model.matrix(smformula, data)` for each candidate value.
  Alternatively, a passive column can be added to the block, with a smcfcs
  passive expression (no `~`) such as `smmethod = c(x1sq = "x1^2")`. Then it
  is updated within the cycle. A mice passive method (`"~I(x1^2)"`) outside
  the block would be stale inside the cycle. Variables outside the block are
  imputed by ordinary mice methods. Inside a cycle they are treated as fixed
  at their current values.
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

### 4.2 Sampler: give block methods their current imputations

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

### 4.3 Differences from calling `smcfcs()` directly

- mice initialises the block with its own starting values (random draws from
  the observed values). smcfcs initialises a missing outcome with an improper
  fit of the substantive model. So the two give different results under the
  same seed. Seed equality is checked only for `smcfcs()` itself.
- `smCoefIter` has no counterpart in mice. Convergence is checked with
  mice's chain means and variances, as for any other method.
- When the outcome is missing, the outcome update refits the substantive
  model on rows with an observed outcome only. The covariate updates fit on
  all rows. This is deliberate, and the mice documentation will say so.

## 5. Tests

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
3. The existing test suite passes unchanged.

In mice:

- The block method recovers the true coefficients in the simulation of
  section 2, under the first-version `smtype`s.
- Under the same seed, `mice()` with an SMC block plus `pmm` for auxiliary
  variables runs end-to-end through `with()` and `pool()`.
- An unsupported `smtype` gives a clear error.
- The `miceadds` example still runs.

## 6. Open questions for Jonathan

1. **Outcome imputation variance (lm).** At line 1316 a missing outcome is
   drawn with `sd = sigmasq^0.5`, the residual variance of the fitted model,
   rather than the posterior draw `outcomeModResVar^0.5` computed on line
   1312. Is that intended? If not, fixing it changes results, so it should
   land before the regression references are generated, or in its own
   commit.
2. **`print()` → `message()`** for `smcfcs()` users: is the change of output
   stream acceptable (section 3.5)?
3. **Caching the setup in `cycleState`**: is that acceptable, or would you
   rather `smcfcs.cycle()` accepted a precomputed `setup` argument?
4. **Transformations: `I(x^2)` in `smformula` or a passive column?** The
   vignette uses `I(x^2)` inside `smformula`. `smcovnames` matches term
   labels against column names, so a term like `I(z^2)` in a fully observed
   `z` does not make `z` a fully observed SM covariate unless `z` also
   appears on its own. And a transformed term never enters the covariate
   models of other partially observed variables. Do you recommend
   `I()` terms, or passive columns, for the mice interface?
5. **Naming and argument order** of `smcfcs.cycle()`: any preferences?

## 7. Timeline

- Draft branch for review (this fork, branch `smcfcs-cycle`): by 30 November
  2026.
- Pull request to `jwb133/smcfcs`, with tests: by 7 December 2026.
- mice changes, on a mice branch: in parallel. They depend on the CRAN
  release of the new smcfcs.
