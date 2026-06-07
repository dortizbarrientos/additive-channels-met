# additive-channels-met

Simulation code for the additive-channel diagnostics in
*Additive channels diagnose when crop genomic prediction transfers across
environments* (Ortiz-Barrientos, Messina & Cooper, in prep.).

The code establishes, by simulation against a known truth, that the paper's
estimators recover the channel diagnostics, and it answers the practical
question a multi-environment trial forces: **how far does breeding-value
estimation error bias the diagnostics, and can a reliability-based correction
remove that bias?**

## The idea in brief

A genotype's breeding value is a vector `z ~ N(0, G)` in trait space. When the
fitness gradient and Hessian change linearly along an environmental coordinate
`u` (`b(u)=b0+du*bu`, `H(u)=H0+du*Hu`), the **genetic reaction norm is exactly
linear in the environment**:

```
g_i(u) = alpha_i + du * theta_i
alpha_i = b0' z_i + 0.5 z_i' H0 z_i       # intercept
theta_i = bu' z_i + 0.5 z_i' Hu z_i       # slope
```

Curvature of the response surface is **not** a quadratic response to `u`; it is
hidden inside the intercept and slope as their quadratic-in-`z` parts.

- **Totals** `c00 = Var(alpha)`, `c11 = Var(theta)`, `c10 = Cov(alpha, theta)`
  need only the target trait.
- **Splits** -- the reference additivity index `Ag` and the GEI additivity
  index `AgGEI` -- are the genetic R-squared of intercept and slope regressed
  *linearly* on the multivariate breeding values, and need a **trait panel**.
  No Hessian is ever estimated.

## Errors in variables (the Berkson correction)

Genomic breeding values are BLUPs (posterior means), so the true value equals
the estimate plus an **orthogonal** error: `z = z_hat + delta`, `delta` independent
of `z_hat` (Berkson), *not* `z_hat = z + e` (classical). The two imply
**opposite-sign** corrections, so the distinction is load-bearing.

Under Berkson error with breeding-value reliability `rho`:

- the naive regression **coefficient** is unbiased, but
- the naive linear-explained variance is deflated by `rho`, so the **naive
  `AgGEI` comes in at about `rho` times the truth**;
- the fix evaluates the (unbiased) coefficient under the reconstructed
  `G = Var(z_hat) + PEV` (equivalently, divide by reliability).

If instead one used de-regressed / unshrunken estimates, the classical model
and its opposite correction would apply.

## Reference output

Truth for the bundled 2-trait example: `Ag = 0.7631`, `AgGEI = 0.6559`.
`scripts/02_eiv_correction.R` reproduces (Monte Carlo means; cross-checked
against `verify/verify_numpy.py`):

```
  rho | AgGEI naive  AgGEI corr |  Ag naive   Ag corr
  0.3 |     ~0.197      ~0.656   |   ~0.229     ~0.762
  0.5 |     ~0.328      ~0.656   |   ~0.382     ~0.763
  0.7 |     ~0.460      ~0.657   |   ~0.535     ~0.764
  0.9 |     ~0.590      ~0.656   |   ~0.687     ~0.763
```

Naive tracks `rho * truth`; corrected lands on the truth at every reliability.

## Does the diagnostic predict transfer? (scripts/03)

For an additive reaction-norm predictor (one that tracks the slope but, being
linear in breeding values, cannot see curvature), three results hold:

1. **Accuracy is the square root of additivity.** Transfer Pearson accuracy to a
   held-out environment satisfies `Pearson(u_t)^2 = Ag(u_t)` exactly -- the
   additivity index is not merely correlated with transfer, it equals it.
2. **Transfer to a shifted environment rises with `AgGEI`.** At a far target
   (`u_t = 1`), accuracy climbs from ~0.74 to ~0.90 (Pearson) and ~0.77 to ~0.95
   (Spearman) as `AgGEI` goes 0.1 -> 0.9.
3. **The Pearson/Spearman gap is itself a failure-mode flag.** When additivity
   is high, ranks are robust and `Spearman >= Pearson` (curvature acts as skewed
   noise). When curvature dominates (`Ag(u_t)` small), ranks reverse and
   `Spearman < Pearson`. The sign of the gap flips near `Ag(u_t) ~ 0.4`.

So curvature always erodes transfer, but it induces genuine *rank reversals*
only once it dominates the local linear signal -- a sharper claim than
"curvature breaks ranking."

## Early warning: does r_env fire before failure? (scripts/04)

The environmental rate signal `r_env` is the rate of change of logit-additivity
as the target environment drifts; its sign is the sign of `d Ag/du`, and it is
computed from quantities the current MET already estimates. Three results:

1. **It is the additivity rate.** `r_env` matches the finite-difference rate of
   `log(Vlin/Vquad)` to four decimals -- a validation of the manuscript formula.
2. **Same present, different futures.** Scenarios with identical present
   reliability `Ag(u0)` (so current cross-validation cannot tell them apart) but
   different `r_env` diverge exactly as the sign predicts: drifting to `u=1`,
   accuracy goes 0.46 (eroding), 0.88 (stable), 0.98 (deepening) from a common
   0.87.
3. **It is an early warning.** In the eroding case (gradient collapsing toward an
   optimum, curvature growing) `r_env < 0` fires while accuracy is still 0.87;
   accuracy only crosses 0.70 after the target has drifted ~0.5 units. That gap
   is the lead time the diagnostic buys.

## Program simulation: the null engine (R/engine.R, scripts/05)

Before layering any additive-channel logic onto selection decisions, the
recurrent-selection engine must reproduce textbook results on a flat-linear
surface (constant selection direction, no curvature). It clears four gates:

1. **Breeder's equation** -- realized response equals `i * G b / sqrt(b' G b)`.
2. **Bulmer effect** -- the index variance `b' G b` declines under truncation
   selection to the closed-form equilibrium `b' Ga b / (1 + k)`,
   `k = i(i - x)`, regenerated each cycle by Mendelian sampling (not to zero).
3. **Truncation intensity** -- realized intensity equals `dnorm(x)/prop`.
4. **No spurious plateau** -- once variance equilibrates, per-cycle gain is
   steady at `i*sqrt(v_eq)`; the engine does not manufacture a stall absent
   curvature.

Only after this null passes do the curvature surface and the channel-guided
selection/crossing protocol (P1 vs P2) get built on top -- so that any later
plateau or transfer failure is attributable to curvature, not to the engine.

## Program simulation: the curvature trap (R/engine.R, scripts/06)

With the null cleared, the engine's selection direction is set each cycle to the
local gradient of a log-performance surface (the genetics are unchanged -- this
is standard truncation selection following local additive merit, i.e. steepest
ascent). On a surface with a tempting curved near peak and a better far peak, a
naive breeder climbs the near peak and stalls. Three signatures appear -- and
none of them appear on a genuine linear surface, which is the non-circularity
guarantee that the pathology is curvature, not the engine:

1. **Performance plateaus** below the far peak (the breeder is trapped at a
   local optimum it cannot leave by local ascent).
2. **A_g crashes** as the mean enters the curved peak region. This is the
   diagnostic payload: a response curve alone cannot tell a curvature trap from
   benign variance exhaustion, but A_g can -- it falls only when curvature, not
   data scarcity, is the cause.
3. **The breeder's equation de-calibrates**: cumulative first-order (tangent)
   predicted gain increasingly exceeds realized gain. Because truncation
   selection is direction-based, the breeder keeps spending its full selection
   differential even as the gradient flattens, so realized performance gain
   falls short of what the linear model promised.

This is the disease. The cure -- a channel-guided protocol that uses the
diagnostic to redirect crossing toward the better peak, and the P1-vs-P2
contrast across a selection x crossing grid -- is the next build.

## Layout

```
R/channels.R               core functions (truth, simulators, estimators)
scripts/01_lane_a_recovery.R   recover Ag, AgGEI, c-coefficients; rank-2 check
scripts/02_eiv_correction.R    attenuation vs reliability + correction + figure
scripts/03_transfer_accuracy.R does AgGEI predict transfer? (Pearson vs Spearman)
scripts/04_temporal_warning.R  does r_env anticipate transfer loss before it bites?
R/engine.R                     infinitesimal-model recurrent-selection engine
scripts/05_null_recurrent_selection.R  NULL gate: recover textbook recurrent selection
scripts/06_curvature_trap.R    the curvature trap: naive breeder stalls
tests/test_engine.R            null gates + curvature-trap assertions
tests/test_channels.R          algebraic + recovery checks (assertions)
verify/verify_numpy.py         independent numpy cross-check (recovery + EIV)
verify/verify_q3_numpy.py      independent numpy cross-check (transfer accuracy)
verify/verify_q4_numpy.py      independent numpy cross-check (temporal warning)
verify/verify_engine_numpy.py  independent numpy cross-check (null engine)
verify/verify_trap_numpy.py    independent numpy cross-check (curvature trap)
figures/                       output (created on run)
```

## Run

```bash
Rscript scripts/01_lane_a_recovery.R
Rscript scripts/02_eiv_correction.R     # writes figures/eiv_attenuation.png
Rscript scripts/03_transfer_accuracy.R  # writes figures/transfer_accuracy.png
Rscript scripts/04_temporal_warning.R   # writes figures/temporal_warning.png
Rscript scripts/05_null_recurrent_selection.R  # NULL gate; writes figures/null_recurrent_selection.png
Rscript scripts/06_curvature_trap.R     # the curvature trap; writes figures/curvature_trap.png
Rscript tests/test_engine.R             # null gates + curvature-trap assertions
Rscript tests/test_channels.R
python3 verify/verify_numpy.py          # optional independent cross-check
```

Base R only (no packages). Seed `20260606` for reproducibility.

## Scope and caveats

The EIV experiment isolates **breeding-value-panel** error (the response,
`theta`, is treated as observed). The full two-sided case additionally corrects
the denominator with the slope BLUP's prediction-error variance,
`c11 = Var(theta_hat) + PEV_theta`; this is the natural next extension. The
framework is a one-gradient, local-Gaussian approximation; the rank-2 adequacy
test guards against a mis-specified or blended environmental coordinate.

## License

MIT (see `LICENSE`).
