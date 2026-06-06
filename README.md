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

## Layout

```
R/channels.R               core functions (truth, simulators, estimators)
scripts/01_lane_a_recovery.R   recover Ag, AgGEI, c-coefficients; rank-2 check
scripts/02_eiv_correction.R    attenuation vs reliability + correction + figure
tests/test_channels.R          algebraic + recovery checks (assertions)
verify/verify_numpy.py         independent numpy implementation (cross-check)
figures/                       output (created on run)
```

## Run

```bash
Rscript scripts/01_lane_a_recovery.R
Rscript scripts/02_eiv_correction.R     # writes figures/eiv_attenuation.png
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
