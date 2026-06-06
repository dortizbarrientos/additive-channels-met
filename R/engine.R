# engine.R -- infinitesimal-model recurrent-selection engine.
# Companion code to Ortiz-Barrientos, Messina & Cooper (in prep.).
# Requires rmvn() from R/channels.R (source channels.R first).
#
# NULL-FIRST DISCIPLINE. On a flat-linear surface (constant selection direction,
# no curvature) this engine must reproduce textbook recurrent selection before
# any additive-channel logic is layered on. The four checks it must pass:
#   (1) breeder's equation: Delta-zbar = i * G b / sqrt(b' G b)
#   (2) Bulmer effect: index variance b' G b -> b' Ga b / (1 + k)
#   (3) truncation intensity: realized i = dnorm(x)/prop, x = qnorm(1-prop)
#   (4) no spurious plateau: steady per-cycle gain i*sqrt(v_eq) once equilibrated
#
# Model: each individual carries a multivariate breeding value z. Selection is
# truncation on the index I = z' bsel (true breeding values; accuracy 1 in the
# null). Offspring = midparent + Mendelian sampling ~ MVN(0, 0.5 Ga). The
# constant genic G-matrix Ga encodes the infinitesimal assumption; the Bulmer
# (linkage-disequilibrium) reduction emerges from selection correlating parental
# breeding values, regenerated each cycle by Mendelian sampling.

## truncation-selection constants for saved fraction `prop`
trunc_intensity <- function(prop) {
  x <- qnorm(1 - prop)
  i <- dnorm(x) / prop
  list(x = x, i = i, k = i * (i - x))   # k = variance-reduction coefficient
}

## one cycle: truncation-select top `prop` on I = Z bsel, random-mate the
## selected parents, add Mendelian sampling. Returns offspring + diagnostics.
sel_cycle <- function(Z, bsel, Ga, prop) {
  N   <- nrow(Z)
  I   <- as.vector(Z %*% bsel)
  thr <- quantile(I, 1 - prop)
  sel <- Z[I >= thr, , drop = FALSE]
  m   <- nrow(sel)
  pa  <- sel[sample.int(m, N, replace = TRUE), , drop = FALSE]
  pb  <- sel[sample.int(m, N, replace = TRUE), , drop = FALSE]
  Zo  <- 0.5 * (pa + pb) + rmvn(N, 0.5 * Ga)
  list(Zo = Zo, diff = colMeans(sel) - colMeans(Z),
       i_real = (mean(sel %*% bsel) - mean(I)) / sd(I))
}

## recurrent selection on a flat-linear surface (constant selection direction).
## Returns a per-cycle trajectory data frame.
recurrent_selection <- function(Ga, bsel, N = 8000, prop = 0.20, ncyc = 30) {
  ti <- trunc_intensity(prop)
  Z  <- rmvn(N, Ga)
  rows <- vector("list", ncyc)
  for (t in seq_len(ncyc)) {
    G    <- cov(Z)
    meanI <- mean(Z %*% bsel)
    varI  <- as.numeric(t(bsel) %*% G %*% bsel)
    pred  <- ti$i * as.vector(G %*% bsel) / sqrt(varI)   # breeder's-equation prediction
    cy    <- sel_cycle(Z, bsel, Ga, prop)
    rows[[t]] <- data.frame(
      cycle  = t,
      meanI  = meanI,
      varI   = varI,
      i_real = cy$i_real,
      bs_err = sqrt(sum((cy$diff - pred)^2)) / sqrt(sum(pred^2))  # breeder's-eq rel. error
    )
    Z <- cy$Zo
  }
  out <- do.call(rbind, rows)
  out$gain <- c(NA, diff(out$meanI))
  out
}
