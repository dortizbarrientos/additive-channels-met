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

# ---------------------------------------------------------------------------
# Curved-surface extension (program simulation, step 2: the curvature trap).
# The engine genetics above are UNTOUCHED. The only change is that the
# selection direction each cycle is the LOCAL gradient of a (possibly curved)
# log-performance surface, instead of a fixed direction. A naive breeder doing
# standard truncation selection on local additive merit thus performs steepest
# ascent on the surface -- and, on a curved surface, can climb a tempting near
# peak and stall there while a better far peak stays out of reach.
# Requires trace_sq() from R/channels.R.

## sum-of-Gaussians log-performance surface; returns value/grad/hess closures.
make_surface <- function(centers, heights, widths) {
  list(
    value = function(z)
      sum(mapply(function(c, a, s) a * exp(-sum((z - c)^2) / (2 * s^2)),
                 centers, heights, widths)),
    grad = function(z) {
      g <- rep(0, length(z))
      for (k in seq_along(centers)) {
        d <- z - centers[[k]]; s2 <- widths[k]^2
        gk <- heights[k] * exp(-sum(d^2) / (2 * s2)); g <- g - gk * d / s2
      }
      g
    },
    hess = function(z) {
      p <- length(z); H <- matrix(0, p, p)
      for (k in seq_along(centers)) {
        d <- z - centers[[k]]; s2 <- widths[k]^2
        gk <- heights[k] * exp(-sum(d^2) / (2 * s2))
        H <- H + gk * ((d %o% d) / s2^2 - diag(p) / s2)
      }
      H
    }
  )
}

## genuine linear surface L(z) = b' z (zero curvature) -- the null contrast.
make_linear_surface <- function(b)
  list(value = function(z) sum(b * z),
       grad  = function(z) b,
       hess  = function(z) matrix(0, length(z), length(z)))

## recurrent selection that follows the LOCAL gradient of `surface` each cycle.
## Tracks performance, A_g at the current mean, and the breeder's-equation
## calibration (cumulative first-order predicted vs realized change in logW).
recurrent_selection_surface <- function(surface, Ga, N = 5000, prop = 0.50, ncyc = 40) {
  Z <- rmvn(N, Ga); cumP <- 0; cumR <- 0; rows <- vector("list", ncyc)
  for (t in seq_len(ncyc)) {
    zbar <- colMeans(Z); G <- cov(Z)
    b <- surface$grad(zbar); H <- surface$hess(zbar)
    Vlin <- as.numeric(t(b) %*% G %*% b); Vq <- 0.5 * trace_sq(H %*% G)
    Ag   <- Vlin / (Vlin + Vq)
    cy   <- sel_cycle(Z, b, Ga, prop)        # standard truncation cycle, local direction
    zbn  <- colMeans(cy$Zo)
    cumP <- cumP + sum(b * (zbn - zbar))                       # first-order (tangent) prediction
    cumR <- cumR + (surface$value(zbn) - surface$value(zbar)) # realized change in logW
    rows[[t]] <- data.frame(cycle = t, z1 = zbar[1], z2 = zbar[2],
                            logW = surface$value(zbar), Ag = Ag,
                            cumPred = cumP, cumReal = cumR,
                            step = sqrt(sum((zbn - zbar)^2)))
    Z <- cy$Zo
  }
  do.call(rbind, rows)
}
