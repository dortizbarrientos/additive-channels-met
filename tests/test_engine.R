# test_engine.R -- the four null gates as assertions (base R; no testthat).
# Run from the repository root:  Rscript tests/test_engine.R
source("R/channels.R")
source("R/engine.R")
set.seed(20260606)

ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("pass: ", msg, "\n")
}

Ga   <- matrix(c(1.0, 0.3, 0.3, 1.0), 2)
bsel <- c(1.0, 0.6); prop <- 0.20
ti   <- trunc_intensity(prop)
v_eq <- as.numeric(t(bsel) %*% Ga %*% bsel) / (1 + ti$k)
traj <- recurrent_selection(Ga, bsel, N = 6000, prop = prop, ncyc = 25)
late <- traj[traj$cycle > 17, ]

ok(mean(late$bs_err) < 0.05,
   "GATE 1 breeder's equation recovered (rel. error < 0.05)")
ok(abs(mean(late$varI) - v_eq) / v_eq < 0.05,
   "GATE 2 Bulmer equilibrium index variance == v_a/(1+k)")
ok(abs(mean(late$i_real) - ti$i) < 0.05,
   "GATE 3 realized truncation intensity == dnorm(x)/prop")
ok(abs(mean(late$gain) - ti$i * sqrt(v_eq)) < 0.05 && sd(late$gain) < 0.10,
   "GATE 4 steady per-cycle gain i*sqrt(v_eq), no plateau")

cat("\nAll null-engine checks passed.\n")

## ---- curvature-trap checks (step 2) ------------------------------------
## On a curved surface a naive breeder must misbehave; on a genuine linear
## surface it must not. These guard against reading curvature pathology into
## an engine artifact.
trap <- make_surface(list(c(6, 0), c(15, 15)), c(1.0, 2.5), c(3.5, 3.5))
blin <- trap$grad(c(0, 0)); blin <- blin / sqrt(sum(blin^2)) * 0.30
flat <- make_linear_surface(blin)
cur  <- recurrent_selection_surface(trap, Ga, N = 5000, prop = 0.50, ncyc = 40)
lin  <- recurrent_selection_surface(flat, Ga, N = 5000, prop = 0.50, ncyc = 40)
far  <- trap$value(c(15, 15))

ok(tail(cur$logW, 1) < 0.6 * far,
   "TRAP curved: performance plateaus well below the far peak")
ok(cur$Ag[1] > 0.8 && tail(cur$Ag, 1) < 0.5,
   "TRAP curved: A_g crashes from high to low (curvature diagnosis)")
ok(tail(cur$cumPred, 1) > 1.3 * tail(cur$cumReal, 1),
   "TRAP curved: breeder's equation de-calibrates (predicted >> realized)")
ok(all(abs(lin$Ag - 1) < 1e-6),
   "TRAP linear: A_g == 1 everywhere (no curvature)")
ok(abs(tail(lin$cumPred, 1) - tail(lin$cumReal, 1)) < 0.05 * tail(lin$cumReal, 1),
   "TRAP linear: predicted == realized (breeder's equation calibrated)")
ok(tail(lin$logW, 1) > 3 * tail(cur$logW, 1),
   "TRAP linear: steady rise, no plateau")

cat("\nAll null-engine + curvature-trap checks passed.\n")
