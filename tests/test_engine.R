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
