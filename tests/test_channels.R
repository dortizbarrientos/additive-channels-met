# test_channels.R -- algebraic and recovery checks (base R; no testthat).
# Run from the repository root:  Rscript tests/test_channels.R
source("R/channels.R")
set.seed(1)

ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("pass: ", msg, "\n")
}

G  <- matrix(c(1.0, 0.3, 0.3, 1.0), 2)
b0 <- c(1.00, 0.60); H0 <- 0.70 * diag(2)
bu <- c(0.45, 0.25); Hu <- 0.40 * diag(2)
tr <- true_channels(G, b0, bu, H0, Hu)

ok(abs((tr$VlinGEI + tr$VquadGEI) - tr$c11) < 1e-12, "VlinGEI + VquadGEI == c11")
ok(tr$c00 * tr$c11 - tr$c10^2 > 0,                   "c00*c11 > c10^2 (Phi2 PD)")
M <- H0 %*% G
ok(abs(trace_sq(M) - sum(diag(M %*% M))) < 1e-10,    "trace_sq == tr(M %*% M)")

N <- 300000
s <- simulate_true(N, G, b0, bu, H0, Hu)
AgGEIh <- vlin_explained(s$theta, s$Z, cov(s$Z)) / var(s$theta)
ok(abs(AgGEIh - tr$AgGEI) < 0.01,                    "AgGEI recovered within 0.01 (true z)")

u <- c(-1, -0.5, 0, 0.5, 1)
gmat <- outer(s$alpha, rep(1, 5)) + outer(s$theta, u)
ev <- sort(eigen(cov(gmat), only.values = TRUE)$values, decreasing = TRUE)
ok(ev[3] / ev[1] < 1e-8,                             "third eigenvalue negligible (rank 2)")

rho <- 0.5; an <- ac <- numeric(30)
for (i in 1:30) {
  e  <- simulate_eiv(N, G, b0, bu, H0, Hu, rho)
  Vo <- cov(e$Zhat)
  an[i] <- vlin_explained(e$theta, e$Zhat, Vo)         / var(e$theta)
  ac[i] <- vlin_explained(e$theta, e$Zhat, Vo + e$PEV) / var(e$theta)
}
ok(abs(mean(an) - rho * tr$AgGEI) < 0.01,            "naive ~ rho*truth at rho=0.5")
ok(abs(mean(ac) - tr$AgGEI)       < 0.01,            "corrected ~ truth at rho=0.5")

cat("\nAll checks passed.\n")

## ---- transfer-accuracy additions (scripts/03) ----
ok(abs(local_channels(0, G, b0, bu, H0, Hu)$Ag - tr$Ag) < 1e-12,
   "local_channels at u=0 equals reference Ag")

g <- make_gei(0.6, 0.5, G, dir_bu = c(0.45, 0.25))
ok(abs(as.numeric(t(g$bu) %*% G %*% g$bu) - 0.6 * 0.5) < 1e-10,
   "make_gei: bu' G bu == AgGEI * c11")
ok(abs(0.5 * trace_sq(g$Hu %*% G) - 0.4 * 0.5) < 1e-10,
   "make_gei: 0.5 tr[(Hu G)^2] == (1-AgGEI) * c11")

rr <- transfer_accuracy(300000, G, b0, g$bu, H0, g$Hu, 1.0)
ok(abs(rr[["pearson"]]^2 - rr[["Ag_ut"]]) < 0.01,
   "transfer Pearson^2 == Ag(u_t) within 0.01")

cat("\nAll transfer checks passed.\n")

## ---- environmental rate signal (scripts/04) ----
lam_u <- function(u, bu, Hu) { lc <- local_channels(u, G, b0, bu, H0, Hu); log(lc$Vlin) - log(lc$Vquad) }
bu_e <- -0.5 * b0; Hu_e <- 0.5 * diag(2)
fd_e <- (lam_u(1e-4, bu_e, Hu_e) - lam_u(-1e-4, bu_e, Hu_e)) / (2e-4)
ok(abs(r_env(1, G, b0, bu_e, H0, Hu_e) - fd_e) < 1e-4,
   "r_env equals logit-additivity finite-difference rate")
ok(r_env(1, G, b0, -0.5 * b0, H0, 0.5 * diag(2)) < 0 &&
   r_env(1, G, b0,  0.5 * b0, H0, -0.3 * diag(2)) > 0,
   "r_env sign: eroding < 0, deepening > 0")
cat("\nAll temporal checks passed.\n")
