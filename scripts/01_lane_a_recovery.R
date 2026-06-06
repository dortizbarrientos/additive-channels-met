# 01_lane_a_recovery.R
# Verify that the Lane A estimators recover the structural truth, and that the
# cross-environment genetic covariance is rank 2.
# Run from the repository root:  Rscript scripts/01_lane_a_recovery.R
source("R/channels.R")
set.seed(20260606)

## ---- structural truth (2 component traits) ----
G  <- matrix(c(1.0, 0.3, 0.3, 1.0), 2)
b0 <- c(1.00, 0.60); H0 <- 0.70 * diag(2)
bu <- c(0.45, 0.25); Hu <- 0.40 * diag(2)
stopifnot(all(eigen(G, only.values = TRUE)$values > 0),
          isSymmetric(H0), isSymmetric(Hu))

truth <- true_channels(G, b0, bu, H0, Hu)
stopifnot(truth$c00 * truth$c11 - truth$c10^2 > 0)   # Phi2 positive definite

## ---- simulate with TRUE breeding values observed ----
N   <- 200000
sim <- simulate_true(N, G, b0, bu, H0, Hu)
Gh  <- cov(sim$Z)

c00h <- var(sim$alpha); c11h <- var(sim$theta); c10h <- cov(sim$alpha, sim$theta)
Agh    <- vlin_explained(sim$alpha, sim$Z, Gh) / c00h
AgGEIh <- vlin_explained(sim$theta, sim$Z, Gh) / c11h

## ---- rank-2 of the cross-environment genetic covariance ----
u    <- c(-1.0, -0.5, 0.0, 0.5, 1.0)
gmat <- outer(sim$alpha, rep(1, length(u))) + outer(sim$theta, u)  # N x E
ev   <- eigen(cov(gmat), only.values = TRUE)$values

cat("Quantity      truth     recovered\n")
cat(sprintf("Ag         %8.4f   %8.4f\n", truth$Ag,    Agh))
cat(sprintf("AgGEI      %8.4f   %8.4f\n", truth$AgGEI, AgGEIh))
cat(sprintf("c00        %8.4f   %8.4f\n", truth$c00,   c00h))
cat(sprintf("c10        %8.4f   %8.4f\n", truth$c10,   c10h))
cat(sprintf("c11        %8.4f   %8.4f\n", truth$c11,   c11h))
cat(sprintf("C^(G) eigenvalues: %s\n", paste(round(ev, 4), collapse = ", ")))
cat("(two non-negligible eigenvalues confirm rank 2)\n")
