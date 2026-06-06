# 04_temporal_warning.R
# Simulation question 4: does a negative r_env anticipate transfer loss BEFORE it
# appears in current-environment accuracy?
#
#  (A) the r_env formula equals the rate of change of logit-additivity (validation)
#  (B) scenarios with IDENTICAL present reliability Ag(u0) but different r_env have
#      divergent futures; r_env evaluated at the present orders the future outcomes
#  (C) early warning: r_env fires while present accuracy still looks healthy
#
# Transfer accuracy to environment u is sqrt(Ag(u)) (established in scripts/03).
# Run from the repository root:  Rscript scripts/04_temporal_warning.R
source("R/channels.R")
set.seed(20260606)

G  <- matrix(c(1.0, 0.3, 0.3, 1.0), 2)
b0 <- c(1.00, 0.60); H0 <- 0.70 * diag(2)

acc <- function(u, bu, Hu) sqrt(local_channels(u, G, b0, bu, H0, Hu)$Ag)
lam <- function(u, bu, Hu) { lc <- local_channels(u, G, b0, bu, H0, Hu); log(lc$Vlin) - log(lc$Vquad) }

## three scenarios with IDENTICAL Ag(0) (b0,H0 fixed); different sensitivities
scen <- list(
  eroding   = list(bu = -0.50 * b0,            Hu =  0.50 * diag(2)),  # optimum approach + curvature up
  stable    = list(bu =  0.30 * c(-0.9, 1.18), Hu =  0.00 * diag(2)),  # gradient turns, no curvature
  deepening = list(bu =  0.50 * b0,            Hu = -0.30 * diag(2))   # gradient grows + curvature down
)

## ---- (A) formula vs finite-difference of logit-additivity ----
cat("(A) r_env formula vs finite difference of logit-additivity (udot = 1)\n")
eps <- 1e-4
for (nm in names(scen)) {
  s <- scen[[nm]]
  fd <- (lam(eps, s$bu, s$Hu) - lam(-eps, s$bu, s$Hu)) / (2 * eps)
  cat(sprintf("  %-10s formula = %+.4f   finite-diff = %+.4f\n",
              nm, r_env(1, G, b0, s$bu, H0, s$Hu), fd))
}

## ---- (B) identical present, divergent futures ----
cat("\n(B) identical present reliability, divergent futures (drift to u = 1)\n")
cat(sprintf("  %-10s %8s %7s %8s %7s %8s\n", "scenario", "r_env", "Ag(0)", "acc(0)", "Ag(1)", "acc(1)"))
for (nm in names(scen)) {
  s <- scen[[nm]]
  cat(sprintf("  %-10s %+8.3f %7.3f %8.3f %7.3f %8.3f\n", nm,
              r_env(1, G, b0, s$bu, H0, s$Hu),
              local_channels(0, G, b0, s$bu, H0, s$Hu)$Ag, acc(0, s$bu, s$Hu),
              local_channels(1, G, b0, s$bu, H0, s$Hu)$Ag, acc(1, s$bu, s$Hu)))
}

## ---- (C) lead time: warning fires while present accuracy still healthy ----
s <- scen$eroding
drift <- seq(0, 1.5, length.out = 301)
a <- sapply(drift, function(u) acc(u, s$bu, s$Hu))
ustar <- drift[which.max(a < 0.70)]
cat(sprintf("\n(C) eroding: acc(0) = %.3f (healthy); r_env = %+.3f fires now;\n    transfer accuracy crosses below 0.70 only at drift u* = %.2f\n",
            a[1], r_env(1, G, b0, s$bu, H0, s$Hu), ustar))
sim <- transfer_accuracy(80000, G, b0, s$bu, H0, s$Hu, 1.0)
cat(sprintf("    check: simulated Pearson at u=1 = %.3f vs sqrt(Ag(1)) = %.3f\n",
            sim[["pearson"]], acc(1, s$bu, s$Hu)))

## ---- figure ----
dir.create("figures", showWarnings = FALSE)
png("figures/temporal_warning.png", width = 1600, height = 720, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

cols <- c(eroding = "#c0392b", stable = "#7f8c8d", deepening = "#2980b9")
plot(NA, xlim = c(0, 1.5), ylim = c(0.4, 1.0), xlab = "target drift  (u - u0)",
     ylab = "transfer accuracy  sqrt(Ag(u))", main = "(A) same present, different futures")
for (nm in names(scen)) {
  s <- scen[[nm]]
  lines(drift, sapply(drift, function(u) acc(u, s$bu, s$Hu)), col = cols[nm], lwd = 2)
}
points(0, acc(0, scen$eroding$bu, scen$eroding$Hu), pch = 19)
abline(h = 0.70, lty = 3, col = "grey50")
legend("bottomleft", bty = "n", lwd = 2, col = cols,
       legend = sprintf("%s (r_env=%+.2f)", names(scen),
                        sapply(scen, function(s) r_env(1, G, b0, s$bu, H0, s$Hu))))

ss <- seq(-1, 1, length.out = 41); re <- ch <- numeric(length(ss))
for (i in seq_along(ss)) {
  bu <- ss[i] * b0; Hu <- -ss[i] * 0.4 * diag(2)
  re[i] <- r_env(1, G, b0, bu, H0, Hu)
  ch[i] <- acc(1, bu, Hu) - acc(0, bu, Hu)
}
plot(re, ch, pch = 19, col = "#8e44ad", xlab = "r_env at present",
     ylab = "future accuracy change  acc(1) - acc(0)", main = "(B) r_env forecasts the change")
abline(h = 0, v = 0, lty = 3, col = "grey50")

invisible(dev.off())
cat("\nWrote figures/temporal_warning.png\n")
