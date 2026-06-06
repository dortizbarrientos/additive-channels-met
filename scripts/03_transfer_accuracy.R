# 03_transfer_accuracy.R
# Simulation question 3: does AgGEI predict cross-environment prediction transfer,
# and do Pearson and Spearman accuracy tell the same story?
#
# Three results:
#  (A) Additive transfer accuracy is exactly the square root of local additivity:
#        Pearson(u_t)^2 = Ag(u_t).
#  (B) Transfer to a shifted (far) environment rises with AgGEI for both metrics.
#  (C) Failure-mode crossover: Spearman - Pearson flips sign with local additivity.
#      High Ag(u_t): ranks robust (Spearman >= Pearson); curvature is skewed noise.
#      Low  Ag(u_t): rank reversals (Spearman < Pearson); curvature scrambles order.
#
# Run from the repository root:  Rscript scripts/03_transfer_accuracy.R
source("R/channels.R")
set.seed(20260606)

G  <- matrix(c(1.0, 0.3, 0.3, 1.0), 2)
b0 <- c(1.00, 0.60); H0 <- 0.70 * diag(2)
u_all <- c(-1.0, -0.5, 0.0, 0.5, 1.0)
N <- 60000

## ---- (A) identity: Pearson^2 = Ag(u_t) ----
gei <- make_gei(0.5, 0.5, G, dir_bu = c(0.45, 0.25))
A <- data.frame()
for (ut in u_all) {
  r <- transfer_accuracy(N, G, b0, gei$bu, H0, gei$Hu, ut)
  A <- rbind(A, data.frame(u_t = ut, Pearson = r[["pearson"]],
                           Pearson2 = r[["pearson"]]^2, Ag_ut = r[["Ag_ut"]],
                           Spearman = r[["spearman"]]))
}
cat("(A) Additive transfer accuracy vs local additivity  (AgGEI = 0.5)\n")
print(round(A, 4), row.names = FALSE)

## ---- (B) far-environment transfer vs AgGEI ----
levels <- c(0.1, 0.3, 0.5, 0.7, 0.9); reps <- 25; ut_far <- 1.0
B <- data.frame()
for (a in levels) {
  g <- make_gei(a, 0.5, G, dir_bu = c(0.45, 0.25))
  P <- S <- numeric(reps)
  for (r in seq_len(reps)) {
    rr <- transfer_accuracy(40000, G, b0, g$bu, H0, g$Hu, ut_far)
    P[r] <- rr[["pearson"]]; S[r] <- rr[["spearman"]]
  }
  B <- rbind(B, data.frame(AgGEI = a, Pearson = mean(P), Spearman = mean(S)))
}
cat("\n(B) Transfer to far environment (u_t = 1) vs AgGEI\n")
print(round(B, 4), row.names = FALSE)

## ---- (C) failure-mode crossover: Spearman - Pearson vs local additivity ----
## weak gradient + increasing curvature drives Ag(u_t) low to expose reversals
b0w <- c(0.20, 0.10); buw <- c(0.10, 0.05); H0w <- matrix(0, 2, 2)
C <- data.frame()
for (hh in c(0.1, 0.3, 0.6, 0.9, 1.2)) {
  rr <- transfer_accuracy(120000, G, b0w, buw, H0w, hh * diag(2), 1.0)
  C <- rbind(C, data.frame(Ag_ut = rr[["Ag_ut"]], Pearson = rr[["pearson"]],
                           Spearman = rr[["spearman"]],
                           gap = rr[["spearman"]] - rr[["pearson"]]))
}
C <- C[order(C$Ag_ut), ]
cat("\n(C) Spearman - Pearson vs local additivity (curvature-dominant sweep)\n")
print(round(C, 4), row.names = FALSE)
cat("gap > 0: ranks robust (curvature = skewed noise); gap < 0: rank reversals.\n")

## ---- figure ----
dir.create("figures", showWarnings = FALSE)
png("figures/transfer_accuracy.png", width = 2100, height = 700, res = 150)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

plot(A$Ag_ut, A$Pearson2, pch = 19, col = "#2980b9", xlim = c(0, 1), ylim = c(0, 1),
     xlab = expression(A[g](u[t])), ylab = expression(Pearson^2),
     main = "(A) accuracy = additivity")
abline(0, 1, lty = 2, col = "grey50")

matplot(B$AgGEI, B[, c("Pearson", "Spearman")], type = "b", pch = c(19, 17),
        col = c("#c0392b", "#2980b9"), lty = 1, ylim = c(0.5, 1),
        xlab = expression(A[g]^(GEI)), ylab = "transfer accuracy (u_t = 1)",
        main = "(B) transfer rises with AgGEI")
legend("bottomright", c("Pearson", "Spearman"), pch = c(19, 17),
       col = c("#c0392b", "#2980b9"), bty = "n")

plot(C$Ag_ut, C$gap, type = "b", pch = 19, col = "#8e44ad",
     xlab = expression(A[g](u[t])), ylab = "Spearman - Pearson",
     main = "(C) failure mode")
abline(h = 0, lty = 2, col = "grey50")

invisible(dev.off())
cat("\nWrote figures/transfer_accuracy.png\n")
