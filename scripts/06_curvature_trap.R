# 06_curvature_trap.R
# Program simulation, step 2: the curvature trap (the disease, before the cure).
# A naive breeder doing standard truncation selection on local additive merit
# climbs a tempting curved near-peak and STALLS there -- while a better far peak
# stays out of reach. Three signatures appear on the curved surface and on NONE
# of them on a genuine linear surface (the non-circularity guarantee):
#   (1) performance plateaus below the far peak
#   (2) A_g crashes (the stall is curvature-driven, not variance exhaustion)
#   (3) the breeder's equation de-calibrates (predicted gain >> realized gain)
# Run from the repository root:  Rscript scripts/06_curvature_trap.R
source("R/channels.R")   # trace_sq, rmvn
source("R/engine.R")
set.seed(20260606)

Ga   <- matrix(c(1.0, 0.3, 0.3, 1.0), 2)
trap <- make_surface(list(c(6, 0), c(15, 15)), c(1.0, 2.5), c(3.5, 3.5))  # near (trap) + far (better)
blin <- trap$grad(c(0, 0)); blin <- blin / sqrt(sum(blin^2)) * 0.30        # matched initial push
flat <- make_linear_surface(blin)                                         # genuine linear contrast

cur <- recurrent_selection_surface(trap, Ga, N = 5000, prop = 0.50, ncyc = 40)
lin <- recurrent_selection_surface(flat, Ga, N = 5000, prop = 0.50, ncyc = 40)
far <- trap$value(c(15, 15))

cat(sprintf("near-peak logW = %.2f   far-peak logW = %.2f\n\n", trap$value(c(6, 0)), far))
cat("CURVED (naive steepest ascent):\n")
cat(sprintf("  performance logW : %.3f -> %.3f   (plateau below far peak %.2f)\n",
            cur$logW[1], tail(cur$logW, 1), far))
cat(sprintf("  A_g              : %.3f -> %.3f   (crashes)\n", cur$Ag[1], tail(cur$Ag, 1)))
cat(sprintf("  cumulative dlogW : predicted %.2f vs realized %.2f   (DE-CALIBRATION)\n\n",
            tail(cur$cumPred, 1), tail(cur$cumReal, 1)))
cat("FLAT-LINEAR (genuine linear surface, no curvature):\n")
cat(sprintf("  performance logW : %.3f -> %.3f   (steady rise, no plateau)\n",
            lin$logW[1], tail(lin$logW, 1)))
cat(sprintf("  A_g              : %.3f -> %.3f   (stays 1)\n", lin$Ag[1], tail(lin$Ag, 1)))
cat(sprintf("  cumulative dlogW : predicted %.2f vs realized %.2f   (MATCH)\n",
            tail(lin$cumPred, 1), tail(lin$cumReal, 1)))

## ---- figure ----
dir.create("figures", showWarnings = FALSE)
png("figures/curvature_trap.png", width = 1600, height = 1300, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

gx <- seq(-2, 18, length.out = 90); gy <- seq(-4, 18, length.out = 90)
M  <- outer(gx, gy, Vectorize(function(a, b) trap$value(c(a, b))))
contour(gx, gy, M, nlevels = 14, col = "grey70", xlab = "breeding value 1",
        ylab = "breeding value 2", main = "(A) naive ascent stalls on the near peak")
lines(cur$z1, cur$z2, lwd = 2, col = "#c0392b")
points(cur$z1[1], cur$z2[1], pch = 19); points(6, 0, pch = 8, col = "#c0392b", cex = 1.6)
points(15, 15, pch = 8, col = "#27ae60", cex = 1.6)
legend("bottomright", bty = "n", pch = c(8, 8), col = c("#c0392b", "#27ae60"),
       legend = c("near peak (trap)", "far peak (better)"))

plot(cur$cycle, cur$logW, type = "l", lwd = 2, col = "#c0392b", ylim = c(0, far * 1.05),
     xlab = "cycle", ylab = "mean log-performance", main = "(B) performance: plateau vs steady rise")
lines(lin$cycle, lin$logW, lwd = 2, col = "#2980b9")
abline(h = far, lty = 3, col = "grey50")
legend("right", bty = "n", lwd = 2, col = c("#c0392b", "#2980b9"),
       legend = c("curved (stalls)", "linear (climbs)"))

plot(cur$cycle, cur$Ag, type = "l", lwd = 2, col = "#c0392b", ylim = c(0, 1),
     xlab = "cycle", ylab = expression(A[g]), main = "(C) additivity diagnoses the stall")
lines(lin$cycle, lin$Ag, lwd = 2, col = "#2980b9")
legend("right", bty = "n", lwd = 2, col = c("#c0392b", "#2980b9"),
       legend = c("curved (crashes)", "linear (== 1)"))

plot(cur$cycle, cur$cumPred, type = "l", lwd = 2, lty = 2, col = "#c0392b",
     ylim = range(0, cur$cumPred, cur$cumReal),
     xlab = "cycle", ylab = "cumulative dlogW",
     main = "(D) breeder's equation de-calibrates (curved)")
lines(cur$cycle, cur$cumReal, lwd = 2, col = "#c0392b")
legend("topleft", bty = "n", lwd = 2, lty = c(2, 1), col = c("#c0392b", "#c0392b"),
       legend = c("predicted (tangent)", "realized"))
invisible(dev.off())
cat("\nWrote figures/curvature_trap.png\n")
