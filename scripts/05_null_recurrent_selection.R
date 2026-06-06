# 05_null_recurrent_selection.R
# NULL gate: the recurrent-selection engine must reproduce textbook results on a
# flat-linear surface BEFORE any additive-channel logic is added.
#   (1) breeder's equation   (2) Bulmer equilibrium   (3) truncation intensity
#   (4) steady gain (no spurious plateau)
# Run from the repository root:  Rscript scripts/05_null_recurrent_selection.R
source("R/channels.R")   # rmvn
source("R/engine.R")
set.seed(20260606)

Ga   <- matrix(c(1.0, 0.3, 0.3, 1.0), 2)   # genic G-matrix
bsel <- c(1.0, 0.6)                        # constant selection direction (flat-linear)
prop <- 0.20

ti    <- trunc_intensity(prop)
v_a   <- as.numeric(t(bsel) %*% Ga %*% bsel)
v_eq  <- v_a / (1 + ti$k)                  # Bulmer equilibrium index variance (closed form)
traj  <- recurrent_selection(Ga, bsel, N = 8000, prop = prop, ncyc = 30)
late  <- traj[traj$cycle > 20, ]

cat(sprintf("theory: i=%.4f  k=%.4f  v_a=%.4f  v_eq=v_a/(1+k)=%.4f\n\n",
            ti$i, ti$k, v_a, v_eq))
cat(sprintf("GATE 1  breeder's equation: late-cycle mean relative error = %.4f\n", mean(late$bs_err)))
cat(sprintf("GATE 2  Bulmer: late index variance = %.4f   (theory v_eq = %.4f)\n",
            mean(late$varI), v_eq))
cat(sprintf("GATE 3  truncation intensity: late realized = %.4f   (theory i = %.4f)\n",
            mean(late$i_real), ti$i))
cat(sprintf("GATE 4  per-cycle gain: late = %.4f   (theory i*sqrt(v_eq) = %.4f); sd = %.4f\n",
            mean(late$gain), ti$i * sqrt(v_eq), sd(late$gain)))

## ---- figure ----
dir.create("figures", showWarnings = FALSE)
png("figures/null_recurrent_selection.png", width = 1600, height = 720, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(traj$cycle, traj$varI, type = "b", pch = 19, col = "#c0392b", ylim = c(0, v_a),
     xlab = "cycle", ylab = "index variance  b' G b", main = "(A) Bulmer equilibrium")
abline(h = v_eq, lty = 2, col = "grey40"); abline(h = v_a, lty = 3, col = "grey70")
legend("topright", bty = "n", lty = c(2, 3), col = c("grey40", "grey70"),
       legend = c("v_eq = v_a/(1+k)", "v_a (genic)"))
plot(traj$cycle, traj$meanI, type = "b", pch = 19, col = "#2980b9",
     xlab = "cycle", ylab = "mean index", main = "(B) steady response (no plateau)")
abline(lm(meanI ~ cycle, data = late), lty = 2, col = "grey40")
invisible(dev.off())
cat("\nWrote figures/null_recurrent_selection.png\n")
