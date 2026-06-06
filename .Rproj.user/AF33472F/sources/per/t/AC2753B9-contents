# 02_eiv_correction.R
# Berkson errors-in-variables: how far does breeding-value estimation error
# attenuate Ag and AgGEI, and does a reliability-based correction restore them?
# Run from the repository root:  Rscript scripts/02_eiv_correction.R
source("R/channels.R")
set.seed(20260606)

G  <- matrix(c(1.0, 0.3, 0.3, 1.0), 2)
b0 <- c(1.00, 0.60); H0 <- 0.70 * diag(2)
bu <- c(0.45, 0.25); Hu <- 0.40 * diag(2)
truth <- true_channels(G, b0, bu, H0, Hu)

rhos <- c(0.3, 0.5, 0.7, 0.9)   # reliability of the breeding-value panel
N    <- 100000
reps <- 60

res <- data.frame()
for (rho in rhos) {
  ng <- cg <- na <- ca <- numeric(reps)
  for (r in seq_len(reps)) {
    s    <- simulate_eiv(N, G, b0, bu, H0, Hu, rho)
    Vobs <- cov(s$Zhat)
    Gc   <- Vobs + s$PEV                     # corrected metric = Var(z_hat)+PEV = G
    c11t <- var(s$theta); c00t <- var(s$alpha)
    ng[r] <- vlin_explained(s$theta, s$Zhat, Vobs) / c11t   # naive AgGEI
    cg[r] <- vlin_explained(s$theta, s$Zhat, Gc)   / c11t   # corrected AgGEI
    na[r] <- vlin_explained(s$alpha, s$Zhat, Vobs) / c00t   # naive Ag
    ca[r] <- vlin_explained(s$alpha, s$Zhat, Gc)   / c00t   # corrected Ag
  }
  res <- rbind(res, data.frame(rho = rho,
    AgGEI_naive = mean(ng), AgGEI_corr = mean(cg),
    Ag_naive    = mean(na), Ag_corr    = mean(ca)))
}

cat(sprintf("TRUTH: Ag = %.4f,  AgGEI = %.4f\n\n", truth$Ag, truth$AgGEI))
print(round(res, 4), row.names = FALSE)
cat("\nExpectation: naive ~ rho * truth (attenuation); corrected ~ truth.\n")

## ---- figure: attenuation and correction vs reliability ----
dir.create("figures", showWarnings = FALSE)
png("figures/eiv_attenuation.png", width = 1500, height = 750, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (q in c("AgGEI", "Ag")) {
  tval <- truth[[q]]
  plot(res$rho, res[[paste0(q, "_naive")]], type = "b", pch = 19, ylim = c(0, 1),
       xlab = "breeding-value reliability  rho", ylab = q, main = q, col = "#c0392b")
  lines(res$rho, res[[paste0(q, "_corr")]], type = "b", pch = 17, col = "#2980b9")
  abline(h = tval, lty = 2, col = "grey40")
  lines(res$rho, res$rho * tval, lty = 3, col = "#c0392b")    # rho * truth reference
  legend("topleft", bty = "n", cex = 0.8,
         legend = c("naive", "corrected", "truth", "rho*truth"),
         col = c("#c0392b", "#2980b9", "grey40", "#c0392b"),
         pch = c(19, 17, NA, NA), lty = c(1, 1, 2, 3))
}
invisible(dev.off())
cat("Wrote figures/eiv_attenuation.png\n")
