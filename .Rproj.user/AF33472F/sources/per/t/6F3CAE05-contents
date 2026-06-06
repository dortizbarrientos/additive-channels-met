# channels.R -- core functions for additive-channel MET diagnostics
# Companion code to Ortiz-Barrientos, Messina & Cooper (in prep.).
# Base R only; no external packages.
#
# Notation
#   p      number of component traits (length of a breeding-value vector)
#   G      p x p trait-space additive genetic covariance matrix (the G-matrix)
#   b0,bu  p-vectors: local fitness gradient and its environmental sensitivity
#   H0,Hu  p x p symmetric matrices: local Hessian and its sensitivity
#   z_i ~ N(0,G)  a genotype's multivariate breeding-value deviation
#
# Lane A (see README). When b(u)=b0+du*bu and H(u)=H0+du*Hu, the genetic
# reaction norm is exactly linear in the environment,
#       g_i(u) = alpha_i + du * theta_i,
# with
#       alpha_i = b0' z_i + 0.5 z_i' H0 z_i      (intercept)
#       theta_i = bu' z_i + 0.5 z_i' Hu z_i      (slope)
# Totals (target trait only):  c00=Var(alpha), c11=Var(theta), c10=Cov(alpha,theta).
# The transferable (linear) part of any such variance is the genetic R^2 of the
# coefficient regressed LINEARLY on the multivariate breeding values -- it needs
# the trait panel, not just the target trait.

## tr(M %*% M) without forming the product: tr(M^2) = sum_{i,j} M_ij M_ji
trace_sq <- function(M) sum(M * t(M))

## closed-form ("oracle") channel quantities from the structural parameters
true_channels <- function(G, b0, bu, H0, Hu) {
  Vlin0  <- as.numeric(t(b0) %*% G %*% b0)            # gradient-aligned variance at u0
  Vquad0 <- 0.5 * trace_sq(H0 %*% G)                  # curvature-exposed variance at u0
  c00    <- Vlin0 + Vquad0
  VlinG  <- as.numeric(t(bu) %*% G %*% bu)            # transferable (linear) GEI
  VquadG <- 0.5 * trace_sq(Hu %*% G)                  # curvature-driven GEI
  c11    <- VlinG + VquadG
  c10L   <- as.numeric(t(b0) %*% G %*% bu)            # linear part of c10
  c10Q   <- 0.5 * sum(diag(H0 %*% G %*% Hu %*% G))    # curvature part of c10
  list(Vlin0 = Vlin0, Vquad0 = Vquad0, c00 = c00, Ag = Vlin0 / c00,
       VlinGEI = VlinG, VquadGEI = VquadG, c11 = c11, AgGEI = VlinG / c11,
       c10L = c10L, c10Q = c10Q, c10 = c10L + c10Q)
}

## draw n rows ~ MVN(0, S) via Cholesky.  chol(S) is upper U with U'U = S,
## so X = N(0,I) %*% U has Cov(X) = U'U = S.
rmvn <- function(n, S) matrix(rnorm(n * ncol(S)), n) %*% chol(S)

## reaction-norm coefficient for each genotype: b' z + 0.5 z' H z
## rowSums((Z %*% H) * Z) computes z_i' H z_i for every row i at once.
rn_coeff <- function(Z, b, H) as.vector(Z %*% b) + 0.5 * rowSums((Z %*% H) * Z)

## Transferable (linear) variance of a coefficient.
## The regression coefficient is taken from the OBSERVED breeding-value cloud;
## the variance is then evaluated under the supplied metric Gmetric:
##   naive       Gmetric = Var(Zobs)
##   corrected   Gmetric = Var(Zobs) + PEV        (Berkson reliability correction)
## Under Berkson (BLUP) error the regression coefficient is UNBIASED, so the only
## fix needed is to evaluate its variance under the reconstructed G = Var+PEV.
vlin_explained <- function(coeff, Zobs, Gmetric) {
  Vobs <- cov(Zobs)
  s    <- as.numeric(cov(Zobs, coeff))   # Cov(z_obs, coeff)
  beta <- solve(Vobs, s)                 # best linear coefficient
  as.numeric(t(beta) %*% Gmetric %*% beta)
}

## simulate genotypes with TRUE breeding values observed (no estimation error)
simulate_true <- function(n, G, b0, bu, H0, Hu) {
  Z <- rmvn(n, G)
  list(Z = Z, alpha = rn_coeff(Z, b0, H0), theta = rn_coeff(Z, bu, Hu))
}

## simulate Berkson (BLUP) estimation error at reliability rho in (0,1]:
##   z_hat ~ N(0, rho*G),  delta ~ N(0,(1-rho)*G)  independent,  z = z_hat + delta
## Then z ~ N(0,G), Cov(z_hat,z)=Var(z_hat)=rho*G, and PEV = (1-rho)*G.
## Reaction-norm coefficients are built from the TRUE z (the response side is
## treated as observed here, isolating breeding-value-panel error).
simulate_eiv <- function(n, G, b0, bu, H0, Hu, rho) {
  Zhat  <- rmvn(n, rho * G)
  delta <- rmvn(n, (1 - rho) * G)
  Z     <- Zhat + delta
  list(Zhat = Zhat, Z = Z,
       alpha = rn_coeff(Z, b0, H0), theta = rn_coeff(Z, bu, Hu),
       PEV = (1 - rho) * G)
}

## ---------------------------------------------------------------------------
## Additions for the transfer-accuracy experiment (scripts/03).
## ---------------------------------------------------------------------------

## Local channel quantities at environment coordinate u (reference u0 = 0):
##   b(u)=b0+u*bu, H(u)=H0+u*Hu, and Ag(u)=Vlin(u)/(Vlin(u)+Vquad(u)).
local_channels <- function(u, G, b0, bu, H0, Hu) {
  b <- b0 + u * bu
  H <- H0 + u * Hu
  Vlin  <- as.numeric(t(b) %*% G %*% b)
  Vquad <- 0.5 * trace_sq(H %*% G)
  list(Vlin = Vlin, Vquad = Vquad, c = Vlin + Vquad, Ag = Vlin / (Vlin + Vquad))
}

## Build (bu, Hu) giving a target GEI additivity index AgGEI at a target total
## GEI variance c11, holding the bu direction (dir_bu) and using isotropic Hu.
##   bu' G bu = AgGEI*c11 ;  0.5 tr[(Hu G)^2] = (1-AgGEI)*c11
make_gei <- function(AgGEI, c11, G, dir_bu = c(1, 0)) {
  Vlin  <- AgGEI * c11
  Vquad <- (1 - AgGEI) * c11
  d  <- dir_bu / sqrt(as.numeric(t(dir_bu) %*% G %*% dir_bu))   # unit in G-metric
  bu <- d * sqrt(Vlin)
  h  <- sqrt(2 * Vquad / trace_sq(G))
  list(bu = bu, Hu = h * diag(nrow(G)))
}

## Additive reaction-norm transfer accuracy to a held-out environment u_t.
## Trains the additive (linear-in-z) intercept and slope maps, predicts the TRUE
## genetic value at u_t, and returns Pearson and Spearman transfer accuracy plus
## the local additivity Ag(u_t).  In the noise-free large-n limit the additive
## reaction-norm fit is fold-invariant, so leave-one-environment-out transfer to
## u_t reduces to transfer accuracy to u_t.
transfer_accuracy <- function(n, G, b0, bu, H0, Hu, ut) {
  Z     <- rmvn(n, G)
  alpha <- rn_coeff(Z, b0, H0)
  theta <- rn_coeff(Z, bu, Hu)
  Vz    <- cov(Z)
  bhat0 <- solve(Vz, as.numeric(cov(Z, alpha)))   # additive intercept map
  bhatu <- solve(Vz, as.numeric(cov(Z, theta)))   # additive slope map
  ghat  <- as.vector(Z %*% (bhat0 + ut * bhatu))  # additive prediction at u_t
  gtrue <- alpha + ut * theta                     # true genetic value at u_t
  c(pearson  = cor(ghat, gtrue),
    spearman = cor(ghat, gtrue, method = "spearman"),
    Ag_ut    = local_channels(ut, G, b0, bu, H0, Hu)$Ag)
}

## Environmental rate signal: the rate of change of logit-additivity as the
## target environment drifts at rate udot.  sign(r_env) = sign(d Ag/du).
## r_env < 0 means the additive-prediction channel is ERODING -- future transfer
## will decline -- even when present additivity Ag(u0) still looks healthy.
## It uses only currently estimable quantities (the local gradient/curvature
## sensitivities the current MET provides), so it forecasts before the drift.
r_env <- function(udot, G, b0, bu, H0, Hu) {
  tc <- true_channels(G, b0, bu, H0, Hu)
  udot * (2 * tc$c10L / (tc$c00 * tc$Ag) -
          2 * tc$c10Q / (tc$c00 * (1 - tc$Ag)))
}
