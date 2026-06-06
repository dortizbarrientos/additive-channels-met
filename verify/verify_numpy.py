# verify_numpy.py -- independent numpy implementation used to cross-check the R.
# Two languages, one set of identities: agreement is the check.
#   python3 verify/verify_numpy.py
import numpy as np
rng = np.random.default_rng(20260606)
G = np.array([[1., .3], [.3, 1.]]); b0 = np.array([1., .6]); H0 = .7 * np.eye(2)
bu = np.array([.45, .25]); Hu = .4 * np.eye(2)
trsq = lambda M: np.sum(M * M.T)
VlinG = bu @ G @ bu; c11 = VlinG + .5 * trsq(Hu @ G); AgGEI = VlinG / c11
Vlin0 = b0 @ G @ b0; c00 = Vlin0 + .5 * trsq(H0 @ G); Ag = Vlin0 / c00
mvn = lambda n, S: rng.standard_normal((n, len(S))) @ np.linalg.cholesky(S).T
rn  = lambda Z, b, H: Z @ b + .5 * np.einsum('ij,jk,ik->i', Z, H, Z)
def vlin(coeff, Zobs, Gmetric):
    s = np.array([np.cov(coeff, Zobs[:, k])[0, 1] for k in range(Zobs.shape[1])])
    beta = np.linalg.solve(np.cov(Zobs.T), s)
    return beta @ Gmetric @ beta
print(f"TRUTH  Ag={Ag:.4f}  AgGEI={AgGEI:.4f}\n")
print(f"{'rho':>5} | {'AgGEI naive':>11} {'AgGEI corr':>10} | {'Ag naive':>9} {'Ag corr':>8}")
N, reps = 100000, 60
for rho in (.3, .5, .7, .9):
    ng = cg = na = ca = 0.0
    for _ in range(reps):
        Zh = mvn(N, rho * G); Z = Zh + mvn(N, (1 - rho) * G)
        th = rn(Z, bu, Hu); al = rn(Z, b0, H0)
        Vo = np.cov(Zh.T); Gc = Vo + (1 - rho) * G
        ng += vlin(th, Zh, Vo) / np.var(th, ddof=1); cg += vlin(th, Zh, Gc) / np.var(th, ddof=1)
        na += vlin(al, Zh, Vo) / np.var(al, ddof=1); ca += vlin(al, Zh, Gc) / np.var(al, ddof=1)
    print(f"{rho:>5.1f} | {ng/reps:>11.4f} {cg/reps:>10.4f} | {na/reps:>9.4f} {ca/reps:>8.4f}")
