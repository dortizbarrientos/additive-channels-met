# verify_q3_numpy.py -- independent numpy cross-check of scripts/03.
#   (A) additive transfer Pearson^2 == Ag(u_t)
#   (B) far-environment transfer rises with AgGEI (both metrics)
#   (C) Spearman - Pearson flips sign with local additivity
#   python3 verify/verify_q3_numpy.py
import numpy as np
rng = np.random.default_rng(20260606)
G = np.array([[1., .3], [.3, 1.]]); b0 = np.array([1., .6]); H0 = .7 * np.eye(2)
u_all = np.array([-1., -.5, 0., .5, 1.])
trsq = lambda M: np.sum(M * M.T)
def rn(Z, b, H): return Z @ b + 0.5 * np.einsum('ij,jk,ik->i', Z, H, Z)
def ranks(x):
    r = np.empty(len(x)); r[np.argsort(x)] = np.arange(len(x)); return r
def spear(x, y): return np.corrcoef(ranks(x), ranks(y))[0, 1]
def Ag_at(u, b0, bu, H0, Hu):
    b = b0 + u * bu; H = H0 + u * Hu
    Vl = b @ G @ b; Vq = 0.5 * trsq(H @ G); return Vl / (Vl + Vq)
def make_gei(a, c11, d=np.array([.45, .25])):
    d = d / np.sqrt(d @ G @ d); return d * np.sqrt(a * c11), np.sqrt(2 * (1 - a) * c11 / trsq(G)) * np.eye(2)
def transfer(N, b0, bu, H0, Hu, ut):
    Z = rng.multivariate_normal(np.zeros(2), G, N)
    al = rn(Z, b0, H0); th = rn(Z, bu, Hu); Vz = np.cov(Z.T)
    bh0 = np.linalg.solve(Vz, [np.cov(Z[:, k], al)[0, 1] for k in range(2)])
    bhu = np.linalg.solve(Vz, [np.cov(Z[:, k], th)[0, 1] for k in range(2)])
    gh = Z @ (bh0 + ut * bhu); gt = al + ut * th
    return np.corrcoef(gh, gt)[0, 1], spear(gh, gt)

print("(A) Pearson^2 == Ag(u_t)   (AgGEI=0.5)")
bu, Hu = make_gei(.5, .5)
for ut in u_all:
    rP, rS = transfer(60000, b0, bu, H0, Hu, ut)
    print(f"  u_t={ut:>4.1f}  Pearson^2={rP**2:.4f}  Ag(u_t)={Ag_at(ut,b0,bu,H0,Hu):.4f}  Spearman={rS:.4f}")

print("\n(B) far-environment (u_t=1) transfer vs AgGEI")
for a in (.1, .3, .5, .7, .9):
    bu, Hu = make_gei(a, .5)
    P = S = 0.0
    for _ in range(25):
        rP, rS = transfer(40000, b0, bu, H0, Hu, 1.0); P += rP; S += rS
    print(f"  AgGEI={a:.1f}  Pearson={P/25:.4f}  Spearman={S/25:.4f}")

print("\n(C) Spearman - Pearson vs local additivity (curvature-dominant)")
b0w, buw = np.array([.2, .1]), np.array([.1, .05])
for hh in (.1, .3, .6, .9, 1.2):
    rP, rS = transfer(120000, b0w, buw, np.zeros((2, 2)), hh * np.eye(2), 1.0)
    print(f"  Ag(u_t)={Ag_at(1,b0w,buw,np.zeros((2,2)),hh*np.eye(2)):.3f}  gap(S-P)={rS-rP:+.4f}")
