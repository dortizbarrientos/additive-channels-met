# verify_engine_numpy.py -- independent numpy cross-check of the null engine.
#   python3 verify/verify_engine_numpy.py
import numpy as np
from statistics import NormalDist
rng = np.random.default_rng(20260606); nd = NormalDist()
Ga = np.array([[1.,.3],[.3,1.]]); b = np.array([1.,.6])
N, prop, ncyc = 8000, 0.20, 30
x = nd.inv_cdf(1-prop); i_ = nd.pdf(x)/prop; k = i_*(i_-x)
v_eq = (b@Ga@b)/(1+k)
Z = rng.multivariate_normal([0,0], Ga, N)
vI=[]; mI=[]; ire=[]; err=[]
for t in range(ncyc):
    G = np.cov(Z.T); I = Z@b
    pred = i_*(G@b)/np.sqrt(b@G@b)
    thr = np.quantile(I, 1-prop); sel = Z[I>=thr]
    diff = sel.mean(0)-Z.mean(0)
    err.append(np.linalg.norm(diff-pred)/np.linalg.norm(pred))
    ire.append((sel@b).mean()-I.mean()); ire[-1]/=I.std(ddof=1)
    vI.append(b@G@b); mI.append(I.mean())
    m=sel.shape[0]
    Z = 0.5*(sel[rng.integers(0,m,N)]+sel[rng.integers(0,m,N)]) + rng.multivariate_normal([0,0],0.5*Ga,N)
mI=np.array(mI)
print(f"theory i={i_:.4f} k={k:.4f} v_eq={v_eq:.4f}")
print(f"GATE 1 breeders-eq rel err (late) = {np.mean(err[-10:]):.4f}")
print(f"GATE 2 Bulmer var (late)          = {np.mean(vI[-10:]):.4f}  vs v_eq {v_eq:.4f}")
print(f"GATE 3 intensity (late)           = {np.mean(ire[-10:]):.4f}  vs i {i_:.4f}")
print(f"GATE 4 gain (late)                = {np.mean(np.diff(mI)[-10:]):.4f}  vs i*sqrt(v_eq) {i_*np.sqrt(v_eq):.4f}")
