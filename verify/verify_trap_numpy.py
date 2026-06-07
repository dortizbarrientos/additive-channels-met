# verify_trap_numpy.py -- independent numpy cross-check of the curvature trap.
#   python3 verify/verify_trap_numpy.py
import numpy as np
from statistics import NormalDist
rng = np.random.default_rng(20260606); nd = NormalDist()
Ga = np.array([[1.,.3],[.3,1.]]); trsq = lambda M: np.sum(M*M.T)
centers=[np.array([6.,0.]),np.array([15.,15.])]; heights=[1.,2.5]; widths=[3.5,3.5]
def Lc(z): return sum(a*np.exp(-((z-c)@(z-c))/(2*s*s)) for c,a,s in zip(centers,heights,widths))
def gc(z):
    o=np.zeros(2)
    for c,a,s in zip(centers,heights,widths):
        gk=a*np.exp(-((z-c)@(z-c))/(2*s*s)); o-=gk*(z-c)/s**2
    return o
def Hc(z):
    M=np.zeros((2,2))
    for c,a,s in zip(centers,heights,widths):
        d=z-c; gk=a*np.exp(-(d@d)/(2*s*s)); M+=gk*(np.outer(d,d)/s**4-np.eye(2)/s**2)
    return M
blin=gc(np.zeros(2)); blin=blin/np.linalg.norm(blin)*0.30
def Ll(z): return blin@z
def gl(z): return blin
def Hl(z): return np.zeros((2,2))
def run(L,g,H,prop=0.5,ncyc=40,N=5000):
    Z=rng.multivariate_normal([0,0],Ga,N); cumP=cumR=0; rec=[]
    for t in range(ncyc):
        zbar=Z.mean(0); G=np.cov(Z.T); b=g(zbar); Hh=H(zbar)
        Vlin=b@G@b; Vq=0.5*trsq(Hh@G); Ag=Vlin/(Vlin+Vq)
        I=Z@b; thr=np.quantile(I,1-prop); sel=Z[I>=thr]; m=sel.shape[0]
        Zn=0.5*(sel[rng.integers(0,m,N)]+sel[rng.integers(0,m,N)])+rng.multivariate_normal([0,0],0.5*Ga,N)
        zbn=Zn.mean(0); cumP+=b@(zbn-zbar); cumR+=L(zbn)-L(zbar)
        rec.append((L(zbar),Ag,cumP,cumR)); Z=Zn
    return np.array(rec)
cur=run(Lc,gc,Hc); lin=run(Ll,gl,Hl); far=Lc(centers[1])
print(f"near={Lc(centers[0]):.2f} far={far:.2f}")
print(f"CURVED  logW {cur[0,0]:.3f}->{cur[-1,0]:.3f}  Ag {cur[0,1]:.3f}->{cur[-1,1]:.3f}  pred {cur[-1,2]:.2f} vs real {cur[-1,3]:.2f}")
print(f"LINEAR  logW {lin[0,0]:.3f}->{lin[-1,0]:.3f}  Ag {lin[0,1]:.3f}->{lin[-1,1]:.3f}  pred {lin[-1,2]:.2f} vs real {lin[-1,3]:.2f}")
