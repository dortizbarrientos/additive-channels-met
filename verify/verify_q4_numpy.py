# verify_q4_numpy.py -- independent numpy cross-check of scripts/04.
#   (A) r_env formula == finite-difference logit-additivity rate
#   (B) identical present Ag(0), divergent futures ordered by r_env
#   python3 verify/verify_q4_numpy.py
import numpy as np
G = np.array([[1.,.3],[.3,1.]]); b0 = np.array([1.,.6]); H0 = .7*np.eye(2)
trsq = lambda M: np.sum(M*M.T)
def loc(u,bu,Hu):
    b=b0+u*bu; H=H0+u*Hu; Vl=b@G@b; Vq=0.5*trsq(H@G); return Vl,Vq,Vl/(Vl+Vq)
def r_env(udot,bu,Hu):
    Vl0=b0@G@b0; Vq0=0.5*trsq(H0@G); c00=Vl0+Vq0; Ag=Vl0/c00
    c10L=b0@G@bu; c10Q=0.5*np.trace(H0@G@Hu@G)
    return udot*(2*c10L/(c00*Ag)-2*c10Q/(c00*(1-Ag)))
def lam(u,bu,Hu): Vl,Vq,_=loc(u,bu,Hu); return np.log(Vl)-np.log(Vq)
scen={"eroding":(-.5*b0,.5*np.eye(2)),"stable":(.3*np.array([-.9,1.18]),np.zeros((2,2))),"deepening":(.5*b0,-.3*np.eye(2))}
print("(A) formula vs finite difference (udot=1)")
for nm,(bu,Hu) in scen.items():
    fd=(lam(1e-4,bu,Hu)-lam(-1e-4,bu,Hu))/2e-4
    print(f"  {nm:10} formula={r_env(1,bu,Hu):+.4f}  fd={fd:+.4f}")
print("\n(B) identical present, divergent futures")
print(f"  {'scenario':10} {'r_env':>7} {'Ag(0)':>6} {'acc(0)':>7} {'acc(1)':>7}")
for nm,(bu,Hu) in scen.items():
    _,_,a0=loc(0,bu,Hu); _,_,a1=loc(1,bu,Hu)
    print(f"  {nm:10} {r_env(1,bu,Hu):>+7.3f} {a0:>6.3f} {np.sqrt(a0):>7.3f} {np.sqrt(a1):>7.3f}")
