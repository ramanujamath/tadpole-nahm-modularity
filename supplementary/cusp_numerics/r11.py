import sys, math, itertools, time
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from cusp import cusp_coeffs, Tmat
def conj_pairs(m):
    r=2*m+1; out=[]
    for z in (-0.5,0,0.5):
        for k in range(m+1):
            x=[0]*(m-k)+[1]*k; B=[]
            for xi in x: B+=[2*z-xi, xi-2*z]
            B.append(z)
            C=-r*(2*r+1)/(48*(r+2))+k*(r+1-k)/(2*(r+2))-k*z+(2*r-1)/4*z*z
            out.append((k,z,B,C))
    return out
m=int(sys.argv[1]); r=2*m+1; A=Tmat(r)
worst=0
for k,z,B,C in conj_pairs(m):
    t=time.time()
    alpha,k1,k2,k3=cusp_coeffs(A,B)
    worst=max(worst,abs(k1-C),abs(k2),abs(k3))
    print(f"k={k} z={z:+.1f}: c1-g0={k1:.12f} predicted C={C:.12f} kappa2={k2:.1e} kappa3={k3:.1e} ({time.time()-t:.0f}s)",flush=True)
print('max deviation',worst)
# control: perturbed shifts
import random
random.seed(1)
for _ in range(2):
    B=[random.choice([-1,-0.5,0,0.5,1]) for _ in range(r)]
    a,k1,k2,k3=cusp_coeffs(A,B); print('control B',B,'kappa2=%.3e kappa3=%.3e'%(k2,k3))
