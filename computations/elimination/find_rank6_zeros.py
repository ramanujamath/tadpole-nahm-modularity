import pickle, itertools, numpy as np, math, sys
from fractions import Fraction as Fr
C0,F=pickle.load(open('r6_N3.pkl','rb'))
den=int(sys.argv[1]); box=Fr(sys.argv[2])
vals=np.array([a/den for a in range(int(-box*den),int(box*den)+1)])
Fs=sorted(F,key=lambda f:max(sum(k) for k in f))
def evalp(p,X):
    out=np.zeros(X.shape[1])
    for k,c in p.items():
        t=np.full(X.shape[1],float(c))
        for i,e in enumerate(k):
            if e: t=t*X[i]**e
        out+=t
    return out
def ev(p,B): return sum(c*math.prod(Fr(b)**e for b,e in zip(B,k)) for k,c in p.items())
sols=set(); total=0
for a0 in vals:
  for a1 in vals:
    rest=np.array(np.meshgrid(*([vals]*4),indexing='ij')).reshape(4,-1)
    grid=np.vstack([np.full(rest.shape[1],a0),np.full(rest.shape[1],a1),rest]); total+=grid.shape[1]
    mask=np.ones(grid.shape[1],bool)
    for f in Fs:
        idx=np.where(mask)[0]
        if len(idx)==0: break
        v=evalp(f,grid[:,idx]); scale=sum(abs(float(c)) for c in f.values())
        mask[idx[np.abs(v)>1e-9*scale]]=False
    for B in grid[:,mask].T:
        Bf=tuple(Fr(b).limit_denominator(1000) for b in B)
        if all(ev(f,Bf)==0 for f in F): sols.add(Bf)
print('grid points',total,' exact zeros',len(sols))
for B in sorted(sols): print('  B=('+', '.join(str(b) for b in B)+')  C=',ev(C0,B))
