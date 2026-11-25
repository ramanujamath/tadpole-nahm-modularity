import sys, itertools, math
from fractions import Fraction as Fr
from classify import analyse
def ev(p,B): return sum(c*math.prod(Fr(b)**e for b,e in zip(B,k)) for k,c in p.items())
r=int(sys.argv[1]); N=int(sys.argv[2]); den=int(sys.argv[3]); box=int(sys.argv[4])
R,F,G,S,St=analyse(r,N)
vals=[Fr(a,den) for a in range(-box*den,box*den+1)]
# quick filter with lowest-degree poly first
F_sorted=sorted(F,key=lambda f:max(sum(k) for k in f))
pts=[]
for B in itertools.product(vals,repeat=r):
    if all(ev(f,B)==0 for f in F_sorted):
        pts.append(B)
C0=R['k1'][0]
print('rational zeros found:',len(pts))
for B in pts: print('  B=',tuple(str(b) for b in B),' C=',ev(C0,B))
