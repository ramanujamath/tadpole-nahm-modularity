import itertools, math, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from wron import nahm, poch
def chi0(v):
    p=1.0
    for i in range(len(v)):
        for j in range(i+1,len(v)): p*=(v[i]**2-v[j]**2)
    return p
def rhs(m,rho,rhos,z,q,sign,R=3):
    r=2*m; L=2*m+2; nu=r*(r-1)//2
    nr=sum(x*x for x in rho); tot=0.0
    for n in itertools.product(range(-R,R+1),repeat=m):
        v=[rho[i]+L*n[i] for i in range(m)]
        s=(-1)**sum(n) if sign else 1
        tot+=s*chi0(v)*q**((sum(x*x for x in v)-nr)/(2*L))
    return poch(-q**(0.5+z),q)**r/poch(q,q)**nu*tot/chi0(rhos)
def cases(m):
    half=[k+0.5 for k in range(m+1)]; ints=list(range(m+1))
    out=[]
    # sector (i): z=0, x chain with k ones; guess rho = half minus element
    for k in range(m+1):
        x=[0]*(m-k)+[1]*k; B=[x[0]]
        for xi in x[1:]: B+=[-xi,xi]
        B.append(0)
        out.append(('i',k,B,0,half,half[:m]))
    for j in range(m):
        y=[1]*j+[0]*(m-1-j); B=[]
        for yi in y: B+=[yi,-yi]
        B+=[0,0.5]; out.append(('ii',j,B,0.5,ints,ints[:m]))
    B=[1,-1]*(m-1)+[1,-0.5]; out.append(('iii',m,B,0.5,ints,ints[:m]))
    return out
if __name__!="__main__": sys.argv=[sys.argv[0],""]
q=float(sys.argv[2]) if len(sys.argv)>2 else 0.3
for m in [int(a) for a in sys.argv[1].split(",") if a]:
    for fam,k,B,z,pool,rhos in cases(m):
        l=nahm(B,q)
        best=[]
        for d in pool:
            rho=[p for p in pool if p!=d]
            for sign in (0,1):
                rr=rhs(m,rho,rhos,z,q,sign)
                best.append((abs(l/rr-1) if rr else 9,d,sign,l/rr if rr else 0))
        best.sort()
        e,d,sign,ratio=best[0]
        print(f"m={m} {fam:3} k={k} B={B}: best delete {d} sign={sign} ratio={ratio:.12f} relerr={e:.1e}")
