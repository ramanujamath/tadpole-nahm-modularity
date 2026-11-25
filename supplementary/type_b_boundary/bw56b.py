"""Reconstructed b->infinity limit of Bartlett-Warnaar (5.6b), level 1, letters (xhat^{±},1).
   H^B_n(t;x) = (1/D) sum_{r in Z^n} Delta_B(-x q^r)/Delta_B(-x) prod_i (-1)^{r_i} t^{r_i} q^{(n+1) r_i^2 + (1/2-n) r_i} x_i^{(2n+2) r_i} (-w x_i/t)_{r_i}/(-w t x_i)_{r_i}
   D = (q)^n prod_i (q x_i^{±2})_inf/(-w t x_i^{±1})_inf prod_{i<j} (q x_i^± x_j^±)_inf,  x_n = 1.
"""
import itertools, math
def pinf(a,q,T=1500):
    p=1.0
    for s in range(T): p*=(1-a*q**s)
    return p
def pfin(a,q,r):  # (a;q)_r for r in Z
    if r>=0:
        p=1.0
        for s in range(r): p*=(1-a*q**s)
        return p
    p=1.0
    for s in range(1,-r+1): p/= (1-a*q**(-s))
    return p
def DeltaB(y):
    n=len(y); p=1.0
    for i in range(n): p*=(1-y[i])
    for i in range(n):
        for j in range(i+1,n): p*=(y[i]-y[j])*(y[i]*y[j]-1)
    return p
def HB(t,xh,q,R=7):
    x=list(xh)+[1.0]; n=len(x); w=q**0.5
    D=pinf(q,q)**n
    for i in range(n):
        D*=pinf(q*x[i]**2,q)*pinf(q/x[i]**2,q)/(pinf(-w*t*x[i],q)*pinf(-w*t/x[i],q))
        for j in range(i+1,n):
            for a in (x[i]*x[j],x[i]/x[j],x[j]/x[i],1/(x[i]*x[j])): D*=pinf(q*a,q)
    base=DeltaB([-v for v in x]); tot=0.0
    for r in itertools.product(range(-R,R+1),repeat=n):
        s=DeltaB([-x[i]*q**r[i] for i in range(n)])/base
        for i in range(n):
            ri=r[i]
            s*=(-1)**(ri%2)*t**ri*q**((n+1)*ri*ri+(0.5-n)*ri)*x[i]**((2*n+2)*ri)*pfin(-w*x[i]/t,q,ri)/pfin(-w*t*x[i],q,ri)
        tot+=s
    return tot/D

def HB2(t,xh,q,R=6):
    """Corrected normalisation: the frozen letter x_n=1 counts once."""
    x=list(xh)+[1.0]; n=len(x); w=q**0.5
    D=pinf(q,q)**n/pinf(-w*t,q)
    for i in range(n-1):
        D*=pinf(q*x[i]**2,q)*pinf(q/x[i]**2,q)/(pinf(-w*t*x[i],q)*pinf(-w*t/x[i],q))
    for i in range(n-1):
        for j in range(i+1,n-1):
            for a in (x[i]*x[j],x[i]/x[j],x[j]/x[i],1/(x[i]*x[j])): D*=pinf(q*a,q)
        D*=pinf(q*x[i],q)*pinf(q/x[i],q)          # cross factor with the frozen letter, counted once
    base=DeltaB([-v for v in x]); tot=0.0
    for r in itertools.product(range(-R,R+1),repeat=n):
        s=DeltaB([-x[i]*q**r[i] for i in range(n)])/base
        for i in range(n):
            ri=r[i]
            s*=(-1)**(ri%2)*t**ri*q**((n+1)*ri*ri+(0.5-n)*ri)*x[i]**((2*n+2)*ri)*pfin(-w*x[i]/t,q,ri)/pfin(-w*t*x[i],q,ri)
        tot+=s
    return tot/D
