import sys, itertools, numpy as np
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from bw56b import *
def PhiB(a,x,q,L,f=None,R=8):
    tot=0.0
    for l in range(-R,R+1):
        y=x*q**l
        val=(y**a+y**(-a)) if f is None else f(y,a)
        tot+=(-1)**(l%2)*q**(L*l*l/2)*x**(L*l)*val
    return tot
def pref(t,xh,q):
    m=len(xh)+1; w=q**0.5
    num=pinf(-w*t,q); den=pinf(q,q)**m
    for i,x in enumerate(xh):
        num*=pinf(-w*t*x,q)*pinf(-w*t/x,q)
        den*=pinf(q*x*x,q)*pinf(q/(x*x),q)*pinf(q*x,q)*pinf(q/x,q)
        for y in xh[i+1:]:
            for a in (x*y,x/y,y/x,1/(x*y)): den*=pinf(q*a,q)
    return num/den
def detform(xh,q,avec,f=None):
    x=list(xh)+[1.0]; m=len(x); L=2*m+2
    N=np.linalg.det(np.array([[PhiB(a,xi,q,L,f) for a in avec] for xi in x]))
    Dn=np.linalg.det(np.array([[xi**a+xi**(-a) for a in avec] for xi in x]))
    return N/Dn
if __name__=='__main__':
    q=0.27
    for m in (2,3):
        half=[k+0.5 for k in range(m)]
        for xh in ([1.3],[0.8],[1.15,0.7],[1.4,1.25]):
            if len(xh)!=m-1: continue
            lhs=HB2(1.0,xh,q,R=7); rhs=pref(1.0,xh,q)*detform(xh,q,half)
            print(f"m={m} xhat={xh}: bilateral {lhs:.12f}  pref*det/den {rhs:.12f}  ratio {lhs/rhs:.12f}")
