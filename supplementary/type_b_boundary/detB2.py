import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from detB import *
q=0.27; w=q**0.5
psi=lambda y,a: (y**a+y**(-a))/(1+w*w+w*(y+1/y))
for m in (2,3):
    half=[k+0.5 for k in range(m)]; ints=[float(k) for k in range(m)]
    for xh in ([1.3],[0.8],[1.15,0.7],[1.4,1.25]):
        if len(xh)!=m-1: continue
        hq=HB2(q,xh,q,R=7); hw=HB2(w,xh,q,R=7)
        rq=pref(1.0,xh,q)*detform(xh,q,half,psi)
        rw_int=pref(w,xh,q)*detform(xh,q,ints)
        print(f"m={m} xhat={xh}: t=q ratio {hq/rq:.12f}   t=w (integer indices) ratio {hw/rw_int:.12f}")
