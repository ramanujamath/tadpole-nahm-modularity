import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from bw56b import *
from evenhl import nahm, poch
import itertools
q=0.27; w=q**0.5
def S(t,R=8):  # raw r-sum for n=1
    tot=0.0
    for r in range(-R,R+1):
        tot+=(1+q**r)/2*(-1)**(r%2)*t**r*q**(2*r*r-0.5*r)*pfin(-w/t,q,r)/pfin(-w*t,q,r)
    return tot
for c in (0.0,0.5,1.0,-0.3):
    t=q**c; H=nahm([0.0,c],q)/poch(-q**(c+0.5),q)
    ratio=H/S(t)
    cands={'(-wt)/(q)':pinf(-w*t,q)/pinf(q,q),'(-wt)^2/(q)^3':pinf(-w*t,q)**2/pinf(q,q)**3,'(-wt)/(q)^2':pinf(-w*t,q)/pinf(q,q)**2,
           '(-wt)(-w/t)/(q)^2':pinf(-w*t,q)*pinf(-w/t,q)/pinf(q,q)**2}
    print(f"c={c}: H/S={ratio:.12f} ", {k: round(ratio/v,10) for k,v in cands.items()})
