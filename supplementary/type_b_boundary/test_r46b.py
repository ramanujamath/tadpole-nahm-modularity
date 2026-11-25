import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from bw56b import *
from evenhl import nahm, poch
q=0.27
for m in (2,3):
    for c in (0.0,0.5,1.0,-0.4):
        t=q**c
        lhs=nahm([0.0]*(2*m-1)+[c],q)/poch(-q**(c+0.5),q)
        out=[]
        for d in (4e-3,2e-3,1e-3):
            xh=[1+d*(k+1) for k in range(m-1)]
            out.append(HB2(t,xh,q))
        r1=2*out[1]-out[0]; r2=2*out[2]-out[1]; ext=(4*r2-r1)/3
        print(f"rank {2*m} c={c}: Nahm {lhs:.11f}  limit {ext:.11f}  rel diff {abs(ext/lhs-1):.1e}")
