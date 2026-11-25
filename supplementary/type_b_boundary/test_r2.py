import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from bw56b import *
from evenhl import nahm, poch
q=0.27
for c in (0.0,0.5,1.0,-0.3,0.37):
    t=q**c
    lhs=nahm([0.0,c],q)/poch(-q**(c+0.5),q)
    print(f"rank 2, c={c}: Nahm side {lhs:.13f}  reconstructed (5.6b) {HB(t,[],q):.13f}")
