import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from bw56b import HB2, pinf
from evenhl import nahm, poch
q=0.27; w=q**0.5
targets={'(i) k=2 (1,-1,1,0)':[1,-1,1,0],'(ii) j=1 (1,-1,0,1/2)':[1,-1,0,0.5],'(iii) (1,-1,1,-1/2)':[1,-1,1,-0.5],'(iv) (2,-2,1,-1/2)':[2,-2,1,-0.5]}
T={k:nahm(v,q) for k,v in targets.items()}
for tn,t in [('1',1.0),('q',q),('w',w),('q^{3/2}',q**1.5),('w^-1',1/w)]:
    vals=[]
    for d in (1e-3,5e-4):
        vals.append(HB2(t,[w*(1+d)],q,R=9))
    v=2*vals[1]-vals[0]
    print(f"t={tn}: H^B(t; x1->w) ~ {v:.10g}   ratios to targets:",{k:round(T[k]/v,9) for k in T})
