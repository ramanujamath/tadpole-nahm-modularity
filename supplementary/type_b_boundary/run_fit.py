import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fitprod import *
def Hlim(t_of_q,q):
    w=q**0.5; vals=[HB2(t_of_q(q),[w*(1+d)],q,R=9) for d in (1e-3,5e-4)]
    return 2*vals[1]-vals[0]
targets={'(i) k=2':[1,-1,1,0],'(ii) j=1':[1,-1,0,0.5],'(iii)':[1,-1,1,-0.5],'(iv)':[2,-2,1,-0.5]}
for tn,tf in [('t=1',lambda q:1.0),('t=w',lambda q:q**0.5),('t=q',lambda q:q),('t=1/w',lambda q:q**-0.5)]:
    for name,B in targets.items():
        fit(lambda q,B=B,tf=tf: nahm(B,q)/Hlim(tf,q), f"{name} / H^B({tn}; w):")
