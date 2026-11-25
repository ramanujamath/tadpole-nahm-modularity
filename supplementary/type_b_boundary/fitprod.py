import sys, itertools, math
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from bw56b import HB2, pinf
from evenhl import nahm
def fit(ratio_fn, label):
    best=[]
    for q in (0.27,0.19):
        pass
    qs=(0.27,0.19)
    R=[ratio_fn(q) for q in qs]
    w=lambda q:q**0.5
    basis=lambda q:[pinf(-q,q),pinf(-q**0.5,q),pinf(q,q),pinf(q,q*q),pinf(-q,q*q),pinf(q*q,q*q)]
    names=['(-q;q)','(-w;q)','(q;q)','(q;q^2)','(-q;q^2)','(q^2;q^2)']
    consts={'1':lambda q:1,'2':lambda q:2,'1/2':lambda q:.5,'w':lambda q:q**.5,'1/w':lambda q:q**-.5,'q':lambda q:q,'1/q':lambda q:1/q,'2w':lambda q:2*q**.5,'2/w':lambda q:2*q**-.5,'w/2':lambda q:q**.5/2,'1+w':lambda q:1+q**.5,'q^-3/2':lambda q:q**-1.5}
    B=[basis(q) for q in qs]
    out=[]
    for ex in itertools.product(range(-3,4),repeat=6):
        for cn,cf in consts.items():
            ok=True
            for qi,q in enumerate(qs):
                val=cf(q)*math.prod(b**e for b,e in zip(B[qi],ex))
                if abs(val/R[qi]-1)>1e-8: ok=False; break
            if ok: out.append((cn,ex))
    print(label, 'matches:', [(c,' '.join(f'{n}^{e}' for n,e in zip(names,ex) if e)) for c,ex in out][:4])
