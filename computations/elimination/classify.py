import sys, math, time, itertools, json
from fractions import Fraction as Fr
from exactcusp import cusp_polys
from gbmodp import GB, standard_monomials, top_part

def to_int(p):
    den=1
    for v in p.values(): den=den*v.denominator//math.gcd(den,v.denominator)
    q={k:int(v*den) for k,v in p.items()}
    g=0
    for v in q.values(): g=math.gcd(g,v)
    return {k:v//g for k,v in q.items()}

def obstructions(r,N):
    R=cusp_polys(r,N)
    gens=R['k1'][1:]+R['k2']+R['k3']+(R.get('k4') or [])
    return R,[to_int(g) for g in gens if g]

def analyse(r,N=3,p=1000003):
    t=time.time()
    R,F=obstructions(r,N)
    print(f'rank {r}: field degree {R["K"].d}, {len(F)} obstruction polys, degrees',[max(sum(k) for k in f) for f in F],f'({time.time()-t:.1f}s)')
    gb=GB(p)
    G=gb.groebner(F); S=standard_monomials(G,r)
    T=[top_part(f) for f in F]
    Gt=gb.groebner(T); St=standard_monomials(Gt,r)
    print('  affine quotient dim mod p:',None if S is None else len(S),'  top-part quotient finite:',St is not None,'' if St is None else f'(dim {len(St)})')
    return R,F,G,S,St

if __name__=='__main__':
    r=int(sys.argv[1]); N=int(sys.argv[2]) if len(sys.argv)>2 else 3
    analyse(r,N)
