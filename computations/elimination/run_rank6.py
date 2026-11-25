import sys, time, pickle
from classify import obstructions
from gbmodp import GB, standard_monomials, top_part
N=int(sys.argv[1]); t=time.time()
R,F=obstructions(6,N)
print('rank 6 N=%d: %d polys, degrees %s (%.0fs)'%(N,len(F),[max(sum(k) for k in f) for f in F],time.time()-t),flush=True)
pickle.dump((R['k1'][0],F),open('r6_N%d.pkl'%N,'wb'))
gb=GB(1000003); G=gb.groebner(F); S=standard_monomials(G,6,cap=40)
print('affine quotient dim mod p:',None if S is None else len(S),'(%.0fs)'%(time.time()-t),flush=True)
pickle.dump(G,open('r6_N%d_gb.pkl'%N,'wb'))
Gt=gb.groebner([top_part(f) for f in F]); St=standard_monomials(Gt,6,cap=40)
print('top-part finite:',St is not None, '' if St is None else len(St),'(%.0fs)'%(time.time()-t),flush=True)
