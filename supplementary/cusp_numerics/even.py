import sys
from fractions import Fraction as Fr
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from cusp import cusp_coeffs, Tmat
def even_pairs(m):
    r=2*m; P=[]
    for k in range(m+1):              # sector ending in 0
        x=[0]*(m-k)+[1]*k            # x_0..x_{m-1} ordered
        B=[x[0]]
        for xi in x[1:]: B+=[-xi,xi]
        B.append(0)
        P.append(('Z',k,B,Fr(-r*(2*r+1),48*(r+2))+Fr(k*(r+1-k),2*(r+2))))
    for j in range(m):               # sector ending in 1/2
        y=[1]*j+[0]*(m-1-j)
        B=[]
        for yi in y: B+=[yi,-yi]
        B+=[0,Fr(1,2)]
        P.append(('H',j,B,Fr(-r*(r-4),24*(r+2))+Fr(j*(r-j),2*(r+2))))
    B=[1,-1]*(m-1)+[1,Fr(-1,2)]
    P.append(('N',m,B,Fr(r,12)))
    B2=list(B); B2[0]+=1; B2[1]-=1
    P.append(('N2',m-1,B2,Fr(-r*(r-4),24*(r+2))+Fr((m-1)*(r-m+1),2*(r+2))+1))
    return P
if __name__=='__main__':
    for m in [int(a) for a in sys.argv[1:]]:
        r=2*m; A=Tmat(r); worst=0
        for fam,k,B,C in even_pairs(m):
            a,k1,k2,k3=cusp_coeffs(A,[float(b) for b in B])
            worst=max(worst,abs(k1-float(C)),abs(k2),abs(k3))
            print(f"r={r} {fam:2} k={k} B={[str(b) for b in B]} C={C} |c1-g0-C|={abs(k1-float(C)):.1e} k2={k2:.1e} k3={k3:.1e}",flush=True)
        print('rank',r,'pairs',len(even_pairs(m)),'max deviation',worst,flush=True)
