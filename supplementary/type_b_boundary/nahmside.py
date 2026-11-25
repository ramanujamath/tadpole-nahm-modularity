import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from evenhl import nahm, poch
from inv import A
q=0.27; w=q**0.5
def pfinite(a,q,j):
    p=1.0
    for s in range(j): p*=(1-a*q**s)
    return p
for m in (2,3,4,5):
    r=2*m
    for k in range(1,m):
        B=[0.0]*r
        for i in range(1,k+1):           # last k pairs (r-2i, r-2i+1) in 1-based -> indices r-2i-1, r-2i
            B[r-2*i-1]-=1; B[r-2*i]+=1
        G=nahm(B,q)
        rhs=sum(A(k,j,q)*pfinite(-w,q,j)*nahm([0.0]*(r-1)+[float(j)],q) for j in range(k+1))
        print(f"rank {r}, k={k}, B={[int(b) for b in B]}: q^(k^2/2) G / sum_j A_kj (-w)_j F_j = {q**(k*k/2)*G/rhs:.12f}")
