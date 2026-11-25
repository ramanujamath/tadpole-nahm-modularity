import sys, math; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from detB import *
def qbin(n,k,q):
    if k<0 or k>n: return 0.0
    num=1.0
    for i in range(k): num*=(1-q**(n-i))/(1-q**(i+1))
    return num
def A(k,j,q): return (-1)**j*q**((k-j)*(k-j-1)/2)*qbin(k+j,2*j,q)
if __name__=='__main__':
    q=0.27; w=q**0.5
    for m in (3,4):
        half=[k+0.5 for k in range(m+1)]
        for xh in ([1.15,0.7],[1.4,1.25],[1.2,0.85,1.33]):
            if len(xh)!=m-1: continue
            x=list(xh)+[1.0]; L=2*m+2
            Dn=np.linalg.det(np.array([[xi**a+xi**(-a) for a in half[:m]] for xi in x]))
            P1=pref(1.0,xh,q)
            Hs=[HB2(q**j,xh,q,R=6) for j in range(m+1)]
            for k in range(1,m):
                tau=[h for h in half if h!=m+0.5-k]
                Jk=P1*np.linalg.det(np.array([[PhiB(a,xi,q,L) for a in tau] for xi in x]))/Dn
                lhs=sum(A(k,j,q)*Hs[j] for j in range(k+1))
                print(f"m={m} xhat={xh} k={k}: sum_j A_kj H(q^j) / (q^(k^2/2) J_k) = {lhs/(q**(k*k/2)*Jk):.12f}")
