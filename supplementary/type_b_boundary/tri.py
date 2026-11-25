import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from detB import *
q=0.27; w=q**0.5
for m in (2,3):
    half=[k+0.5 for k in range(m+1)]
    tau1=[h for h in half if h!=m-0.5]
    for xh in ([1.3],[0.8],[1.15,0.7],[1.4,1.25]):
        if len(xh)!=m-1: continue
        diff=HB2(1.0,xh,q,R=7)-HB2(q,xh,q,R=7)
        # J_1 analogue: P1 * det(Phi^+_tau1)/det(x^a+x^-a, a=1/2..m-1/2)
        x=list(xh)+[1.0]; L=2*m+2
        N=np.linalg.det(np.array([[PhiB(a,xi,q,L) for a in tau1] for xi in x]))
        Dn=np.linalg.det(np.array([[xi**a+xi**(-a) for a in half[:m]] for xi in x]))
        J1=pref(1.0,xh,q)*N/Dn
        print(f"m={m} xhat={xh}: [H(1)-H(q)]/J1 = {diff/J1:.12f}   (w={w:.12f})")
