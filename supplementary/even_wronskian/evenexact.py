import sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from prodexp import series
def chi0(v):
    p=Fr(1)
    for i in range(len(v)):
        for j in range(i+1,len(v)): p*=(v[i]**2-v[j]**2)
    return p
def smul(a,b,N):
    r=[0]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j in range(N+1-i): r[i+j]+=x*b[j]
    return r
def inv(a,N):
    r=[Fr(0)]*(N+1); r[0]=Fr(1)/a[0]
    for n in range(1,N+1): r[n]=-sum(a[k]*r[n-k] for k in range(1,n+1))/a[0]
    return r
def rhs(m,tau,taus,z,c,N):  # series in t=q^(1/2), N+1 terms
    r=2*m; L=2*m+2; nu=r*(r-1)//2
    th=[Fr(0)]*(N+1); c0=chi0(taus)
    R=1
    while L*R*R-2*max(tau)*R<=N: R+=1
    for n in itertools.product(range(-R,R+1),repeat=m):
        e=sum(2*tau[i]*n[i]+L*n[i]*n[i] for i in range(m))
        if 0<=e<=N:
            v=[tau[i]+L*n[i] for i in range(m)]
            th[int(e)]+=(1 if sum(n)%2==0 else -1)*chi0(v)/c0
    num=[Fr(1)]+[Fr(0)]*N
    s=0
    while True:
        e=int(1+2*z+2*s)
        if e>N: break
        f=[0]*(N+1); f[0]=1; f[e]+=1
        for _ in range(r): num=smul(num,f,N)
        s+=1
    den=[Fr(1)]+[Fr(0)]*N
    for s in range(1,N//2+1):
        f=[0]*(N+1); f[0]=1; f[2*s]=-1
        g=inv(f,N)
        for _ in range(nu): den=smul(den,g,N)
    return [c*x for x in smul(smul(num,den,N),th,N)]
def fam(m):
    half=[Fr(2*k+1,2) for k in range(m+1)]; ints=[Fr(k) for k in range(m+1)]
    out=[]
    for k in range(m+1):
        x=[0]*(m-k)+[1]*k; B=[x[0]]
        for xi in x[1:]: B+=[-xi,xi]
        B.append(0)
        out.append(('i',k,B,[h for h in half if h!=Fr(2*m+1,2)-k],half[:m],0,1))
    for j in range(m+1):
        if j<m:
            y=[1]*j+[0]*(m-1-j); B=[]
            for yi in y: B+=[yi,-yi]
            B+=[0,Fr(1,2)]; c=1
        else:
            B=[1,-1]*(m-1)+[1,Fr(-1,2)]; c=2
        out.append(('ii' if j<m else 'iii',j,B,[t for t in ints if t!=m-j],ints[:m],Fr(1,2),c))
    return out
N=int(sys.argv[1])
for m in [int(a) for a in sys.argv[2:]]:
    for name,k,B,tau,taus,z,c in fam(m):
        den,mn,lhs=series([Fr(b) for b in B],N)
        assert den==2 and mn==0, (den,mn)
        rr=rhs(m,tau,taus,z,c,N)
        ok=all(Fr(a)==b for a,b in zip(lhs,rr))
        print(f"rank {2*m} sector {name:3} index {k}: identity holds to q^{N//2}: {ok}   leading coeff {lhs[0]}",flush=True)
