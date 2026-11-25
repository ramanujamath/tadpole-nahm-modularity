# Numerical test of the cusp-expansion procedure of Section 2.1 against the Nahm sum itself.
import math, itertools, sys
from fractions import Fraction
from decimal import Decimal as D, getcontext
import numpy as np
getcontext().prec=50
from collections import defaultdict

def Tmat(r):
    A=[[0]*r for _ in range(r)]
    for i in range(r):
        A[i][i]=2 if i<r-1 else 1
        if i+1<r: A[i][i+1]=A[i+1][i]=-1
    return A
def saddle(A):
    r=len(A); Q=[0.3]*r
    for _ in range(500):
        Fv=np.array([1-Q[i]-math.prod(Q[j]**A[i][j] for j in range(r)) for i in range(r)])
        J=np.zeros((r,r))
        for i in range(r):
            p=math.prod(Q[j]**A[i][j] for j in range(r))
            for j in range(r):
                J[i,j]=(-1 if i==j else 0)-p*A[i][j]/Q[j]
        step=np.linalg.solve(J,Fv); lam=1.0
        while True:
            Qn=list(np.array(Q)-lam*step)
            if all(0<x<1 for x in Qn): break
            lam/=2
        Q=Qn
        if max(abs(Fv))<1e-15: break
    return Q
def Li2(x): return sum(x**n/n**2 for n in range(1,400))
def Lrog(x): return Li2(x)+0.5*math.log(x)*math.log(1-x)

# polynomials in X: dict tuple->float
def padd(a,b,c=1.0):
    r=defaultdict(float,a)
    for k,v in b.items(): r[k]+=c*v
    return dict(r)
def pmul(a,b):
    r=defaultdict(float)
    for k1,v1 in a.items():
        for k2,v2 in b.items(): r[tuple(x+y for x,y in zip(k1,k2))]+=v1*v2
    return dict(r)
def bern(n):
    Bn=[Fraction(1)]
    for m in range(1,n+1):
        Bn.append(-sum(math.comb(m+1,k)*Bn[k] for k in range(m))/Fraction(m+1))
    return Bn  # B1=-1/2 convention
def Ppoly(n):  # coefficients of P_n(x)
    P={1:Fraction(1)}
    for _ in range(n-1):
        d={e-1:c*e for e,c in P.items() if e>0}
        new=defaultdict(Fraction)
        for e,c in d.items(): new[e+1]+=c; new[e+2]+=c
        P=dict(new)
    return P
def evalP(P,x): return sum(float(c)*x**e for e,c in P.items())

def cusp_coeffs(A,B,N=3):
    r=len(A); Q=saddle(A); xi=[q/(1-q) for q in Q]
    M=np.linalg.inv(np.array(A,float)+np.diag(xi))
    Bn=bern(2*N+4)
    def unit(i,p):
        t=[0]*r; t[i]=p; return tuple(t)
    H={}
    for k in range(1,2*N+1):
        h={}
        if k==1:
            for i in range(r): h=padd(h,{unit(i,1):B[i]+xi[i]/2})
        for i in range(r):
            for g in range(3,2*N+6):
                j=k+2-g
                if 0<=j<=g:
                    c=float(Bn[j])/(math.factorial(j)*math.factorial(g-j))*evalP(Ppoly(g-1),xi[i])
                    h=padd(h,{unit(i,g-j):-c})
        H[k]=h
    Dl={0:{tuple([0]*r):1.0}}
    for l in range(1,2*N+1):
        acc={}
        for k in range(1,l+1): acc=padd(acc,pmul(H[k],Dl[l-k]),k/l)
        Dl[l]=acc
    memo={}
    def mom(a):
        if sum(a)%2: return 0.0
        if sum(a)==0: return 1.0
        if a in memo: return memo[a]
        i=next(t for t in range(r) if a[t]>0)
        s=0.0
        for j in range(r):
            pre=a[j]-(1 if i==j else 0)
            if pre:
                b=list(a); b[i]-=1; b[j]-=1; s+=pre*M[i,j]*mom(tuple(b))
        memo[a]=s; return s
    c=[sum(v*mom(k) for k,v in Dl[2*j].items()) for j in range(1,N+1)]
    gamma0=r/24+sum(xi)/12
    alpha=sum(math.pi**2/6-Lrog(q) for q in Q)
    k1=c[0]-gamma0; k2=c[1]-c[0]**2/2; k3=c[2]-c[0]*c[1]+c[0]**3/3
    return alpha,k1,k2,k3

def logF(A,B,eps,N):
    # T_r tadpole only: exponent sum (n_i-n_{i-1})^2/2 + B.n ; Decimal DP
    q=D(-eps).exp() if False else (D(-1)*D(repr(eps))).exp()
    qp=[D(1)]
    for n in range(1,N+1): qp.append(qp[-1]*(1-q**n))
    r=len(A)
    sq=[q**(D(d*d)/2) for d in range(N+1)]
    f=[q**(D(n*n)/2+D(repr(B[0]))*n)/qp[n] for n in range(N+1)]
    for i in range(1,r):
        bi=D(repr(B[i]))
        f=[sum(f[k]*sq[abs(n-k)] for k in range(N+1))*q**(bi*n)/qp[n] for n in range(N+1)]
    return float(sum(f).ln())

if __name__=='__main__':
    r=int(sys.argv[1]); B=[float(x) for x in sys.argv[2].split(',')]
    A=Tmat(r)
    alpha,k1,k2,k3=cusp_coeffs(A,B)
    print('predicted: coeff eps^1 =',-0+k1,' (C=0)  kappa2=',k2,' kappa3=',k3)
    eps=[0.05,0.06,0.07,0.08,0.09,0.10,0.12,0.14,0.16,0.18,0.20,0.22]
    g=[]
    for e in eps:
        Nn=int(12/e)+20
        g.append(logF(A,B,e,Nn)-alpha/e)
    X=np.vander(np.array(eps),7,increasing=True)
    coef=np.linalg.lstsq(X,np.array(g),rcond=None)[0]
    print('fitted   : eps^1 =',coef[1],' eps^2 =',coef[2],' eps^3 =',coef[3])
