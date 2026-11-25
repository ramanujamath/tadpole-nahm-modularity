"""Exact implementation of the cusp-coefficient procedure (Section 2.1) for tadpole matrices.

Field: K = Q(x), x = Q_r = 1 - cos(pi/(r+2)).  All arithmetic is exact (Fractions).
Output: rational component polynomials of c1-gamma0, kappa2, kappa3 in the shift B.
"""
from fractions import Fraction as Fr
import math, itertools, sys
from collections import defaultdict

# ---------- univariate polynomials over Q (list of Fr, low degree first) ----------
def ptrim(a):
    a=list(a)
    while a and a[-1]==0: a.pop()
    return a
def padd(a,b):
    n=max(len(a),len(b)); return ptrim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(n)])
def psub(a,b): return padd(a,[-x for x in b])
def pmul(a,b):
    if not a or not b: return []
    r=[Fr(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b): r[i+j]+=x*y
    return ptrim(r)
def pdivmod(a,b):
    a=ptrim(a); b=ptrim(b); q=[Fr(0)]*max(len(a)-len(b)+1,1)
    while len(a)>=len(b) and a:
        c=a[-1]/b[-1]; k=len(a)-len(b); q[k]=c
        a=psub(a,[Fr(0)]*k+[c*x for x in b])
    return ptrim(q),a

def minpoly_x(n):
    """Minimal polynomial of x = 1 - cos(pi/n) over Q (monic)."""
    import decimal
    ks=[k for k in range(1,n) if math.gcd(k,2*n)==1]
    # integer minimal polynomial of Y=2cos(pi/n): product (Y - 2cos(k pi/n))
    coeffs=[1.0]
    for k in ks:
        c=2*math.cos(k*math.pi/n)
        new=[0.0]*(len(coeffs)+1)
        for i,a in enumerate(coeffs):
            new[i+1]+=a; new[i]-=c*a
        coeffs=new
    Y=[Fr(round(c)) for c in coeffs]  # low degree first, integer
    assert all(abs(c-round(c))<1e-6 for c in coeffs)
    # substitute Y = 2 - 2x
    res=[]; pw=[Fr(1)]
    for i,c in enumerate(Y):
        res=padd(res,[c*v for v in pw]); pw=pmul(pw,[Fr(2),Fr(-2)])
    lead=res[-1]; return [c/lead for c in res]

class Field:
    def __init__(self,f,root):
        self.f=f; self.d=len(f)-1; self.root=root
    def el(self,a):
        a=ptrim(a)
        if len(a)>self.d: a=pdivmod(a,self.f)[1]
        return tuple(a+[Fr(0)]*(self.d-len(a)))
    def const(self,c): return self.el([Fr(c)])
    def gen(self): return self.el([Fr(0),Fr(1)])
    def add(self,a,b): return tuple(x+y for x,y in zip(a,b))
    def sub(self,a,b): return tuple(x-y for x,y in zip(a,b))
    def neg(self,a): return tuple(-x for x in a)
    def mul(self,a,b): return self.el(pmul(list(a),list(b)))
    def scal(self,c,a): return tuple(c*x for x in a)
    def inv(self,a):
        # extended Euclid
        r0,r1=list(self.f),ptrim(list(a)); s0,s1=[],[Fr(1)]
        while r1:
            q,r=pdivmod(r0,r1)
            r0,r1=r1,r; s0,s1=s1,psub(s0,pmul(q,s1))
        assert len(r0)==1
        return self.el([c/r0[0] for c in s0])
    def div(self,a,b): return self.mul(a,self.inv(b))
    def iszero(self,a): return all(x==0 for x in a)
    def num(self,a): return sum(float(c)*self.root**i for i,c in enumerate(a))

def tadpole(r):
    A=[[0]*r for _ in range(r)]
    for i in range(r):
        A[i][i]=2 if i<r-1 else 1
        if i+1<r: A[i][i+1]=A[i+1][i]=-1
    return A

def setup(r):
    n=r+2
    f=minpoly_x(n); root=1-math.cos(math.pi/n)
    K=Field(f,root)
    x=K.gen(); one=K.const(1)
    Q=[None]*(r+1)  # 1-indexed
    Q[r]=x
    Q[r-1]=K.div(x,K.sub(one,x)) if r>=2 else None
    for i in range(r-1,1,-1):
        Q[i-1]=K.div(K.mul(Q[i],Q[i]),K.mul(K.sub(one,Q[i]),Q[i+1]))
    Q=Q[1:]
    A=tadpole(r)
    # verify saddle equations exactly: (1-Q_i) * prod Q_j^{-A_ij} = 1
    for i in range(r):
        lhs=K.sub(one,Q[i]); rhs=one
        for j in range(r):
            e=A[i][j]
            for _ in range(abs(e)):
                rhs=K.mul(rhs,Q[j]) if e>0 else K.div(rhs,Q[j])
        assert K.iszero(K.sub(lhs,rhs)), 'saddle equation fails'
    for q in Q: assert 0<K.num(q)<1
    return K,A,Q

def bernoulli(n):
    B=[Fr(1)]
    for m in range(1,n+1):
        B.append(-sum(math.comb(m+1,k)*B[k] for k in range(m))/Fr(m+1))
    return B
def Ppoly(n):
    P={1:Fr(1)}
    for _ in range(n-1):
        d={e-1:c*e for e,c in P.items() if e>0}
        new=defaultdict(Fr)
        for e,c in d.items(): new[e+1]+=c; new[e+2]+=c
        P=dict(new)
    return P

def cusp_polys(r,N=3,verbose=False):
    K,A,Q=setup(r)
    one=K.const(1); zero=K.const(0)
    xi=[K.div(q,K.sub(one,q)) for q in Q]
    # M = (A + diag xi)^{-1} over K by Gauss-Jordan
    Mat=[[K.add(K.const(A[i][j]), xi[i] if i==j else zero) for j in range(r)]+[one if i==j else zero for j in range(r)] for i in range(r)]
    for c in range(r):
        p=next(i for i in range(c,r) if not K.iszero(Mat[i][c])); Mat[c],Mat[p]=Mat[p],Mat[c]
        iv=K.inv(Mat[c][c]); Mat[c]=[K.mul(iv,v) for v in Mat[c]]
        for i in range(r):
            if i!=c and not K.iszero(Mat[i][c]):
                fct=Mat[i][c]; Mat[i]=[K.sub(a,K.mul(fct,b)) for a,b in zip(Mat[i],Mat[c])]
    M=[row[r:] for row in Mat]
    Bn=bernoulli(2*N+6)
    # polynomials in (X, b): key = (Xexp tuple, bexp tuple) -> K element
    def P_add(a,b,c=None):
        res=dict(a)
        for k,v in b.items():
            v2=v if c is None else K.mul(c,v)
            res[k]=K.add(res[k],v2) if k in res else v2
            if K.iszero(res[k]): del res[k]
        return res
    def P_mul(a,b):
        res={}
        for (x1,b1),v1 in a.items():
            for (x2,b2),v2 in b.items():
                k=(tuple(i+j for i,j in zip(x1,x2)),tuple(i+j for i,j in zip(b1,b2)))
                v=K.mul(v1,v2)
                res[k]=K.add(res[k],v) if k in res else v
        return {k:v for k,v in res.items() if not K.iszero(v)}
    z0=tuple([0]*r)
    def ux(i,p):
        t=[0]*r; t[i]=p; return tuple(t)
    H={}
    for k in range(1,2*N+1):
        h={}
        if k==1:
            for i in range(r):
                h=P_add(h,{(ux(i,1),ux(i,1)):one})
                h=P_add(h,{(ux(i,1),z0):K.scal(Fr(1,2),xi[i])})
        for i in range(r):
            for g in range(3,2*N+8):
                j=k+2-g
                if 0<=j<=g and Bn[j]!=0:
                    Pp=Ppoly(g-1)
                    val=zero; xp=one
                    for e in range(0,max(Pp)+1):
                        if e in Pp: val=K.add(val,K.scal(Pp[e],xp))
                        xp=K.mul(xp,xi[i])
                    c=-Bn[j]/(math.factorial(j)*math.factorial(g-j))
                    h=P_add(h,{(ux(i,g-j),z0):K.scal(c,val)})
        H[k]=h
    Dl={0:{(z0,z0):one}}
    for l in range(1,2*N+1):
        acc={}
        for k in range(1,l+1):
            acc=P_add(acc,P_mul(H[k],Dl[l-k]),K.const(Fr(k,l)))
        Dl[l]=acc
        if verbose: print('D',l,len(acc),file=sys.stderr)
    memo={}
    def mom(a):
        s=sum(a)
        if s%2: return zero
        if s==0: return one
        if a in memo: return memo[a]
        i=next(t for t in range(r) if a[t]>0)
        tot=zero
        for j in range(r):
            pre=a[j]-(1 if i==j else 0)
            if pre:
                b=list(a); b[i]-=1; b[j]-=1
                tot=K.add(tot,K.scal(Fr(pre),K.mul(M[i][j],mom(tuple(b)))))
        memo[a]=tot; return tot
    def expect(poly):
        res={}
        for (xe,be),v in poly.items():
            m=mom(xe)
            if not K.iszero(m):
                t=K.mul(v,m)
                res[be]=K.add(res[be],t) if be in res else t
        return {k:v for k,v in res.items() if not K.iszero(v)}
    c=[expect(Dl[2*j]) for j in range(1,N+1)]
    # polynomial ops in b over K
    def badd(a,b,s=Fr(1)):
        res=dict(a)
        for k,v in b.items():
            v=K.scal(s,v); res[k]=K.add(res[k],v) if k in res else v
        return {k:v for k,v in res.items() if not K.iszero(v)}
    def bmul(a,b):
        res={}
        for k1,v1 in a.items():
            for k2,v2 in b.items():
                k=tuple(i+j for i,j in zip(k1,k2)); v=K.mul(v1,v2)
                res[k]=K.add(res[k],v) if k in res else v
        return {k:v for k,v in res.items() if not K.iszero(v)}
    gamma0=K.add(K.const(Fr(r,24)),K.scal(Fr(1,12),xi[0]))
    for v in xi[1:]: gamma0=K.add(gamma0,K.scal(Fr(1,12),v))
    c1,c2,c3=c[:3]
    k1=badd(c1,{z0:gamma0},Fr(-1))
    k2=badd(c2,bmul(c1,c1),Fr(-1,2))
    k3=badd(badd(c3,bmul(c1,c2),Fr(-1)),bmul(c1,bmul(c1,c1)),Fr(1,3))
    if N>=4:
        c4=c[3]
        k4=badd(badd(badd(badd(c4,bmul(c1,c3),Fr(-1)),bmul(c2,c2),Fr(-1,2)),bmul(bmul(c1,c1),c2)),bmul(bmul(c1,c1),bmul(c1,c1)),Fr(-1,4))
    else: k4=None
    def comps(p):
        out=[]
        for t in range(K.d):
            out.append({k:v[t] for k,v in p.items() if v[t]!=0})
        return out
    res={'K':K,'k1':comps(k1),'k2':comps(k2),'k3':comps(k3)}
    if k4 is not None: res['k4']=comps(k4)
    return res

def num_eval(p,B):
    return sum(float(c)*math.prod(b**e for b,e in zip(B,k)) for k,c in p.items())
