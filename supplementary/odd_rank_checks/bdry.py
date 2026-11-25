# polynomials in w (q=w^2), Laurent allowed: dict exp->coef
from collections import defaultdict
def add(*ps):
    r=defaultdict(int)
    for p in ps:
        for k,v in p.items(): r[k]+=v
    return {k:v for k,v in r.items() if v}
def mul(a,b):
    r=defaultdict(int)
    for i,x in a.items():
        for j,y in b.items(): r[i+j]+=x*y
    return {k:v for k,v in r.items() if v}
def sc(c,p): return {k:c*v for k,v in p.items() if c*v}
def W(e): return {e:1}
ONE={0:1}
def qpoch(n):  # (q;q)_n
    r=ONE
    for j in range(1,n+1): r=mul(r,add(ONE,{2*j:-1}))
    return r
def div_exact(a,b):  # polynomial division in w (nonneg exps)
    a=dict(a); q={}
    db=max(b); lb=b[db]
    while a:
        da=max(a)
        if da<db: raise Exception('not exact')
        c=a[da]//lb; assert c*lb==a[da]
        q[da-db]=c
        a=add(a,sc(-c,mul({da-db:1},b)))
    return q
def gb(n,k):
    if k<0 or k>n or n<0: return {}
    return div_exact(qpoch(n),mul(qpoch(k),qpoch(n-k)))
def A(k,j): return sc((-1)**j,mul(W(2*((k-j)*(k-j-1)//2)),gb(k+j,2*j)))
def Ud(d,t):
    M=[[{} for _ in range(d)] for _ in range(d)]
    for i in range(d):
        M[i][i]=add(ONE,mul(t,t)) if i<d-1 else ONE
        if i+1<d: M[i][i+1]=t; M[i+1][i]=t
    return M
def mm(X,Y):
    d=len(X); return [[add(*[mul(X[i][k],Y[k][j]) for k in range(d)]) for j in range(d)] for i in range(d)]
ok=True
for m in (1,2,3,4):
    d=m+1; V=[[ONE if i==j else {} for j in range(d)] for i in range(d)]; Vs=[V]
    for j in range(1,m+4):
        V=mm(V,Ud(d,W(2*(j-1)+1))); Vs.append(V)
    # Prop coeffs
    for j in range(0,m+1):
        for h in range(0,j+1):
            D=add(gb(2*j,j-h), sc(-1,mul(W(2*(2*h+1)),gb(2*j,j-h-1))))
            lhs=Vs[j][d-h-1][d-1]
            ok&= lhs==mul(W(h*h),D)
    # C_{j,l}
    def C(j,l):
        tot={}
        for h in range(l,m+1):
            t=Vs[j][d-h-1][d-1]
            assert all(e>=h*h for e in t) if t else True
            t={e-h*h:c for e,c in t.items()}
            tot=add(tot,sc((-1)**h,mul(t,A(h,l))))
        return tot
    for j in range(m+1):
        for l in range(m+1): ok&= C(j,l)==(ONE if j==l else {})
    for l in range(m+1):
        c=C(m+1,l); ok&= c==sc((-1)**m,A(m+1,l)); ok&= all(e%2==0 and e>=0 for e in c)
    for j in range(m+4):
        for l in range(m+1): ok&= all(e%2==0 and e>=0 for e in C(j,l))
print('Prop coeffs, C_{j,l} in Z[q], C_{jj}=delta, C_{m+1,l}=(-1)^m A:',ok)
