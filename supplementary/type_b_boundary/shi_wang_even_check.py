import itertools, math
def poch(a,q,T=2000):
    p=1.0
    for s in range(T): p*=(1-a*q**s)
    return p
def qp(n,q):
    p=1.0
    for j in range(1,n+1): p*=(1-q**j)
    return p
def chi_lhs(xs,q,N=40):
    r=len(xs); tot=0.0
    # DP over coordinates with quadratic sum (n_i-n_{i-1})^2/2
    f={n: q**(n*n/2)*xs[0]**n/qp(n,q) for n in range(N)}
    for i in range(1,r):
        f={n: sum(v*q**((n-k)**2/2) for k,v in f.items())*xs[i]**n/qp(n,q) for n in range(N)}
    return sum(f.values())
def chi_rhs(xs,x0,q,N=22,L=12):
    R=len(xs)//2; x=[x0]+list(xs)   # x[0..2R]
    pref=poch(-q**0.5*x[2*R],q)/poch(q,q)**(R-1)
    tot=0.0
    for n in itertools.product(range(N),repeat=R+1):
        e=n[0]**2/2+(n[0]+n[1])**2/2+sum((n[i-1]-n[i])**2 for i in range(2,R+1))/4  # n_i indices 1..R+1 -> n[0..R]
        if e>60: continue
        t=q**e/math.prod(qp(k,q) for k in n)*x[0]**(-n[0])
        for j in range(1,R+1): t*=(x[2*j]/x[2*j-2])**n[j]
        for k in range(2,R+1):
            d=(n[k-1]-n[k])/2
            t*=sum(q**((l-d)**2)*x[2*k-2]**(-l) for l in range(-L,L+1))
        tot+=t
    return pref*tot
q=0.3
for xs,x0 in [([1.3,1.7,1/1.7,0.9],1/1.3),([0.8,1.2,1/1.2,1.1],1/0.8)]:
    print('rank4: lhs',chi_lhs(xs,q),' rhs',chi_rhs(xs,x0,q))
