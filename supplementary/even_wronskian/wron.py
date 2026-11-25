import itertools, math
def nahm(B,q,N=90):
    r=len(B); qp=[1.0]
    for n in range(1,N+1): qp.append(qp[-1]*(1-q**n))
    f=[q**(n*n/2+B[0]*n)/qp[n] for n in range(N+1)]
    for i in range(1,r):
        g=[]
        for n in range(N+1):
            s=sum(f[k]*q**((n-k)**2/2) for k in range(N+1))
            g.append(s*q**(B[i]*n)/qp[n])
        f=g
    return sum(f)
def poch(a,q,T=3000):
    p=1.0
    for s in range(T): p*=(1-a*q**s)
    return p
def chi(v):
    p=1.0
    for x in v: p*=x
    for i in range(len(v)):
        for j in range(i+1,len(v)): p*=(v[i]**2-v[j]**2)
    return p
def rhs(m,k,z,q,R=None):
    L=2*m+3; r=2*m+1; p=m*(2*m+1)
    if z==0:
        rho=[x for x in range(1,m+2) if x!=m+1-k]; rs=list(range(1,m+1))
    else:
        rho=[x+0.5 for x in range(0,m+1) if x!=k]; rs=[x+0.5 for x in range(0,m)]
    R = R or (4 if m<=4 else 2)
    nr=sum(x*x for x in rho); tot=0.0
    for n in itertools.product(range(-R,R+1),repeat=m):
        v=[rho[i]+L*n[i] for i in range(m)]
        tot+=chi(v)*q**((sum(x*x for x in v)-nr)/(2*L))
    return poch(-q**(0.5+z),q)**r/poch(q,q)**p*tot/chi(rs)
def Bvec(m,k,z):
    x=[0]*(m-k)+[1]*k; B=[]
    for xi in x: B+= [2*z-xi, xi-2*z]
    return B+[z]
if __name__=='__main__':
    for q in (0.2,0.35):
      for m in (1,2,3,4,8):
        for k in (range(m+1) if m<8 else (6,7,8)):
          for z in (0,0.5):
            l=nahm(Bvec(m,k,z),q); rr=rhs(m,k,z,q)
            print(f"q={q} m={m} k={k} z={z}: rel.err={abs(l-rr)/abs(l):.1e}")
