# Exact q-series check of identity (12), in t = q^(1/2), to order t^NT.
import sys
from fractions import Fraction as Fr
NT=int(sys.argv[1]) if len(sys.argv)>1 else 24
def smul(a,b):
    r=[0]*(NT+1)
    for i,x in enumerate(a):
        if x:
            for j in range(NT+1-i): r[i+j]+=x*b[j]
    return r
def inv(a):  # a[0]=1
    r=[0]*(NT+1); r[0]=1
    for n in range(1,NT+1): r[n]=-sum(a[k]*r[n-k] for k in range(1,n+1))
    return r
def poch_t(step0,stepq,count=None):  # prod_{s>=0} (1 + sign t^{step0+2s}) generic helper not used
    pass
# 1/(q;q)_n in t
qinv=[None]
cur=[1]+[0]*NT; qinvs=[cur[:]]
for n in range(1,NT+40):
    f=[0]*(NT+1); f[0]=1
    if 2*n<=NT: f[2*n]=-1
    cur=smul(cur,inv(f)); qinvs.append(cur[:])
def nahm(m,k,z):
    x=[0]*(m-k)+[1]*k; r=2*m+1; OFF=2*m+2; W=NT+OFF+2*m+2
    # states: n_prev -> series indexed shifted by OFF (allows negative exponents)
    st={0:[1 if i==OFF else 0 for i in range(W+1)]}
    for i in range(1,r+1):
        new={}
        for npv,ser in st.items():
            amin=-npv
            a=amin
            while True:
                # exponent in t: a^2 ; even step i=2j : -2*c*a  (c=2z-x_j)
                e=a*a
                if i%2==0 and i<r:
                    c=2*z-x[i//2-1]; e-= int(round(2*c*a))
                if i==r: e+= int(round(2*z*(npv+a)))
                if a>0 and a*a-2*abs(a) > W+NT: break
                n=npv+a
                if e<=W and n<len(qinvs):
                    tgt=new.setdefault(n,[0]*(W+1))
                    qi=qinvs[n]
                    for p,cf in enumerate(ser):
                        if cf and p+e<=W and p+e>=0:
                            # multiply by t^e and 1/(q)_n
                            base=p+e
                            for j in range(0,W+1-base):
                                if qi[j] if j<=NT else 0: tgt[base+j]+=cf*qi[j]
                        elif cf and p+e<0: raise Exception('neg')
                a+=1
        st=new
    tot=[0]*(W+1)
    for s in st.values():
        for i,c in enumerate(s): tot[i]+=c
    assert all(c==0 for c in tot[:OFF]), "negative powers"
    return tot[OFF:OFF+NT+1]
def chi(v):
    p=1
    for a in v: p*=a
    for i in range(len(v)):
        for j in range(i+1,len(v)): p*=(v[i]**2-v[j]**2)
    return p
def rhs(m,k,z):
    L=2*m+3; r=2*m+1; p=m*(2*m+1)
    if z==0:
        rho=[2*x for x in range(1,m+2) if x!=m+1-k]; rs=[2*x for x in range(1,m+1)]
    else:
        rho=[2*x+1 for x in range(0,m+1) if x!=k]; rs=[2*x+1 for x in range(0,m)]
    # doubled coordinates; exponent in t: (|v|^2-|rho|^2)/L with v doubled -> /(4L)*2
    th=[Fr(0)]*(NT+1); c0=chi(rs)
    def rec(i,n,expo):
        if i==m:
            v=[rho[j]+2*L*n[j] for j in range(m)]
            th[expo]+=Fr(chi(v),c0); return
        nn=0
        for s in (1,-1):
            nn=0 if s==1 else -1
            while True:
                e=rho[i]*nn+L*nn*nn   # t-exponent contribution: 2*(rho/2)*n + L n^2
                if e>NT- (0): 
                    if (s==1 and nn>0) or (s==-1): break
                if expo+e<=NT and expo+e>=0: rec(i+1,n+[nn],expo+e)
                elif expo+e>NT and ((s==1 and nn>0) or s==-1): break
                nn+=s
    rec(0,[],0)
    # products
    num=[1]+[0]*NT
    for s in range(NT+1):
        e=(1+2*z)+2*s  # (-q^{1/2+z} q^s) -> t^{1+2z+2s}
        e=int(round(e))
        if e>NT: break
        f=[0]*(NT+1); f[0]=1; f[e]=1
        for _ in range(r): num=smul(num,f)
    den=[1]+[0]*NT
    for s in range(1,NT+1):
        if 2*s>NT: break
        f=[0]*(NT+1); f[0]=1; f[2*s]=-1
        g=inv(f)
        for _ in range(p): den=smul(den,g)
    return smul(smul(num,den),th)
cases=[(int(a),int(b),Fr(c)) for a,b,c in (s.split(',') for s in sys.argv[2:])]
for m,k,z in cases:
    l=nahm(m,k,z); rr=rhs(m,k,z)
    ok=all(Fr(a)==b for a,b in zip(l,rr))
    print(f"m={m} k={k} z={z}: match to t^{NT}: {ok}; first coeffs {l[:8]}")
