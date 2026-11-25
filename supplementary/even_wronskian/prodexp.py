# exact q-series of F_{r,B,0} (t = q^{1/2} if needed) and its Euler-product exponents
import sys
from fractions import Fraction as Fr
def series(B,NT,MARG=None):
    """Exact coefficients of F_{r,B,0} in t=q^(1/den), from the lowest power, NT+1 terms.
    MARG: truncation margin (t-units) to absorb later negative contributions."""
    import math
    r=len(B); den=2
    for b in B: den=den*Fr(b).denominator//math.gcd(den,Fr(b).denominator)
    MARG=MARG if MARG is not None else 2*NT
    TOP=NT+MARG
    def qinv(n):
        cur=[1]+[0]*TOP
        for j in range(1,n+1):
            st=j*den
            for i in range(st,TOP+1): cur[i]+=cur[i-st]
        return cur
    QI={}
    st={0:{0:1}}
    for i in range(r):
        new={}
        for npv,ser in st.items():
            n=0
            while True:
                a=n-npv; e=Fr(a*a,2)+Fr(B[i])*n; ei=int(e*den)
                if n>npv and ei>TOP+MARG: break
                if n>4*TOP: break
                if ei<=TOP:
                    tgt=new.setdefault(n,{})
                    for k,v in ser.items():
                        kk=k+ei
                        if -MARG<=kk<=TOP: tgt[kk]=tgt.get(kk,0)+v
                n+=1
        st={}
        for n,ser in new.items():
            if n not in QI: QI[n]=qinv(n)
            q=QI[n]; out={}
            for k,v in ser.items():
                for j in range(0,min(TOP-k,TOP)+1):
                    c=q[j]
                    if c: out[k+j]=out.get(k+j,0)+v*c
            st[n]=out
    tot={}
    for ser in st.values():
        for k,v in ser.items(): tot[k]=tot.get(k,0)+v
    mn=min(k for k,v in tot.items() if v)
    return den,mn,[tot.get(mn+i,0) for i in range(NT+1)]
def euler_exponents(c,M):
    # c[0]!=0 ; normalise and write c/c0 = prod_{n>=1} (1-x^n)^{-a_n}
    c0=c[0]; f=[Fr(x,c0) for x in c[:M+1]]
    # log derivative: x f'/f = sum n a_n x^n/(1-x^n)
    # compute g = x f'/f coefficients
    g=[Fr(0)]*(M+1); inv=[Fr(0)]*(M+1); inv[0]=Fr(1)
    for n in range(1,M+1): inv[n]=-sum(f[k]*inv[n-k] for k in range(1,n+1))
    for n in range(1,M+1): g[n]=sum(k*f[k]*inv[n-k] for k in range(1,n+1))
    # g_n = sum_{d|n} d a_d
    a=[Fr(0)]*(M+1)
    for n in range(1,M+1):
        a[n]=(g[n]-sum(d*a[d] for d in range(1,n) if n%d==0))/n
    return a
if __name__=='__main__':
    B=[Fr(x) for x in sys.argv[1].split(',')]
    den,mn,c=series(B,int(sys.argv[2]))
    a=euler_exponents(c,len(c)-1)
    print('B=',sys.argv[1],'t=q^(1/%d)'%den,'min exp',Fr(mn,den),'c0=',c[0],' first coeffs',c[:12])
    print('  product exponents a_n (in t):',[str(x) for x in a[1:61]])
