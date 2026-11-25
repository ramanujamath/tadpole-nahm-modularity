from collections import defaultdict
V=['a','b','c','d','e','w','u']
def mono(**kw):
    t=[0]*7
    for k,v in kw.items(): t[V.index(k)]=v
    return {tuple(t):1}
def add(*ps):
    r=defaultdict(int)
    for p in ps:
        for k,v in p.items(): r[k]+=v
    return {k:v for k,v in r.items() if v}
def mul(*ps):
    r={(0,)*7:1}
    for p in ps:
        s=defaultdict(int)
        for k1,v1 in r.items():
            for k2,v2 in p.items():
                s[tuple(x+y for x,y in zip(k1,k2))]+=v1*v2
        r={k:v for k,v in s.items() if v}
    return r
def sc(c,p): return {k:c*v for k,v in p.items()}
def neg(p): return sc(-1,p)
def sig(var,p):
    i=V.index(var); r={}
    for k,v in p.items():
        t=list(k); t[5]+=2*k[i]; r[tuple(t)]=v
    return r
one=mono()
a,b,c,d,e,w,u=[mono(**{x:1}) for x in V]
def pw(x,n): return mul(*([x]*n)) if n>0 else one
ui=mono(u=-1); wi=mono(w=-1)
K=add(one,mul(w,ui)); R=add(mul(K,e),neg(one)); M=mul(a,b,c,d,e)
def A2(P): return add(mul(a,c,add(one,neg(b)),P), neg(mul(pw(u,2),pw(b,2),sig('b',P))))
def A3(P): return add(mul(b,d,add(one,neg(c)),P), neg(mul(mono(u=-2),pw(c,2),sig('c',P))))
def A4(P): return add(mul(c,e,add(one,neg(d)),P), neg(mul(pw(u,2),pw(d,2),sig('d',P))))
def A5(P): return add(mul(d,add(one,neg(e)),P), neg(mul(mono(w=-1,u=-1),e,sig('e',P))))
w2,w3,w4,w5=[pw(w,n) for n in (2,3,4,5)]
P2=add(mul(d,e,add(mul(K,pw(e,2)), neg(mul(add(mul(K,add(one,w4)),w2),e)), w2)),
       mul(w3,ui,pw(e,2),R), neg(mul(w3,u,pw(d,2),e,R)))
P3=mul(w3,u,a,d,e,R)
P4=mul(a,add(mul(w3,ui,b,e,R), neg(mul(w4,mono(u=-2),c,e)), neg(mul(w4,K,e)), mul(w3,u,K,d,add(one,neg(e)))))
P5=add(mul(a,c,e,add(mul(K,e), neg(mul(w2,add(w2,mul(w,u),one))), neg(mul(w3,ui,c)))),
       mul(pw(u,2),pw(b,2),e,add(mul(w3,u,K,d),w2,neg(mul(K,e)))),
       mul(w5,pw(u,3),K,a,pw(d,2)))
lhs=add(A2(P2),A3(P3),A4(P4),A5(P5))
rhs=add(mul(w4,pw(u,2),pw(b,2),pw(d,2),e),
        neg(mul(M, add(mul(K,add(one,mul(w3,ui)),pw(e,2)), neg(mul(add(one,w2,w4),K,e)), w2))))
diff=add(lhs,neg(rhs))
print('terms lhs',len(lhs),'identity (41) holds:',not diff)
if diff: print(list(diff.items())[:10])
