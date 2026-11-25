"""Buchberger algorithm over F_p with graded reverse lexicographic order.
Used with Lemma 2.3 (finite specialization): we need
  d  = dim F_p[X]/(f_i mod p)            (affine ideal)
  and finiteness of F_p[X]/(top parts mod p).
"""
from itertools import combinations
import itertools

def key(e):  # grevlex key
    return (sum(e), tuple(-x for x in reversed(e)))
def lt(p): return max(p,key=key)
def divides(a,b): return all(x<=y for x,y in zip(a,b))
class GB:
    def __init__(self,p): self.p=p
    def norm(self,f):
        f={k:v%self.p for k,v in f.items() if v%self.p}
        if not f: return f
        c=pow(f[lt(f)],-1,self.p); return {k:v*c%self.p for k,v in f.items()}
    def reduce(self,f,G,LT):
        p=self.p; f=dict(f); r={}
        while f:
            m=lt(f); c=f[m]
            for g,lg in zip(G,LT):
                if divides(lg,m):
                    sh=tuple(y-x for x,y in zip(lg,m))
                    for k,v in g.items():
                        kk=tuple(a+b for a,b in zip(k,sh))
                        nv=(f.get(kk,0)-c*v)%p
                        if nv: f[kk]=nv
                        else: f.pop(kk,None)
                    break
            else:
                r[m]=c; del f[m]
        return r
    def spoly(self,f,g,lf,lg):
        p=self.p; L=tuple(max(a,b) for a,b in zip(lf,lg))
        s1=tuple(a-b for a,b in zip(L,lf)); s2=tuple(a-b for a,b in zip(L,lg))
        res={}
        for k,v in f.items():
            kk=tuple(a+b for a,b in zip(k,s1)); res[kk]=(res.get(kk,0)+v)%p
        for k,v in g.items():
            kk=tuple(a+b for a,b in zip(k,s2)); res[kk]=(res.get(kk,0)-v)%p
        return {k:v for k,v in res.items() if v}
    def groebner(self,F):
        G=[]; LT=[]
        for f in F:
            f=self.norm(f)
            if f:
                f=self.norm(self.reduce(f,G,LT)) if G else f
                if f: G.append(f); LT.append(lt(f))
        pairs=[(i,j) for i in range(len(G)) for j in range(i)]
        while pairs:
            # select pair with smallest lcm degree (normal strategy)
            pairs.sort(key=lambda ij: -sum(max(a,b) for a,b in zip(LT[ij[0]],LT[ij[1]])))
            i,j=pairs.pop()
            lf,lg=LT[i],LT[j]
            if all(a==0 or b==0 for a,b in zip(lf,lg)): continue
            h=self.reduce(self.spoly(G[i],G[j],lf,lg),G,LT)
            if h:
                h=self.norm(h); G.append(h); LT.append(lt(h)); n=len(G)-1
                pairs+=[(n,k) for k in range(n)]
        # minimal basis
        keep=[i for i in range(len(G)) if not any(j!=i and divides(LT[j],LT[i]) and (LT[j]!=LT[i] or j<i) for j in range(len(G)))]
        G=[G[i] for i in keep]; LT=[LT[i] for i in keep]
        G=[self.norm(self.reduce(g,[h for h in G if h is not g],[l for h,l in zip(G,LT) if h is not g])|{}) for g in G]
        # reduce tails properly
        out=[]
        for idx,g in enumerate(G):
            others=[h for k,h in enumerate(G) if k!=idx]; olt=[LT[k] for k in range(len(G)) if k!=idx]
            lg=LT[idx]; tail={k:v for k,v in g.items() if k!=lg}
            tail=self.reduce(tail,others,olt)
            out.append(self.norm({lg:1,**tail}))
        return out
def standard_monomials(G,n,cap=60):
    L=[lt(g) for g in G]
    out=[]; frontier=[tuple([0]*n)]; seen=set(frontier)
    while frontier:
        e=frontier.pop()
        if any(divides(l,e) for l in L): continue
        out.append(e)
        if len(out)>100000: return None
        if sum(e)>cap: return None
        for i in range(n):
            f=list(e); f[i]+=1; f=tuple(f)
            if f not in seen: seen.add(f); frontier.append(f)
    return out
def top_part(f):
    d=max(sum(k) for k in f); return {k:v for k,v in f.items() if sum(k)==d}
