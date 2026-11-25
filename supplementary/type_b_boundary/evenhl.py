"""Experimental search for even-rank boundary determinant formulas.

E(v,c) = F_{2m,(2v,-2v,...,2v,c),0}(q).
Candidate: E(v,c) = K * Pref(t,u) * W(Phi_tau)(u) / W(delta_tau*)(u),
where W(.) is the Wronskian in the operator x d/dx of m functions evaluated at x=u=q^v
(this is the confluent limit of det(f_j(x_i))/det(g_j(x_i)) at x_i -> u).
"""
import math, itertools, sys
import numpy as np

def qp(n,q):
    p=1.0
    for j in range(1,n+1): p*=(1-q**j)
    return p
def poch(a,q,T=3000):
    p=1.0
    for s in range(T): p*=(1-a*q**s)
    return p
def nahm(B,q,N=120):
    r=len(B)
    f=np.array([q**(n*n/2+B[0]*n)/qp(n,q) for n in range(N)])
    n=np.arange(N)
    for i in range(1,r):
        K=q**((n[:,None]-n[None,:])**2/2)   # K[n,k]
        f=(K@f)*np.array([q**(B[i]*k)/qp(k,q) for k in range(N)])
    return f.sum()

def phi_derivs(a,u,q,L,m,sign=-1,kind='cos',R=8):
    """(x d/dx)^d Phi_a(x) at x=u for d=0..m-1, Phi_a = sum_l s^l q^{L l^2/2} (q^{a l} x^{L l + a} +/- q^{-a l} x^{L l - a})."""
    out=[]
    for d in range(m):
        tot=0.0
        for l in range(-R,R+1):
            base=(sign**l)*q**(L*l*l/2)
            e1=L*l+a; e2=L*l-a
            t1=q**(a*l)*u**e1*e1**d
            t2=q**(-a*l)*u**e2*e2**d
            tot+=base*(t1+t2 if kind=='cos' else t1-t2)
        out.append(tot)
    return out
def mono_derivs(a,u,m,kind='cos'):
    out=[]
    for d in range(m):
        out.append(a**d*u**a+((-a)**d*u**(-a) if kind=='cos' else -((-a)**d)*u**(-a)))
    return out
def wron(cols): return np.linalg.det(np.array(cols).T)

def E(m,v,c,q): return nahm([2*v if i%2==0 else -2*v for i in range(2*m-1)]+[c],q)

def candidate(m,v,c,q,tau,taus,kind='cos',sign=-1,pref='A'):
    L=2*m+2; u=q**v; w=q**0.5; t=q**(c+v)
    Wn=wron([phi_derivs(a,u,q,L,m,sign,kind) for a in tau])
    Wd=wron([mono_derivs(a,u,m,kind) for a in taus])
    # type-D affine denominator at x_i = u (all equal)
    D=poch(q,q)**m*(poch(q*u*u,q)*poch(q,q)**2*poch(q/(u*u),q))**(m*(m-1)//2)
    num=(poch(-w*t*u,q)*poch(-w*t/u,q))**m
    return num/D*Wn/Wd
