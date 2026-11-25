import math
def T(r):
    A=[[0]*r for _ in range(r)]
    for i in range(r):
        A[i][i]=2 if i<r-1 else 1
        if i+1<r: A[i][i+1]=A[i+1][i]=-1
    return A
def inv(M):
    n=len(M); A=[row[:]+[1.0 if i==j else 0.0 for j in range(n)] for i,row in enumerate(M)]
    for c in range(n):
        p=max(range(c,n),key=lambda r:abs(A[r][c])); A[c],A[p]=A[p],A[c]
        pv=A[c][c]; A[c]=[x/pv for x in A[c]]
        for r in range(n):
            if r!=c:
                f=A[r][c]; A[r]=[x-f*y for x,y in zip(A[r],A[c])]
    return [row[n:] for row in A]
def check(r,Q):
    A=T(r)
    res=[1-Q[i]-math.prod(Q[j]**A[i][j] for j in range(r)) for i in range(r)]
    return max(abs(x) for x in res)
s=math.sqrt(5)
Q3=[(3-s)/2,s-2,(3-s)/4]; print('r3 saddle resid',check(3,Q3))
xi=[q/(1-q) for q in Q3]; print('xi',xi,[(s-1)/2,(s-1)/4,s-2])
A=T(3); M=inv([[A[i][j]+(xi[i] if i==j else 0) for j in range(3)] for i in range(3)])
Mc=[[-1.5+0.9*s,-1+0.6*s,0.5-s/10],[-1+0.6*s,0.4*s,0.5+s/10],[0.5-s/10,0.5+s/10,0.5+0.4*s]]
print('M diff',max(abs(M[i][j]-Mc[i][j]) for i in range(3) for j in range(3)))
s3=math.sqrt(3); print('r4',check(4,[(s3-1)/2,(3-s3)/6,(2*s3-3)/3,(2-s3)/2]))
t=1-math.cos(math.pi/7); print('r5 poly',8*t**3-20*t**2+12*t-1, check(5,[-4*t*t+4*t,2*t,12*t*t-20*t+2,-8*t*t+12*t-1,t]))
t=1-math.cos(math.pi/9); print('r7 poly',8*t**3-24*t**2+18*t-1,t)
t=1-math.cos(math.pi/11); print('r9 poly',32*t**5-144*t**4+224*t**3-140*t**2+30*t-1,t)
