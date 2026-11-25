import sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from cusp import cusp_coeffs, Tmat
m=5; r=11; A=Tmat(r)
for x in ([1,0,0,0,0],[0,1,0,0,0],[1,0,1,0,1],[0,0,0,1,0],[1,1,0,1,1]):
    for z in (0,0.5):
        B=[]
        for xi in x: B+=[2*z-xi, xi-2*z]
        B.append(z)
        a,k1,k2,k3=cusp_coeffs(A,B)
        print('x=',x,'z=',z,' kappa2=%.3e kappa3=%.3e'%(k2,k3),flush=True)
