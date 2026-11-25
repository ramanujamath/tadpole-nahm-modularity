import sys
from fractions import Fraction as Fr
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from prodexp import series
N=int(sys.argv[1])
for m in range(1,int(sys.argv[2])+1):
    Biii=[1,-1]*(m-1)+[1,Fr(-1,2)]; Biv=list(Biii); Biv[0]+=1; Biv[1]-=1
    y=[1]*(m-1); Bii=[]
    for yi in y: Bii+=[yi,-yi]
    Bii+=[0,Fr(1,2)]
    d1,mn1,a=series(Biv,N); d2,mn2,b=series(Bii,N)
    print(f"rank {2*m}: F_(iv) starts at q^{Fr(mn1,d1)}, F_(ii,j=m-1) at q^{Fr(mn2,d2)}; F_(iv) = 2*q^-1*F_(ii) through q^{N//2}:", Fr(mn1,d1)==-1 and mn2==0 and all(x==2*y for x,y in zip(a,b)))
