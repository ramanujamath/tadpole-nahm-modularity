import sys, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from evenw import rhs, nahm
for q in (0.2,0.3):
  for m in (1,2,3):
    ints=list(range(m+1)); B=[1,-1]*(m-1)+[1,-0.5]
    l=nahm(B,q)
    for z in (0.5,-0.5):
        rho=[p for p in ints if p!=0]
        rr=rhs(m,rho,ints[:m],z,q,1)
        print(f"q={q} m={m} (iii) z={z}: ratio LHS/RHS = {l/rr:.12f}")
