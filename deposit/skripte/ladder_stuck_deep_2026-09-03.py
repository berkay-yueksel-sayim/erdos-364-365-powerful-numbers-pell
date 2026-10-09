# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
#!/usr/bin/env python3
# ladder_stuck_deep_2026-09-03.py
# Deeper witness pass (primes 20000 <= p < 800000) for the 14 stuck ladder positions (m, k, stage).
# Per position: finds the fundamental solution T_1 + U_1*sqrt(m) of x^2 - m*y^2 = 1 by continued fractions, computes T_k mod
#   p^2
# for each prime p not dividing m, and records p as a witness if p | T_k but p^2 does not divide T_k (v_p(T_k) = 1).
# Stops after 2 witnesses. With witnesses it prints the next ladder index k' = k * (product of the witnesses) and the
#   approximate
# digit count of T_k'; without any it reports that no witness exists below 800000.
# Reads: nothing (the positions are hard-coded in `STUCK`). Writes: nothing (stdout only).
# Usage: python ladder_stuck_deep_2026-09-03.py   (no arguments). Controls: none.
import math, sys, io
from sympy import primerange
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding='utf-8',errors='replace')
# `STUCK` = stuck positions as tuples (kernel m, index k, ladder stage `Stufe`).
STUCK=[(159,17260434899,2),(247,240196912525063,2),(271,122537257,2),(719,3219255554629,2),(767,40651,2),(815,2388074418715,2),(951,951,1),(967,1042560160613,2),(1031,12756563,2),(1111,31825388597299,2),(1327,894720461,2),(1423,866821873,2),(1535,4150333,2),(1671,557,1)]
# `fund_pos` = fundamental solution: smallest positive (p, q) with p^2 - D*q^2 = 1,
# found among the continued-fraction convergents of sqrt(D).
def fund_pos(D):
    a0=math.isqrt(D); m,d,a=0,1,a0; pp,p=1,a0; qq,q=0,1
    if p*p-D*q*q==1: return p,q
    while True:
        m=d*a-m; d=(D-m*m)//d; a=(a0+m)//d
        pp,p=p,a*p+pp; qq,q=q,a*q+qq
        if p*p-D*q*q==1: return p,q
# `T_mod` = T_k mod M, the rational part (coefficient of 1) of (T1 + U1*sqrt(m))^k, computed in Z[sqrt(m)]/M by binary
#   exponentiation.
def T_mod(T1,U1,m,k,M):
    ra,rb=1%M,0; ba,bb=T1%M,U1%M
    while k:
        if k&1: ra,rb=(ra*ba+m*rb*bb)%M,(ra*bb+rb*ba)%M
        ba,bb=(ba*ba+m*bb*bb)%M,(2*ba*bb)%M; k>>=1
    return ra
# Candidate witness primes.
PR=list(primerange(20000,800000))
res=[]
for m,k,stage in STUCK:
    # `le` = log10(2*T1): T_k has about k*le decimal digits. `w` = witnesses found for this position.
    T1,U1=fund_pos(m); le=math.log10(2*T1); w=[]
    for p in PR:
        if m%p==0: continue
        a=T_mod(T1,U1,m,k,p*p)
        if a%p==0 and a!=0:
            w.append(p)
            if len(w)>=2: break
    # Witnesses found: next ladder index `kn` = k * product of the witnesses (printed as `GELOEST` = solved).
    if w:
        kn=k
        for p in w: kn*=p
        res.append((m,k,stage,w,kn*le)); print(f"m={m:5} Stufe {stage} k={k}: Zeugen {w} -> naechstes k = {kn}  -> T hat ~10^{math.log10(kn*le):.1f} Stellen   GELOEST")
    # No witness below 800000: position stays open (printed as `OFFEN`).
    else:
        res.append((m,k,stage,[],None)); print(f"m={m:5} Stufe {stage} k={k}: KEIN Zeuge < 800000  (Stellenzahl von T_k: ~{k*le:.3g})   OFFEN")
print(f"\ngeloest: {sum(1 for r in res if r[3])} von {len(res)}")
