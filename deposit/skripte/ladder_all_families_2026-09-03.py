# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# ladder_all_families_2026-09-03.py
# Purpose: for every squarefree kernel m ≡ 7 (mod 8) with m <= MMAX, build the witness ladder of the Lucas sequence T_k and
#   print a lower bound for the digit count of the first possible triple middle with kernel m.
# Method: base k0 = m' (product of the primes p | m with p not dividing U1); the witnesses of T_k0 are the primes p with
#   p || T_k0 (exponent exactly 1, found from T_k0 mod p²); k1 = k0 * (product of the witnesses); the witnesses of T_k1
#   give k2. The bound holds rigorously under the lifting-the-exponent statement for Lucas (V) sequences (checked numerically
#   elsewhere in the series).
# Reads: nothing. Writes: nothing (standard output only).
# Usage: python ladder_all_families_2026-09-03.py   (no arguments; MMAX = 2000)
# Controls: none inside this script.
import math, sys, io
from sympy import factorint, primerange
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding='utf-8',errors='replace')
MMAX=2000; PR=list(primerange(3,20000)); NW=3   # `PR` = primes searched for witnesses; `NW` = max. witnesses kept per step
def fund_pos(D):   # fundamental solution (p, q) of p^2 - D q^2 = 1 by the continued fraction of sqrt(D)
    a0=math.isqrt(D); m,d,a=0,1,a0; pp,p=1,a0; qq,q=0,1
    if p*p-D*q*q==1: return p,q
    while True:
        m=d*a-m; d=(D-m*m)//d; a=(a0+m)//d
        pp,p=p,a*p+pp; qq,q=q,a*q+qq
        if p*p-D*q*q==1: return p,q
def T_mod(T1,U1,m,k,M):   # T_k modulo M, by binary exponentiation of T1 + U1*sqrt(m)
    ra,rb=1%M,0; ba,bb=T1%M,U1%M
    while k:
        if k&1: ra,rb=(ra*ba+m*rb*bb)%M,(ra*bb+rb*ba)%M
        ba,bb=(ba*ba+m*bb*bb)%M,(2*ba*bb)%M; k>>=1
    return ra
def witnesses(T1,U1,m,k,nmax):   # first `nmax` primes p in `PR` with p || T_k (exponent exactly 1), p not dividing m
    out=[]
    for p in PR:
        if m%p==0: continue
        a=T_mod(T1,U1,m,k,p*p)
        if a%p==0 and a!=0:
            out.append(p)
            if len(out)>=nmax: break
    return out
rows=[]; dead=0; stuck=[]   # `dead` = families with T1 odd (no candidate); `stuck` = (m, k, level) where no witness was found
for m in range(7,MMAX+1,8):
    if any(m%(p*p)==0 for p in range(2,int(m**0.5)+1)): continue
    T1,U1=fund_pos(m); le=math.log10(2*T1)
    if T1%2==1: dead+=1; continue
    k0=1
    for p in factorint(m):
        if U1%p: k0*=p
    w1=witnesses(T1,U1,m,k0,NW)
    if not w1: stuck.append((m,k0,1)); rows.append((m,k0,k0*le,[],None,None,[],None)); continue
    k1=k0
    for p in w1: k1*=p
    w2=witnesses(T1,U1,m,k1,NW)
    if not w2: stuck.append((m,k1,2)); rows.append((m,k0,k0*le,w1,k1,k1*le,[],None)); continue
    k2=k1
    for p in w2: k2*=p
    rows.append((m,k0,k0*le,w1,k1,k1*le,w2,k2*le))
print(f"Familien m<={MMAX} (m=7 mod 8, quadratfrei): {len(rows)+dead} | 2-adisch tot (T1 ungerade): {dead} | Leiter gebaut: {len(rows)-len(stuck)} | steckengeblieben: {len(stuck)}")
print()
print(" m    k0   Stellen(T_k0)  Zeugen1        k1            Stellen(T_k1)   Zeugen2         log10(Stellen T_k2)")
for r in rows[:25]:
    m,k0,d0,w1,k1,d1,w2,d2=r
    print(f"{m:5} {k0:5} {d0:12.0f}   {str(w1):14} {str(k1):14} {('%.3g'%d1) if d1 else '-':>14}   {str(w2):15} {('%.1f'%math.log10(d2)) if d2 else '-':>8}")
if len(rows)>25: print(f"   ... ({len(rows)-25} weitere Familien)")
ok=[r for r in rows if r[7]]
print()
print(f"Stufe-2-Turm erreicht fuer {len(ok)} Familien. Kleinste Stufe-2-Schranke: {min(r[7] for r in ok):.3g} Stellen  (m={min(ok,key=lambda r:r[7])[0]})")
print(f"Kleinste Stufe-1-Schranke: {min(r[5] for r in rows if r[5]):.3g} Stellen  (m={min((r for r in rows if r[5]),key=lambda r:r[5])[0]})")
if stuck: print("STECKENGEBLIEBEN (kein Zeuge < 20000):", stuck)
