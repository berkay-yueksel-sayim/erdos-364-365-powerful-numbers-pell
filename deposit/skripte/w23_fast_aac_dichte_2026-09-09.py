# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w23_fast_aac_dichte_2026-09-09.py
# Purpose: density of the "fast-AAC" kernels, the measurement that decides Conjecture 2.6. For the squarefree kernels m of one
#   residue class mod 8 it computes m' (the product of the primes p | m with p not dividing U1), the contribution
#   w(m) = 1/(2*m'*log10(eps_m)) to c_M, the counts N_j(X) = #{m <= X : m' = j} per decade, and the share of the fast-AAC
#   kernels (m' < m) in c_M. m/m' is the "AAC fraction"; a kernel with m' = 1 is an AAC kernel.
# Reads: nothing (the kernels are generated up to MMAX).
# Writes: w23_fastaac_c<KLASSE>_M<MMAX>_result.json in the current directory.
# Usage: python w23_fast_aac_dichte_2026-09-09.py [MMAX] [KLASSE]   (defaults 1000000 and 7; `KLASSE` = class of m mod 8, 7 or 3)
# Controls (abort on failure): PK-A, PK-B, PK-C, PK-D, described below.
#
# Question. Conj. 2.6: a kernel with m' = j contributes w(m) = 1/(2*m'*log10(eps_m)) to the pair count.
#   If N_j(X) = #{m <= X : m' = j} grows like log X, c_M diverges -> case (B).
#   If the set is finite or thinner than log X / (log log X)^(1+d), c_M converges -> case (A).
#   So the script measures N_j(X) per decade, plus the contribution of the fast-AAC kernels to c_M.
#
# Method. As in w9: continued fraction computed modularly (`fund_fast`), no big integers. The computation is done mod 2m so that
#   both the PARITY of T1 (living family <=> T1 even) and U1 mod m can be read off.
#   log10(eps) = log10(T1 + U1*sqrt(m)) = log10(2*T1) up to O(1/T1^2), since T1^2 - m U1^2 = 1.
#
# Controls:
#   PK-A (cross-check against w13, an independent route): the values per decade printed in Conjecture 2.6 must be reproduced
#     (class 7: 0.059, 0.134, 0.166, 0.176, 0.176, 0.179, 0.181 for M = 10^1..10^7; class 3: 0.291, 0.367, 0.398, 0.403, 0.406;
#     tolerance 0.0016 in the code). The printed values are the HEIGHT-RESTRICTED sum: only kernels whose first tower rung lies
#     below 10^2000, i.e. m' * log10(eps_m) <= 2000. Comparing them with the UNRESTRICTED sum fails (0.1708 vs 0.166 at X = 10^3).
#     Both quantities are therefore kept apart and printed: the control runs against the restricted one, the question of
#     Conj. 2.6 depends on the unrestricted one.
#   PK-B: for class 7 and MMAX >= 4099215 the AAC kernel 4099215 (the only AAC kernel with even T1 below 10^8) must appear with
#     m' = 1; for smaller MMAX no AAC kernel with even T1 may appear. 209991 is also AAC but has odd T1 (see below).
#   PK-C (negative): m = 7 must have m' = 7, i.e. it is NOT fast-AAC.
#   PK-D (tool): `fund_fast` must agree with the exact big-integer computation at 200 sample kernels.

import sys, math, json, time
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
from math import isqrt

MMAX   = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000
KLASSE = int(sys.argv[2]) if len(sys.argv) > 2 else 7  # 7 or 3 (mod 8)
assert KLASSE in (3, 7), 'KLASSE muss 3 oder 7 sein'

# ---- smallest prime factor, for fast factoring (`spf_sieve`, `faktor` = factorization into {prime: exponent}) ----
def spf_sieve(n):
    s = list(range(n+1))
    i = 2
    while i*i <= n:
        if s[i] == i:
            for j in range(i*i, n+1, i):
                if s[j] == j: s[j] = i
        i += 1
    return s

def faktor(n, spf):
    f = {}
    while n > 1:
        p = spf[n]; e = 0
        while n % p == 0: n //= p; e += 1
        f[p] = e
    return f

# `fund_exakt` = exact fundamental solution (T1, U1) of x^2 - m*y^2 = 1 by continued fraction (big integers)
def fund_exakt(m):
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0

def fund_fast(m):
    """(T1 mod 2m, U1 mod 2m, log10(T1)) without big integers."""
    M = 2*m
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    hm1, hm0 = 1 % M, a0 % M
    km1, km0 = 0, 1 % M
    hf1, hf0 = 1.0, float(a0)
    s = 0; i = 0
    while True:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        i += 1
        if Q == 1 and i % 2 == 0:
            return hm0, km0, math.log10(hf0) + s
        hm1, hm0 = hm0, (a*hm0 + hm1) % M
        km1, km0 = km0, (a*km0 + km1) % M
        hf1, hf0 = hf0, a*hf0 + hf1
        if hf0 > 1e100:
            hf0 /= 1e100; hf1 /= 1e100; s += 100

print(f'W23 -- Fast-AAC-Dichte, Klasse {KLASSE} (mod 8), Kerne bis {MMAX}')
t0 = time.time()
spf = spf_sieve(MMAX)
print(f'  Sieb fertig ({time.time()-t0:.0f}s)')

# ---- PK-D: tool against the exact computation (`proben` = samples) ----
proben = 0
for m in range(KLASSE, min(MMAX, 60000)+1, 8):
    f = faktor(m, spf)
    if any(e > 1 for e in f.values()): continue
    proben += 1
    if proben > 200: break
    tm, um, l10 = fund_fast(m)
    T1, U1 = fund_exakt(m)
    if (T1 % (2*m), U1 % (2*m)) != (tm, um) or abs(math.log10(T1) - l10) > 1e-6:
        sys.exit(f'ABBRUCH PK-D: fund_fast weicht bei m={m} ab: {(tm,um,l10)} vs {(T1%(2*m), U1%(2*m), math.log10(T1))}')
print(f'  PK-D ok: fund_fast = exakte Rechnung an {min(proben,200)} Kernen (Werte UND log10).')

# ---- main run ----
# `DEK` = decade bounds 10^1, 10^2, ...
DEK = [10**e for e in range(1, 9) if 10**e <= MMAX*10]
H_TURM = 2000.0  # height bound of the w9/w13 runs (middle < 10^2000)
LOG2   = math.log10(2.0)

def log10_eps(l10):
    """log10(T1 + U1*sqrt(m)) from log10(T1). Since T1^2 - m U1^2 = 1, eps = T1 + sqrt(T1^2 - 1),
    so log10(eps) = l10 + log10(1 + sqrt(1 - T1^-2)). For large T1 this is l10 + log10(2);
    for SMALL T1 it is not -- the error is about 5 % at m = 3 and would shift class 3 by 0.015
    (caught by the exact route)."""
    if l10 < 15.0:
        T1 = round(10.0**l10)
        if T1 >= 2:
            return math.log10(T1 + math.sqrt(T1*T1 - 1.0))
    return l10 + LOG2

Nj      = {}  # `Nj` = N_j: m' -> list of the fast-AAC kernels with that m'
c_dek   = {}  # decade -> cumulative c_M, UNRESTRICTED (the quantity of Conj. 2.6)
cR_dek  = {}  # decade -> cumulative c_M, RESTRICTED to m'*log10(eps) <= 2000
fastaac_dek = {}  # decade -> (# fast-AAC, their unrestricted c contribution)
# `c`, `cR` = running sums (unrestricted, restricted); `fa_n`, `fa_c` = number and c contribution of the fast-AAC kernels
c = 0.0; cR = 0.0; fa_c = 0.0; fa_n = 0
# `fam` = squarefree kernels seen, `lebend` = living families (T1 even), `aac_kerne` = kernels with m' = 1
# (`aac_kerne` = AAC kernels), `di` = index into `DEK`
fam = 0; lebend = 0
aac_kerne = []
di = 0
for m in range(KLASSE, MMAX+1, 8):
    while di < len(DEK) and m > DEK[di]:
        c_dek[DEK[di]] = c; cR_dek[DEK[di]] = cR; fastaac_dek[DEK[di]] = (fa_n, fa_c); di += 1
    f = faktor(m, spf)
    if any(e > 1 for e in f.values()): continue  # squarefree
    fam += 1
    tm, um, l10 = fund_fast(m)
    if tm % 2: continue  # only living families: T1 even
    lebend += 1
    mp = 1
    for p in f:
        if um % p: mp *= p  # p does not divide U1  ->  p enters m'
    leps = log10_eps(l10)
    w = 1.0 / (2.0 * mp * leps)
    c += w
    # height filter exactly as in w9 (kmax from log10(T1), not from log10(eps), with integer rounding up):
    # `cR` collects the restricted sum
    if mp <= int((H_TURM + 0.31) / l10) + 1: cR += w
    if mp < m:  # at least one p | m divides U1  ->  fast-AAC
        fa_n += 1; fa_c += w
        Nj.setdefault(mp, []).append(m)
        if mp == 1: aac_kerne.append(m)
while di < len(DEK):
    c_dek[DEK[di]] = c; cR_dek[DEK[di]] = cR; fastaac_dek[DEK[di]] = (fa_n, fa_c); di += 1

print(f'  Lauf fertig ({time.time()-t0:.0f}s): {fam} quadratfreie Kerne, davon {lebend} lebend (T1 gerade).')

# ---- PK-A: cross-check against the c_M printed in Conj. 2.6 (`SOLL` = expected values per class and decade X) ----
SOLL = {7: {10:0.059, 100:0.134, 1000:0.166, 10000:0.176, 100000:0.176, 1000000:0.179, 10000000:0.181},
        3: {10:0.291, 100:0.367, 1000:0.398, 10000:0.403, 100000:0.406}}
geprueft = 0
for X, soll in SOLL[KLASSE].items():
    if X > MMAX: continue
    ist = cR_dek.get(X)  # against the RESTRICTED sum, see the header
    if ist is None: continue
    geprueft += 1
    if abs(ist - soll) > 0.0016:
        sys.exit(f'ABBRUCH PK-A: beschraenktes c_M bei X={X} ist {ist:.4f}, Conj. 2.6 druckt {soll:.3f}. '
                 f'Entweder W13 oder W23 rechnet falsch -- nicht weiterrechnen.')
print(f'  PK-A ok: {geprueft} Dekaden-Werte der BESCHRAENKTEN Summe stimmen mit Conjecture 2.6 ueberein.')

# ---- PK-B / PK-C ----
# 209991 (class 7) and 1752299 (class 3) are AAC kernels, but with ODD T1: they are parity-dead and are sorted out
# before m' is computed, so they must not be expected in the m' = 1 list.
# The only AAC kernel with even T1 below 10^8 is 4099215 (class 7).
if KLASSE == 7:
    if MMAX >= 4099215:
        if 4099215 not in aac_kerne:
            sys.exit('ABBRUCH PK-B: der AAC-Kern 4099215 (T1 gerade) erscheint nicht mit m\' = 1.')
        print('  PK-B ok: AAC-Kern 4099215 mit m\' = 1 gefunden.')
    else:
        if aac_kerne:
            sys.exit(f'ABBRUCH PK-B: unterhalb {MMAX} darf es keinen AAC-Kern mit geradem T1 geben, '
                     f'gefunden: {aac_kerne}. Widerspricht Prop. 6.2/6.4 -- pruefen, nicht ignorieren.')
        print(f'  PK-B ok (Abwesenheit erwartet): kein AAC-Kern mit geradem T1 unterhalb {MMAX}; '
              f'209991 ist AAC, aber T1 ungerade (paritaetstot) und faellt vor der m\'-Rechnung heraus.')
if KLASSE == 7:
    tm7, um7, _ = fund_fast(7)
    mp7 = 7 if um7 % 7 else 1
    if mp7 != 7: sys.exit('ABBRUCH PK-C: m = 7 wird als fast-AAC gefuehrt.')
    print('  PK-C ok (negativ): m = 7 hat m\' = 7, ist also nicht fast-AAC.')

# ---- output ----
print(f'\n=== c_M je Dekade (Klasse {KLASSE}) ===')
print(f'  {"X":>6}  {"c_M unbeschraenkt":>18} {"Zuwachs":>9}   {"c_M beschraenkt":>16} {"Zuwachs":>9}   fast-AAC')
pu = pr = 0.0
for X in sorted(c_dek):
    n_fa, c_fa = fastaac_dek[X]
    anteil = 100*c_fa/c_dek[X] if c_dek[X] else 0
    print(f'  10^{len(str(X))-1:<3} {c_dek[X]:>18.4f} {c_dek[X]-pu:>+9.4f}   {cR_dek[X]:>16.4f} {cR_dek[X]-pr:>+9.4f}   {n_fa:5d} Kerne, {c_fa:.4f} ({anteil:.1f} %)')
    pu, pr = c_dek[X], cR_dek[X]
print('  🔑 Lesart: die Frage von Conjecture 2.6 haengt an der UNBESCHRAENKTEN Spalte.')
print('     Zuwachs je Dekade faellt weiter  -> Fall (A), c_M konvergiert.')
print('     Zuwachs je Dekade wird konstant  -> Fall (B), c_M ~ const * log M, divergiert.')
print('     Die BESCHRAENKTE Spalte (= die in Conj. 2.6 gedruckte) flacht durch die Hoehenschranke')
print('     zwangslaeufig ab und darf fuer diese Frage NICHT gelesen werden.')

print(f'\n=== N_j(X): Anzahl Kerne mit m\' = j (nur fast-AAC, also m\' < m) ===')
print(f'{"j":>8}  {"Anzahl":>7}   kleinste Kerne')
for j in sorted(Nj)[:25]:
    ms = Nj[j]
    print(f'{j:>8}  {len(ms):>7}   {ms[:5]}')
print(f'  ... insgesamt {len(Nj)} verschiedene Werte von m\', {fa_n} fast-AAC-Kerne.')

print(f'\n=== Wachstum von #fast-AAC je Dekade (der Kern der Frage) ===')
prev = 0
for X in sorted(fastaac_dek):
    n = fastaac_dek[X][0]
    print(f'  bis 10^{len(str(X))-1:<2}: {n:6d}   neu in dieser Dekade: {n-prev:6d}')
    prev = n
print('  Lesart: bleibt "neu je Dekade" ungefaehr KONSTANT, waechst die Zahl wie log X  -> Fall (B).')
print('          faellt sie je Dekade,                                                   -> Fall (A).')

res = {'MMAX': MMAX, 'klasse': KLASSE, 'familien': fam, 'lebend': lebend,
       'c_dekade': {str(k): v for k, v in c_dek.items()},
       'fastaac_dekade': {str(k): v for k, v in fastaac_dek.items()},
       'Nj': {str(j): {'anzahl': len(v), 'kleinste': v[:10]} for j, v in sorted(Nj.items())[:60]},
       'aac_kerne': aac_kerne, 'sekunden': round(time.time()-t0, 1)}
fn = f'w23_fastaac_c{KLASSE}_M{MMAX}_result.json'
json.dump(res, open(fn, 'w'), indent=1)
print(f'\nErgebnis in {fn}  ({time.time()-t0:.0f}s)')
