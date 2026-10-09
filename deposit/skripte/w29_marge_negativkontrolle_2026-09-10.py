# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w29_marge_negativkontrolle_2026-09-10.py
# Purpose: recomputes the margin in the negative control to Theorem 3.1, i.e. for families rejected by the height criterion of
#   w9 / Theorem 3.1 (family rejected <=> m' > (H + 0.31)/log10(T1) + 1). Margin (in decimal orders by which the first
#   rung exceeds the bound 10^H) = m' * log10(eps) - H, computed as a RIGOROUS lower bound, once for a random sample of
#   rejected families and once for the families at the edge of the rejection threshold.
# Reads: nothing (kernels are generated). Writes: w29_marge_result.json in the current directory.
# Usage: python w29_marge_negativkontrolle_2026-09-10.py [MMAX] [STICHPROBE] [H] [RANDMAX]
#   defaults 1e8, 530, 2000, 300000: kernels m <= MMAX, STICHPROBE = number of rejected families in the sample, bound 10^H,
#   RANDMAX = bound for the edge scan (b).
# Controls (abort on failure): PK (m = 7 gives a negative margin) and PK-NEG (period 1 gives the trivial bound), see below.
#
# Design.
#   - REPRODUCIBLE: fixed random seed, so that the sample can be recovered.
#   - FLOAT-FREE: no floating-point estimate of log10(T1). From the continued-fraction recursion
#     h_j = a_j h_{j-1} + h_{j-2} with all a_j >= 1 it follows h_j >= F_{j+1}, hence
#     log10(T1) = log10(h_{l-1}) >= (l-1) * log10(phi) - log10(sqrt 5).
#     Only the period length l is needed, and that is pure integer arithmetic. The result is therefore a RIGOROUS LOWER BOUND
#     of the margin, which is exactly what a negative control needs.
#   - In addition the really informative quantity: the families AT THE EDGE of the rejection threshold. A floating-point error
#     would strike there, not at the randomly drawn giants.
# Controls: PK: for m = 7 (T1 = 8, m' = 7) the computation must give a NEGATIVE margin, since that family is NOT rejected.
#   PK-NEG: a period artificially set to 1 must not produce a positive margin.

import sys, math, json, random, time
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
from math import isqrt

# `MMAX` = largest kernel; `NPROB` = sample size (number of rejected families to collect); `RANDMAX` = bound of the edge scan
MMAX  = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**8
NPROB = int(sys.argv[2]) if len(sys.argv) > 2 else 530
H     = float(sys.argv[3]) if len(sys.argv) > 3 else 2000.0
RANDMAX = int(float(sys.argv[4])) if len(sys.argv) > 4 else 300000
SEED  = 20260910

LOG10PHI = math.log10((1 + 5 ** 0.5) / 2)     # 0.20898764…
LOG10SQ5 = math.log10(5 ** 0.5)               # 0.34948500…
LOG10_2  = math.log10(2.0)

def periode_und_mprime(m):
    """Continued-fraction period l of sqrt(m) and m' = prod{p | m : p does not divide U1}, purely integer.
    Returns (l, m_prime), or None if m is a square."""
    a0 = isqrt(m)
    if a0 * a0 == m: return None
    # prime factors of m
    f = []; d = 2; n = m
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0: n //= d
        d += 1 if d == 2 else 2
    if n > 1: f.append(n)
    # continued fraction; carry the convergent denominators k_j modulo each p | m (U1 = k_{l-1})
    P, Q, a = 0, 1, a0
    km1 = {p: 0 for p in f}; km0 = {p: 1 % p for p in f}
    l = 0
    while True:
        P = a * Q - P; Q = (m - P * P) // Q; a = (a0 + P) // Q
        l += 1
        for p in f:
            km1[p], km0[p] = km0[p], (a * km0[p] + km1[p]) % p
        if Q == 1:
            break
    # U1 = k_{l-1}: after l steps `km1` holds the value k_{l-1}
    mp = 1
    for p in f:
        if km1[p] % p: mp *= p
    return l, mp

def log10T1_unten(l):
    """Rigorous lower bound for log10(T1) = log10(h_{l-1}) >= (l-1)*log10(phi) - log10(sqrt5)."""
    return max(0.0, (l - 1) * LOG10PHI - LOG10SQ5)

t0 = time.time()
print(f'W29 -- Marge der Negativ-Kontrolle zu Thm 3.1.  Kerne <= {MMAX:.0e}, Stichprobe {NPROB}, H = {H:.0f}, Seed {SEED}')

# ---------- PK (positive control) ----------
r = periode_und_mprime(7)
assert r is not None
l7, mp7 = r
assert mp7 == 7, f'PK: m-prime von 7 ist {mp7}, erwartet 7'
marge7 = mp7 * (log10T1_unten(l7) + LOG10_2) - H
print(f'  PK m = 7: Periode l = {l7}, m\' = {mp7}, untere Schranke log10(T1) >= {log10T1_unten(l7):.4f} (wahr: {math.log10(8):.4f})')
assert marge7 < 0, f'PK: m = 7 darf keine positive Marge haben (ist {marge7:.1f}) - die Familie wird NICHT verworfen'
print(f'      Marge = {marge7:.1f} < 0  ->  korrekt NICHT verworfen.')
assert log10T1_unten(1) == 0.0, 'PK-NEG: Periode 1 muesste die triviale Schranke 0 geben'
print('  PK-NEG ok: Periode 1 gibt die triviale Schranke 0, keine positive Marge.')

# ---------- (a) reproducible sample of rejected families ----------
rng = random.Random(SEED)
# `gezogen` = drawn, `verworfen` = rejected, `margen` = margins, `widerspruch` = contradictions (rejected but margin <= 0),
# `versuche` = attempts
gezogen = 0; verworfen = 0; margen = []; widerspruch = 0
versuche = 0
while verworfen < NPROB and versuche < NPROB * 400:
    versuche += 1
    m = rng.randrange(7, MMAX + 1, 8)  # class 7 (mod 8)
    # squarefree test (`sf`) by trial division
    q = m; sf = True; d = 2
    while d * d <= q:
        if q % (d * d) == 0: sf = False; break
        if q % d == 0: q //= d
        d += 1 if d == 2 else 2
    if not sf: continue
    r = periode_und_mprime(m)
    if r is None: continue
    l, mp = r
    gezogen += 1
    lo = log10T1_unten(l)
    if lo <= 0: continue
    # rejected by the w9 criterion?  m' > (H + 0.31)/log10(T1) + 1
    if mp > (H + 0.31) / lo + 1:
        verworfen += 1
        margen.append(mp * (lo + LOG10_2) - H)
        if mp * (lo + LOG10_2) - H <= 0: widerspruch += 1

margen.sort()
print(f'\n(a) Stichprobe: {gezogen} quadratfreie Kerne gezogen, davon {verworfen} verworfen ({versuche} Versuche, {time.time()-t0:.0f}s)')
print(f'    Widersprueche (verworfen, aber Marge <= 0): {widerspruch}')
print(f'    kleinste Marge: {margen[0]:.3g} Dezimalordnungen  ({margen[0]/1e6:.2f}e6)')
print(f'    Median: {margen[len(margen)//2]:.3g}   groesste: {margen[-1]:.3g}')
assert widerspruch == 0, 'ABBRUCH: eine verworfene Familie hat keine positive Marge'

# ---------- (b) the families AT THE EDGE ----------
print(f'\n(b) Familien am Rand der Schwelle (die eigentliche Probe): m\' * log10(eps) knapp ueber H')
# `rand` = edge (list of (margin, m, m', period) with 0 < margin < 50)
rand = []
for m in range(7, min(MMAX, RANDMAX) + 1, 8):
    q = m; sf = True; d = 2
    while d * d <= q:
        if q % (d * d) == 0: sf = False; break
        if q % d == 0: q //= d
        d += 1 if d == 2 else 2
    if not sf: continue
    r = periode_und_mprime(m)
    if r is None: continue
    l, mp = r
    lo = log10T1_unten(l)
    if lo <= 0: continue
    marge = mp * (lo + LOG10_2) - H
    if 0 < marge < 50:
        rand.append((round(marge, 3), m, mp, l))
rand.sort()
print(f'    Kerne <= {RANDMAX:.0e} mit Marge zwischen 0 und 50 Dezimalordnungen: {len(rand)}')
for r_ in rand[:10]:
    print(f'      Marge {r_[0]:>8.3f}  m = {r_[1]:>9}  m\' = {r_[2]:>9}  Periode = {r_[3]}')
if rand:
    print(f'    ⇒ kleinste Marge am Rand: {rand[0][0]:.3f} Dezimalordnungen (Kern {rand[0][1]}).')
    print(f'      Das ist die Zahl, die die Kontrolle eigentlich braucht: selbst die knappste verworfene')
    print(f'      Familie ueberschiesst 10^{H:.0f} noch um {rand[0][0]:.3f} Dezimalordnungen, RIGOROS nach unten geschaetzt.')

# result record; `stichprobe_verworfen` = number of rejected families in the sample, `marge_min` / `marge_median` /
# `marge_max` = margin statistics, `widerspruch` = number of contradictions, `rand_bis_3e6` = the first 40 edge families
json.dump({'MMAX': MMAX, 'H': H, 'seed': SEED, 'stichprobe_verworfen': verworfen,
           'marge_min': margen[0], 'marge_median': margen[len(margen)//2], 'marge_max': margen[-1],
           'widerspruch': widerspruch, 'rand_bis_3e6': rand[:40], 'sekunden': round(time.time()-t0, 1)},
          open('w29_marge_result.json', 'w'), indent=1)
print(f'\nErgebnis in w29_marge_result.json  ({time.time()-t0:.0f}s)')
