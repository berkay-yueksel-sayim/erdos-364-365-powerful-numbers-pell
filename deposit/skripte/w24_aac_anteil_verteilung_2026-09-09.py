# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w24_aac_anteil_verteilung_2026-09-09.py
# Purpose: distribution of the AAC fraction b = m/m'; makes the heuristic layer of Question 2.6 quantitative.
# Background: m' = product of the primes p | m with p not dividing U1. The rest b = m/m' is the product of the primes that DO
#   divide U1; they drop out of the index and make the tower faster: w(m) = 1/(2 m' log10 eps_m) = b/(2 m log10 eps_m).
#   A large b is therefore a "fast tower".
# Three measurements, each with a prediction:
#   (1) For each small prime p: how often does p also divide U1, given p | m? Naive model: 1/p.
#       This is the only assumption behind the whole heuristic; here it is MEASURED, not assumed.
#   (2) Distribution of b itself (how often b = 1, 3, 5, 7, ...) and the mean of b.
#       Model prediction from (1): E[b] = prod over p | m of ((1 - 1/p) * 1 + (1/p) * p) = prod (2 - 1/p) ~ 2^omega(m).
#   (3) The contribution-weighted tail: how much of c_M comes from kernels with b >= B? This decides whether the rare fast towers
#       can carry the sum.
# Reads: nothing (the squarefree kernels of one class mod 8 are generated up to MMAX).
# Writes: w24_aacanteil_c<KLASSE>_M<MMAX>_result.json in the current directory.
# Usage: python w24_aac_anteil_verteilung_2026-09-09.py [MMAX] [KLASSE]   (defaults 1000000 and 7;
#   `KLASSE` = class of m mod 8, 3 or 7)
# Controls (all three abort the run on failure):
#   PK-A: m = 87 = 3*29 has U1 = 3, so 3 | U1 and 29 does not divide it -> b = 3, m' = 29. Must come out exactly so.
#   PK-B (negative): m = 7 has U1 = 3 and 7 does not divide 3 -> b = 1. Must not be counted as a fast tower.
#   PK-C: the sum of the w(m) must reproduce the c_M values of w23 (same quantity, computed independently again).
#
# Own code.
import sys, math, json, time
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
from math import isqrt

MMAX   = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000   # `MMAX` = kernel limit
KLASSE = int(sys.argv[2]) if len(sys.argv) > 2 else 7
assert KLASSE in (3, 7)

# `spf_sieve(n)` = smallest-prime-factor sieve; `faktor(n, spf)` = factorization of n as {prime: exponent}
def spf_sieve(n):
    s = list(range(n + 1)); i = 2
    while i * i <= n:
        if s[i] == i:
            for j in range(i * i, n + 1, i):
                if s[j] == j: s[j] = i
        i += 1
    return s

def faktor(n, spf):
    d = {}
    while n > 1:
        p = spf[n]; e = 0
        while n % p == 0: n //= p; e += 1
        d[p] = e
    return d

def fund_fast(m):
    """(T1 mod 2m, U1 mod 2m, log10 T1) without big integers."""
    M = 2 * m; a0 = isqrt(m); P, Q, a = 0, 1, a0
    hm1, hm0 = 1 % M, a0 % M; km1, km0 = 0, 1 % M
    hf1, hf0 = 1.0, float(a0); s = 0; i = 0
    while True:
        P = a * Q - P; Q = (m - P * P) // Q; a = (a0 + P) // Q; i += 1
        if Q == 1 and i % 2 == 0: return hm0, km0, math.log10(hf0) + s
        hm1, hm0 = hm0, (a * hm0 + hm1) % M; km1, km0 = km0, (a * km0 + km1) % M
        hf1, hf0 = hf0, a * hf0 + hf1
        if hf0 > 1e100: hf0 /= 1e100; hf1 /= 1e100; s += 100

LOG2 = math.log10(2.0)
# `log10_eps(l10)` = log10 of eps = T1 + U1*sqrt(m) from l10 = log10 T1 (exact for small T1, else l10 + log10 2)
def log10_eps(l10):
    if l10 < 15.0:
        T1 = round(10.0 ** l10)
        if T1 >= 2: return math.log10(T1 + math.sqrt(T1 * T1 - 1.0))
    return l10 + LOG2

def fund_exakt(m):                              # `fund_exakt` = exact fundamental solution (T1, U1) with big integers
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0 * h0 - m * k0 * k0 != 1:
        P = a * Q - P; Q = (m - P * P) // Q; a = (a0 + P) // Q
        h1, h0 = h0, a * h0 + h1; k1, k0 = k0, a * k0 + k1
    return h0, k0

print(f'W24 -- Verteilung des AAC-Anteils b = m/m\', Klasse {KLASSE} (mod 8), Kerne bis {MMAX}')
t0 = time.time(); spf = spf_sieve(MMAX); print(f'  Sieb fertig ({time.time()-t0:.0f}s)')

# ---- PK-A / PK-B on the reference kernels ----
# `mref` = reference kernel, `b_soll`/`mp_soll` = expected b and m'; `T1e`/`U1e` = exact fundamental solution
for mref, b_soll, mp_soll in ((87, 3, 29), (7, 1, 7)):
    T1e, U1e = fund_exakt(mref)
    f = faktor(mref, spf); mp = 1
    for p in f:
        if U1e % p: mp *= p
    b = mref // mp
    if (b, mp) != (b_soll, mp_soll):
        sys.exit(f'ABBRUCH PK: m={mref} ergibt b={b}, m\'={mp}; erwartet b={b_soll}, m\'={mp_soll}')
print(f'  PK-A ok: m = 87 -> b = 3, m\' = 29 (3 teilt U1 = 3).   PK-B ok (negativ): m = 7 -> b = 1 (kein schneller Turm).')

# ---- main run ----
p_treffer = {}      # `p_treffer`: p -> [kernels with p | m, of these p | U1]
b_verteilung = {}   # `b_verteilung`: b -> count of kernels
b_gewicht = {}      # `b_gewicht`: b -> contribution to c_M
# `lebend` = living kernels (T1 even); `summe_b`, `summe_Eb` = running sums of b and of the model value E[b]
c = 0.0; lebend = 0; summe_b = 0; summe_Eb = 0.0
for m in range(KLASSE, MMAX + 1, 8):
    f = faktor(m, spf)
    if any(e > 1 for e in f.values()): continue
    tm, um, l10 = fund_fast(m)
    if tm % 2: continue
    lebend += 1
    mp = 1
    for p in f:
        t = p_treffer.setdefault(p, [0, 0]); t[0] += 1
        if um % p: mp *= p
        else: t[1] += 1
    b = m // mp
    # `leps` = log10 eps, `w` = contribution w(m) to c_M
    leps = log10_eps(l10); w = 1.0 / (2.0 * mp * leps)
    c += w
    b_verteilung[b] = b_verteilung.get(b, 0) + 1
    b_gewicht[b] = b_gewicht.get(b, 0.0) + w
    summe_b += b
    summe_Eb += math.prod(2.0 - 1.0 / p for p in f)      # model expectation for this m

print(f'  Lauf fertig ({time.time()-t0:.0f}s): {lebend} lebende Kerne, c_M = {c:.4f}')

# ---- PK-C: c_M against w23 ----
# `SOLL_C`: (class, MMAX) -> c_M value found by w23
SOLL_C = {(7, 1000000): 0.2151, (7, 10000000): 0.2217, (3, 1000000): 0.4513, (3, 10000000): 0.4566}
soll = SOLL_C.get((KLASSE, MMAX))
if soll is not None:
    if abs(c - soll) > 0.0005:
        sys.exit(f'ABBRUCH PK-C: c_M = {c:.4f}, W23 ergab {soll:.4f} -- zwei Rechnungen derselben Groesse weichen ab.')
    print(f'  PK-C ok: c_M stimmt mit W23 ueberein ({c:.4f} vs {soll:.4f}).')

# ---- (1) does p also divide U1? ----
print(f'\n=== (1) Wie oft teilt p auch U1, gegeben p | m?  (naives Modell: 1/p) ===')
print(f'{"p":>7} {"Kerne mit p|m":>14} {"davon p|U1":>11} {"gemessen":>10} {"Modell 1/p":>11} {"Verhaeltnis":>12}')
for p in sorted(p_treffer)[:14]:
    n_ges, n_tr = p_treffer[p]   # kernels with p | m / of these with p | U1; `gem` = measured fraction
    if n_ges < 50: continue
    gem = n_tr / n_ges
    print(f'{p:>7} {n_ges:>14} {n_tr:>11} {gem:>10.4f} {1/p:>11.4f} {gem*p:>12.2f}')

# ---- (2) distribution of b ----
print(f'\n=== (2) Verteilung des AAC-Anteils b = m/m\' ===')
print(f'{"b":>9} {"Kerne":>9} {"Anteil":>9} {"Beitrag zu c_M":>16} {"Anteil an c_M":>15}')
for b in sorted(b_verteilung)[:15]:
    n = b_verteilung[b]
    print(f'{b:>9} {n:>9} {100*n/lebend:>8.2f}% {b_gewicht[b]:>16.5f} {100*b_gewicht[b]/c:>14.2f}%')
print(f'  verschiedene Werte von b: {len(b_verteilung)}; groesstes b: {max(b_verteilung)}')
print(f'  Mittelwert von b: gemessen {summe_b/lebend:.4f}  |  Modell E[b] = prod(2 - 1/p): {summe_Eb/lebend:.4f}')

# ---- (3) contribution-weighted tail ----
print(f'\n=== (3) Wieviel von c_M kommt von Kernen mit b >= B? ===')
print(f'{"B":>9} {"Kerne":>9} {"Beitrag":>11} {"Anteil an c_M":>15}')
for B in (1, 3, 5, 7, 11, 21, 51, 101, 1001):
    n = sum(v for k, v in b_verteilung.items() if k >= B)
    g = sum(v for k, v in b_gewicht.items() if k >= B)
    print(f'{B:>9} {n:>9} {g:>11.5f} {100*g/c:>14.2f}%')
print('  Lesart: bleibt der Anteil bei wachsendem B schnell klein, tragen die seltenen schnellen Tuerme die Summe NICHT')
print('          -> Fall (A). Waechst er, koennen sie es doch -> Fall (B) moeglich.')

# result keys: `lebend` = living kernels, `b_mittel_gemessen` / `b_mittel_modell` = measured / model mean of b, `sekunden` =
# seconds
res = {'MMAX': MMAX, 'klasse': KLASSE, 'lebend': lebend, 'c_M': c,
       'b_mittel_gemessen': summe_b/lebend, 'b_mittel_modell': summe_Eb/lebend,
       'p_treffer': {str(k): v for k, v in sorted(p_treffer.items())[:40]},
       'b_verteilung': {str(k): v for k, v in sorted(b_verteilung.items())[:60]},
       'b_gewicht': {str(k): v for k, v in sorted(b_gewicht.items())[:60]},
       'sekunden': round(time.time()-t0, 1)}
fn = f'w24_aacanteil_c{KLASSE}_M{MMAX}_result.json'
json.dump(res, open(fn, 'w'), indent=1)
print(f'\nErgebnis in {fn}')
