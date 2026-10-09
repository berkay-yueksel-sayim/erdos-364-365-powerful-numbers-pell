# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w107_paar_korrektur_2026-09-25.py
# Purpose: neighbor correction for pairs of consecutive powerful numbers. Compares the observed number of pairs (n, n + 1) of
#   powerful numbers, up to 2^52, with the naive expectation that treats n and n + 1 as independent and with a corrected one.
# Observation (script w100): 24 pairs n, n + 1 up to 2⁴⁶, i.e. 0.52 per binary section; the naive calculation that treats n and
#   n + 1 as independent gives (c/2)²·ln 2 = 0.818 per section (c = ζ(3/2)/ζ(3)).
# Model (as for the Hardy–Littlewood constant of twin primes): the independent baseline is corrected, prime by prime, by
#   the factor C_p = P(v_p(n) ≠ 1 and v_p(n+1) ≠ 1) / P(v_p(n) ≠ 1)².  With a := P(v_p(n) = 1) = (p − 1)/p², and since n and
#   n + 1 cannot both be divisible by p, the two events "v_p = 1" are disjoint, so C_p = (1 − 2a) / (1 − a)² and C = ∏_p C_p.
#   Expected pairs in [A, B]: C · ∫ ρ(x)² dx, with the density ρ = F' from the counting formula F(x) = c₁√x + c₂x^(1/3)
#   (the second term is a small-x correction). This is a HEURISTIC (locally independent model), not a theorem; here it is
#   measured against the data.
# Data: complete enumeration of all powerful numbers, section by section up to 2^52 (each section as a²b³ with numpy; checksum
#   against the counting formula of w100); pairs = neighbors at distance 1. Pairs across a section boundary (2^k − 1, 2^k): by
#   w103, 2^k − 1 is never powerful for 2 ≤ k ≤ 63, so there are none.
# Reads:  ergebnisse/w100_brille_binaer_result.json (keys `abschnitte` = sections, `paare` = pairs).
# Writes: ergebnisse/w107_paar_korrektur_result.json and ergebnisse/w107_paar_korrektur_output.txt.
# Usage:  python w107_paar_korrektur_2026-09-25.py   (no arguments)
# Expectations, stated before the run:
#   E1  C ≈ 0.75 (rough estimate: 0.889 · 0.918 · 0.964 · 0.980 · …); asymptotically ≈ 0.61 pairs per section.
#   E2  cumulatively, the corrected expectation is closer to the measured 10 / 14 / 18 / 24 (up to 10⁸ / 10¹⁰ / 10¹² / 2⁴⁶) than
#       the naive one [assumption]; for such small counts Poisson noise (±√N) decides.
#   E3  in the sections 46 … 51 together about 3–4 new pairs (6 × 0.61).
# Controls (asserted): the enumeration equals the counting formula in every section (w100 values), and the 24 known pairs up to
#   2⁴⁶ are found again.

import sys, json, math, time, pathlib
from math import isqrt
import numpy as np
import mpmath
sys.stdout.reconfigure(encoding='utf-8')
ERG = pathlib.Path(__file__).resolve().parent.parent/'ergebnisse'   # `ERG` = results folder
aus = []   # `aus` = output lines, written to the output file at the end
def sag(s=''):   # `sag` = "say": print a line and keep it for the output file
    print(s, flush=True); aus.append(s)
t0 = time.time()
sag('='*100); sag('w107 — NACHBARSCHAFTS-KORREKTUR FUER POWERFUL PAARE   ' + time.strftime('%Y-%m-%d %H:%M')); sag('='*100)
# `c1`, `c2` = constants of the counting formula F(x) = c₁√x + c₂x^(1/3): c₁ = ζ(3/2)/ζ(3), c₂ = ζ(2/3)/ζ(2)
c1 = float(mpmath.zeta(1.5)/mpmath.zeta(3)); c2 = float(mpmath.zeta(mpmath.mpf(2)/3)/mpmath.zeta(2))

# ---------------- (1) the constant C
N_P = 10**7   # `N_P` = the product C = ∏ C_p runs over the primes up to N_P
sieb = bytearray([1]) * (N_P + 1); sieb[0] = sieb[1] = 0   # `sieb` = sieve of Eratosthenes up to N_P
for p in range(2, isqrt(N_P) + 1):
    if sieb[p]: sieb[p*p::p] = bytearray(len(range(p*p, N_P + 1, p)))
logC, erste = 0.0, []   # `logC` = running sum of log C_p; `erste` = the first values C_p for p <= 13
for p in range(2, N_P + 1):
    if not sieb[p]: continue
    a = (p - 1) / (p * p); Cp = (1 - 2*a) / (1 - a)**2
    logC += math.log(Cp)
    if p <= 13: erste.append((p, round(Cp, 4)))
C = math.exp(logC)
# `schwanz` = tail beyond 10⁷: |Σ log C_p| ≲ Σ 1/p² < 10⁻⁷, negligible (times 0)
schwanz = sum(1/(p*p) for p in (10**7 + 1,)) * 0
naiv = c1*c1/4 * math.log(2)   # `naiv` = naive expected number of pairs per binary section, (c₁/2)²·ln 2
sag(f'(1) C_p fuer p = 2 … 13: ' + ' · '.join(f'{p}: {c}' for p, c in erste))
sag(f'    C = ∏_(p ≤ 10⁷) C_p = {C:.5f}   ⇒ asymptotisch {naiv:.4f} · C = {naiv*C:.4f} Paare je Binaerabschnitt (naiv {naiv:.4f})')

# ---------------- (2) enumeration per section up to 2^52
w100 = json.load(open(ERG/'w100_brille_binaer_result.json', encoding='utf-8'))
PW = {a['k']: a['powerful'] for a in w100['abschnitte']}   # `PW` = number of powerful numbers in section k (from w100)
bekannt = set(int(s) for s in w100['paare'])   # `bekannt` = the 24 known pairs (smaller member)
KMAX = 52   # `KMAX` = enumerate the sections k = 1 … KMAX − 1, i.e. all numbers below 2^KMAX
# Every powerful number is a²b³ with b squarefree; b needs b³ ≤ 2^KMAX, hence b ≤ 2^(KMAX/3). (A bound of 2^(KMAX/4) is too small:
# the count check against the counting formula fails at section 39.)
BMAX = int(round(2 ** (KMAX / 3))) + 3
sfb = bytearray([1]) * (BMAX + 1); sfb[0] = 0   # `sfb` = squarefree flags for b up to BMAX
for q in range(2, isqrt(BMAX) + 1):
    sfb[q*q::q*q] = bytearray(len(range(q*q, BMAX + 1, q*q)))
B_SF = [b for b in range(1, BMAX + 1) if sfb[b]]   # `B_SF` = the squarefree b
paare, je_k = [], {}   # `paare` = all pairs found (smaller member); `je_k` = number of pairs per section
for k in range(1, KMAX):
    lo, hi = 2**k, 2**(k+1) - 1
    teile = []   # `teile` = parts: arrays of the numbers a²·b³ in this section, one array per b
    for b in B_SF:
        b3 = b**3
        if b3 > hi: break
        a0 = isqrt((lo - 1) // b3) + 1; a1 = isqrt(hi // b3)
        if a1 >= a0:
            a = np.arange(a0, a1 + 1, dtype=np.int64); teile.append(a * a * np.int64(b3))
    if not teile: je_k[k] = 0; continue
    v = np.sort(np.concatenate(teile)); del teile
    assert v.size == PW[k], ('PK VERLETZT: Anzahl ≠ Zaehlformel', k, v.size, PW[k])
    idx = np.nonzero(np.diff(v) == 1)[0]
    neu = [int(v[i]) for i in idx]
    paare.extend(neu); je_k[k] = len(neu)
    if k >= 44: sag(f'    Abschnitt {k}: {v.size:,} powerful · {len(neu)} Paare  ({time.time()-t0:.0f}s)')
    del v
bis46 = {n for n in paare if n < 2**46}
assert bis46 == bekannt, ('PK VERLETZT: die 24 bekannten Paare', sorted(bis46 ^ bekannt))
sag(f'PK ✅  jede Abschnittszaehlung = Zaehlformel (w100); die 24 Paare bis 2⁴⁶ exakt wiedergefunden.')
sag(f'(2) Paare bis 2^{KMAX}: {len(paare)} — neu in den Abschnitten 46 … {KMAX-1}: ' + ', '.join(str(n) for n in paare if n >= 2**46))

# ---------------- (3) expectation versus measurement
def rho(x): return c1/(2*math.sqrt(x)) + c2/(3*x**(2/3))   # `rho` = density ρ = F'(x)
def int_rho2(A, B, n=4000):                     # ∫_A^B ρ(x)² dx on a log scale (midpoint rule)
    la, lb = math.log(A), math.log(B); h = (lb - la)/n; s = 0.0
    for i in range(n):
        x = math.exp(la + (i + 0.5)*h); r = max(rho(x), 0.0); s += r*r*x*h
    return s
# `START`: the first pair is (8, 9); below it the density approximation is meaningless
START = 8
zeilen = []   # `zeilen` = rows of the comparison table (stored in the result file)
for X in (10**8, 10**10, 10**12, 2**46, 2**48, 2**50, 2**52):
    # `gem` = measured number of pairs up to X; `e_naiv`, `e_korr` = naive / corrected expectation
    gem = sum(1 for n in paare if n + 1 <= X)
    e_naiv = int_rho2(START, X); e_korr = C * e_naiv
    z_naiv = (gem - e_naiv) / math.sqrt(e_naiv); z_korr = (gem - e_korr) / math.sqrt(e_korr)
    zeilen.append(dict(X=X if X < 2**40 else f'2^{round(math.log2(X))}', gemessen=gem, naiv=round(e_naiv, 2), korrigiert=round(e_korr, 2),
                       z_naiv=round(z_naiv, 2), z_korr=round(z_korr, 2)))
    sag(f'    bis {("10^" + str(round(math.log10(X)))) if X < 2**40 else "2^" + str(round(math.log2(X))):<6}: gemessen {gem:>3} · naiv {e_naiv:6.2f} (z {z_naiv:+.2f}) · '
        f'korrigiert {e_korr:6.2f} (z {z_korr:+.2f})')
# `hinten` = the last sections 30 … KMAX − 1, compared separately (`gem_h` measured, `e_h` expected)
hinten = [k for k in range(30, KMAX)]
gem_h = sum(je_k[k] for k in hinten); e_h = C * sum(int_rho2(2**k, 2**(k+1)) for k in hinten)
sag(f'    Abschnitte 30 … {KMAX-1} allein: gemessen {gem_h}, korrigiert erwartet {e_h:.2f} (z {(gem_h - e_h)/math.sqrt(e_h):+.2f}), '
    f'naiv {e_h/C:.2f} (z {(gem_h - e_h/C)/math.sqrt(e_h/C):+.2f})')
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), C=C, C_p_erste=erste, naiv_je_abschnitt=naiv,
           korrigiert_je_abschnitt=naiv*C, paare=[str(n) for n in paare], je_abschnitt=je_k, vergleich=zeilen,
           hinten=dict(von=30, bis=KMAX-1, gemessen=gem_h, erwartet_korr=round(e_h, 3), erwartet_naiv=round(e_h/C, 3)), laufzeit_s=round(time.time()-t0))
(ERG/'w107_paar_korrektur_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w107_paar_korrektur_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
sag(f'Laufzeit {time.time()-t0:.0f}s')
