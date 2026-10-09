# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w129_x4plus1_aac_2026-09-28.py
# Purpose: for every prime q ≡ 1 (mod 8) up to GRENZE, find the smallest solution (v0, u0) of V² − q·U² = −1 (solvable, since
#   q ≡ 1 (mod 4) is prime) and test whether q | u0. This is the computation behind the remark that x⁴ + 1 = q³a² has no solution
#   for prime q (no q ≡ 1 (mod 8) below 40 000 with q | u0).
# Why this is the Ankeny–Artin–Chowla (AAC) question: for q ≡ 1 (mod 8), q | u0 iff q | u for the fundamental unit (t + u√q)/2.
#   If that unit is half-integral, u0 = u(3t² + qu²)/8, and 3t² + qu² ≡ 3t² ≡ −12 ≢ 0 (mod q), since t² ≡ −4 (mod q).
# Method: continued fraction of √q in integers (m, d, a); the convergents p_n/q_n are carried only modulo q. For odd period
#   length L, (p_{L−1}, q_{L−1}) is the smallest solution of x² − q·y² = −1, and u0 ≡ q_{L−1} (mod q).
# Expectation: 0 hits up to 10^5; a hit would be a counterexample to the AAC conjecture.
#   The run checks our own claim; it is not a new result.
# Reads: nothing. Writes: ergebnisse/w129_x4plus1_aac_result.json and ergebnisse/w129_x4plus1_aac_output.txt.
# Usage: python w129_x4plus1_aac_2026-09-28.py   (no arguments; limit GRENZE = 100 000)
# Controls (PK): (1) for q < 2000, exact integer arithmetic gives v0² − q·u0² = −1 and u0 mod q equals the carried residue;
#   (2) the detector must find a known divisor: for q = 41, (v0, u0) = (32, 5), so "5 | u0" must be reported;
#   (3) the period length is odd for all q.
import sys, json, time, pathlib
from math import isqrt
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'   # `HIER` = this folder; `ERG` = results folder
T0 = time.time(); zeilen = []   # `zeilen` = output lines, written to the output file at the end
def sag(s=''): print(s, flush=True); zeilen.append(s)   # `sag` = "say": print a line and keep it for the output file
GRENZE = 100_000   # `GRENZE` = limit for q

def primzahlen(n):   # sieve of Eratosthenes: all primes <= n
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def minus_loesung(q, modul=None, exakt=False):   # `minus_loesung` = solution of the "minus" equation x² − q·y² = −1
    """Continued fraction of √q; returns (L, p_{L-1}, q_{L-1}), either modulo `modul` or exactly (`exakt`=True)."""
    a0 = isqrt(q); m, d, a = 0, 1, a0
    p_alt, p = 1, a0; r_alt, r = 0, 1          # p_{-1}, p_0 ; q_{-1}, q_0
    L = 0
    while True:
        m = d * a - m; d = (q - m * m) // d; a = (a0 + m) // d; L += 1
        if a == 2 * a0:                        # end of the period: the convergent with index L−1 is (p, r)
            return L, p, r
        if exakt:
            p_alt, p = p, a * p + p_alt; r_alt, r = r, a * r + r_alt
        else:
            p_alt, p = p, (a * p + p_alt) % modul; r_alt, r = r, (a * r + r_alt) % modul

ok = True
# PK 1: exact versus modular computation, q < 2000
for q in [p for p in primzahlen(2000) if p % 8 == 1]:
    L, v, u = minus_loesung(q, exakt=True)
    Lm, _, um = minus_loesung(q, modul=q)
    ok &= (L % 2 == 1) and (v * v - q * u * u == -1) and (u % q == um) and (L == Lm)
sag(f'PK 1 (q < 2000, exakt: v₀² − q·u₀² = −1, Periode ungerade, u₀ mod q stimmt): {"✅" if ok else "🔴"}')
# PK 2: known divisor
L, v, u = minus_loesung(41, exakt=True)
_, _, u5 = minus_loesung(41, modul=5)
pk2 = (v, u) == (32, 5) and u5 == 0
sag(f'PK 2 (q = 41: (v₀, u₀) = ({v}, {u}); Detektor meldet 5 | u₀: {u5 == 0}): {"✅" if pk2 else "🔴"}'); ok &= pk2

# `treffer` = hits (q with q | u0); `ungerade_fehlt` = q whose period length is even (should stay empty);
# n = number of primes tested; maxL = longest period length seen
treffer, ungerade_fehlt, n, maxL = [], [], 0, 0
for q in primzahlen(GRENZE):
    if q % 8 != 1: continue
    n += 1
    L, _, u = minus_loesung(q, modul=q)
    maxL = max(maxL, L)
    if L % 2 == 0: ungerade_fehlt.append(q)
    if u == 0: treffer.append(q)
bis40k = [q for q in treffer if q < 40_000]
sag(f'Primzahlen q ≡ 1 (mod 8) bis {GRENZE:,}: {n:,} · groesste Periodenlaenge {maxL:,} · gerade Periode (duerfte nicht vorkommen): {len(ungerade_fehlt)}')
sag(f'q | u₀: {len(treffer)} Treffer bis {GRENZE:,} · davon unter 40 000: {len(bis40k)} — Behauptung der Synthese (kein Treffer < 40 000): {"bestaetigt" if not bis40k else "WIDERLEGT"}')
sag(f'PK gesamt: {"GRUEN" if ok and not ungerade_fehlt else "ROT"} · Laufzeit {time.time() - T0:.0f} s')
ERG.joinpath('w129_x4plus1_aac_result.json').write_text(json.dumps(dict(skript=pathlib.Path(__file__).name, grenze=GRENZE, anzahl_q=n, treffer=treffer,
    treffer_unter_40000=bis40k, groesste_periode=maxL, gerade_perioden=ungerade_fehlt, pk=bool(ok), laufzeit_s=round(time.time() - T0, 1)),
    ensure_ascii=False, indent=1), encoding='utf-8')
ERG.joinpath('w129_x4plus1_aac_output.txt').write_text('\n'.join(zeilen) + '\n', encoding='utf-8')
