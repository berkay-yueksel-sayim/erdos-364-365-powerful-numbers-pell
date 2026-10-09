# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w127_dichte_mitten_2026-09-28.py
# Which middles can carry triples: the share of squares, cubes and double squares among the powerful numbers n with 4 | n.
#
# Reads:    nothing (enumerates the powerful numbers up to N = 10^10, 10^12, 10^14 directly).
# Writes:   ../ergebnisse/w127_dichte_mitten_result.json and ../ergebnisse/w127_dichte_mitten_output.txt
# Usage:    python w127_dichte_mitten_2026-09-28.py   (no arguments; uses numpy)
# Controls: PK for every N: the exact counts of squares, cubes and double squares with 4 | n up to N (formulas below) must
#           match the enumeration, and no number may be both a square and a double square.
#
# MOTIVATION: the naive estimate "the share of squares among the powerful numbers tends to about 46 %" (= 1/2.173) counts ALL
#   powerful numbers, but the middle of a triple must be divisible by 4 (Lemma 1.1). Cube middles have density zero.
# THEORY: #{powerful <= X} ~ (zeta(3/2)/zeta(3))*sqrt(X) = 2.1732*sqrt(X) (Bateman-Grosswald). The local factor at 2 is
#   sum over e in {0,2,3,...} of 2^(-e/2) = 1 + (1/2)/(1 - 2^(-1/2)) = 2.70711; the part with e >= 2 (i.e. 4 | n) is 1.70711
#   => fraction 0.63060 => #{powerful <= X, 4 | n} ~ 1.37046*sqrt(X). Squares with 4 | n are (2t)^2 => sqrt(X)/2
#   => SHARE -> 0.5/1.37046 = 0.36484. Cubes with 4 | n are (2t)^3 => X^(1/3)/2 => share -> 0.
#   Expected: the measured square shares at 10^10, 10^12, 10^14 approach 0.3648 (not 0.46); the cube share falls toward 0.
# Theorem A (Thm 5.3) covers middles that are a square, a DOUBLE square or a fourth power, so double squares are measured as well.
#   Double squares with 4 | n are 8y^2 => sqrt(X)/(2*sqrt(2)) => share -> (1/(2*sqrt(2)))/1.37046 = 0.25798; together with the
#   squares (disjoint) -> 0.62282. Expected: the measured share of squares plus double squares approaches 0.6228.
# Own code (enumeration a^2*b^3 with b squarefree, which is unique, after Golomb).
import sys, json, time, pathlib
from math import isqrt
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'   # `HIER` = this directory, `ERG` = results directory
T0 = time.time(); zeilen = []   # `zeilen` = lines of the output file
def sag(s=''): print(s, flush=True); zeilen.append(s)   # `sag` = print a line and record it for the output file
def quadratfrei(n):   # `quadratfrei` = squarefree test by trial division
    p = 2
    while p * p <= n:
        if n % (p * p) == 0: return False
        p += 1
    return True
def powerful_bis(N):   # `powerful_bis` = sorted array of all powerful numbers <= N, as a^2*b^3 with b squarefree
    t = []; b = 1
    while b ** 3 <= N:
        if quadratfrei(b):
            t.append(np.arange(1, isqrt(N // b ** 3) + 1, dtype=np.int64) ** 2 * (b ** 3))
        b += 1
    a = np.concatenate(t); a.sort(); return a
def ist_quadrat(a): r = np.floor(np.sqrt(a.astype(np.float64))).astype(np.int64); return (r * r == a) | ((r + 1) * (r + 1) == a)   # `ist_quadrat` = square test
def ist_kubus(a):   # `ist_kubus` = cube test (checks the neighbors of the rounded floating-point root)
    r = np.round(np.cbrt(a.astype(np.float64))).astype(np.int64)
    return ((r - 1) ** 3 == a) | (r ** 3 == a) | ((r + 1) ** 3 == a)
GRENZ_Q = 0.5 / (2.173243 * (1.70711 / 2.70711))   # `GRENZ_Q` = limit of the square share (theory)
assert abs(GRENZ_Q - 0.36484) < 2e-4, GRENZ_Q
GRENZ_2Q = (1 / (2 * 2 ** 0.5)) / (2.173243 * (1.70711 / 2.70711))   # `GRENZ_2Q` = limit of the double-square share (theory)
def ist_doppelquadrat(a): return (a % 2 == 0) & ist_quadrat(a // 2)   # `ist_doppelquadrat` = double-square test (n = 2*y^2)
erg = {}   # `erg` = results per N
for N in (10**10, 10**12, 10**14):
    a = powerful_bis(N); m = a[a % 4 == 0]
    q = ist_quadrat(m); k = ist_kubus(m); d = ist_doppelquadrat(m)
    # PK: the number of squares with 4|n <= N is floor(sqrt(N)/2), of cubes floor(N^(1/3)/2) (even bases), of double squares
    #   8y^2 <= N: floor(sqrt(N/8))
    pk_q = int(q.sum()) == isqrt(N) // 2
    kb = int(round(N ** (1 / 3))); kb = kb - 1 if kb ** 3 > N else kb
    pk_k = int(k.sum()) == kb // 2
    pk_d = int(d.sum()) == isqrt(N // 8) and not bool((q & d).any())
    # keys of `e`: `mit_4` = count of powerful n with 4 | n; `quadrate`/`kuben`/`doppelquadrate` = squares/cubes/double squares
    #   among them; `beides` = both square and cube; `anteil_*` = share (`weder` = neither square nor cube, `satz_a` = square or
    #   double square); `verhaeltnis_4_zu_sqrtN` = ratio of `mit_4` to sqrt(N); `pk` = all three positive controls passed
    e = dict(powerful=int(len(a)), mit_4=int(len(m)), quadrate=int(q.sum()), kuben=int(k.sum()), beides=int((q & k).sum()),
             doppelquadrate=int(d.sum()), anteil_quadrat=float(q.mean()), anteil_kubus=float(k.mean()), anteil_weder=float((~q & ~k).mean()),
             anteil_doppelquadrat=float(d.mean()), anteil_satz_a=float((q | d).mean()),
             verhaeltnis_4_zu_sqrtN=len(m) / N ** 0.5, pk=bool(pk_q and pk_k and pk_d))
    erg[str(N)] = e
    sag(f'N = 10^{len(str(N)) - 1}: {len(m):,} powerful mit 4|n ({e["verhaeltnis_4_zu_sqrtN"]:.4f}·√N; Theorie 1,3705) · Quadrate {e["anteil_quadrat"]:.4%} '
        f'· Kuben {e["anteil_kubus"]:.4%} · weder {e["anteil_weder"]:.4%} · doppelte Quadrate {e["anteil_doppelquadrat"]:.4%} '
        f'· Bereich Satz A {e["anteil_satz_a"]:.4%} · PK {"✅" if e["pk"] else "🔴"}')
sag(f'Grenzwert Quadrat-Anteil (Theorie): {GRENZ_Q:.4%} — die „~46 %" der Analyse vom 30.08. rechnete ohne die Bedingung 4 | n.')
sag(f'Grenzwert Bereich von Satz A (Quadrate + doppelte Quadrate): {GRENZ_Q + GRENZ_2Q:.4%} (doppelte Quadrate allein {GRENZ_2Q:.4%}).')
ERG.joinpath('w127_dichte_mitten_result.json').write_text(json.dumps(dict(skript=pathlib.Path(__file__).name, ergebnis=erg, grenzwert_quadrat=GRENZ_Q,
    grenzwert_doppelquadrat=GRENZ_2Q, grenzwert_satz_a=GRENZ_Q + GRENZ_2Q, laufzeit_s=round(time.time() - T0, 1)), ensure_ascii=False, indent=1), encoding='utf-8')
ERG.joinpath('w127_dichte_mitten_output.txt').write_text('\n'.join(zeilen) + '\n', encoding='utf-8')
sag(f'Laufzeit {time.time() - T0:.0f} s')
