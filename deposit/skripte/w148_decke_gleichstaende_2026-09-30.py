# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w148_decke_gleichstaende_2026-09-30.py
# The ceiling table of Part III (remark rem:decke), made tie-proof: correction of a tie artifact in w49.
#
# Reads:    ../ergebnisse/w49_mechanismen_und_10920_result.json (key `teil_c_gcd`, for the cross-check PK2);
#           otherwise the powerful numbers up to 10^14 are enumerated directly (numpy).
# Writes:   ../ergebnisse/w148_decke_gleichstaende_result.json and ../ergebnisse/w148_decke_gleichstaende_output.txt
# Usage:    python w148_decke_gleichstaende_2026-09-30.py   (no arguments)
# Controls: PK1, PK2, PK3 (positive) and NK inside `anteil` (negative); a failed control aborts with an assertion.
#
# BACKGROUND: w49 takes the "8 non-powerful comparison values" with [...][:8] after sorting by "frequency descending, value
#   ascending". With 89 occurrences, however, FOUR non-powerful values tie, and w49 keeps two of them. The same artifact was
#   already removed in w125 for the head of the table. This script is the same measurement on the same data as w49, only without
#   the arbitrary order among equal frequencies.
# WHAT: up to 10^14, (1) the group G200 = all gap values with occ >= occ of the 200th value (tie-proof, as in w125) and its
#   non-powerful members; per value: sqfull(g), occ, the share of occurrences with gcd(n, g) = sqfull(g), and the rank of
#   sqfull(g) as a SPAN [1 + #(occ > c), #(occ >= c)]. (2) the same quantity for the group G20 of the most frequent values (occ
#   >= occ of the 20th), as a min-max span of the share.
# EXPECTATIONS (set before the run): PK1 the inventory is the same as in w125 (21 663 503 powerful numbers, 8 235 300 gap values,
#   top value 44100 with 163 occurrences). PK2 the 8 values printed by w49 have the same occ and the same share as in w49
#   (`teil_c_gcd`: `gcd_max` / `anz`). PK3 G200 has 201 values, 10 of them non-powerful (w125).
#   NK: every gcd(n, g) is powerful and divides sqfull(g) (Lemma III.1.1); a violation would be an error in the tool.
# Own code, same convention as w125 (1 is powerful; a^2*b^3 with b squarefree).
import sys, json, time, pathlib, math
from math import isqrt, gcd
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'   # `HIER` = this directory, `ERG` = results directory
T0 = time.time(); zeilen = []   # `zeilen` = lines of the output file
def sag(s=''): print(s, flush=True); zeilen.append(s)   # `sag` = print a line and record it for the output file

def faktoren(n):   # `faktoren` = prime factorization of n as {prime: exponent} (trial division)
    f = {}; p = 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1 if p == 2 else 2
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
def ist_powerful(n): return n > 0 and all(e >= 2 for e in faktoren(n).values())   # `ist_powerful` = powerful test
def quadratfrei(n): return all(e == 1 for e in faktoren(n).values())   # `quadratfrei` = squarefree test
def sqfull(n): return math.prod(p ** e for p, e in faktoren(n).items() if e >= 2)   # `sqfull` = square-full part of n

def powerful_bis(N):   # `powerful_bis` = sorted array of all powerful numbers <= N, as a^2*b^3 with b squarefree
    teile = []; b = 1
    while b ** 3 <= N:
        if quadratfrei(b):
            amax = isqrt(N // b ** 3)
            teile.append(np.arange(1, amax + 1, dtype=np.int64) ** 2 * (b ** 3))
        b += 1
    x = np.unique(np.concatenate(teile)); return x[x <= N]

N = 10 ** 14
# `P` = powerful numbers <= N; `L` = gaps between consecutive ones; `vals`/`cnt` = distinct gap values and
#   their number of occurrences
P = powerful_bis(N); L = np.diff(P); vals, cnt = np.unique(L, return_counts=True)
top = int(vals[np.argmax(cnt)]); pk1 = (len(P), len(vals), top, int(cnt.max())) == (21663503, 8235300, 44100, 163)
sag(f'PK1 {"✅" if pk1 else "❌"}  {len(P)} powerful, {len(vals)} Werte, Top 1 = {top} mit {int(cnt.max())}'); assert pk1
ordn = np.argsort(-cnt, kind='stable'); c_sorted = cnt[ordn]   # `c_sorted` = occurrence counts in descending order
# `gruppe` = returns (c, idx): c = count of the k-th most frequent value, idx = indices of all values with count >= c
def gruppe(k):
    c = int(c_sorted[k - 1]); idx = np.nonzero(cnt >= c)[0]; return c, idx
occ = dict(zip(vals.tolist(), cnt.tolist()))   # `occ` = occurrences: gap value -> number of occurrences
def rangspanne(g):   # `rangspanne` = rank span of the gap value g: (1 + #(occ > c), #(occ >= c)) with c = occ[g]
    c = occ.get(g, 0); return int((cnt > c).sum()) + 1, int((cnt >= c).sum())
def anteil(g):
    # `anteil` = share: over all positions i with gap g, count how often gcd(P[i], g) = sqfull(g); returns (hits, positions).
    # NK: each gcd must be powerful and divide sqfull(g) (Lemma III.1.1).
    pos = np.nonzero(L == g)[0]; s = sqfull(g); treffer = 0
    for i in pos.tolist():
        d = gcd(int(P[i]), g)
        assert ist_powerful(d) and s % d == 0, ('NK Lemma verletzt', g, int(P[i]), d)
        treffer += d == s
    return treffer, len(pos)

c200, g200 = gruppe(200); c20, g20 = gruppe(20)
# `npw` = non-powerful gap values of G200, sorted by occurrence (descending), then by value
npw = sorted([int(vals[i]) for i in g200 if not ist_powerful(int(vals[i]))], key=lambda g: (-occ[g], g))
pk3 = (len(g200), len(npw)) == (201, 10)
sag(f'PK3 {"✅" if pk3 else "❌"}  G200: {len(g200)} Werte (occ ≥ {c200}), davon nicht powerful {len(npw)}'); assert pk3
w49 = json.loads((ERG/'w49_mechanismen_und_10920_result.json').read_text(encoding='utf-8'))['teil_c_gcd']
zeilen_npw = []   # `zeilen_npw` = table rows for the non-powerful values of G200
for g in npw:
    t, n = anteil(g); s = sqfull(g); r = rangspanne(s)
    zeilen_npw.append(dict(g=g, sqfull=s, occ=n, gcd_max=t, anteil=round(100 * t / n, 1), rang_sqfull=list(r), in_w49=str(g) in w49))
    sag(f'   g = {g:>7}  sqfull = {s:>6}  occ = {n:>4}  gcd = sqfull: {t:>3} ({100*t/n:.1f} %)  Rang sqfull {r[0]}–{r[1]}  in w49: {str(g) in w49}')
pk2 = all(z['occ'] == w49[str(z['g'])]['anz'] and z['gcd_max'] == w49[str(z['g'])]['gcd_max'] for z in zeilen_npw if z['in_w49'])
sag(f'PK2 {"✅" if pk2 else "❌"}  die {sum(z["in_w49"] for z in zeilen_npw)} Werte aus w49: occ und gcd-Anteil identisch'); assert pk2
zeilen_20 = []   # `zeilen_20` = table rows for the values of G20
for i in g20:
    g = int(vals[i]); t, n = anteil(g); zeilen_20.append(dict(g=g, occ=n, gcd_max=t, anteil=round(100 * t / n, 1), powerful=ist_powerful(g)))
# shares (in percent) for G20 and for the non-powerful values
a20 = [z['anteil'] for z in zeilen_20]; an = [z['anteil'] for z in zeilen_npw]
sag(f'G20: {len(g20)} Werte (occ ≥ {c20}), alle powerful: {all(z["powerful"] for z in zeilen_20)}, Anteil {min(a20):.1f}–{max(a20):.1f} %')
sag(f'nicht-powerful in G200: Anteil {min(an):.1f}–{max(an):.1f} %')
erg = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), N=N, pk=dict(pk1=pk1, pk2=pk2, pk3=pk3),
           g200=dict(schwelle=c200, groesse=len(g200), nicht_powerful=zeilen_npw),
           g20=dict(schwelle=c20, groesse=len(g20), werte=zeilen_20, anteil_min=min(a20), anteil_max=max(a20)),
           anteil_npw_min=min(an), anteil_npw_max=max(an), laufzeit_s=round(time.time() - T0, 1))
(ERG/'w148_decke_gleichstaende_result.json').write_text(json.dumps(erg, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w148_decke_gleichstaende_output.txt').write_text('\n'.join(zeilen) + '\n', encoding='utf-8')
sag(f'Ergebnis: w148_decke_gleichstaende_result.json ({time.time() - T0:.0f} s)')
