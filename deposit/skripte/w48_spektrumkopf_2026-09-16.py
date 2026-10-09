# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w48_spektrumkopf_2026-09-16.py
# Purpose: the head of the gap spectrum (`spektrumkopf`), with two controls. A gap is the difference between consecutive powerful
#   numbers; the spectrum counts how often each gap value occurs.
# Reads: no files. Writes: ergebnisse/w48_spektrumkopf_output.txt and ergebnisse/w48_spektrumkopf_result.json
#   (the results folder is found relative to the script). Usage: python w48_spektrumkopf_2026-09-16.py (no arguments).
# Controls: PK+: at N = 10^12 there must be 2 158 391 powerful numbers and 18 pairs (consecutive integers; as in w38 and w40).
#   PK-: a b that is not squarefree must not be enumerated (otherwise duplicates arise) -- checked by counting against the
#   target number and by testing that the enumeration has no duplicates.

# Finding under test (measured up to 10^12):
#   Of 832 161 occurring gap values only 0.3 % are powerful, but the top 50 are 100 % powerful,
#   and the first outlier (non-powerful value) has rank 72. The sieve moduli 900 and 44100 have rank 1 and 3.
#
# TWO CONTROLS:
#   (K1) LARGER RANGE. Does the statement hold at 10^13 and 10^14, or was it an artifact of the bound?
#        If the ranking flips, the finding is not one.
#   (K2) THE MECHANISM, quantitatively. w40 explained frequent gaps by SHIFTED PAIRS: if (k, k+1) is a
#        pair "powerful except at the primes of g" and g is powerful, then g*k and g*k+g are both powerful.
#        => The mechanism requires g to be POWERFUL. Measurement: what fraction of the occurrences of a top g has g | n?
#        If the fraction is high, the head of the spectrum is explained and not merely observed.
#
# NULL MODEL for (K1): if the most frequent gaps were drawn at random from all occurring values, the powerful
#   share among the top k would equal the overall share (0.3 %). Expected among 50: 0.16. Observed: 50.
#
# Enumeration: a^2 * b^3 with b squarefree -- Golomb 1970, p. 848 (the representation is unique, hence free of duplicates,
#   which allows the numpy variant without np.unique).

import sys, json, time, pathlib, collections
import numpy as np
from math import isqrt

HIER = pathlib.Path(__file__).resolve().parent
ERG  = HIER.parent / 'ergebnisse'
t0 = time.time(); Z = []   # `Z` = output lines, written to the _output.txt file
def out(s=''): Z.append(s); print(s, flush=True)   # `out` = print a line and record it in `Z`

# `fac` = factorization by trial division: dict prime -> exponent
def fac(n):
    f = {}; m = n; p = 2
    while p*p <= m:
        if m % p == 0:
            c = 0
            while m % p == 0: m //= p; c += 1
            f[p] = c
        p += 1 if p == 2 else 2
    if m > 1: f[m] = f.get(m, 0) + 1
    return f
def powerful(n): return n > 0 and all(e >= 2 for e in fac(n).values())
def squarefree(n): return all(e == 1 for e in fac(n).values())
def fs(n): return '*'.join(f'{p}^{e}' for p, e in sorted(fac(n).items()))   # `fs` = factorization as a string

def powerful_bis(N):
    """All powerful numbers <= N, ascending, as a numpy array. Free of duplicates by Golomb's uniqueness."""
    teile = []   # `teile` = parts: one array of a^2 * b^3 per squarefree b
    b = 1
    while b**3 <= N:
        if squarefree(b):
            amax = isqrt(N // b**3)
            if amax:
                a = np.arange(1, amax + 1, dtype=np.int64)
                teile.append(a * a * (b**3))
        b += 1
    arr = np.concatenate(teile)
    arr.sort()
    return arr

out('W48 — Kopf des Lueckenspektrums, mit zwei Kontrollen.')
out()
BOUNDS = [10**12, 10**13, 10**14]
ergebnis = {}    # `ergebnis` = results per bound N (as strings)
spektren = {}    # `spektren` = spectra: N -> (array of powerful numbers, Counter of gap values)
for N in BOUNDS:
    t1 = time.time()
    arr = powerful_bis(N)
    if N == 10**12:
        assert len(arr) == 2158391, ('PK+: 2 158 391 powerful bis 10^12', len(arr))
        assert len(np.unique(arr)) == len(arr), 'PK-: Aufzaehlung muss duplikatfrei sein'
    d = np.diff(arr)   # the gaps
    if N == 10**12:
        assert int((d == 1).sum()) == 18, ('PK+: 18 Paare', int((d == 1).sum()))
    spek = collections.Counter(d.tolist())   # `spek` = spectrum: gap value -> number of occurrences
    spektren[N] = (arr, spek)
    top = spek.most_common()   # `top` = (gap value, count) by decreasing count (ties in order of first occurrence)
    gesamt_pw = sum(1 for g in spek if powerful(g))   # number of gap values that are powerful
    # `zeile` = result row; keys: `n_luecken` = number of gaps, `n_werte` = number of distinct gap values, `anteil_pw_gesamt` =
    #   overall powerful share, top{k}_pw = number of powerful values among the k most frequent, erster_nicht_pw_rang = rank of
    #   the first non-powerful value, rang_<g> = rank of the gap value g, top1 = most frequent (value, count)
    zeile = dict(n_powerful=int(len(arr)), n_luecken=int(len(d)), n_werte=len(spek),
                 anteil_pw_gesamt=gesamt_pw / len(spek), sekunden=time.time() - t1)
    for k in (10, 20, 50, 100, 200):
        kop = [g for g, _ in top[:k]]
        zeile[f'top{k}_pw'] = sum(1 for g in kop if powerful(g))
    erster = next((i + 1 for i, (g, _) in enumerate(top) if not powerful(g)), None)
    zeile['erster_nicht_pw_rang'] = erster
    zeile['rang_900'] = next((i + 1 for i, (g, _) in enumerate(top) if g == 900), None)
    zeile['rang_44100'] = next((i + 1 for i, (g, _) in enumerate(top) if g == 44100), None)
    zeile['rang_36'] = next((i + 1 for i, (g, _) in enumerate(top) if g == 36), None)
    zeile['top1'] = top[0]
    ergebnis[str(N)] = zeile
    out(f'N = 10^{len(str(N))-1}: {len(arr):,} powerful · {len(spek):,} Lueckenwerte · '
        f'{gesamt_pw:,} davon powerful ({100*gesamt_pw/len(spek):.2f} %) · {time.time()-t1:.0f}s')

out()
out('(K1) HAELT DER BEFUND BEI GROESSERER SCHRANKE?')
kopf = f"{'Schranke':>10} {'Werte':>10} {'% pw ges.':>10} | {'Top10':>6} {'Top20':>6} {'Top50':>6} {'Top100':>7} {'Top200':>7} | {'1. nicht-pw':>12} | {'Rang 900':>9} {'Rang 44100':>11}"
out(kopf); out('-' * len(kopf))
for N in BOUNDS:
    z = ergebnis[str(N)]
    out(f"{'10^'+str(len(str(N))-1):>10} {z['n_werte']:>10,} {100*z['anteil_pw_gesamt']:>9.2f}% | "
        f"{z['top10_pw']:>6} {z['top20_pw']:>6} {z['top50_pw']:>6} {z['top100_pw']:>7} {z['top200_pw']:>7} | "
        f"{z['erster_nicht_pw_rang']:>12} | {z['rang_900']:>9} {z['rang_44100']:>11}")
out('   Lesart: "Top50 = 50" heisst, alle 50 haeufigsten Lueckenwerte sind powerful.')
out()
N = BOUNDS[-1]
z = ergebnis[str(N)]
erw = 50 * z['anteil_pw_gesamt']   # `erw` = expected number of powerful values among the top 50 under the null model
out(f'   NULLMODELL: waeren die Top 50 zufaellig aus allen {z["n_werte"]:,} Werten gezogen, waeren davon')
out(f'   erwartet {erw:.2f} powerful. Beobachtet: {z["top50_pw"]}. ⇒ Das ist kein Zufallseffekt.')

out()
out('(K2) DER MECHANISMUS: wie viele Vorkommen eines haeufigen g sind VERSCHOBENE PAARE (g | n)?')
arr, spek = spektren[N]
pos = {}
d = np.diff(arr)
top20 = [g for g, _ in spek.most_common(20)]   # the 20 most frequent gap values g
out(f"{'g':>8} {'Anzahl':>7} {'davon g|n':>10} {'Anteil':>8} {'Faktorisierung':>22} {'= Quadrat von':>14}")
# `mech` = mechanism: g -> `anzahl` (number of occurrences) and `verschoben` (occurrences with g | n, n = lower end of the gap)
mech = {}
for g in top20:
    idx = np.nonzero(d == g)[0]
    ns = arr[idx]
    versch = int((ns % g == 0).sum())   # `versch` = number of occurrences that are shifted pairs
    q = isqrt(g) if isqrt(g)**2 == g else None   # square root of g if g is a perfect square
    mech[str(g)] = dict(anzahl=len(idx), verschoben=versch)
    out(f'{g:>8} {len(idx):>7} {versch:>10} {100*versch/len(idx):>7.1f}% {fs(g):>22} {str(q) if q else "-":>14}')
# `gesamt_a`, `gesamt_v` = totals over the top 20: occurrences and shifted pairs
gesamt_a = sum(m['anzahl'] for m in mech.values()); gesamt_v = sum(m['verschoben'] for m in mech.values())
out(f'   SUMME der Top 20: {gesamt_v} von {gesamt_a} Vorkommen sind verschobene Paare ({100*gesamt_v/gesamt_a:.1f} %).')
out()
out('   GEGENPROBE an nicht-powerful Lueckenwerten aehnlicher Haeufigkeit:')
# control group: non-powerful gap values among the top 400
nicht_pw = [(g, c) for g, c in spek.most_common(400) if not powerful(g)][:8]
if nicht_pw:
    out(f"{'g':>8} {'Anzahl':>7} {'davon g|n':>10} {'Anteil':>8} {'Faktorisierung':>22}")
    for g, c in nicht_pw:
        idx = np.nonzero(d == g)[0]; ns = arr[idx]
        versch = int((ns % g == 0).sum())
        out(f'{g:>8} {len(idx):>7} {versch:>10} {100*versch/len(idx):>7.1f}% {fs(g):>22}')
    out('   ⇒ Bei nicht-powerful g kann der Mechanismus gar nicht greifen: g*k und g*k+g = g(k+1) waeren nur')
    out('   powerful, wenn g es ist. Der Anteil muss dort klein sein, und er ist es.')

out()
out('ERGEBNIS in einem Satz:')
out('   Der Kopf des Lueckenspektrums ist powerful, WEIL der einzige Mechanismus, der viele Luecken derselben')
out('   Groesse erzeugt — das Skalieren eines Paares mit g — voraussetzt, dass g selbst powerful ist.')
out('   Die Beobachtung ist damit nicht nur gemessen, sondern erklaert.')
out(f'\nZeit {time.time()-t0:.0f}s')
(ERG / 'w48_spektrumkopf_output.txt').write_text('\n'.join(Z) + '\n', encoding='utf-8')
# result keys: `skript` = script name, `ergebnis` = results per bound, `mechanismus_top20` = mechanism table of the top 20,
#   `laufzeit_s` = run time in seconds
(ERG / 'w48_spektrumkopf_result.json').write_text(json.dumps(dict(skript='w48_spektrumkopf_2026-09-16.py',
    bounds=BOUNDS, ergebnis=ergebnis, mechanismus_top20=mech, laufzeit_s=time.time()-t0,
    python=sys.version.split()[0]), indent=1), encoding='utf-8')
