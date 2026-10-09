# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# -*- coding: utf-8 -*-
# w121_zuwachs_zerlegung_2026-09-27.py
# Where does the growth of the sum over families of 1/ln F come from? Decomposes the increment per decade of the bound B.
# Background: the lower bound liminf N(x)/ln x >= sum 1/ln F (1.6822 at B = 10^6) grows per factor 10 in B by +0.119, +0.088,
#   +0.058 (families with a square) and +0.060, +0.045, +0.031 (families without a square). A tail model would be a heuristic
#   and does not belong in the paper; instead this script counts which families carry the increment per decade, split by period
#   (m_b, g, resp. M) and by regulator ln eps. By Walker (1976, end of Section 3) EVERY pair lies in exactly one of these
#   families. Also recorded: the run time per decade, as a basis for a cost estimate for B = 10^7.
# Single source: `pell` and `ln_eps` are taken from w113; `ist_quadrat`, `primfaktoren_qf`, `zerlegungen`, `restklasse` from
#   w115 (loaded via `ast`, because both scripts start their run on import). Only the family selection is duplicated (the loop
#   bodies of w113 and w115, without counting or controls); hence PK1: the partial sums must match the canonical values of w113
#   (B = 10^6) and w115 (D = 10^6) to within 1e-9.
# Unlike w113/w115 this script keeps no list of families but accumulates the counts while computing (B = 10^7 would otherwise
#   need several GB).
# Usage: python w121_zuwachs_zerlegung_2026-09-27.py [B]   (B = bound, default 1 000 000; another B adds the suffix _B<B> to
#   the output file names).
# Reads: the scripts w113_familien_weit_2026-09-25.py and w115_walker_familien_2026-09-26.py (parsed, not imported) and, for B
#   >= 10^6, ergebnisse/w113_familien_weit_B1000000_result.json and ergebnisse/w115_walker_familien_D1000000_result.json.
# Writes: ergebnisse/w121_zuwachs_zerlegung<suffix>_result.json and ..._output.txt.
# Controls: PK1 (partial sums against w113 and w115, tolerance 1e-9), PK2 (the share of periods <= 3 in the increment of the
#   square families above 10^4 must reproduce the 18.1 % of w113). Both are asserted; there is no negative control.
# Expectation, fixed in advance: the earlier expectation in w113 (increment carried by few b with small m_b) was not confirmed
#   (measured share 18.1 %); the increment is expected to spread over many families with medium period; which regulator range
#   carries it is [open].
import ast, sys, json, math, time, pathlib
from array import array
from math import isqrt, gcd, log
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'
BMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
TAG = '' if BMAX == 1_000_000 else f'_B{BMAX}'
# `aus` = collected output lines; `sag` (= say) prints a line and records it for the output file.
aus = []
def sag(s=''):
    print(s, flush=True); aus.append(s)
t0 = time.time()
sag('='*100); sag(f'w121 — ZUWACHS VON Σ 1/ln F, ZERLEGT JE DEKADE (beide Familienarten), B = {BMAX:,}   ' + time.strftime('%Y-%m-%d %H:%M')); sag('='*100)

# `lade` = load: executes only the named function definitions of a script (parsed with `ast`) and returns their namespace.
def lade(datei, namen):
    ns = {'isqrt': isqrt, 'gcd': gcd, 'log': log, 'math': math}
    for k in ast.parse((HIER/datei).read_text(encoding='utf-8')).body:
        if isinstance(k, ast.FunctionDef) and k.name in namen:
            exec(compile(ast.Module([k], []), datei, 'exec'), ns)
    assert all(n in ns for n in namen), (datei, [n for n in namen if n not in ns])
    return ns
n113 = lade('w113_familien_weit_2026-09-25.py', ('pell', 'ln_eps'))
n115 = lade('w115_walker_familien_2026-09-26.py', ('ist_quadrat', 'primfaktoren_qf', 'zerlegungen', 'restklasse'))
pell, ln_eps = n113['pell'], n113['ln_eps']
ist_quadrat, primfaktoren_qf, zerlegungen, restklasse = (n115[n] for n in ('ist_quadrat', 'primfaktoren_qf', 'zerlegungen', 'restklasse'))
spf = array('i', range(BMAX + 1))               # `spf` = smallest prime factor (w115 uses a list; an array saves memory at 10^7)
for p in range(2, isqrt(BMAX) + 1):
    if spf[p] == p:
        for q in range(p*p, BMAX + 1, p):
            if spf[q] == q: spf[q] = p
n115['spf'] = spf                               # `primfaktoren_qf` reads `spf` as a global variable of its namespace

PER = [(1, 1), (2, 3), (4, 10), (11, 100), (101, 1000), (1001, 10**18)]           # `PER` = period classes (inclusive bounds)
LE = [(0, 5), (5, 20), (20, 100), (100, 1000), (1000, 1e18)]                     # `LE` = regulator classes of ln eps
# `NDEK` = number of decades up to B; `dekade(key)` = decade index j of key; `stat[art][j]` = accumulator per family kind
# `art` ('Q' = families with a square, 'W' = families without a square, Walker) and decade j; `leer` builds an empty one:
# `n` = number of families, `summe` = sum of 1/ln F, `max` / `arg` = largest single contribution and its family,
# `per` / `le` = sums per period class / regulator class, `beitraege` = all contributions.
NDEK = len(str(BMAX)) - 1 + (0 if BMAX == 10**(len(str(BMAX)) - 1) else 1)
def dekade(key): return len(str(key - 1)) - 1   # key in (10^j, 10^(j+1)] ⇒ j
def leer():
    return dict(n=0, summe=0.0, max=0.0, arg=None, per=[0.0]*len(PER), le=[0.0]*len(LE), beitraege=array('d'))
stat = {art: [leer() for _ in range(NDEK)] for art in ('Q', 'W')}
# `eintragen` = record: adds the contribution 1/ln F of one family (`wer` = label of the family) to its decade and classes.
def eintragen(art, key, periode, le, lnF, wer):
    c = 1/lnF; s = stat[art][dekade(key)]
    s['n'] += 1; s['summe'] += c; s['beitraege'].append(c)
    if c > s['max']: s['max'], s['arg'] = c, dict(wer=wer, periode=periode, ln_eps=round(le, 3))
    s['per'][next(i for i, (a, b) in enumerate(PER) if a <= periode <= b)] += c
    s['le'][next(i for i, (a, b) in enumerate(LE) if a <= le < b)] += c

# `zeit` = elapsed seconds when each power of 10 is passed.
zeit = {}; naechste = 10
for m in range(2, BMAX + 1):
    if m > naechste: zeit[str(naechste)] = round(time.time() - t0, 1); naechste *= 10
    # Progress line; it sits before the `continue` below because 100 000*k is never squarefree (after it, it would never
    #   print).
    if m % 100_000 == 0: sag(f'   … B bis {m:,} ({time.time() - t0:.0f} s)')
    # `ps` = prime factors of m, or None if m is not squarefree.
    ps = primfaktoren_qf(m)
    if ps is None: continue
    (x1, y1), neg = pell(m); le = ln_eps(x1, y1, m)
    # family (b,1) with b = m, as in w113
    mb = m // gcd(m, y1 % m); eintragen('Q', m, mb, le, 2*mb*le, f'(b,1) b={m}')
    # family (1,d) with d = m, as in w113 (only if the negative Pell equation is solvable)
    if neg:
        y0, x0 = neg; g = m // gcd(m, x0 % m)
        if g % 2 == 1: eintragen('Q', m, g, le, 2*g*le, f'(1,d) d={m}')
    # families without a square (Walker), as in w115: m = b*d with d*X0^2 - b*Y0^2 = 1
    if len(ps) >= 2:
        for b, d in zerlegungen(ps):
            if (x1 + 1) % (2*d) or (x1 - 1) % (2*b): continue
            X0 = ist_quadrat((x1 + 1)//(2*d)); Y0 = ist_quadrat((x1 - 1)//(2*b))
            if not X0 or Y0 is None or d*X0*X0 - b*Y0*Y0 != 1: continue
            rk = restklasse(b, d, x1, y1, X0, Y0)
            if rk is None: continue
            eintragen('W', m, rk[1], le, 2*rk[1]*le, f'(b,d)=({b},{d})')
zeit[str(BMAX)] = round(time.time() - t0, 1)

# ---- PK1: cumulative partial sums against the canonical values
def kum(art, B):
    j = len(str(B)) - 2; return sum(stat[art][i]['summe'] for i in range(j + 1))
pk1 = []
if BMAX >= 10**6:
    r113 = json.load(open(ERG/'w113_familien_weit_B1000000_result.json', encoding='utf-8'))['teilsummen']
    r115 = json.load(open(ERG/'w115_walker_familien_D1000000_result.json', encoding='utf-8'))['summe_1_durch_lnF']
    for B in (10, 100, 1000, 10_000, 100_000, 1_000_000):
        pk1.append(('Q', B, kum('Q', B), r113[str(B)]))
    for B in (1000, 10_000, 100_000, 1_000_000):
        pk1.append(('W', B, kum('W', B), r115[str(B)]))
    ok1 = all(abs(a - b) < 1e-9 for _, _, a, b in pk1)
    sag(f'PK1 Teilsummen gegen w113 (B = 10⁶) und w115 (D = 10⁶), {len(pk1)} Werte: ' + ('✅ alle auf 1e-9 gleich' if ok1 else '🔴 ' + str([x for x in pk1 if abs(x[2] - x[3]) >= 1e-9])))
    assert ok1
    zu = sum(stat['Q'][j]['summe'] for j in range(4, NDEK)); k3 = sum(sum(stat['Q'][j]['per'][:2]) for j in range(4, NDEK))
    sag(f'PK2 Anteil der Perioden ≤ 3 am Zuwachs der Quadrat-Familien ueber 10⁴: {k3/zu:.1%} (w113: 18,1 %) ' + ('✅' if abs(k3/zu - 0.181) < 0.0006 else '🔴'))
    assert abs(k3/zu - 0.181) < 0.0006
sag()

# ---- decomposition per decade
# `anteil_oben(beitr, q)` = share of the total carried by the largest fraction q of the contributions.
def anteil_oben(beitr, q):
    v = sorted(beitr, reverse=True); k = max(1, int(len(v)*q)); t = sum(v); return (sum(v[:k]) / t) if t else 0.0
erg = {}
for art, name in (('Q', 'MIT Quadrat (b,1)/(1,d)'), ('W', 'OHNE Quadrat (Walker, D = b·d)')):
    sag(f'Familien {name}:')
    sag(f'   {"Dekade":>16} {"Familien":>9} {"Σ 1/ln F":>9} {"groesster":>9} {"Top 1 %":>7} {"Top 10":>7} | Periode 1 · 2–3 · 4–10 · 11–100 · 101–1000 · >1000 | ln ε <5 · 5–20 · 20–100 · 100–1000 · >1000')
    erg[art] = []
    for j in range(NDEK):
        s = stat[art][j]
        if not s['n']: continue
        v = sorted(s['beitraege'], reverse=True); top10 = sum(v[:10]) / s['summe']
        z = dict(dekade=f'(1e{j}, 1e{j+1}]', familien=s['n'], summe=s['summe'], groesster=s['max'], groesster_wer=s['arg'],
                 top1pct=anteil_oben(v, 0.01), top10=top10,
                 periode={f'{a}-{b if b < 10**17 else "inf"}': x for (a, b), x in zip(PER, s['per'])},
                 ln_eps={f'{a}-{b if b < 1e17 else "inf"}': x for (a, b), x in zip(LE, s['le'])})
        erg[art].append(z)
        pp = ' · '.join(f'{x/s["summe"]:5.1%}' for x in s['per']); ll = ' · '.join(f'{x/s["summe"]:5.1%}' for x in s['le'])
        sag(f'   {z["dekade"]:>16} {s["n"]:>9,} {s["summe"]:9.4f} {s["max"]:9.5f} {z["top1pct"]:7.1%} {top10:7.1%} | {pp} | {ll}')
    sag()
# Run time up to B, and a rough extrapolation of the local exponent to B = 10^7.
sag('Laufzeit bis B (s): ' + ' · '.join(f'{k}: {v}' for k, v in zeit.items()))
ks = sorted((int(k), v) for k, v in zeit.items() if v > 5)
if len(ks) >= 2:
    (b1, v1), (b2, v2) = ks[-2], ks[-1]; ex = math.log(v2/v1)/math.log(b2/b1)
    sag(f'   lokaler Exponent der Laufzeit (B {b1:,} → {b2:,}): {ex:.2f} ⇒ B = 10⁷ grob {v2*(10**7/b2)**ex/3600:.1f} h [Hochrechnung, nicht gemessen]')
json.dump(dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), BMAX=BMAX,
               periodenklassen=PER[:-1] + [[1001, 'inf']], lnepsklassen=LE[:-1] + [[1000, 'inf']],
               pk1=[dict(art=a, B=B, w121=x, kanonisch=y) for a, B, x, y in pk1], zerlegung=erg, laufzeit_bis_B=zeit),
          open(ERG/f'w121_zuwachs_zerlegung{TAG}_result.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
open(ERG/f'w121_zuwachs_zerlegung{TAG}_output.txt', 'w', encoding='utf-8').write('\n'.join(aus) + '\n')
sag(f'Ergebnis: w121_zuwachs_zerlegung{TAG}_result.json · Laufzeit {time.time() - t0:.0f} s')
