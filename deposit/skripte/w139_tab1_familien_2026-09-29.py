# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w139_tab1_familien_2026-09-29.py
# -*- coding: utf-8 -*-
# Purpose: Table 1 of Part II: the families of the 39 pairs (n, n+1) of consecutive powerful numbers up to 3.888*10^21
#   (the 39 terms of OEIS A060355 in this range), each row cross-checked against THREE sources. A row: family (b, d), period,
#   log10 n_1, log10 F, number of pairs. No OEIS values appear in the table.
# Reads ONLY result files in ergebnisse/: w111_paar_familien_result.json (which family carries how many of the 39 pairs,
#   taken from the data; OEIS terms 27-39 only as factorizations), w112_quadrat_familien_result.json (families with a square:
#   period, ln F, birth t = ln n_1, PREDICTED count, own iteration), w115_walker_familien_D1000000_result.json (families
#   without a square: period M, ln F, t; the self-generated pairs up to 3.89*10^21 as exact numbers).
# Writes: ergebnisse/w139_tab1_familien_result.json. Usage: python w139_tab1_familien_2026-09-29.py   (no arguments)
# Controls (all abort on failure): E1 to E6, described below.
#
# Expectations and controls (`quadrat` = square):
#   E1  Every one of the 16 families of w111 appears in w112 (with a square, b = 1 or d = 1) or in w115 (without a square).
#   E2  With a square: the count predicted in w112 equals the count found in the data (w111), for EVERY family.
#   E3  Without a square: from the exact pairs of w115 (sqfree part of n and n+1, own trial division) the same family and the
#       same count as in w111.
#   E4  Sums: 39 pairs, of which 32 with a square and 7 without, as in w111.
#   E5  NK: the order (b, d) carries information (see below).
#   E6  The period column is Walker's quantity (see below).

import json, math, pathlib, sys
from collections import Counter
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'
# `J(n)` = load the result file `n` from `ERG`;  `sqfree` = squarefree part of n
def J(n): return json.load(open(ERG/n, encoding='utf-8'))
def sqfree(n):
    s, q = 1, 2
    while q*q <= n:
        e = 0
        while n % q == 0: n //= q; e += 1
        if e % 2: s *= q
        q += 1 if q == 2 else 2
    return s * n if n > 1 else s  # the remainder > 1 is a prime with exponent 1

# `mitQ` = families with a square (index of w112, key (b, d)); `ohneQ` = families without a square (index of w115);
# `fam` = the 16 families of w111; `fehlt` = missing ones
w111, w112, w115 = J('w111_paar_familien_result.json'), J('w112_quadrat_familien_result.json'), J('w115_walker_familien_D1000000_result.json')
mitQ = {(x['b'], x['d']): x for x in w112['geburten_erste']}
ohneQ = {(x['b'], x['d']): x for x in w115['frueheste']}
fam = w111['familien']
fehlt = [(f['b'], f['d']) for f in fam if (f['b'], f['d']) not in (mitQ if f['quadrat'] else ohneQ)]
print(f'PK E1 {"✅" if not fehlt else "❌"}  {len(fam)} Familien, alle in w112/w115 gefunden' + (f' — FEHLT: {fehlt}' if fehlt else '.'))
assert not fehlt and len(fam) == 16
# `e2` = families whose predicted and counted numbers differ
e2 = [(f['b'], f['d'], f['anzahl'], mitQ[(f['b'], f['d'])]['anzahl']) for f in fam if f['quadrat'] and f['anzahl'] != mitQ[(f['b'], f['d'])]['anzahl']]
print(f'PK E2 {"✅" if not e2 else "❌"}  mit Quadrat: vorhergesagt (w112) = gezaehlt (w111) fuer alle {sum(1 for f in fam if f["quadrat"])}' + (f' — ABWEICHUNG {e2}' if e2 else '.'))
assert not e2
# E3: sqfree(n) and sqfree(n+1) of the self-generated pairs (<= 22 digits). Trial division up to the root of the
# remainder suffices, because n = b*Y^2 with b | Y: the remainder after the small factors is small or a square.
# (`eigen` = own count, `soll` = expected count from w111)
eigen = Counter()
for s in w115['paare_ohne_quadrat_bis_XMAX']:
    n = int(s); b, d = sqfree(n), sqfree(n + 1)
    assert (n % b == 0 and math.isqrt(n // b)**2 == n // b and (n + 1) % d == 0 and math.isqrt((n + 1) // d)**2 == (n + 1) // d), n
    eigen[(b, d)] += 1
soll = {(f['b'], f['d']): f['anzahl'] for f in fam if not f['quadrat']}
print(f'PK E3 {"✅" if dict(eigen) == soll else "❌"}  ohne Quadrat: Familien und Anzahlen aus den exakten Paaren von w115 = w111: {dict(eigen)}')
assert dict(eigen) == soll
nQ = sum(f['anzahl'] for f in fam if f['quadrat']); nO = sum(f['anzahl'] for f in fam if not f['quadrat'])
e4 = (nQ + nO == len(w111['paare']) == 39 and nQ == w111['quadrat_paare'] == 32 and nO == w111['nicht_quadrat_paare'] == 7)
print(f'PK E4 {"✅" if e4 else "❌"}  {nQ} + {nO} = {nQ + nO} Paare.')
assert e4
# E5 (negative control): (1, 5) and (5, 1) are BOTH genuine families (n + 1 = 5*X^2 with n a square, resp.
# n = 5*Y^2 with n + 1 a square), so (d, b) must not be required to be absent. Where (d, b) is also a family,
# the two entries must carry DIFFERENT births t: an assignment that loses the order would be caught exactly here.
alle_fam = {**mitQ, **ohneQ}
paare_beide = [(f['b'], f['d']) for f in fam if (f['d'], f['b']) in alle_fam and f['b'] != f['d']]
e5 = bool(paare_beide) and all(alle_fam[(b, d)]['t'] != alle_fam[(d, b)]['t'] for b, d in paare_beide)
print(f'NK E5 {"✅" if e5 else "❌"}  Familien, deren Vertauschung auch eine ist: {paare_beide} — Geburten verschieden: {e5}')
assert e5

# E6: the column 'period' is meant to be Walker's quantity: b' (type (b,1)), g = d/gcd(d, x0) (type (1,d)),
# l = product of the odd primes of b*d that do not divide x0*y0 (type II, Walker Thm 3.2(3)).
# w112/w115 compute their period by their OWN route (iteration, resp. residue class via CRT). Here Walker's
# definition is applied directly to the smallest solution, for each of the 16 rows; a deviation aborts the run
# (the table label would then be wrong).
# `pell` = smallest solution (x, y) of x^2 - D*y^2 = 1 by continued fraction; `neg_pell` = smallest solution of
# y^2 - d*x^2 = -1, or None; `primteiler` = list of the distinct prime factors;
# `walker_periode` = the period by Walker's definition
def pell(D):
    a0 = math.isqrt(D); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - D*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (D - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
def neg_pell(d):
    # smallest solution of y^2 - d*x^2 = -1 via the continued fraction (convergents up to the end of the period), else None
    a0 = math.isqrt(d); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    for _ in range(10**6):
        if h0*h0 - d*k0*k0 == -1: return h0, k0
        if h0*h0 - d*k0*k0 == 1: return None
        Pp = a*Q - Pp; Q = (d - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
def primteiler(n):
    r, q = [], 2
    while q*q <= n:
        if n % q == 0:
            r.append(q)
            while n % q == 0: n //= q
        q += 1
    return r + ([n] if n > 1 else [])
def walker_periode(b, d):
    if d == 1:
        T1, U1 = pell(b); return math.prod(p for p in primteiler(b) if U1 % p)
    if b == 1:
        y0, x0 = neg_pell(d); return d // math.gcd(d, x0)
    x1, y1 = pell(b*d)
    X0, Y0 = math.isqrt((x1 + 1)//(2*d)), math.isqrt((x1 - 1)//(2*b))
    assert d*X0*X0 - b*Y0*Y0 == 1, (b, d)
    return math.prod(p for p in primteiler(b*d) if p % 2 and (X0*Y0) % p)
# `abw6` = deviations between Walker's period and the period stored in w112 (`m_period`) resp. w115 (`M`)
abw6 = [(f['b'], f['d'], walker_periode(f['b'], f['d']), (mitQ if f['quadrat'] else ohneQ)[(f['b'], f['d'])]['m_period' if f['quadrat'] else 'M'])
        for f in fam]
abw6 = [x for x in abw6 if x[2] != x[3]]
print(f'PK E6 {"✅" if not abw6 else "❌"}  Periode nach Walkers Definition = Periode aus w112/w115 fuer alle {len(fam)} Familien' + (f' — ABW. {abw6}' if abw6 else '.'))
assert not abw6

LN10 = math.log(10)
# `zeilen` = the table rows, sorted by birth t
zeilen = []
for f in sorted(fam, key=lambda f: ((mitQ if f['quadrat'] else ohneQ)[(f['b'], f['d'])]['t'])):
    q = (mitQ if f['quadrat'] else ohneQ)[(f['b'], f['d'])]
    zeilen.append(dict(b=f['b'], d=f['d'], quadrat=f['quadrat'], periode=q['m_period'] if f['quadrat'] else q['M'],
                       log10_n1=q['t'] / LN10, log10_F=q['lnF'] / LN10, anzahl=f['anzahl']))
for z in zeilen:
    print(f"  ({z['b']:>4}, {z['d']:>2})  Periode {z['periode']:>2}  log10 n1 {z['log10_n1']:6.2f}  log10 F {z['log10_F']:6.2f}  Paare {z['anzahl']:>2}")
res = dict(skript=pathlib.Path(__file__).name, datum='2026-09-29', familien=len(zeilen), paare=nQ + nO, mit_quadrat=nQ, ohne_quadrat=nO,
           grenze='n + 1 <= 3.888e21 (39 Glieder A060355, Johnson 2011 unzertifiziert vollstaendig)', zeilen=zeilen)
(ERG/'w139_tab1_familien_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
print('Ergebnis: w139_tab1_familien_result.json')
