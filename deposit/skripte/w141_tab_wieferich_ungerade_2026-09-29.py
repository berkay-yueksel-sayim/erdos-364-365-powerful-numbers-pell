# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w141_tab_wieferich_ungerade_2026-09-29.py
# -*- coding: utf-8 -*-
# Part I, Prop. 6.11: table of the Wieferich events of odd rank d >= 5 on units with T_1 even.
# Purpose: list EVERY event (m, d, p) in the searched range (the hypothesis of 6.11 holds only for T_1 even) together with
#   the status of the block M_d: a witness exists (a prime of rank d with v = 1 below 3*10^6), the block is completely
#   factored, or it is open.
# Reads (result files only): ergebnisse/w137_wieferich_79kerne_result.json (events for the 79 kernels, witnesses with v = 1 below
#   3*10^6 per rank) and ergebnisse/w131_wieferich_zweite_route_result.json (complete factorization of M_5(319)).
# Writes: ergebnisse/w141_tab_wieferich_result.json. Usage: python w141_tab_wieferich_ungerade_2026-09-29.py (no arguments).
# The parity of T_1 is computed here from the Pell solution itself (own continued fraction), not taken from labels in other files.
# Controls (an assert aborts the run if one fails):
#   E1 (PK)  On the 25 old kernels exactly the five triples (m, d, p): (15,45,181), (31,39,157), (79,253,4049), (143,13,311),
#            (319,5,20161); (327,417,1667) drops out because T_1 is odd.
#   E2       Every row has p = +-1 (mod 4d) and v = 2 (w137); no row with v >= 3 (the three cases v = 3 occur at d = 1).
#   E3 (PK)  Status `zeuge` (witness): at least one prime of rank d with v = 1 below 3*10^6 (w137, uncapped). Status `vollzerlegt`
#            (completely factored): only (319, 5); w131 must give 20161^2 * 1375502561 * 158585313121 there (besides the
#            intrinsic 5).
#            Otherwise `offen` (open).
#   E4 (NK)  A kernel with T_1 odd (327) is recognized as odd; m = 7 (T_1 = 8) as even.
import json, math, pathlib, sys
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'   # `ERG` = results folder (`ergebnisse`)
# `fund(m)` = T_1 (first coordinate) of the fundamental solution of x^2 - m*y^2 = 1, from the continued fraction of sqrt(m):
# (Pp, Q, a) = state of the expansion, h = convergent numerators, k = convergent denominators.
def fund(m):
    a0 = math.isqrt(m); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (m - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0
e4 = fund(327) % 2 == 1 and fund(7) == 8
print(f'NK E4 {"✅" if e4 else "❌"}  T₁(327) ungerade, T₁(7) = 8.')
assert e4
# `w137`/`w131` = the two result files; `ereignisse` = events, `bloecke_reichweite` = blocks in range
w137 = json.load(open(ERG/'w137_wieferich_79kerne_result.json', encoding='utf-8'))
w131 = json.load(open(ERG/'w131_wieferich_zweite_route_result.json', encoding='utf-8'))
# `ev` = events of odd rank d >= 5; `t1` = kernel m -> True if T_1(m) is even; `zeilen` = rows of the table
ev = [e for e in w137['ereignisse'] if e['d'] >= 5 and e['d'] % 2 == 1]
t1 = {m: fund(m) % 2 == 0 for m in {e['m'] for e in ev}}
zeilen = []
for e in sorted(ev, key=lambda e: (e['m'], e['d'])):
    if not t1[e['m']]: continue
    assert e['p'] % (4*e['d']) in (1, 4*e['d'] - 1) and e['v'] == 2, e
    # status: `zeuge` = witness below 3e6 exists, `voll` = block completely factored (`weg` = route in w131), `offen` = open
    if e['zeugen_unter_3e6'] > 0: status = 'zeuge'
    else:
        b = [x for x in w131['bloecke_reichweite'] if (x['m'], x['d']) == (e['m'], e['d'])]
        status = 'voll' if b and b[0]['weg'].startswith('vollstaendige Zerlegung') else 'offen'
    # row: kernel `m`, rank `d`, prime `p`, `kern` = kernel origin ('alt' = old kernels), `zeugen` = number of witnesses below
    # 3e6, `status`, `frage` = field `in_frage_population` copied from w137
    zeilen.append(dict(m=e['m'], d=e['d'], p=e['p'], kern=e['kern'], zeugen=e['zeugen_unter_3e6'], status=status,
                       frage=e['in_frage_population']))
# `kern` = 'alt' marks the 25 old kernels; `alt5` = their events (m, d, p)
alt5 = {(z['m'], z['d'], z['p']) for z in zeilen if z['kern'] == 'alt'}
e1 = alt5 == {(15, 45, 181), (31, 39, 157), (79, 253, 4049), (143, 13, 311), (319, 5, 20161)}
print(f'PK E1 {"✅" if e1 else "❌"}  alte Kerne: {sorted(alt5)}')
assert e1
# `b319` = factorization of the block M_5(319) from w131 (`zerlegung`: prime -> exponent)
b319 = [x for x in w131['bloecke_reichweite'] if (x['m'], x['d']) == (319, 5)][0]['zerlegung']
e3 = ({int(k): v for k, v in b319.items()} == {5: 1, 20161: 2, 1375502561: 1, 158585313121: 1}
      and [z['status'] for z in zeilen if z['status'] == 'voll'] == ['voll'] and all(z['zeugen'] > 0 for z in zeilen if z['status'] == 'zeuge'))
print(f'PK E3 {"✅" if e3 else "❌"}  M₅(319) = 20161² · 1375502561 · 158585313121 (w131), einziger Fall „vollzerlegt".')
assert e3
from collections import Counter
# `st` = number of rows per status
st = Counter(z['status'] for z in zeilen)
print(f'   {len(zeilen)} Zeilen: Zeuge {st["zeuge"]} · vollzerlegt {st["voll"]} · offen {st["offen"]} — offen: '
      f'{[(z["m"], z["d"], z["p"]) for z in zeilen if z["status"] == "offen"]}')
# JSON keys: `skript` = script, `datum` = date, `ereignisse` = number of events, `verschiedene_p` = number of distinct primes p,
# `ungerade_rang` = events of odd rank, `t1_gerade` = those with T_1 even, `zeuge` / `voll` / `offen` = rows per status,
# `zeilen` = rows, `m319_zerlegung` = factorization of the block M_5(319)
res = dict(skript=pathlib.Path(__file__).name, datum='2026-09-29', ereignisse=len(w137['ereignisse']),
           verschiedene_p=len({e['p'] for e in w137['ereignisse']}), ungerade_rang=len(ev), t1_gerade=len(zeilen),
           zeuge=st['zeuge'], voll=st['voll'], offen=st['offen'], zeilen=zeilen, m319_zerlegung=b319)
(ERG/'w141_tab_wieferich_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
print('Ergebnis: w141_tab_wieferich_result.json')
