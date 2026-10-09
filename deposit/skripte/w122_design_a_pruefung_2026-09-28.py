# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w122_design_a_pruefung_2026-09-28.py
# -*- coding: utf-8 -*-
# Check of "design A" before the search run: do constructed kernels with a small fundamental unit carry candidates?
# Idea under test: kernels m = t^2 + r with small |r| have tiny fundamental units, so that kernels far beyond 10^15 would become
#   reachable in the same height window. Objection: a triple with kernel m requires m' | k with m' = m / gcd(m, U_1); a small unit
#   eps = T_1 + U_1*sqrt(m) > 2*U_1*sqrt(m) means a small U_1, hence a large m', and Theorem 2.1 gives
#   log10 T_k >= m'*(1.5*log10 m - log10 m') - log10 2. The lever "small T_1" only acts when m' = 1.
# Check on the two families of kernels m = 7 (mod 8): m = t^2 - 2 (t odd) and m = t^2 - 1 (t = 0 mod 4), squarefree, t < 400:
#   (1) the fundamental unit is the expected one ((t^2-1) + t*sqrt(m) resp. t + sqrt(m)), computed by continued fraction
#       (function `pell` taken from w113);
#   (2) m' = m / gcd(m, U_1), expected m' = m;
#   (3) the first triple candidate T_{m'}: exact digit count (for m' <= 3000) against the Theorem 2.1 bound;
#   (4) PK: the kernels 7, 23, 47, 79, 119 (= t^2 - 2) have continued-fraction period 7, 23, 47, 79, 119 in w113 (B),
#       and m' = m here.
# Question answered: up to which m does the first candidate lie below 10^2000 (the height window of the search, kernels <= 10^8)?
# Reads: w113_familien_weit_2026-09-25.py (only the function `pell`, extracted from the syntax tree) and
#   ergebnisse/w113_familien_weit_B1000000_result.json.
# Writes: ergebnisse/w122_design_a_pruefung_result.json and ergebnisse/w122_design_a_pruefung_output.txt. No arguments.
# Controls (an assert aborts the run if one fails): checks (1) to (4) above.
import ast, json, math, sys, pathlib, time
from math import isqrt, gcd, log, log10
sys.stdout.reconfigure(encoding='utf-8')
sys.set_int_max_str_digits(0)                   # T_(m') has up to ~11,000 digits; the default str conversion limit is 4300
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'   # `ERG` = `ergebnisse` (results)
# load only the function `pell` of w113 (returns the fundamental solution (T_1, U_1) first) into the namespace `ns`
ns = {'isqrt': isqrt, 'gcd': gcd, 'log': log, 'math': math}
for k in ast.parse((HIER/'w113_familien_weit_2026-09-25.py').read_text(encoding='utf-8')).body:
    if isinstance(k, ast.FunctionDef) and k.name == 'pell': exec(compile(ast.Module([k], []), 'w113', 'exec'), ns)
pell = ns['pell']
aus = []                                        # `aus` = output lines (also written to the output file)
def sag(s=''):                                  # `sag` = say: print a line and record it
    print(s, flush=True); aus.append(s)
def quadratfrei(n):                             # `quadratfrei` = squarefree test
    return all(n % (p*p) for p in range(2, isqrt(n) + 1))
# `T_potenz` = T_e of the unit x + y*sqrt(m), exact, by square-and-multiply on the pair (X, Y)
def T_potenz(x, y, m, e):
    X, Y, bx, by = 1, 0, x, y
    while e:
        if e & 1: X, Y = X*bx + m*Y*by, X*by + Y*bx
        bx, by = bx*bx + m*by*by, 2*bx*by; e >>= 1
    return X
sag('='*100); sag('w122 — DESIGN A VOR DEM LAUF: Kerne mit kleiner Einheit, Rang m′ und erste Tripel-Kandidatin   ' + time.strftime('%Y-%m-%d %H:%M')); sag('='*100)
# `zeilen` = table rows; `FENSTER` = height window in digits; `art` = family; `t_iter` = range of t;
# `erwart` = expected unit (T_1, U_1) as a function of t
zeilen = []; FENSTER = 2000
for art, t_iter, erwart in (('t²−2', range(3, 400, 2), lambda t: (t*t - 1, t)), ('t²−1', range(4, 400, 4), lambda t: (t, 1))):
    for t in t_iter:
        m = t*t - 2 if art == 't²−2' else t*t - 1
        if m % 8 != 7 or not quadratfrei(m): continue
        (x1, y1), _ = pell(m)
        ok_einheit = (x1, y1) == erwart(t)
        # `mp` = m'; `schranke` = Theorem 2.1 bound on log10 T_{m'}; `exakt` = exact digit count of T_{m'};
        # `naeherung` = approximation of log10 T_{m'} from the unit (T_{m'} is about eps^m' / 2)
        mp = m // gcd(m, y1)
        schranke = mp*(1.5*log10(m) - log10(mp)) - log10(2)
        exakt = len(str(T_potenz(x1, y1, m, mp))) if mp <= 3000 else None
        naeherung = mp*log10(x1 + y1*math.sqrt(m)) - log10(2)
        # row keys: `einheit_wie_erwartet` = unit as expected, `m_strich` = m', `stellen_T_mstrich` = digits of T_{m'},
        # `log10_naeherung` = approximation of log10 T_{m'}, `thm21` = Theorem 2.1 bound
        zeilen.append(dict(art=art, t=t, m=m, T1=x1, U1=y1, einheit_wie_erwartet=ok_einheit, m_strich=mp,
                           stellen_T_mstrich=exakt, log10_naeherung=round(naeherung, 1), thm21=round(schranke, 1)))
# `ok_e`, `ok_m`, `ok_s` = checks (1), (2), (3) passed for all rows
ok_e = all(z['einheit_wie_erwartet'] for z in zeilen); ok_m = all(z['m_strich'] == z['m'] for z in zeilen)
ok_s = all(z['stellen_T_mstrich'] is None or z['stellen_T_mstrich'] - 1 >= z['thm21'] for z in zeilen)
sag(f'(1) Grundeinheit wie erwartet an allen {len(zeilen)} Kernen: {"✅" if ok_e else "🔴"}')
sag(f'(2) m′ = m an allen: {"✅" if ok_m else "🔴 " + str([z["m"] for z in zeilen if z["m_strich"] != z["m"]][:8])}')
sag(f'(3) exakte Stellenzahl von T_(m′) ≥ Theorem-2.1-Schranke, wo exakt gerechnet: {"✅" if ok_s else "🔴"}')
for art in ('t²−2', 't²−1'):
    zs = [z for z in zeilen if z['art'] == art]
    im = [z for z in zs if z['log10_naeherung'] < FENSTER]   # `im` = rows whose first candidate lies inside the window
    sag(f'   {art}: {len(zs)} Kerne (t < 400, quadratfrei, ≡ 7 mod 8); erste Kandidatin unter 10^{FENSTER} nur bis m = {max((z["m"] for z in im), default=0):,} '
        f'({len(im)} Kerne) — groesster Kern hier m = {zs[-1]["m"]:,}: erste Kandidatin ≈ 10^{zs[-1]["log10_naeherung"]:,.0f}')
    sag('     ' + ' · '.join(f'm={z["m"]}: 10^{z["log10_naeherung"]:.0f}' for z in zs[:8]))
w113 = json.load(open(ERG/'w113_familien_weit_B1000000_result.json', encoding='utf-8'))
per = {k['m']: k['periode'] for k in w113['kerne']}   # `per` = period of the continued fraction per kernel (check (4))
pk = {m: (per.get(m), next((z['m_strich'] for z in zeilen if z['m'] == m), None)) for m in (7, 23, 47, 79, 119)}
ok_pk = all(a == b == m for m, (a, b) in pk.items())
sag(f'(4) PK gegen w113 (Periode der Kerne 7, 23, 47, 79, 119): {pk} {"✅" if ok_pk else "🔴"}')
sag()
# `mmax` = largest kernel whose first candidate lies inside the window
mmax = max(z['m'] for z in zeilen if z['log10_naeherung'] < FENSTER)
sag(f'FOLGERUNG (nur fuer diese Familien, gemessen): ihre Kandidaten unter 10^{FENSTER} liegen alle bei m ≤ {mmax:,}, also in der Karte (Kerne ≤ 10⁸);')
sag('   darueber beginnt die erste Kandidatin weit jenseits jedes Rechenfensters. Kleine Einheit bringt hier keine neuen Kerne ins Fenster.')
assert ok_e and ok_m and ok_s and ok_pk
# JSON keys: `skript` = script, `datum` = date, `fenster_stellen` = window in digits, `kerne` = kernels (rows),
# `pk_w113` = PK against w113
json.dump(dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), fenster_stellen=FENSTER, kerne=zeilen,
               pk_w113=pk), open(ERG/'w122_design_a_pruefung_result.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
open(ERG/'w122_design_a_pruefung_output.txt', 'w', encoding='utf-8').write('\n'.join(aus) + '\n')
