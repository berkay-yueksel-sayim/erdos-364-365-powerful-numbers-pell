# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w113_familien_weit_2026-09-25.py
# -*- coding: utf-8 -*-
# Purpose: the square families up to b = 2*10^5 with the rank formula (partial sums of 1/ln F, counts of births), and the 25
#   kernels of the main problem as pair families (m, 1). Two types: (b, 1): n + 1 = X^2, n = b*Y^2; and (1, d): n = Y^2,
#   n + 1 = d*X^2.
# Reads: ergebnisse/w112_quadrat_familien_result.json, ergebnisse/w67_beweismenge_result.json, and the module
#   w114_a060355_quelle_2026-09-26.py (the list of terms, function `glieder`).
# Writes: ergebnisse/w113_familien_weit{TAG}_result.json and ..._output.txt (TAG empty for B2 = 200000, else _B<B2>).
# Usage: python w113_familien_weit_2026-09-25.py [B2]   (B2 = upper bound for b and d, default 200000)
# Controls: PK1 and PK2 against the stored results of w112 (below); the run aborts on a mismatch.
#
# Tool. w112 determined the period by iteration modulo b (effort ~b per b), too slow for b up to 2*10^5. Formulas used here:
#   (b, 1): n + 1 = X^2, n = b*Y^2:  m_b = b / gcd(b, Y_1)  (confirmed directly in w112 PK1 for all 6082 squarefree b <= 10^4).
#   (1, d): n = Y^2, n + 1 = d*X^2: negative Pell unit eta = y0 + x0*sqrt(d) (if it exists), solutions eta^n with n odd;
#           X_n = x0*U_n (Lucas sequence P = 2*y0, Q = -1). For odd p | d, p divides the discriminant 4*d*x0^2, so
#           U_n = n*y0^(n-1) (mod p), hence p | X_n <=> p | x0 or p | n; for p = 2, U_n is odd for odd n.
#           => g = d / gcd(d, x0); g even => NO family; otherwise n = g, 3g, 5g, ...
#           => birth t = g*ln(eps) - ln 4, growth ln F = 2g*ln(eps).
#   This (1, d) formula is derived here, hence the control PK1 against w112.
# PK1: the first births stored in w112 (direct iteration, b, d <= 10^4) are found again exactly, with type, birth, ln F, period
#   and count.
# PK2: the sum of 1/ln F at B = 10^4 equals 1.2235 as in w112.
# Questions. (A) How far does b <= 2*10^5 raise the lower bound liminf N(x)/ln x >= sum of 1/ln F?
#   (C) How do the partial sums keep growing per factor 10 in B (so far +0.29, +0.21, +0.12)?
#   (B) The 25 kernels of the main problem as pair families (m, 1).
# Expectations, stated in advance: the increment from 10^4 to 10^5 is about +0.07 [assumption: the increments shrink roughly
#   geometrically by a factor 0.6], carried by few b with small m_b (b | Y_1 or nearly). Births up to T: with more families
#   closer to the true rate; the local exponent lies between 2/3 and 1 [open].

import sys, json, math, time, pathlib
from math import isqrt, gcd, log
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'
# `aus` = output lines (also written to the output file); `sag` = say: print a line and keep it
aus = []
def sag(s=''):
    print(s, flush=True); aus.append(s)
t0 = time.time()
B2 = int(sys.argv[1]) if len(sys.argv) > 1 else 200_000
TAG = '' if B2 == 200_000 else f'_B{B2}'  # other bound -> separate files; the canonical ones (B2 = 2*10^5) stay untouched
sag('='*100); sag(f'w113 — QUADRAT-FAMILIEN BIS b = {B2:,} (Rang-Formel) · DIE 25 KERNE   ' + time.strftime('%Y-%m-%d %H:%M')); sag('='*100)
import importlib.util  # the terms of the sequence come from w114 (read from the b-file), not hard-coded in this script
_s = importlib.util.spec_from_file_location('w114', pathlib.Path(__file__).resolve().parent/'w114_a060355_quelle_2026-09-26.py')
w114 = importlib.util.module_from_spec(_s); _s.loader.exec_module(w114)
# `XMAX` = last listed term + 1, `LX` = ln XMAX
XMAX = w114.glieder()[-1] + 1; LX = log(XMAX)
# `pell` = returns ((x1, y1), neg): the fundamental solution of x^2 - m*y^2 = 1 and the first solution with value -1
# seen along the continued-fraction iteration (None if there is none)
def pell(m):
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1; neg = None
    while True:
        v = h0*h0 - m*k0*k0
        if v == -1 and neg is None: neg = (h0, k0)
        if v == 1: return (h0, k0), neg
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
def ln_eps(x, y, m):  # ln(x + y*sqrt(m)); for large x, y*sqrt(m) = sqrt(x^2 -+ 1) ~ x, hence ln(2x) (error < 1e-24)
    return log(x) + (log(2) if x > 10**12 else log(1 + y*math.sqrt(m)/x))
# `zaehle` = count: the exact number of terms with n + 1 <= XMAX; start/step = (X, Y) pairs used as units
def zaehle(start, schritt, m, art):
    (Xa, Ya), (Xp, Yp), anz = start, schritt, 0
    while (Xa*Xa if art == 'b1' else Xa*Xa + 1) <= XMAX:
        anz += 1; Xa, Ya = Xa*Xp + m*Ya*Yp, Xa*Yp + Ya*Xp
    return anz
# `potenz` = (X, Y) with (x + y*sqrt(m))^e = X + Y*sqrt(m)
def potenz(x, y, m, e):
    X, Y, bx, by = 1, 0, x, y
    while e:
        if e & 1: X, Y = X*bx + m*Y*by, X*by + Y*bx
        bx, by = bx*bx + m*by*by, 2*bx*by; e >>= 1
    return X, Y
# `sf` = sieve of squarefree numbers up to B2; `fam` = list of families, `melde` = next progress report
sf = bytearray([1]) * (B2 + 1)
for q in range(2, isqrt(B2) + 1): sf[q*q::q*q] = bytearray(len(range(q*q, B2 + 1, q*q)))
fam = []; melde = 25_000
for m in range(2, B2 + 1):
    if not sf[m]: continue
    # (b, 1) family for every squarefree b = m: `mb` = m_b (period), `lnF` = ln F, `t` = birth, `anz` = number of terms
    # up to XMAX (counted exactly when the birth lies in the range of the list)
    (x1, y1), neg = pell(m); le = ln_eps(x1, y1, m)
    mb = m // gcd(m, y1 % m); lnF = 2*mb*le; t = lnF - log(4); anz = 0
    if t <= LX + 5:
        s = potenz(x1, y1, m, mb); assert (s[0]*s[0] - 1) == m*s[1]*s[1] and s[1] % m == 0
        t = log(s[0]*s[0] - 1); anz = zaehle(s, s, m, 'b1')
    fam.append(dict(typ='(b,1)', b=m, d=1, m_period=mb, ln_eps=le, lnF=lnF, t=t, anzahl=anz))
    # (1, d) family for d = m if a negative Pell solution exists and g is odd
    if neg:
        y0, x0 = neg; g = m // gcd(m, x0 % m)
        if g % 2 == 1:
            t1 = g*le - log(4); lnF1 = 2*g*le; anz = 0
            if t1 <= LX + 5:
                Ya, Xa = potenz(y0, x0, m, g)  # eta^g = Y + X*sqrt(m)
                assert Ya*Ya - m*Xa*Xa == -1 and Xa % m == 0
                t1 = log(Ya*Ya); anz = zaehle((Ya, Xa), potenz(x1, y1, m, g), m, '1d')
            fam.append(dict(typ='(1,d)', b=1, d=m, m_period=g, j0=(g - 1)//2, ln_eps=le, lnF=lnF1, t=t1, anzahl=anz))
    if m >= melde: sag(f'   … b bis {m:,} gerechnet ({time.time()-t0:.0f}s)'); melde += 25_000
# ---- PK1/PK2 against w112 (`abw` = deviations)
w112 = json.load(open(ERG/'w112_quadrat_familien_result.json', encoding='utf-8'))
idx = {(f['typ'], f['b'], f['d']): f for f in fam}
abw = []
for g_ in w112['geburten_erste']:
    f = idx.get((g_['typ'], g_['b'], g_['d']))
    if not f or abs(f['t'] - g_['t']) > 1e-3 or abs(f['lnF'] - g_['lnF']) > 1e-3 or f['anzahl'] != g_['anzahl'] or f['m_period'] != g_['m_period']:
        abw.append((g_['typ'], g_['b'], g_['d'], g_['t'], f and round(f['t'], 4)))
sag(f'PK1 die {len(w112["geburten_erste"])} ersten Geburten aus w112 (direkte Iteration): ' + ('✅ alle exakt wieder (Typ, Geburt, ln F, Periode, Anzahl)' if not abw else f'🔴 {abw[:5]}'))
def summe(B): return sum(1/f['lnF'] for f in fam if max(f['b'], f['d']) <= B)
s4 = summe(10_000); ref = w112['summe_1_durch_lnF']['10000']
sag(f'PK2 Σ 1/ln F bei B = 10⁴: {s4:.6f} · w112: {ref:.6f} ' + ('✅' if abs(s4 - ref) < 1e-9 else '🔴'))
assert not abw and abs(s4 - ref) < 1e-9
sag()
# ---- (A)/(C) partial sums (`summe(B)` = sum of 1/ln F over families with b or d <= B)
Bs = [b for b in (10, 100, 1000, 10_000, 30_000, 100_000, 200_000, 1_000_000) if b <= B2]
ts = [summe(B) for B in Bs]
sag('(A)/(C) Σ 1/ln F ueber die Quadrat-Familien mit b bzw. d ≤ B  (= exakte Untergrenze fuer liminf N(x)/ln x, soweit gerechnet):')
sag('   ' + ' · '.join(f'B = {B:,}: {s:.4f}' for B, s in zip(Bs, ts)))
dek = [(Bs[i], Bs[i+1], ts[i+1] - ts[i]) for i in range(len(Bs) - 1)]
sag('   Zuwachs: ' + ' · '.join(f'{a:,} → {b:,}: +{d:.4f}' for a, b, d in dek))
# `neu` = the 12 largest contributions with b or d > 10^4
neu = sorted((f for f in fam if max(f['b'], f['d']) > 10_000), key=lambda f: -1/f['lnF'])[:12]
sag(f'   die groessten Beitraege mit b bzw. d > 10⁴ (1/ln F · Periode · ln ε):')
for f in neu:
    sag(f'     {f["typ"]} b = {f["b"]:>7,}, d = {f["d"]:>7,}: 1/ln F = {1/f["lnF"]:.5f} · Periode {f["m_period"]} · ln ε = {f["ln_eps"]:.2f} · Geburt t = {f["t"]:.1f}')
anteil_m1 = sum(1/f['lnF'] for f in fam if max(f['b'], f['d']) > 10_000 and f['m_period'] <= 3) / max(1e-12, ts[-1] - s4)
sag(f'   Anteil der Perioden ≤ 3 am Zuwachs ueber 10⁴: {anteil_m1:.1%}')
Ts = [25, 50, 75, 100, 150, 200, 300, 500, 1000, 2000]
nb = [sum(1 for f in fam if f['t'] <= T) for T in Ts]
sag('Geburten bis T (b bzw. d ≤ B2 — bei grossem T fehlen Familien mit groesserem b):  ' + ' · '.join(f'{T}: {n}' for T, n in zip(Ts, nb)))
sag('   lokale Exponenten: ' + ' · '.join(f'{(log(nb[i+1]) - log(nb[i]))/(log(Ts[i+1]) - log(Ts[i])):.2f}' for i in range(len(Ts) - 1) if nb[i] > 0))
# ---- (B) the 25 kernels (`KERNE`)
KERNE = json.load(open(ERG/'w67_beweismenge_result.json', encoding='utf-8'))['kerne']
sag()
sag('(B) die 25 Kerne aus #364 als Paar-Familien (m, 1): n + 1 = T_k², n = m·U_k², powerful ⟺ m | U_k ⟺ k ≡ 0 mod m/ggT(m, U₁):')
kr = []
for m in KERNE:
    f = idx.get(('(b,1)', m, 1))
    if f is None:  # kernels beyond B2 (209991, 4099215): treated individually with the same formula
        (x1, y1), _ = pell(m); le = ln_eps(x1, y1, m); mb = m // gcd(m, y1 % m)
        f = dict(m_period=mb, ln_eps=le, lnF=2*mb*le, t=2*mb*le - log(4), anzahl=0)
        assert f['t'] > LX + 5, 'Kern mit Paar im Listenbereich muesste exakt gezaehlt werden'
    kr.append(dict(m=m, periode=f['m_period'], ln_eps=round(f['ln_eps'], 4), geburt_t=round(f['t'], 3), geburt_log10=round(f['t']/log(10), 1), anzahl_bis_liste=f['anzahl']))
    sag(f'   m = {m:>6}: Periode {f["m_period"]:>6} · ln ε = {f["ln_eps"]:9.3f} · erstes Paar bei n ≈ 10^{f["t"]/log(10):,.1f}' + (f' · {f["anzahl"]} Paar(e) bis 3,9·10²¹' if f['anzahl'] else ''))
# result record: `teilsummen` = partial sums, `groesste_ueber_1e4` = largest contributions above 10^4, `geburten_bis_T` =
# births up to T, `kerne` = per kernel: `periode` = period, `geburt_t` = birth (ln scale), `geburt_log10` = birth (log10),
# `anzahl_bis_liste` = number of terms up to the end of the list
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), b2=B2, familien=len(fam),
           teilsummen={str(B): s for B, s in zip(Bs, ts)}, groesste_ueber_1e4=[{k: (round(v, 6) if isinstance(v, float) else v) for k, v in f.items()} for f in neu],
           geburten_bis_T=dict(zip(map(str, Ts), nb)), kerne=kr, laufzeit_s=round(time.time() - t0))
(ERG/f'w113_familien_weit{TAG}_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/f'w113_familien_weit{TAG}_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
sag(f'Laufzeit {time.time()-t0:.0f}s')
