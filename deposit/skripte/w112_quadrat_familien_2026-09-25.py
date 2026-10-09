# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w112_quadrat_familien_2026-09-25.py
# -*- coding: utf-8 -*-
# The SQUARE FAMILIES from the Pell equation: when is each family born, and how fast does it grow?
#
# Reads:    ../ergebnisse/w111_paar_familien_result.json (pair families of the data),
#   ../ergebnisse/w107_paar_korrektur_result.json
#           (key `C`, a correction factor of the random model), and w114_a060355_quelle_2026-09-26.py (imported as a module;
#           `glieder()` returns the terms of the list read from the b-file). Needs mpmath.
# Writes:   ../ergebnisse/w112_quadrat_familien_result.json and ../ergebnisse/w112_quadrat_familien_output.txt
# Usage:    python w112_quadrat_familien_2026-09-25.py   (no arguments)
# Controls: PK1 and PK2 (below), and an assertion that the random-model formula E(L) reproduces the value 50.44 at 10^21.
#
# CLAIM TO BE TESTED: in the family (b, 1) (n + 1 = X², n = b·Y², b squarefree), n is powerful if and only if b | Y; with
#   ε_b = X₁ + Y₁√b (fundamental solution of x² − b·y² = 1) one has b | Y_k ⟺ m_b | k, m_b = b / gcd(b, Y₁).
#   Then the growth factor per term is F_b = ε_b^(2·m_b), and the family is born at t_b = ln n₁ ≈ 2·m_b·ln ε_b − ln 4.
#   This is a CLAIM: here m_b is determined DIRECTLY by iteration modulo b and compared with the formula.
# REVERSE DIRECTION (1, d): n = Y² (a square), n + 1 = d·X², d | X  ⇔  Y² − d·X² = −1 (negative Pell equation) with d | X;
#   the solutions are the odd powers η^(2j+1) of the fundamental solution η of x² − d·y² = −1 (if it exists); sought: the
#   smallest j with d | X, and the period.
# CHECKS:
#   PK1  m_b determined directly = formula, for all b ≤ 10⁴.
#   PK2  every square family from w111 is refound with EXACTLY its count up to 3.9·10²¹ (for b, d ≤ 10⁴), and conversely: every
#        family predicted here with terms ≤ 3.9·10²¹ is also in the data. (This is at the same time a PARTIAL check of the
#        completeness of the OEIS list: all square pairs with squarefree part ≤ 10⁴.)
# QUESTION: does the birth rate fall (births up to T ~ T^(2/3)) or stay constant? And: does Σ 1/ln F converge?
#   Prediction, noted before the run: next births (b, 1): b = 15 at t ≈ 60.5; b = 21 at ≈ 64.4; b = 11 at ≈ 64.5.

import sys, json, math, time, pathlib
from math import isqrt, gcd, log
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'   # `HIER` = this directory, `ERG` = results directory
aus = []   # `aus` = lines of the output file
def sag(s=''):   # `sag` = print a line and record it for the output file
    print(s, flush=True); aus.append(s)
t0 = time.time()
sag('='*100); sag('w112 — QUADRAT-FAMILIEN: Geburt, Wachstum, Vergleich mit den 39 Paaren   ' + time.strftime('%Y-%m-%d %H:%M')); sag('='*100)
BMAX = 10_000   # bound for b and d
import importlib.util                      # the terms are obtained via w114 from the b-file, not hard-coded here
_s = importlib.util.spec_from_file_location('w114', pathlib.Path(__file__).resolve().parent/'w114_a060355_quelle_2026-09-26.py')
w114 = importlib.util.module_from_spec(_s); _s.loader.exec_module(w114)
XMAX = w114.glieder()[-1] + 1              # n + 1 ≤ XMAX (last term of the list + 1)
LX = log(XMAX)   # `LX` = ln XMAX
def pell(m):
    """smallest solutions of x² − m·y² = +1 and (if it exists) = −1, from the continued fraction of √m."""
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1; neg = None
    while True:
        v = h0*h0 - m*k0*k0
        if v == -1 and neg is None: neg = (h0, k0)
        if v == 1: return (h0, k0), neg
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
def lneps(x, y, m): return log(x) + log(1 + y*math.sqrt(m)/x)       # ln(x + y√m), stable for large x
sf = bytearray([1]) * (BMAX + 1)   # `sf` = squarefree sieve: sf[q] = 1 if q is squarefree (q <= BMAX)
for q in range(2, isqrt(BMAX) + 1): sf[q*q::q*q] = bytearray(len(range(q*q, BMAX + 1, q*q)))
w111 = json.load(open(ERG/'w111_paar_familien_result.json', encoding='utf-8'))
daten = {(f['b'], f['d']): f['anzahl'] for f in w111['familien']}   # `daten` = data: family (b, d) -> number of pairs (`anzahl`)
fam, pk1_fehler = [], []   # `fam` = predicted families (one dict each); `pk1_fehler` = failures of PK1
for m in range(2, BMAX + 1):
    if not sf[m]: continue
    (x1, y1), neg = pell(m); le = lneps(x1, y1, m)   # `le` = ln ε_m
    # ---- (m, 1): n = m·Y², n + 1 = X²; m | Y_k. Iteration modulo m.
    x, y, k = x1 % m, y1 % m, 1
    while y % m:
        x, y = (x*x1 + m*y*y1) % m, (x*y1 + y*x1) % m; k += 1
        if k > 2*m + 2: raise RuntimeError(('keine Periode', m))
    mb = k   # `mb` = m_b: least k with m | Y_k
    if mb != m // gcd(m, y1): pk1_fehler.append((m, mb, m // gcd(m, y1)))
    lnF = 2*mb*le; tb = lnF - log(4)   # `lnF` = ln F (growth factor per term), `tb` = birth time t_b
    anz = 0   # `anz` = exact number of terms up to XMAX
    if tb <= LX + 5:                         # count exactly: terms ε^(j·m_b), j = 1, 2, …
        Xs, Ys = 1, 0; Xm, Ym = 1, 0
        for _ in range(mb): Xm, Ym = Xm*x1 + m*Ym*y1, Xm*y1 + Ym*x1          # ε^(m_b)
        Xs, Ys = Xm, Ym
        n1 = Xs*Xs - 1; assert n1 == m*Ys*Ys and Ys % m == 0; tb = log(n1)
        while Xs*Xs <= XMAX:
            anz += 1; Xs, Ys = Xs*Xm + m*Ys*Ym, Xs*Ym + Ys*Xm
    fam.append(dict(typ='(b,1)', b=m, d=1, m_period=mb, ln_eps=le, lnF=lnF, t=tb, anzahl=anz))
    # ---- (1, m): n = Y², n + 1 = m·X², m | X; solutions η·ε^j of Y² − m·X² = −1.
    if neg:
        yn, xn = neg                          # yn² − m·xn² = −1
        Y, X, j, treffer = yn % m, xn % m, 0, []   # `treffer` = the first two exponents j with m | X
        while len(treffer) < 2 and j <= 2*m + 2:
            if X % m == 0: treffer.append(j)
            Y, X = (Y*x1 + m*X*y1) % m, (Y*y1 + X*x1) % m; j += 1
        if treffer:
            j0 = treffer[0]; per = treffer[1] - treffer[0] if len(treffer) > 1 else None   # `j0` = first exponent, `per` = period
            lnet = lneps(yn, xn, m)                           # ln η = ½ ln ε
            t1 = 2*(lnet + j0*le) - log(4); lnF1 = 2*per*le if per else float('inf')
            anz = 0
            if t1 <= LX + 5:
                Ya, Xa = yn, xn
                for _ in range(j0): Ya, Xa = Ya*x1 + m*Xa*y1, Ya*y1 + Xa*x1
                Xp, Yp = 1, 0
                for _ in range(per): Xp, Yp = Xp*x1 + m*Yp*y1, Xp*y1 + Yp*x1              # ε^per
                assert Ya*Ya - m*Xa*Xa == -1 and Xa % m == 0; t1 = log(Ya*Ya)
                while Ya*Ya + 1 <= XMAX:
                    anz += 1; Ya, Xa = Ya*Xp + m*Xa*Yp, Ya*Yp + Xa*Xp
            fam.append(dict(typ='(1,d)', b=1, d=m, m_period=per, j0=j0, ln_eps=le, lnF=lnF1, t=t1, anzahl=anz))
sag(f'PK1 m_b direkt (Iteration mod b) gegen die Formel b/ggT(b, Y₁), {sum(1 for f in fam if f["typ"] == "(b,1)")} quadratfreie b ≤ {BMAX}: '
    + ('✅ ueberall gleich' if not pk1_fehler else f'🔴 {len(pk1_fehler)} Abweichungen, z. B. {pk1_fehler[:5]}'))
# ---- PK2: comparison with the data
# `vorher` = predicted counts per family with at least one term
vorher = {(f['b'], f['d']): f['anzahl'] for f in fam if f['anzahl'] > 0}
quad_daten = {k: v for k, v in daten.items() if 1 in k}   # `quad_daten` = the square families (b = 1 or d = 1) in the data
fehlt_in_daten = {k: v for k, v in vorher.items() if quad_daten.get(k, 0) != v}   # predicted, but the data differ
fehlt_hier = {k: v for k, v in quad_daten.items() if max(k) <= BMAX and vorher.get(k, 0) != v}   # in the data, but the prediction differs
sag(f'PK2 Quadrat-Familien bis 3,9·10²¹: vorhergesagt {sum(vorher.values())} Paare in {len(vorher)} Familien · in den Daten {sum(quad_daten.values())} in '
    f'{len(quad_daten)} · Abweichungen: ' + ('keine ✅' if not fehlt_in_daten and not fehlt_hier else f'🔴 hier {fehlt_in_daten} · Daten {fehlt_hier}'))
sag()
# ---- births
geb = sorted(fam, key=lambda f: f['t'])   # `geb` = families sorted by birth time t
sag('Die ersten 24 Geburten (t = ln n₁; F = Wachstumsfaktor je Glied, ln F; Anzahl bis 3,9·10²¹; Periode m):')
for f in geb[:24]:
    sag(f'   {f["typ"]} b = {f["b"]:>4}, d = {f["d"]:>4}:  t = {f["t"]:7.2f} · ln F = {f["lnF"]:8.2f} · m = {f["m_period"]} · bis 3,9·10²¹: {f["anzahl"]}'
        + ('   ← jenseits der Liste' if f['t'] > LX else ''))
for bb, tp in ((15, 60.5), (21, 64.4), (11, 64.5)):
    f = next(f for f in fam if f['typ'] == '(b,1)' and f['b'] == bb)
    sag(f'   Vorhersage b = {bb}: t ≈ {tp} · hier {f["t"]:.2f} ' + ('✅' if abs(f['t'] - tp) < 0.3 else '🔴'))
Ts = [25, 50, 75, 100, 150, 200, 300, 500, 1000]   # `Ts` = thresholds T for the cumulative number of births
nb = [sum(1 for f in fam if f['t'] <= T) for T in Ts]   # `nb` = number of births up to each T
sag('Geburten bis T (Quadrat-Familien, b bzw. d ≤ 10⁴ — ab t ≈ 1000 fehlen Familien mit groesserem b, daher nur bis T = 1000 gezeigt):')
sag('   ' + ' · '.join(f'T = {T}: {n}' for T, n in zip(Ts, nb)))
st = [(log(nb[i+1]) - log(nb[i]))/(log(Ts[i+1]) - log(Ts[i])) for i in range(len(Ts) - 1) if nb[i] > 0]   # `st` = local exponents d ln N / d ln T
sag('   lokale Exponenten d ln N / d ln T: ' + ' · '.join(f'{s:.2f}' for s in st) + '  (konstante Rate = 1, Vorhersage ~0,67)')
# ---- Σ 1/ln F
teil = {}   # `teil` = partial sums of 1/ln F by bound B
for B in (10, 100, 1000, 10_000):
    teil[B] = sum(1/f['lnF'] for f in fam if max(f['b'], f['d']) <= B and f['lnF'] < float('inf'))
sag('Σ 1/ln F ueber die Quadrat-Familien mit b bzw. d ≤ B (= asymptotische Steigung dieser Familien je Einheit ln x):')
sag('   ' + ' · '.join(f'B = {B}: {v:.4f}' for B, v in teil.items()) + f'   (Zufallsmodell gesamt 0,8835)')
# ---- EXACT LOWER BOUND against the random models. Already at B = 100 the sum Σ 1/ln F exceeds 0.8835 (the random-model total),
#   which motivates this exact count. Every family f has exactly max(0, ⌊(L − t_f)/ln F_f⌋ + 1) terms up to x = e^L (up to
#   rounding at the boundary): NO heuristic. The square families with b, d ≤ 10⁴ are only a PART of all pairs ⇒ N_Q(L) is a LOWER
#   BOUND of the true count. Random model in closed form: E(L) = ∫_8^(e^L) ρ² dx = (c₁²/4)(L − ln 8) + 2c₁c₂(8^(−1/6) − e^(−L/6))
#   + (c₂²/3)(8^(−1/3) − e^(−L/3)).
import mpmath
c1 = float(mpmath.zeta(1.5)/mpmath.zeta(3)); c2 = float(mpmath.zeta(mpmath.mpf(2)/3)/mpmath.zeta(2))   # `c1` = ζ(3/2)/ζ(3), `c2` = ζ(2/3)/ζ(2)
C = json.load(open(ERG/'w107_paar_korrektur_result.json', encoding='utf-8'))['C']   # `C` = correction factor of the random model (from w107)
def E(L): return c1*c1/4*(L - log(8)) + 2*c1*c2*(8**(-1/6) - math.exp(-L/6)) + c2*c2/3*(8**(-1/3) - math.exp(-L/3))
assert abs(E(log(10**21)) - 50.44) < 0.01                  # PK: equals w110 at 10²¹ (numerically integrated there)
endl = [f for f in fam if f['lnF'] < float('inf')]   # `endl` = families with a finite growth factor
def NQ(L): return sum(max(0, math.floor((L - f['t'])/f['lnF']) + 1) for f in endl if f['t'] <= L)   # `NQ` = exact number of square pairs up to e^L
# `steig` = momentary slope Σ 1/ln F over the families born by L
def steig(L): return sum(1/f['lnF'] for f in endl if f['t'] <= L)
assert NQ(LX) == 32
sag('Exakte Anzahl der Quadrat-Paare (b, d ≤ 10⁴) bis x = e^L — eine UNTERGRENZE aller Paare — gegen die Zufallsmodelle fuer ALLE Paare:')
zeil = []   # `zeil` = table rows for the result file
for L in (49.7, 75, 100, 150, 200, 300, 500, 1000):
    zeil.append(dict(L=L, NQ=NQ(L), korr=round(C*E(L), 2), naiv=round(E(L), 2), steigung=round(steig(L), 4)))
    sag(f'   L = {L:>6} (x ≈ 10^{L/log(10):.0f}): Quadrat-Paare {NQ(L):>5} · Zufallsmodell korrigiert {C*E(L):8.2f} · naiv {E(L):8.2f} · '
        f'momentane Steigung der Quadrat-Familien {steig(L):.3f}')
# `kreuz_k` / `kreuz_n` = first L (in steps of 0.1) at which the square pairs alone exceed the corrected / the naive random model
kreuz_k = next((L/10 for L in range(497, 30001) if NQ(L/10) > C*E(L/10)), None)
kreuz_n = next((L/10 for L in range(497, 30001) if NQ(L/10) > E(L/10)), None)
sag(f'   ⇒ die Quadrat-Paare ALLEIN ueberholen die korrigierte Vorhersage fuer alle Paare ab L ≈ {kreuz_k} (x ≈ 10^{kreuz_k/log(10):.0f})'
    if kreuz_k else '   ⇒ kein Ueberholen bis L = 3000')
sag(f'   ⇒ … und die NAIVE Vorhersage ab L ≈ {kreuz_n} (x ≈ 10^{kreuz_n/log(10):.0f})' if kreuz_n else '   ⇒ die naive Vorhersage wird bis L = 3000 nicht ueberholt')
res_unter = dict(zeilen=zeil, ueberholt_korrigiert_ab_L=kreuz_k, ueberholt_naiv_ab_L=kreuz_n)   # `res_unter` = lower-bound block of the result file
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), bmax=BMAX, pk1_abweichungen=pk1_fehler, untergrenze=res_unter,
           pk2=dict(vorhergesagt=sum(vorher.values()), daten=sum(quad_daten.values()), abw_hier=str(fehlt_in_daten), abw_daten=str(fehlt_hier)),
           geburten_erste=[{k: (round(v, 4) if isinstance(v, float) else v) for k, v in f.items()} for f in geb[:60]],
           geburten_bis_T=dict(zip(map(str, Ts), nb)), summe_1_durch_lnF={str(k): v for k, v in teil.items()}, laufzeit_s=round(time.time() - t0))
(ERG/'w112_quadrat_familien_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w112_quadrat_familien_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
sag(f'Laufzeit {time.time()-t0:.0f}s')
