# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# -*- coding: utf-8 -*-
# w110_paare_bis_1e21_2026-09-25.py
# Compares the pair heuristic of w107 with the known pairs n, n+1 of consecutive powerful numbers up to 3.9*10^21.
# Background: w107 counts pairs up to 2^52 (26 pairs; expected 26.87 with the correction factor C, 35.92 naively). The
#   counting question is Erdős problem #365 (is the number of n <= x with n, n+1 both powerful bounded by (log x)^O(1)?). The
#   heuristic predicts more sharply: count ~ C*(c1^2/4)*ln x = 0.8835*ln x (0.6124 per binary section, i.e. per doubling of x).
# Data (external, outputs comparison only): OEIS A060355, b-file, terms 1 ... 39. The OEIS data are CC BY-SA 4.0, so they are
#   only cited: the values are read at run time via w114 from the reader's own copy of the b-file. Terms 1 ... 26 were confirmed
#   independently by our own enumeration (w107); terms 27 ... 39 are taken over unchecked [assumption: the b-file is complete up
#   to term 39].
# Reads: ergebnisse/w107_paar_korrektur_result.json and, via import, the sibling script w114_a060355_quelle_2026-09-26.py
#   (function `glieder()` returns the terms, `QUELLE` names the source).
# Writes: ergebnisse/w110_paare_bis_1e21_result.json and ergebnisse/w110_paare_bis_1e21_output.txt. No command-line arguments.
# Expectation (rough estimate from the w107 numbers): up to 3.9*10^21 (= 2^71.7) corrected ~ 26.87 + 19.7*0.6124 ~ 39, naive ~
#   52. Whether the measured count is closer to the corrected value is the test [assumption]; Poisson noise +-sqrt(39) ~ +-6.
# Controls (asserted): PK: our 26 pairs from w107 equal terms 1 ... 26; for terms 27 ... 39, n and n+1 are both powerful
#   (trial division). Two cross-checks follow the table (see below). There is no negative control.

import sys, json, math, time, pathlib
import mpmath
sys.stdout.reconfigure(encoding='utf-8')
ERG = pathlib.Path(__file__).resolve().parent.parent/'ergebnisse'
# `aus` = collected output lines; `sag` (= say) prints a line and records it for the _output.txt file.
aus = []
def sag(s=''):
    print(s, flush=True); aus.append(s)
sag('='*100); sag('w110 — PAAR-HEURISTIK GEGEN OEIS A060355 BIS 3,9·10²¹   ' + time.strftime('%Y-%m-%d %H:%M')); sag('='*100)
import importlib.util                      # load the terms via w114 from the b-file (they are not hard-coded here)
_s = importlib.util.spec_from_file_location('w114', pathlib.Path(__file__).resolve().parent/'w114_a060355_quelle_2026-09-26.py')
w114 = importlib.util.module_from_spec(_s); _s.loader.exec_module(w114)
A060355 = w114.glieder()
w107 = json.load(open(ERG/'w107_paar_korrektur_result.json', encoding='utf-8'))
# `eigene` = our own pairs (terms 1 ... 26) from w107.
eigene = [int(s) for s in w107['paare']]
assert A060355[:len(eigene)] == eigene, 'Glieder 1 … 26 weichen von w107 ab'
# spot check: every term and its successor must be powerful (trial division up to n^(1/3))
def ist_powerful(n):
    m, p = n, 2
    while p*p*p <= m:
        if m % p == 0:
            e = 0
            while m % p == 0: m //= p; e += 1
            if e == 1: return False
        p += 1 if p == 2 else 2
    # the remainder m has no prime factor <= m^(1/3): m = 1, q^2, q or q*r; only m = 1 and q^2 are powerful
    r = mpmath.mpf(m)
    s = int(mpmath.sqrt(r)); return m == 1 or s*s == m
# `pruef` = checks: (term, n powerful, n+1 powerful) for the terms 27 ... 39.
pruef = [(n, ist_powerful(n), ist_powerful(n + 1)) for n in A060355[26:]]
sag(f'PK  die 26 eigenen Paare = Glieder 1 … 26 ✅ · Glieder 27 … 39: n und n+1 powerful: ' + ' '.join('✅' if a and b else '🔴' for _, a, b in pruef))
assert all(a and b for _, a, b in pruef)
# c1, c2: constants of the count c1*sqrt(x) + c2*x^(1/3) of powerful numbers <= x; `rho` = its derivative (density of
# powerful numbers near x); `int_rho2(A, B)` = integral of rho^2 over [A, B] (midpoint rule in log x, n steps) =
# naive expected number of pairs; `C` = correction factor from w107.
c1 = float(mpmath.zeta(1.5)/mpmath.zeta(3)); c2 = float(mpmath.zeta(mpmath.mpf(2)/3)/mpmath.zeta(2)); C = w107['C']
def rho(x): return c1/(2*math.sqrt(x)) + c2/(3*x**(2/3))
def int_rho2(A, B, n=8000):
    la, lb = math.log(A), math.log(B); h = (lb - la)/n; s = 0.0
    for i in range(n):
        x = math.exp(la + (i + 0.5)*h); r = max(rho(x), 0.0); s += r*r*x*h
    return s
# `zeilen` = table rows per limit X: `gemessen` = measured pairs with n+1 <= X, `naiv` = integral of rho^2 from 8 to X,
# `korrigiert` = C * naiv, `z_naiv` / `z_korr` = deviations in Poisson standard deviations.
zeilen = []
for X in [10**k for k in range(8, 22)] + [A060355[-1] + 1]:
    gem = sum(1 for n in A060355 if n + 1 <= X); e = int_rho2(8, X); ek = C*e
    zeilen.append(dict(X=(str(X) if X <= 10**21 else f'{X:.4e}'), gemessen=gem,   # last limit only rounded (term 39 + 1)
                        # stored with 6 places; the final rounding is left to the consumer
                        naiv=round(e, 6), korrigiert=round(ek, 6), z_naiv=round((gem - e)/math.sqrt(e), 4),
                       # (rounding earlier would turn 16.125 into 16.13)
                       z_korr=round((gem - ek)/math.sqrt(ek), 4)))
    sag(f'   bis {("10^" + str(round(math.log10(X)))) if X in [10**k for k in range(8, 22)] else "3,888·10²¹":<11}: gemessen {gem:>2} · naiv {e:6.2f} '
        f'(z {(gem - e)/math.sqrt(e):+.2f}) · korrigiert {ek:6.2f} (z {(gem - ek)/math.sqrt(ek):+.2f})')
g = zeilen[-1]
# Cross-checks that count nothing twice: (a) only the part the heuristic has never seen, terms 27 ... 39, i.e. (2^52,
#   3.9*10^21]; (b) slope WITH a free intercept from 10^8 on (a fit through the origin ignores the intercept of the prediction
#   (about -5) and gives 0.73 instead of ~0.88, a spurious contradiction).
X0, X1 = 2**52, A060355[-1] + 1
neu = sum(1 for n in A060355 if X0 < n + 1 <= X1); e_n = int_rho2(X0, X1); ek_n = C*e_n
sag(f'Gegenprobe (a) NEUER Teil (2⁵², 3,9·10²¹]: gemessen {neu} · korrigiert {ek_n:.2f} (z {(neu - ek_n)/math.sqrt(ek_n):+.2f}) · '
    f'naiv {e_n:.2f} (z {(neu - e_n)/math.sqrt(e_n):+.2f})')
# `pts` = (ln(n+1), index) for the terms with n+1 >= 10^8; `a_fit` = least-squares slope of index against ln x.
pts = [(math.log(n + 1), k) for k, n in enumerate(A060355, 1) if n + 1 >= 10**8]
mx = sum(l for l, _ in pts)/len(pts); my = sum(k for _, k in pts)/len(pts)
a_fit = sum((l - mx)*(k - my) for l, k in pts)/sum((l - mx)**2 for l, _ in pts)
sag(f'Gegenprobe (b) Steigung mit freiem Achsenabschnitt ab 10⁸ ({len(pts)} Glieder): a = {a_fit:.4f} je Einheit ln x · Heuristik C·c₁²/4 = '
    f'{C*c1*c1/4:.4f} · naiv c₁²/4 = {c1*c1/4:.4f}  [kumulative Punkte sind korreliert — nur grob]')
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), quelle=w114.QUELLE,
           C=C, vergleich=zeilen, neuer_teil=dict(gemessen=neu, korrigiert=round(ek_n, 6), naiv=round(e_n, 6)),
           steigung_frei_ab_1e8=a_fit, steigung_heuristik=C*c1*c1/4, steigung_naiv=c1*c1/4)
(ERG/'w110_paare_bis_1e21_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w110_paare_bis_1e21_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
