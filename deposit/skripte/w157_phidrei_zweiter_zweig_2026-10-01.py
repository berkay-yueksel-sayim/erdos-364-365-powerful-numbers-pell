# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w157_phidrei_zweiter_zweig_2026-10-01.py
# -*- coding: utf-8 -*-
# Purpose: Part III, Prop. III.2.2 (`phidrei` = Phi_3): the equation y^2 - 28m^2 = -3 has a SECOND BRANCH that the recursion does
#   not hit. Context: Lavi (Zenodo 2026, DOI 10.5281/zenodo.23043300) reports a PRIME x with x^2 + x + 1 = 7^3 * a^2
#   (both numbers are taken from the description of the record). Prop. III.2.2 builds the family from (y_0, m_0) = (37, 7) with
#   (y, m) -> (127y + 672m, 24y + 127m). The text says that the reduction of De Koninck-Doyon-Luca leads for D = 7 exactly to
#   y^2 - 28m^2 = -3 with 7 | m, but the recursion hits only ONE part of these solutions.
# Reads: nothing. Writes: ergebnisse/w157_phidrei_zweiter_zweig_result.json and ..._output.txt.
# Usage: python w157_phidrei_zweiter_zweig_2026-10-01.py [JMAX=200]
# Controls: PK, NK and the completeness check below; the run ends with an assertion that the total number of deviations is 0.
#
# What it does. (1) The branch through (5, 1) with the same map: y^2 - 28m^2 = -3, y odd, m_j = j + 1 (mod 7), so 7 | m_j exactly
#   for j = 6 (mod 7); for these j, x = (y - 1)/2 has Phi_3(x) = 343 (m_j/7)^2 powerful and Phi_6(x + 1) = Phi_3(x).
#   (2) The two branches are different: no y of one occurs in the other (j <= 200).
#   (3) Lavi's x is the term j = 13 of the branch through (5, 1), not a term of the family of Prop. III.2.2.
#   (4) For context: how many x of the two branches with 7 | m_j are probably prime for j <= JMAX (SymPy, BPSW; no primality
#   proof, since Lavi's primality proof is his own and is cited).
# Expectations: (1) holds for all j <= 200; (2) empty intersection; (3) hit at j = 13.
#   PK: family III.2.2 at j = 0 gives x = 18, Phi_3 = 343.
#   NK: with the wrong starting point (5, 2) (y^2 - 28m^2 = -87) the norm check fires; it runs through the same check `pruefe`.

import sys, json, time, pathlib, math
from sympy import isprime
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
sys.set_int_max_str_digits(0)
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent / 'ergebnisse'
JMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
# `aus` = output lines (also written to the output file); `sag` = say: print a line and keep it
aus = []
def sag(s=''): print(s, flush=True); aus.append(s)
LAVI_X = 47116128896261596331524418282573  # from the description of the record
LAVI_A = 2544031832644333879196802002161

# `zweig` = branch: the terms (j, y, m) for j = 0..J, generated from (y, m) by the map above
def zweig(y, m, J):
    z = []
    for j in range(J + 1):
        z.append((j, y, m)); y, m = 127 * y + 672 * m, 24 * y + 127 * m
    return z

# `pruefe` = check: counts the violations `verl` (norm, parity of y, residue of j mod 7 against `rest7`, Phi relations)
# and returns them with the hits `treffer` = (j, x) for the terms with 7 | m
def pruefe(z, rest7):
    verl = 0; treffer = []
    for j, y, m in z:
        if y * y - 28 * m * m != -3 or y % 2 == 0: verl += 1
        if (m % 7 == 0) != (j % 7 == rest7): verl += 1
        if m % 7 == 0:
            x = (y - 1) // 2; f3 = x * x + x + 1
            if f3 != 343 * (m // 7) ** 2 or (x + 1) ** 2 - (x + 1) + 1 != f3: verl += 1
            treffer.append((j, x))
    return verl, treffer

sag('=' * 100); sag(f'w157 — Prop. III.2.2: zweiter Zweig von y² − 28m² = −3 · JMAX = {JMAX} · {time.strftime("%Y-%m-%d %H:%M")}'); sag('=' * 100)
A = zweig(37, 7, JMAX); B = zweig(5, 1, JMAX)
assert (A[0][1] - 1) // 2 == 18 and 18 * 18 + 18 + 1 == 343
sag('PK ✅ Familie III.2.2, j = 0: x = 18, Φ_3(18) = 343')
# the NK runs through the SAME check `pruefe`, not through a formula of its own
nkv, _ = pruefe(zweig(5, 2, 20), 6); assert nkv > 0, 'NK: falscher Startpunkt nicht gemeldet'
sag(f'NK ✅ Startpunkt (5, 2) (y² − 28m² = −87): pruefe() meldet {nkv} Verletzungen')
# `vA`, `vB` = violations of the two branches; `tA`, `tB` = their hits
vA, tA = pruefe(A, 0); vB, tB = pruefe(B, 6)
sag(f'Zweig (37, 7) = Prop. III.2.2: Verletzungen {vA}; Glieder mit 7 | m_j: {len(tA)} (j ≡ 0 mod 7)')
sag(f'Zweig (5, 1):                 Verletzungen {vB}; Glieder mit 7 | m_j: {len(tB)} (j ≡ 6 mod 7); erstes: j = {tB[0][0]}, x mit {len(str(tB[0][1]))} Stellen')
# `schnitt` = intersection of the y values of the two branches
schnitt = {y for _, y, _ in A} & {y for _, y, _ in B}
sag(f'Schnitt der beiden Zweige (y-Werte, j ≤ {JMAX}): {len(schnitt)}')
# `lavi` / `lavi_A` = indices j at which Lavi's x occurs in the branch (5, 1) / in the family III.2.2
lavi = [j for j, x in tB if x == LAVI_X]; lavi_A = [j for j, x in tA if x == LAVI_X]
assert LAVI_X ** 2 + LAVI_X + 1 == 343 * LAVI_A ** 2
sag(f'Lavis x (32 Stellen): x² + x + 1 = 7³·a² ✅; im Zweig (5, 1) bei j = {lavi}; in der Familie III.2.2: {lavi_A or "nicht enthalten"}')
# `pA`, `pB` = indices j whose x is probably prime (BPSW test)
pA = [j for j, x in tA if isprime(x)]; pB = [j for j, x in tB if isprime(x)]
sag(f'Wahrscheinlich prime x (BPSW, kein Beweis), j ≤ {JMAX}: Zweig (37, 7): {pA} · Zweig (5, 1): {pB}')
# COMPLETENESS, second route. The proof in the text: Z[sqrt 7] has class number 1, eps = 8 + 3*sqrt 7 has norm +1, so all
# solutions of y^2 - 7z^2 = -3 have the form +-(2 +- sqrt 7)*eps^n, and z = 2m is even exactly for odd n: exactly the two
# branches. Independently here: ALL (y, m) with 1 <= m <= MBRUTE and y^2 = 28m^2 - 3 by direct search, and each must lie on
# one of the two branches. (`gefunden` = found, `ausserhalb` = outside the two branches)
MBRUTE = 10 ** 6
ya = {y for _, y, _ in A}; yb = {y for _, y, _ in B}
gefunden = [(math.isqrt(28 * m * m - 3), m) for m in range(1, MBRUTE + 1) if math.isqrt(28 * m * m - 3) ** 2 == 28 * m * m - 3]
ausserhalb = [(y, m) for y, m in gefunden if y not in ya and y not in yb]
assert (37, 7) in gefunden and (5, 1) in gefunden, 'PK: die Startpunkte fehlen in der direkten Suche'
sag(f'Vollstaendigkeit (direkte Suche, alle m ≤ ' + f'{MBRUTE:,}'.replace(',', ' ') + f'): {len(gefunden)} Loesungen, davon ausserhalb der zwei Zweige: {len(ausserhalb)} {ausserhalb}')
# `gesamt` = total number of deviations (must be 0)
gesamt = vA + vB + len(schnitt) + (0 if lavi == [13] else 1) + len(ausserhalb)
sag(f'\nGESAMT: Abweichungen {gesamt}')
(ERG / 'w157_phidrei_zweiter_zweig_result.json').write_text(json.dumps(dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'),
    jmax=JMAX, verletzungen=[vA, vB], glieder=[len(tA), len(tB)], schnitt=len(schnitt), lavi_j=lavi, bpsw_prim_A=pA, bpsw_prim_B=pB,
    mbrute=MBRUTE, direkt_gefunden=len(gefunden), direkt_ausserhalb=len(ausserhalb), gesamt=gesamt),
    indent=1, ensure_ascii=False), encoding='utf-8')
(ERG / 'w157_phidrei_zweiter_zweig_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
sag('Ergebnis: w157_phidrei_zweiter_zweig_result.json')
assert gesamt == 0
