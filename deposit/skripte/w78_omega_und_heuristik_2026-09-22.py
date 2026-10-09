# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w78_omega_und_heuristik_2026-09-22.py
# Purpose: (④) distribution of ω(Φ_d(m)), the number of distinct prime factors of the cyclotomic value Φ_d(m), over the
#   completely factored cases (Φ_d(m) written Φ_d below);
#   (⑤) the heuristic expected number of powerful Φ_d, as a number; (③) a test of whether a theorem of Juricevic transfers.
# Reads:  ergebnisse/w30_vollfaktorisierung_S55_N25_result.json (158 completely factored Φ_d, threshold 55) and
#         ergebnisse/w70_zeitlimit_faelle_result.json (5 further cases): 163 cases in total, keyed by (kernel m, rank d).
# Writes: ergebnisse/w78_omega_heuristik_result.json and ergebnisse/w78_omega_heuristik_output.txt.
# Usage:  python w78_omega_und_heuristik_2026-09-22.py   (no arguments)
# ④ Distribution of ω(Φ_d). How often is Φ_d prime? A prime power?
# ⑤ Heuristic: P(Φ_d powerful) ≈ Π_{p | Φ_d} 1/p (each prime factor is Wieferich independently with probability 1/p), so the
#   expected number of powerful cases among the 163 is Σ Π 1/p. This is labeled a heuristic throughout, not a proof.
# ③ Test: Juricevic 2008, Thm 2.1 gives at least 2 primitive divisors of u_n for real pairs and n ≥ 31 (n/(ηκ) odd). Our rank
#   d is the U-rank 2d, so a transfer would require ω(Φ_d) ≥ 2 for 2d ≥ 31, i.e. d ≥ 16. If some case has ω = 1 with d ≥ 16,
#   the transfer does not apply (conditions (1.12), η, κ exclude our pair, or the rank correspondence is wrong). This is
#   measured, not assumed.
# Controls (expectations stated before the run):
#   E1  exactly 163 factored cases, none of them powerful (otherwise abort loudly: counterexample); asserted.
#   E2  there are ω = 1 cases (known: Φ_7(143), Φ_19(95)); which ranks d they have is open and measured.
#   E3  Σ Π 1/p ≪ 1; expected order of magnitude 10⁻⁴ or smaller, dominated by the smallest prime factors.

import json, sys, math, pathlib
from collections import Counter
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'
aus = []   # `aus` = output lines, written to the output file at the end
def sag(s=''):   # `sag` = "say": print a line and keep it for the output file
    print(s, flush=True); aus.append(s)
sag('='*96); sag('w78 — ω(Φ_d), HEURISTIK, PRUEFSTEIN FUER JURICEVIC'); sag('='*96)

s55 = json.load(open(ERG/'w30_vollfaktorisierung_S55_N25_result.json', encoding='utf-8'))
w70 = json.load(open(ERG/'w70_zeitlimit_faelle_result.json', encoding='utf-8'))
F = {}   # `F` = complete factorizations: (kernel m, rank d) -> {prime: exponent} of Φ_d
for x in s55['teil_B']['raenge']:
    if x['vollstaendig']: F[(x['m'], x['d'])] = {int(p): int(e) for p, e in x['phi_faktoren'].items()}
for x in w70['faelle']:
    F[(x['m'], x['d'])] = {int(p): int(e) for p, e in x['phi'].items()}
assert len(F) == 163, ('E1: erwartet 163', len(F))
assert not any(f and all(e >= 2 for e in f.values()) for f in F.values()), 'E1: POWERFUL GEFUNDEN — GEGENBEISPIEL!'
sag(f'PK E1 ✅  {len(F)} vollstaendig faktorisierte Φ_d, keines powerful.')
sag()

# ---------------------------------------------------------------- ④ distribution of ω
om = Counter(len(f) for f in F.values())   # `om` = number of cases per value of ω
sag('--- ④ ω(Φ_d) = Zahl der verschiedenen Primteiler:')
for k in sorted(om): sag(f'    ω = {k}: {om[k]:>3} Faelle  ({100*om[k]/len(F):.0f} %)')
# `prim` = cases where Φ_d is prime; `pp` = cases where Φ_d is a prime power with exponent e ≥ 2
prim = sorted([(m, d) for (m, d), f in F.items() if len(f) == 1 and list(f.values())[0] == 1], key=lambda x: x[1])
pp   = sorted([(m, d) for (m, d), f in F.items() if len(f) == 1 and list(f.values())[0] >= 2])
sag(f'    Φ_d ist PRIM in {len(prim)} Faellen: ' + ', '.join(f'Φ_{d}({m})' for m, d in prim))
sag(f'    Φ_d ist Primzahlpotenz (ω = 1, e ≥ 2): {len(pp)}  {pp if pp else ""}')
sag()

# ---------------------------------------------------------------- ③ test of the transfer
sag('--- ③ PRUEFSTEIN: Juricevic-Uebertragung wuerde ω ≥ 2 fuer d ≥ 16 verlangen')
verletzt = [(m, d) for m, d in prim if d >= 16]   # `verletzt` = violations: ω = 1 although d ≥ 16
klein    = [(m, d) for m, d in prim if d < 16]    # `klein` = ω = 1 with d < 16 (outside the range of the theorem)
sag(f'    ω = 1 mit d < 16  (ausserhalb der Reichweite des Satzes): {klein}')
sag(f'    ω = 1 mit d ≥ 16  (WIDERSPRAECHE der Uebertragung):      {verletzt if verletzt else "keine"}')
if verletzt:
    sag('    ⇒ ⚠️ Die naive Uebertragung greift NICHT fuer alle d ≥ 16. Entweder schliessen (1.12)/η/κ unser Paar')
    sag('       fuer diese d aus, oder die Rang-Korrespondenz (Rang d = U-Rang 2d) braucht eine Zusatzbedingung.')
    sag('       Das ist der Punkt, an dem Jones 2011 („still a mystery") recht behaelt. Bedingungen am Volltext pruefen.')
else:
    sag('    ⇒ vertraeglich mit der Uebertragung (kein Gegenbeleg unter den 163 — das ist Vertraeglichkeit, kein Beweis).')
d_max_om1 = max((d for m, d in prim), default=None)
sag(f'    groesstes d mit ω = 1: {d_max_om1}')
sag()

# ---------------------------------------------------------------- ⑤ heuristic
sag('--- ⑤ HEURISTIK (Etikett: Modell, kein Beweis): P(Φ_d powerful) ≈ Π_{p|Φ_d} 1/p')
# `beitr` = contributions: Π 1/p per case, sorted in descending order; `E` = their sum (expected number of powerful cases)
beitr = sorted(((math.prod(1/p for p in f), m, d) for (m, d), f in F.items()), reverse=True)
E = sum(b for b, _, _ in beitr)
sag(f'    E[# powerful Φ_d unter den 163] = Σ Π 1/p = {E:.3e}')
sag('    die 5 groessten Beitraege:')
for b, m, d in beitr[:5]:
    f = F[(m, d)]
    sag(f'      Φ_{d}({m}) = ' + ' · '.join(f'{p}' + (f'^{e}' if e > 1 else '') for p, e in sorted(f.items())) + f'   →  {b:.2e}')
sag(f'    ⇒ ein Gegenbeispiel unter den 163 waere ein Ereignis der Groessenordnung 10^{math.log10(E):.0f}.')
sag('    ⚠️ Das ist die INNERE Heuristik (jeder Primteiler unabhaengig, W\'keit 1/p). Die 16 realen Wieferich-')
sag('       Ereignisse aus w77 zeigen, dass 1/p die richtige Groessenordnung hat (16 beobachtet gegen 20,03 erwartet).')

res = dict(skript=pathlib.Path(__file__).name, datum='2026-09-22', n=len(F),
           omega_verteilung={str(k): v for k, v in sorted(om.items())},
           phi_prim=[list(x) for x in prim], phi_primpotenz=[list(x) for x in pp],
           pruefstein_verletzt=[list(x) for x in verletzt], d_max_omega1=d_max_om1,
           heuristik_E=E, top5=[[m, d, b] for b, m, d in beitr[:5]])
(ERG/'w78_omega_heuristik_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w78_omega_heuristik_output.txt').write_text('\n'.join(aus)+'\n', encoding='utf-8')
print('\nErgebnis: w78_omega_heuristik_result.json')
