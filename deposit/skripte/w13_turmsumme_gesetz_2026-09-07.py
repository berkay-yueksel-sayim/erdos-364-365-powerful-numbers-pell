# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w13_turmsumme_gesetz_2026-09-07.py
# Purpose: the pair count is an EXACT tower sum (`turmsumme`). For every kernel m with T1 even, the pairs (T_k - 1, T_k + 1) are
#   exactly the odd multiples k of m' (Lemma 1.3' and Lemma L). With T_k = (eps^k + eps^-k)/2 and eps = T1 + U1*sqrt(m),
#   digits(T_k) <= h holds iff k*log10(eps) - log10(2) <= h (up to the negligible term eps^-k). Hence
#   P_2(10^h; kernels <= M) = sum over m of #{ k = m', 3m', 5m', ... : k <= (h + log10 2)/log10 eps_m }. The script evaluates this
#   sum over the candidate families of the two w9 runs and compares it with the COUNTED points of the same runs at six heights.
# Reads: w9_v31_M100000000_H2000_result.json and w9_v32_c3_M10000000_H2000_result.json (current directory). Writes: nothing.
# Usage: python w13_turmsumme_gesetz_2026-09-07.py
# Controls: positive: exact agreement of the sum with the counted pairs at every height (abort otherwise).
#   Negative: the naive formula with log10 T1 in place of log10 eps must be off (it was: 462 instead of 371).
import json, math, sys
# (result file, residue class of m mod 8)
RUNS = (('w9_v31_M100000000_H2000_result.json', 7), ('w9_v32_c3_M10000000_H2000_result.json', 3))
HEIGHTS = (100, 250, 500, 1000, 1500, 2000)

def log10eps(l10T1):
    if l10T1 >= 15: return l10T1 + math.log10(2.0)
    T1 = 10.0 ** l10T1
    return math.log10(T1 + math.sqrt(T1 * T1 - 1.0))

# `ev` = list of (m, m', log10 eps, log10 T1), one entry per candidate family with T1 even
def tower_sum(ev, h, use_T1=False):
    tot = 0
    for m, mp, le, l10T1 in ev:
        step = l10T1 if use_T1 else le
        kmax = (h + math.log10(2.0)) / step
        if kmax >= mp: tot += int((kmax / mp + 1) // 2)   # odd multiples of m' up to kmax
    return tot

total_slope = 0.0
for fn, cls in RUNS:
    # `cands` = candidate families, `pts` = lattice points, `paare` = the counted points that are pairs (status not 'parity')
    d = json.load(open(fn)); cands = d['cands']; pts = d['points']
    paare = [p for p in pts if p[4] != 'parity']
    ev = [(c[0], c[1], log10eps(c[3]), c[3]) for c in cands if c[6] == 'T1even']
    slope = sum(1.0 / (2 * mp * le) for m, mp, le, _ in ev); total_slope += slope
    print(f"Klasse {cls} (mod 8), Kerne <= {d['MMAX']}: {len(cands)} Kandidaten-Familien, {len(ev)} mit T1 gerade; gezaehlte Paare bis 10^{d['H']}: {len(paare)}")
    print(f"   lineare Steigung Summe 1/(2 m' log10 eps) = {slope:.4f} Paare pro Dezimalstelle")
    top = sorted(ev, key=lambda t: -1.0 / (2 * t[1] * t[2]))[:6]
    print("   groesste Beitraege (m, m', log10 eps, Anteil): " + "; ".join(f"({m}, {mp}, {le:.3f}, {1/(2*mp*le)/slope:.3f})" for m, mp, le, _ in top))
    ok = True
    for h in HEIGHTS:
        g = sum(1 for p in paare if p[3] <= h); t = tower_sum(ev, h); n = tower_sum(ev, h, use_T1=True)
        flag = 'OK' if g == t else 'MISMATCH'; ok &= (g == t)
        print(f"   bis 10^{h:4d}: gezaehlt {g:4d} | Turm-Summe {t:4d} {flag} | naiv (log10 T1) {n:4d} | linear {slope*h:7.1f}")
    assert ok, "PK FEHLGESCHLAGEN: Turm-Summe trifft die gezaehlten Paare nicht"
    assert tower_sum(ev, 2000, use_T1=True) != len(paare), "Negativ-Kontrolle: naive Formel darf nicht zufaellig passen"
print(f"Beide Klassen zusammen: Steigung {total_slope:.4f} Paare pro Dezimalstelle (Kerne <= 1e8 bzw. 1e7).")
print("PK ok: exakte Uebereinstimmung an allen Hoehen in beiden Klassen; naive log10-T1-Formel liegt daneben (Negativ-Kontrolle).")
