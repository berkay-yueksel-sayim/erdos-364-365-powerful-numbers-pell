# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w125_spektrumkopf_gleichstaende_2026-09-28.py
# Purpose: head of the gap spectrum (`spektrumkopf`) and ties (`gleichstaende`): the spectrum counts how often each gap between
#   consecutive powerful numbers occurs. Which statements about its head (for example: the top 50 gap values are powerful, the
#   first non-powerful value has rank 53) are TIE-PROOF, i.e. do not depend on how values of equal frequency are ordered?
# Reads: nothing. Writes (ergebnisse/, found relative to the script): w125_spektrumkopf_gleichstaende_result.json and
#   w125_spektrumkopf_gleichstaende_output.txt.
# Usage: python w125_spektrumkopf_gleichstaende_2026-09-28.py (no arguments).
# Controls: positive: counts of powerful numbers, counts of distinct gap values and the most frequent gap (when unique) must
#   equal the w48 values for N = 10^12, 10^13, 10^14. Negative-type check: the enumeration must be free of duplicates.

# Background: w48 assigns ranks with Counter.most_common(); for equal frequencies the order of FIRST OCCURRENCE in the gap
#   sequence decides there, an implementation artifact and not a mathematical rank. Even the top 20 at 10^14 contain ties
#   (noted in the form "134 x3, 121 x4").
# Question: which statements are tie-proof, i.e. independent of how ties are ordered?
#   (T1) Frequency threshold c_k := frequency of the k-th value; the GROUP of all values with occ >= c_k (ties included): how many
#       of them are powerful?
#   (T2) Largest threshold c* above which ALL values are powerful: all gap values with occ > c* are powerful, and at occ = c*
#        there is a non-powerful value. Tie-proof wording: "the R* most frequent values (all with occ > c*)".
#   (T3) The tie group containing the first non-powerful value: the range of ranks it can take depending on the ordering.
# Expectation (stated before the run): the numbers of w48 (2 158 391 / 6 840 384 / 21 663 503 powerful numbers; 832 161 /
#   2 593 517 / 8 235 300 distinct gap values; most frequent gap 900 resp. 44100) are reproduced exactly (positive control).
#   Whether "top 50" is tie-proof is OPEN; no expectation about the direction.
# Convention as in w48: 1 is powerful (a = b = 1); enumeration as a^2 * b^3 with b squarefree (Golomb; unique representation).
import sys, json, time, pathlib
from math import isqrt
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent   # `HIER` = this script's folder
ERG = HIER.parent/'ergebnisse'                   # `ERG` = results folder
T0 = time.time(); zeilen = []                    # `zeilen` = output lines, written to the _output.txt file
def sag(s=''): print(s, flush=True); zeilen.append(s)   # `sag` = say: print and record a line

# `faktoren` = factorization by trial division: dict prime -> exponent
def faktoren(n):
    f = {}; p = 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1 if p == 2 else 2
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
def ist_powerful(n): return n > 0 and all(e >= 2 for e in faktoren(n).values())   # `ist_powerful` = is powerful
def quadratfrei(n): return all(e == 1 for e in faktoren(n).values())               # `quadratfrei` = is squarefree

# `powerful_bis` = all powerful numbers <= N as a sorted array: a^2 * b^3 with b squarefree (`teile` = parts, one per b)
def powerful_bis(N):
    teile = []; b = 1
    while b ** 3 <= N:
        if quadratfrei(b):
            amax = isqrt(N // b ** 3)
            teile.append(np.arange(1, amax + 1, dtype=np.int64) ** 2 * (b ** 3))
        b += 1
    a = np.concatenate(teile); a.sort(); return a

# `SOLL` = target values from w48: N -> (`n_soll` = number of powerful numbers <= N, `w_soll` = number of distinct gap values,
# `top_soll` = most frequent gap)
SOLL = {10**12: (2158391, 832161, 900), 10**13: (6840384, 2593517, 44100), 10**14: (21663503, 8235300, 44100)}
ergebnis = {}   # `ergebnis` = results per N
for N, (n_soll, w_soll, top_soll) in SOLL.items():
    t = time.time()
    arr = powerful_bis(N)
    assert len(np.unique(arr)) == len(arr), 'PK-: duplikatfrei'
    d = np.diff(arr)   # gaps between consecutive powerful numbers
    werte, occ = np.unique(d, return_counts=True)   # `werte` = distinct gap values, `occ` = their frequencies
    # Ordering 1: frequency descending, ties by value ascending (mathematical order, reproducible); W = values, C = frequencies
    o = np.lexsort((werte, -occ)); W = werte[o]; C = occ[o]
    pk = (len(arr) == n_soll, len(werte) == w_soll, int(W[0]) == top_soll if C[0] > C[1] else True)   # positive control
    # `kopf` = head: the first 300 entries as (value, frequency, value is powerful)
    kopf = [(int(W[i]), int(C[i]), ist_powerful(int(W[i]))) for i in range(300)]
    # (T1) groups at the thresholds k = 10, 20, 50, 100, 200; `gruppe` = all entries with frequency >= c_k
    # keys: `schwelle` = threshold, `gruppengroesse` = group size, `davon_powerful` = how many are powerful,
    #   `gleichstand_ueber_k` = number of entries whose frequency equals the threshold
    t1 = {}
    for k in (10, 20, 50, 100, 200):
        c = int(C[k - 1]); gruppe = [x for x in kopf if x[1] >= c]
        t1[k] = dict(schwelle=c, gruppengroesse=len(gruppe), davon_powerful=sum(x[2] for x in gruppe),
                     gleichstand_ueber_k=sum(1 for x in kopf if x[1] == c))
    # (T2, T3) first non-powerful value (index `i_np`, frequency `c_np`) and its tie group;
    #   `rang_min`, `rang_max` = smallest and largest rank it can take; `nicht_pw_in_gruppe` = non-powerful values in the group
    i_np = next(i for i, x in enumerate(kopf) if not x[2]); c_np = kopf[i_np][1]
    rang_min = sum(1 for x in kopf if x[1] > c_np) + 1; rang_max = sum(1 for x in kopf if x[1] >= c_np)
    nicht_pw_in_gruppe = [x[0] for x in kopf if x[1] == c_np and not x[2]]
    alle_ueber = sum(1 for x in kopf if x[1] > c_np)          # number of values with occ > c_np; ALL of them are powerful
    # result row: `n_powerful`, `n_werte` = number of distinct gap values, `pk` = positive control passed, `top1` = most frequent
    #   gap, `top1_eindeutig` = it is unique, `schwellen` = thresholds (T1), `erster_nicht_pw` = first non-powerful value (`wert`
    #   = value, `occ`, `rang_spanne` = rank range, `nicht_pw_in_gruppe`), `alle_powerful_bis_rang` = all values up to this rank
    #   are powerful, `kopf60` = first 60 entries, `sekunden` = seconds
    ergebnis[str(N)] = dict(n_powerful=int(len(arr)), n_werte=int(len(werte)), pk=all(pk),
                            top1=[kopf[0][0], kopf[0][1]], top1_eindeutig=bool(C[0] > C[1]),
                            schwellen=t1, erster_nicht_pw=dict(wert=kopf[i_np][0], occ=c_np, rang_spanne=[rang_min, rang_max],
                                                               nicht_pw_in_gruppe=nicht_pw_in_gruppe),
                            alle_powerful_bis_rang=alle_ueber, kopf60=[[x[0], x[1]] for x in kopf[:60]],
                            sekunden=round(time.time() - t, 1))
    e = ergebnis[str(N)]
    sag(f'N = 10^{len(str(N)) - 1}: {e["n_powerful"]:,} powerful · {e["n_werte"]:,} Werte · PK {"✅" if e["pk"] else "🔴"} · Top 1 = {e["top1"]} '
        f'({"eindeutig" if e["top1_eindeutig"] else "GLEICHSTAND"}) · {e["sekunden"]} s')
    for k, v in t1.items():
        sag(f'   Schwelle Rang {k:>3}: occ ≥ {v["schwelle"]:>4} → Gruppe {v["gruppengroesse"]:>4} Werte, davon powerful {v["davon_powerful"]:>4}'
            f'  (am Rand {v["gleichstand_ueber_k"]} gleich haeufige)')
    enp = e['erster_nicht_pw']
    sag(f'   erster nicht-powerful Wert {enp["wert"]} (occ {enp["occ"]}): Rang je nach Ordnung {enp["rang_spanne"][0]}–{enp["rang_spanne"][1]}; '
        f'nicht-powerful in dieser Gruppe: {enp["nicht_pw_in_gruppe"]}')
    sag(f'   ⇒ TIE-FEST: die {alle_ueber} haeufigsten Werte (alle mit occ > {enp["occ"]}) sind ausnahmslos powerful.')
# result file: `skript` = script name, `ergebnis` = results per N, `laufzeit_s` = run time in seconds
ERG.joinpath('w125_spektrumkopf_gleichstaende_result.json').write_text(json.dumps(dict(skript=pathlib.Path(__file__).name, ergebnis=ergebnis,
    laufzeit_s=round(time.time() - T0, 1)), ensure_ascii=False, indent=1), encoding='utf-8')
ERG.joinpath('w125_spektrumkopf_gleichstaende_output.txt').write_text('\n'.join(zeilen) + '\n', encoding='utf-8')
sag(f'Laufzeit {time.time() - T0:.0f} s · Ergebnis: w125_spektrumkopf_gleichstaende_result.json')
