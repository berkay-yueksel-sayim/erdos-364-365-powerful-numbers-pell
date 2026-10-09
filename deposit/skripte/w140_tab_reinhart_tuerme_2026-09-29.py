# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# -*- coding: utf-8 -*-
# w140_tab_reinhart_tuerme_2026-09-29.py
# Part I, Prop. 6.3: table of ladders for the known OEIS A135735 kernels m = 7 (mod 8) with T_1 even.
# Reads (ergebnisse/) only result files of w11: the runs with witness bound P = 3000 and with P = 2*10^5 (per ladder stage:
#   the witnesses p <= P with v_p(T_k) = 1, the first possible index k, and the minimum digit count of the middle n). For each
#   kernel the run with the larger bound P in which it occurs is used.
# Writes: ergebnisse/w140_tab_reinhart_result.json (one row per kernel with T_1 even). Usage: no command-line arguments.
# Controls (all asserted):
#   E1 (PK)  The P = 3000 run holds four kernels (4 099 215, 117 477 414 815, 39 028 039 587 479, 209 991), three of them with
#            T_1 even; the P = 2*10^5 run holds the first two.
#   E2 (PK)  Per kernel: stages are numbered 1, 2, ...; the first possible index k exceeds the index of the last stage if that
#            stage has witnesses; one bound P per kernel; the run with the larger P has at least as many stage-1 witnesses as
#            the smaller one.
#   E3 (PK)  Witnesses 701 (m = 4 099 215) and 5 (m = 39 028 039 587 479) occur in stage 1, the same as in Reinhart (2024), p.
#            6 (literature; used only as an expectation on our own output, no value is taken over).
#   E4 (NK)  Kernel 209 991 (T_1 odd) must not appear in the table.
import json, re, pathlib, sys
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'

# `lies` = read: parses one w11 output file into a dict kernel m -> block (`bloecke`) with the parity and digit count of T_1,
# the list of stages (`stufen`), the first possible index `k_min` and the digit bound `mitte_stellen` of the middle n.
def lies(name):
    t = (ERG/name).read_text(encoding='utf-8').replace('\r\n', '\n')
    bloecke = {}
    for kopf in re.finditer(r'^m = (\d+) \(T_1 (gerade|UNGERADE)[^\n]*\((\d+) Stellen\)', t, flags=re.M):
        m, par, st = int(kopf.group(1)), kopf.group(2), int(kopf.group(3))
        rest = t[kopf.end():]; ende = re.search(r'^\s+==> erster moeglicher Index k >= (\d+)[^\n]*Mitte n hat >= (\S+) Stellen', rest, flags=re.M)
        teil = rest[:ende.start()]
        stufen = []
        for s in re.finditer(r'^\s+Stufe (\d+): k = (\d+)\s+-> (?:Zeugen v_p = 1 \(p <= (\d+)\): \[([^\]]*)\]|kein Zeuge <= (\d+))', teil, flags=re.M):
            zeugen = [int(x) for x in s.group(4).split(',')] if s.group(4) else []
            stufen.append(dict(nr=int(s.group(1)), k=int(s.group(2)), P=int(s.group(3) or s.group(5)), zeugen=zeugen))
        bloecke[m] = dict(m=m, t1_gerade=(par == 'gerade'), t1_stellen=st, stufen=stufen, k_min=int(ende.group(1)), mitte_stellen=float(ende.group(2)))
    return bloecke

# `klein` / `gross` = blocks per kernel from the P = 3000 run / from the P = 2*10^5 run.
# Per stage: `zeugen` = witnesses (primes p <= P with v_p(T_k) = 1), `t1_gerade` = T_1 is even.
klein = lies('w11_a135735_tuerme_2026-09-07_output.txt'); gross = lies('w11_a135735_tuerme_P200000_2026-09-07_output.txt')
e1 = (sorted(klein) == [209991, 4099215, 117477414815, 39028039587479] and sorted(gross) == [4099215, 117477414815]
      and sum(1 for b in klein.values() if b['t1_gerade']) == 3)
print(f'PK E1 {"✅" if e1 else "❌"}  P = 3000: {sorted(klein)} · P = 2·10⁵: {sorted(gross)}')
assert e1
# `fehler` = list of detected inconsistencies (E2).
fehler = []
for quelle in (klein, gross):
    for b in quelle.values():
        nr = [s['nr'] for s in b['stufen']]
        if nr != list(range(1, len(nr) + 1)): fehler.append((b['m'], 'Nummern'))
        if b['k_min'] <= b['stufen'][-1]['k'] and b['stufen'][-1]['zeugen']: fehler.append((b['m'], 'k_min'))
        if len({s['P'] for s in b['stufen']}) != 1: fehler.append((b['m'], 'P'))
for m in gross:
    if len(gross[m]['stufen'][0]['zeugen']) < len(klein[m]['stufen'][0]['zeugen']): fehler.append((m, 'P gross hat weniger Zeugen'))
print(f'PK E2 {"✅" if not fehler else "❌"}  Stufen, Schranken und Zeugenzahlen stimmig' + (f' — {fehler}' if fehler else '.'))
assert not fehler
# `wahl` = chosen block per kernel: the P = 2*10^5 block if present, else the P = 3000 block.
wahl = {m: (gross[m] if m in gross else klein[m]) for m in klein}
e3 = 701 in wahl[4099215]['stufen'][0]['zeugen'] and 5 in wahl[39028039587479]['stufen'][0]['zeugen']
print(f'PK E3 {"✅" if e3 else "❌"}  Zeugen 701 und 5 in Stufe 1.')
assert e3
# Table rows: one per kernel with T_1 even; `mit` = the stages that have witnesses, `stufen` = their number,
# `stufen_zeugen` = the witness lists per stage, `t1_stellen` = digit count of T_1.
zeilen = []
for m in sorted(wahl):
    b = wahl[m]
    if not b['t1_gerade']: continue
    mit = [s for s in b['stufen'] if s['zeugen']]
    zeilen.append(dict(m=m, t1_stellen=b['t1_stellen'], P=b['stufen'][0]['P'], stufen_zeugen=[s['zeugen'] for s in mit],
                       stufen=len(mit), k_min=b['k_min'], mitte_stellen=b['mitte_stellen']))
e4 = 209991 not in {z['m'] for z in zeilen} and len(zeilen) == 3
print(f'NK E4 {"✅" if e4 else "❌"}  209 991 (T₁ ungerade) nicht in der Tabelle; {len(zeilen)} Zeilen.')
assert e4
for z in zeilen: print(f"  m = {z['m']:>15}  T1 {z['t1_stellen']:>8} St.  P = {z['P']:>6}  Stufen {z['stufen']}  {z['stufen_zeugen']}  k ≥ {z['k_min']}  Mitte ≥ {z['mitte_stellen']:.3g} St.")
# `paritaet` (key `paritaetstot`) = the kernel excluded by parity (T_1 odd), recorded for the negative control E4.
paritaet = dict(m=209991, t1_stellen=klein[209991]['t1_stellen'])
(ERG/'w140_tab_reinhart_result.json').write_text(json.dumps(dict(skript=pathlib.Path(__file__).name, datum='2026-09-29', zeilen=zeilen,
    paritaetstot=paritaet), indent=1, ensure_ascii=False), encoding='utf-8')
print('Ergebnis: w140_tab_reinhart_result.json')
