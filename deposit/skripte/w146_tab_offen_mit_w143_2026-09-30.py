# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# -*- coding: utf-8 -*-
# w146_tab_offen_mit_w143_2026-09-30.py
# Part I, table of open blocks (Prop. I.6.8): the ECM effort spent on each open block, including the overnight run w143.
# Successor of w142, which for the three new open blocks read only the curves of w135 (300 at B1 = 5*10^4, 72 at B1 = 10^6).
#   w143 ran further ECM curves on exactly these three blocks (no factor found). This script adds the curve counts of both runs
#   and, like w142, asks GMP-ECM for each number and each B1 (with ONE curve) for the chosen stage-2 bound B2 and for the
#   expected number of curves for 35- and 40-digit factors (w120: above ~680 digits GMP-ECM chooses a smaller B2). The older
#   open blocks come from w92 as in w142; blocks that w158 closes (algebraic parts probable prime) are dropped.
# Reads (ergebnisse/): w92_karte_811_daten.json, w135_stufe2_ecm_result.json, w135_stufe2_ecm_checkpoint.json,
#   w143_nacht_ecm_checkpoint.json, w143_nacht_ecm_result.json, w142_tab_offen_result.json,
#   w158_teile_wahrscheinlich_prim_result.json; imports w142_tab_offene_raenge_2026-09-29.py (function `frage`), which calls the
#   GMP-ECM program.
# Writes: ergebnisse/w146_tab_offen_result.json and ergebnisse/w146_tab_offen_output.txt. No command-line arguments.
# Controls (expectations fixed in advance; all asserted):
#   E1       No ECM block is running (no lock file), so the query curves do not compete with a block for the processors.
#   E2 (PK)  As in w142: B1 = 10^6 on a 151-digit number gives B2 = 1 045 563 762 and 1071 curves for 35 digits; B1 = 5*10^4
#            gives B2 = 12 746 592 and 88265 curves.
#   E3 (PK)  Two routes for the curve counts of the new blocks: (a) the w143 checkpoint (`kurven_w135` + `kurven`), (b) the
#            w135 result + the w143 result. Both must agree per block and B1; the remaining open number must be the same in
#            w135 and w143.
#   E4       Number of old open blocks = 6 minus those closed by w158 (after the factor found for (95, 475) in a later ECM
#            block, w92 holds 6 open old blocks, previously 7); 3 new blocks, none of them closed by w158.
#   E5       t35 and t40 can only grow compared with w142 (more curves, same B2 per B1); a decrease would be an error.
import json, re, subprocess, sys, time, pathlib, importlib.util
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'; BOX = HIER.parent.parent
# `BOX` = directory above the data folder; the lock file `ECM_LAEUFT.lock` (exists while an ECM block runs) is looked up
#   there.
# `aus` = collected output lines; `sag` (= say) prints a line and records it for the output file.
aus = []
def sag(s=''):
    print(s, flush=True); aus.append(s)

_spec = importlib.util.spec_from_file_location('w142', HIER/'w142_tab_offene_raenge_2026-09-29.py')
# only `frage` (= ask) is used; w142's main() does not run (__main__ guard)
w142 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(w142)

def main():
    assert not (BOX/'ECM_LAEUFT.lock').exists(), 'E1: ein ECM-Block laeuft (Sperrdatei) — spaeter starten'
    sag('E1 ✅  kein ECM-Block aktiv.')
    # `probe` = 151-digit test number for E2; `frage(N, B1)` returns B2, `k35` / `k40` (expected curves for 35- / 40-digit
    # factors), `faktor` (a factor found by the query curve, if any) and `sek` (seconds).
    probe = 10**150 + 7
    p1, p2 = w142.frage(probe, 1_000_000), w142.frage(probe, 50_000)
    e2 = (p1['B2'], p1['k35'], p2['B2'], p2['k35']) == (1_045_563_762, 1071.0, 12_746_592, 88265.0)
    sag(f"PK E2 {'✅' if e2 else '❌'}  B1 = 10⁶: B2 = {p1['B2']}, 35 St. {p1['k35']} · B1 = 5·10⁴: B2 = {p2['B2']}, 35 St. {p2['k35']}")
    assert e2
    # `w135` = earlier ECM result with the open blocks; `ck135` / `ck143` = checkpoints of w135 / w143; `r143` = w143 result
    # per block (key `m_d`); `w92` = block map; `alt142` = rows of w142 by (m, d), used in E5.
    w135 = json.loads((ERG/'w135_stufe2_ecm_result.json').read_text(encoding='utf-8'))
    ck135 = json.loads((ERG/'w135_stufe2_ecm_checkpoint.json').read_text(encoding='utf-8'))
    ck143 = json.loads((ERG/'w143_nacht_ecm_checkpoint.json').read_text(encoding='utf-8'))
    r143 = {f"{p['m']}_{p['d']}": p for p in json.loads((ERG/'w143_nacht_ecm_result.json').read_text(encoding='utf-8'))['paare']}
    w92 = json.loads((ERG/'w92_karte_811_daten.json').read_text(encoding='utf-8'))
    alt142 = {(z['m'], z['d']): z for z in json.loads((ERG/'w142_tab_offen_result.json').read_text(encoding='utf-8'))['zeilen']}
    # `zeilen` = rows of the table; `alt` = old open blocks (route `offen` in w92).
    zeilen = []
    alt = [r for r in w92['raenge'] if r['route'] == 'offen']
    # w158 closes some blocks via their algebraic parts (classes `beide_prp`, `zeuge_prp`); they are no longer open.
    k158 = json.loads((ERG/'w158_teile_wahrscheinlich_prim_result.json').read_text(encoding='utf-8'))['klassen']
    zu158 = {tuple(x) for x in k158['beide_prp'] + k158['zeuge_prp']}
    alt_vorher = len(alt)
    alt = [r for r in alt if (r['m'], r['d']) not in zu158]
    assert alt_vorher - len(alt) == len(zu158), ('w158 schliesst einen Block, der hier nicht offen war', sorted(zu158))
    sag(f'w158: {sorted(zu158)} geschlossen (wahrscheinlich prime algebraische Teile) — nicht mehr in der Tabelle')
    for r in alt:
        w = r['w89']
        zeilen.append(dict(m=r['m'], d=r['d'], kern='alt', stellen=int(r['lg']) + 1, kurven={int(b): n for b, n in w['kurven'].items()},
                           t35=w['t35'], t40=w['t40']))
    # `neu` = new open blocks from w135; each entry is (m, d, ..., curve counts per B1).
    neu = w135['offen']
    for o in neu:
        m, d = o[0], o[1]; k = f'{m}_{d}'; c = ck143[k]
        # `route_a` / `route_b` = curve counts per B1 by route (a) / (b) of E3; `rest_gleich` = same remaining open number.
        route_a = {int(b): c['kurven_w135'].get(b, 0) + c['kurven'].get(b, 0) for b in set(c['kurven_w135']) | set(c['kurven'])}
        r = r143[k]; b135 = {str(b): n for b, n in o[3].items()}
        route_b = {int(b): b135.get(b, 0) + r['kurven_w143'].get(b, 0) for b in set(b135) | set(r['kurven_w143'])}
        rest_gleich = int(c['rest']) == int(ck135[k]['rest'])
        e3 = route_a == route_b and rest_gleich
        sag(f"PK E3 {'✅' if e3 else '❌'}  ({m}, {d}): Checkpoint {dict(sorted(route_a.items()))} · Ergebnisse {dict(sorted(route_b.items()))} · Rest gleich {rest_gleich}")
        assert e3
        # `N` = the remaining open number of the block; `kurven` = nonzero curve counts per B1.
        N = int(c['rest']); kurven = {b: n for b, n in route_a.items() if n}
        # t35 = sum over B1 of (curves done / expected curves for a 35-digit factor), i.e. the effort in units of one
        # complete 35-digit search; t40 likewise; `b2s` = the B2 chosen by GMP-ECM per B1.
        t35 = t40 = 0.0; b2s = {}
        for B1, n in sorted(kurven.items()):
            f = w142.frage(N, B1)
            if f['faktor']: sag(f'  ‼ FAKTOR bei der Frage gefunden ({m}, {d}): {f["faktor"]} — NICHT gebucht; auswerten wie w143')
            assert f['B2'] and f['k35'] and f['k40'], ('keine Tabelle', m, d, B1, f)
            t35 += n / f['k35']; t40 += n / f['k40']; b2s[B1] = f['B2']
            sag(f"  ({m}, {d}) {len(str(N))} St.  B1 = {B1:>9,}: {n:>5} Kurven, B2 = {f['B2']:,}, erwartet 35/40 St.: {f['k35']}/{f['k40']} ({f['sek']} s)")
        z = dict(m=m, d=d, kern='neu', stellen=len(str(N)), kurven=kurven, t35=round(t35, 3), t40=round(t40, 3), b2=b2s)
        a = alt142[(m, d)]; e5 = z['t35'] >= a['t35'] and z['t40'] >= a['t40']
        sag(f"PK E5 {'✅' if e5 else '❌'}  ({m}, {d}): t35 {a['t35']} → {z['t35']}, t40 {a['t40']} → {z['t40']}")
        assert e5
        zeilen.append(z)
    # E4: after the factor found for (95, 475), w92 holds only 6 open old blocks (previously 7).
    e4 = len(alt) == 6 - len(zu158) and len(neu) == 3 and not any((o[0], o[1]) in zu158 for o in neu)
    sag(f"PK E4 {'✅' if e4 else '❌'}  {len(alt)} alte + {len(neu)} neue = {len(zeilen)} offene Bloecke"); assert e4
    zeilen.sort(key=lambda z: (z['stellen'], z['m']))
    for z in zeilen:
        sag(f"  {z['kern']}  ({z['m']}, {z['d']})  {z['stellen']} St.  Kurven {z['kurven']}  t35 {z['t35']}  t40 {z['t40']}")
    (ERG/'w146_tab_offen_result.json').write_text(json.dumps(dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'),
        quellen=['w92_karte_811_daten.json', 'w135_stufe2_ecm_result.json', 'w135_stufe2_ecm_checkpoint.json', 'w143_nacht_ecm_checkpoint.json',
                 'w143_nacht_ecm_result.json'],
        pk=dict(b1_1e6=p1, b1_5e4=p2), zeilen=[{**z, 'kurven': {str(k): v for k, v in z['kurven'].items()},
        **({'b2': {str(k): v for k, v in z['b2'].items()}} if 'b2' in z else {})} for z in zeilen]), indent=1, ensure_ascii=False), encoding='utf-8')
    sag('Ergebnis: w146_tab_offen_result.json')
    (ERG/'w146_tab_offen_output.txt').write_text(f'w146 · {time.strftime("%Y-%m-%d %H:%M")}\n' + '\n'.join(aus) + '\n', encoding='utf-8')

if __name__ == '__main__':
    main()
