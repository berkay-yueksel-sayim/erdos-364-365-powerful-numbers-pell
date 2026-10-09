# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w160_anteil_powerful_haeufige_luecken_2026-10-02.py
# -*- coding: utf-8 -*-
# Purpose: Part III, section 1: the share of powerful numbers among the frequent gap values (a comparison base). The statement
#   "only 0.10 % of the gap values are powerful" (w48) refers to ALL 8 235 300 distinct values up to 10^14; the figure is
#   dominated by large values that occur only once. Here the share is measured on a base of comparable frequency: the values
#   with at least 20 occurrences, on which the counting model w59 (Rem. III.1.4) also measures.
# Reads: nothing (the powerful numbers up to 10^14 are generated).
# Writes: ergebnisse/w160_anteil_powerful_haeufige_luecken_result.json and ..._output.txt.
# Usage: python w160_anteil_powerful_haeufige_luecken_2026-10-02.py   (no arguments)
# Controls (abort on failure): E1 to E6 below.
#
# Method (own code): all powerful numbers up to 10^14 as a^2 b^3 (b squarefree), sorted; gaps = differences; number of occurrences
#   per gap value; powerful test of a value by factorization (trial division up to the root; the values are small enough).
# Expectations and controls:
#   E1  PK inventory: 21 663 503 powerful numbers up to 10^14 (as in w59) and 8 235 300 distinct gap values (as in w125,
#       key `lswerte`).
#   E2  PK share: the share of powerful numbers among ALL values, rounded to two decimals, is 0.10 % (w48, key `lsanteil`).
#   E3  The values with at least 20 occurrences number 62 118 (w59, key `zmn`); their powerful share is measured.
#   E4  Null model as in w48 for the 50 most frequent values: expected number of powerful ones if they were drawn at random from
#       this base = 50 * share.
#   E5  NK: the powerful test rejects 63900 = 2^2*3^2*5^2*71 and accepts 44100 = 2^2*3^2*5^2*7^2.
#   E6  Second route: the vectorized sieve test over all values agrees with the single test on a sample (3000 random values plus
#       the first 500 powerful ones).
#   In addition (for context only, not used in the text): the share for the other minimum counts 2, 5, 100.

import sys, json, time, pathlib
from math import isqrt
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'
# `aus` = output lines (also written to the output file); `sag` = say: print a line and keep it
aus = []
def sag(s=''):
    print(s, flush=True); aus.append(s)

# `quadratfrei` = squarefree test; `powerful_bis(N)` = sorted array of all powerful numbers <= N as a^2 * b^3
# (`teile` = parts, one array per b); `ist_powerful(g)` = single powerful test by trial division
def quadratfrei(b):
    p = 2
    while p*p <= b:
        if b % (p*p) == 0: return False
        p += 1
    return True
def powerful_bis(N):
    teile, b = [], 1
    while b**3 <= N:
        if quadratfrei(b):
            amax = isqrt(N // b**3)
            if amax: teile.append(np.arange(1, amax + 1, dtype=np.int64)**2 * (b**3))
        b += 1
    a = np.concatenate(teile); a.sort(); return a
def ist_powerful(g):
    g = int(g); p = 2
    while p*p <= g:
        if g % p == 0:
            e = 0
            while g % p == 0: g //= p; e += 1
            if e == 1: return False
        p += 1 if p == 2 else 2
    return g == 1

e5 = ist_powerful(44100) and not ist_powerful(63900) and ist_powerful(1) and not ist_powerful(2)
sag(f'NK/PK E5 {"✅" if e5 else "❌"}  powerful-Test: 44100 ja, 63900 nein'); assert e5

# `arr` = the powerful numbers; `werte`, `anz` = the distinct gap values and their numbers of occurrences
t0 = time.time(); arr = powerful_bis(10**14)
d = np.diff(arr); werte, anz = np.unique(d, return_counts=True)
e1 = len(arr) == 21663503 and len(werte) == 8235300
sag(f'PK E1 {"✅" if e1 else "❌"}  {len(arr):,} powerful bis 10^14, {len(werte):,} verschiedene Lueckenwerte ({time.time() - t0:.0f} s)'); assert e1

# powerful test for all values, vectorized: every prime p up to the root of the largest value is divided out;
# exponent 1 => not powerful; a remainder > 1 at the end is a prime with exponent 1 => not powerful.
# Second route: the single test on a sample.
gmax = int(werte.max()); wurzel = isqrt(gmax)
sieb = np.ones(wurzel + 1, dtype=bool); sieb[:2] = False
for i in range(2, isqrt(wurzel) + 1):
    if sieb[i]: sieb[i*i::i] = False
# `pw` = Boolean array: the value is powerful (computed by the sieve below); `rest` = remainder after dividing out primes
rest = werte.copy(); pw = np.ones(len(werte), dtype=bool)
for p in np.nonzero(sieb)[0]:
    p = int(p); idx = np.nonzero(rest % p == 0)[0]
    if not len(idx): continue
    r = rest[idx]; e = np.zeros(len(idx), dtype=np.int64)
    while True:
        mm = r % p == 0
        if not mm.any(): break
        r[mm] //= p; e[mm] += 1
    pw[idx[e == 1]] = False; rest[idx] = r
pw &= rest == 1
rng = np.random.default_rng(20261002); stich = rng.choice(len(werte), 3000, replace=False)
stich = np.concatenate([stich, np.nonzero(pw)[0][:500]])  # also powerful values in the sample
abw = [int(werte[i]) for i in stich if ist_powerful(werte[i]) != bool(pw[i])]
sag(f'PK E6 {"✅" if not abw else "❌"}  Sieb gegen Einzeltest auf {len(stich)} Werten (davon 500 powerful): {len(abw)} Abweichungen; groesster Wert {gmax:,}')
assert not abw, abw[:10]
anteil_alle = 100 * pw.mean()
e2 = f'{anteil_alle:.2f}' == '0.10'
sag(f'PK E2 {"✅" if e2 else "❌"}  Anteil powerful unter allen Werten: {anteil_alle:.4f} % (gerundet {anteil_alle:.2f} %), {int(pw.sum()):,} Werte'); assert e2

# `ergebnis` = result per minimum number of occurrences (`mind`): number of values, number powerful, share in percent
ergebnis = {}
for mind in (2, 5, 20, 100):
    m = anz >= mind
    ergebnis[mind] = dict(werte=int(m.sum()), powerful=int(pw[m].sum()), anteil_prozent=round(100 * pw[m].mean(), 4))
    sag(f'          mindestens {mind:>3} Vorkommen: {int(m.sum()):>9,} Werte, davon {int(pw[m].sum()):>7,} powerful = {100 * pw[m].mean():.2f} %')
e3 = ergebnis[20]['werte'] == 62118
sag(f'E3 {"✅" if e3 else "❌"}  Basis mit mindestens 20 Vorkommen: {ergebnis[20]["werte"]:,} Werte (w59: 62 118), Anteil powerful {ergebnis[20]["anteil_prozent"]:.2f} %'); assert e3
# `erw50` = expected number of powerful values among 50 drawn at random from the base (E4)
erw50 = 50 * ergebnis[20]['powerful'] / ergebnis[20]['werte']
sag(f'E4     Nullmodell: von 50 zufaellig aus dieser Basis gezogenen Werten waeren etwa {erw50:.1f} powerful (aus allen Werten: {50 * pw.mean():.2f})')

res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), n_powerful=int(len(arr)), n_werte=int(len(werte)),
           powerful_werte=int(pw.sum()), anteil_alle_prozent=round(anteil_alle, 4), je_mindestzahl={str(k): v for k, v in ergebnis.items()},
           basis_mind=20, basis_werte=ergebnis[20]['werte'], basis_powerful=ergebnis[20]['powerful'], basis_anteil_prozent=ergebnis[20]['anteil_prozent'],
           nullmodell_50_basis=round(erw50, 2), nullmodell_50_alle=round(50 * pw.mean(), 4))
(ERG/'w160_anteil_powerful_haeufige_luecken_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w160_anteil_powerful_haeufige_luecken_output.txt').write_text(f'w160 · {time.strftime("%Y-%m-%d %H:%M")}\n' + '\n'.join(aus) + '\n', encoding='utf-8')
sag('Ergebnis: w160_anteil_powerful_haeufige_luecken_result.json')
