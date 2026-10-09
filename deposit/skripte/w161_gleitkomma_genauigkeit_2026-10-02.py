# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# -*- coding: utf-8 -*-
# w161_gleitkomma_genauigkeit_2026-10-02.py
# Part I, Theorem I.3.1 ("The computation"): how accurate is the floating-point value of log10 T_1 used in the kernel search?
# Motivation: the text states "to a relative accuracy of about 10^{-12}", but this had never been measured. The search
#   w9_hoehenbox_v31 compares floating-point and exact values only for the surviving kernels (|delta log10 T_1| < 10^-6), not
#   for the discarded ones, where the accuracy matters (required: relative error < 3*10^-4, Part I).
# Method: `fund_fast` is repeated here verbatim from w9_hoehenbox_v31_2026-09-06.py (that script cannot be imported, because
#   it starts the whole search on import). Exact reference: the smallest solution (T_1, U_1) of x^2 - m*y^2 = 1 by continued
#   fractions in integers (`fund`), log10 by math.log10 of the integer. Sample: ALL squarefree m = 7 (mod 8) up to 2*10^4 and
#   1500 random ones up to 10^8 (fixed seed 20261002).
# Reads: nothing. Writes: ergebnisse/w161_gleitkomma_genauigkeit_result.json and the matching _output.txt.
# Usage: no command-line arguments.
# Controls (expectations fixed in advance; E2, E3 and E4 are asserted):
#   E1 (PK)  fund_fast(7) gives residues (1, 3) and log10 8 (T_1(7) = 8, U_1(7) = 3), like the control in w9.
#   E2 (PK)  the residues (T_1 mod m, U_1 mod m) returned by `fund_fast` equal the exact ones for every kernel of the sample.
#   E3       measured: largest relative error |delta| / log10 T_1 over the sample; must be below 3*10^-4 (expected far below).
#   E4 (NK)  the same computation in single precision (every intermediate rounded to float32), on the first 400 kernels, gives
#            an error more than 100 times larger; otherwise E3 would not measure the floating-point error.
import sys, json, time, math, random, pathlib, struct
from math import isqrt
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'
# `aus` = collected output lines; `sag` (= say) prints a line and records it for the output file.
aus = []
def sag(s=''):
    print(s, flush=True); aus.append(s)

# `fund_fast` returns (T_1 mod m, U_1 mod m, floating-point log10 T_1, index of the last convergent). The residues come from
# the convergents of sqrt(m) taken mod m; the size is tracked in floating point, rescaled by 10^100 (10^30 for float32) with
#   the
# shift counted in `s`. Verbatim from w9_hoehenbox_v31_2026-09-06.py; `f32` is used only for the negative control E4.
def fund_fast(m, f32=False):
    r = (lambda x: struct.unpack('f', struct.pack('f', x))[0]) if f32 else (lambda x: x)
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    hm1, hm0 = 1 % m, a0 % m
    km1, km0 = 0, 1 % m
    hf1, hf0 = 1.0, float(a0)
    s = 0; i = 0
    while True:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        i += 1
        if Q == 1 and i % 2 == 0:
            return hm0, km0, math.log10(hf0) + s, i - 1
        hm1, hm0 = hm0, (a*hm0 + hm1) % m
        km1, km0 = km0, (a*km0 + km1) % m
        hf1, hf0 = hf0, r(a*hf0 + hf1)
        if hf0 > (1e30 if f32 else 1e100):
            d = 1e30 if f32 else 1e100; hf0 = r(hf0 / d); hf1 = r(hf1 / d); s += 30 if f32 else 100
# `fund` = exact fundamental solution (T_1, U_1) of x^2 - m*y^2 = 1 in integer arithmetic (continued fractions).
def fund(m):
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
# `quadratfrei` = squarefree test by trial division.
def quadratfrei(m):
    p = 2
    while p*p <= m:
        if m % (p*p) == 0: return False
        p += 1
    return True

tm, um, l10, per = fund_fast(7)
e1 = (tm, um) == (1, 3) and abs(l10 - math.log10(8)) < 1e-12
sag(f'PK E1 {"✅" if e1 else "❌"}  fund_fast(7): log10 T_1 = {l10:.15f} gegen log10 8 = {math.log10(8):.15f}'); assert e1

# `klein` = all squarefree m = 7 (mod 8) up to 2*10^4; `gross` = 1500 random ones in (2*10^4, 10^8); `stich` = the sample.
klein = [m for m in range(7, 20001, 8) if quadratfrei(m)]
rng = random.Random(20261002); gross = set()
while len(gross) < 1500:
    m = rng.randrange(7, 10**8, 8)
    if m > 20000 and quadratfrei(m): gross.add(m)
stich = klein + sorted(gross)
# `fmax` = (largest relative error, kernel where it occurs); `f32max` = largest float32 error; `laengste` = longest period.
t0 = time.time(); fmax = (0.0, None); f32max = 0.0; laengste = 0
for m in stich:
    T1, U1 = fund(m); tm, um, l10, per = fund_fast(m)
    assert (T1 % m, U1 % m) == (tm, um), ('E2 verletzt', m)
    exakt = math.log10(T1); rel = abs(l10 - exakt) / exakt
    if rel > fmax[0]: fmax = (rel, m)
    laengste = max(laengste, per)
for m in stich[:400]:                        # negative control on part of the sample (the float32 simulation is slow)
    T1, _ = fund(m); l32 = fund_fast(m, f32=True)[2]; exakt = math.log10(T1)
    f32max = max(f32max, abs(l32 - exakt) / exakt)
sag(f'PK E2 ✅  Reste mod m stimmen fuer alle {len(stich)} Kerne ({len(klein)} bis 2·10^4, {len(gross)} zufaellig bis 10^8); laengste Periode {laengste} ({time.time() - t0:.0f} s)')
sag(f'E3     groesster relativer Fehler von log10 T_1 (doppelte Genauigkeit): {fmax[0]:.2e} bei m = {fmax[1]} — benoetigt < 3·10^-4')
e4 = f32max > 100 * max(fmax[0], 1e-16)
sag(f'NK E4 {"✅" if e4 else "❌"}  einfache Genauigkeit auf 400 Kernen: groesster relativer Fehler {f32max:.2e} (deutlich groesser)'); assert e4
assert fmax[0] < 3e-4
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), n_kerne=len(stich), n_bis_2e4=len(klein), n_zufall_bis_1e8=len(gross),
           max_rel_fehler=fmax[0], bei_m=fmax[1], laengste_periode=laengste, f32_max_rel_fehler=f32max, seed=20261002)
(ERG/'w161_gleitkomma_genauigkeit_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w161_gleitkomma_genauigkeit_output.txt').write_text(f'w161 · {time.strftime("%Y-%m-%d %H:%M")}\n' + '\n'.join(aus) + '\n', encoding='utf-8')
sag('Ergebnis: w161_gleitkomma_genauigkeit_result.json')
