# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w159_alekseyev_klasse3_bis_1e8_2026-10-02.py
# Purpose: Part I, Prop. I.6.2(b): comparison with Alekseyev's table (OEIS A076445), now with the class-3 kernels (m ≡ 3 mod 8)
#   up to 10^8 instead of 10^6. The proposition used class 3 only up to 10^6 (key `klassedreimax` of w116), whereas Part II
#   (Prop. II.6.1) has class 3 complete up to 10^8 (w9_v32_c3_M100000000); this script redoes the comparison with that range.
# Method (as in w116, own code rewritten here): the OEIS entry is licensed CC BY-SA and is only cited;
#   Alekseyev's values are NOT read in.
#   Instead the SHA-256 hash of our own values (the smaller member a = n − 1 of each pair, ascending, one "a\n" per value, up to
#   67 digits) is compared with the hash that w116 stored for Alekseyev's 33 terms. The hash does not reveal the values.
# Reads:  ergebnisse/w9_v31_M100000000_H2000_result.json (class 7 up to 10^8),
#         ergebnisse/w9_v32_c3_M1000000_H2000_result.json and ergebnisse/w9_v32_c3_M100000000_H2000_result.json (class 3).
# Writes: ergebnisse/w159_alekseyev_klasse3_1e8_result.json and ergebnisse/w159_alekseyev_klasse3_1e8_output.txt.
# Usage:  python w159_alekseyev_klasse3_bis_1e8_2026-10-02.py   (no arguments)
# Controls (E1 and E2 are positive controls, E4 is a negative control; E1, E2, E4 are asserted):
#   E1  class 3 up to 10^6 together with class 7 up to 10^8 gives exactly the hash of w116 (33 values).
#   E2  the class-3 list up to 10^8, restricted to m ≤ 10^6, equals the list up to 10^6 (same hits up to 70 digits).
#   E3  with class 3 up to 10^8: number of values up to 67 digits and their hash. Equal hash: the comparison also holds with 10^8.
#       Unequal: the additional values (only count, digits, kernel) are reported and nothing is adjusted.
#   E4  dropping one value must give a different hash.
import sys, json, time, pathlib, hashlib
from math import isqrt
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'
SHA_ALEKSEYEV_33 = 'a4af793c863e8a4d6114c7bf446cb66ad5fda3eea03be2936cff3426af0e7a9a'   # from w116: hash of the 33 terms
aus = []   # `aus` = output lines, written to the output file at the end
def sag(s=''):   # `sag` = "say": print a line and keep it for the output file
    print(s, flush=True); aus.append(s)
def sha(werte): return hashlib.sha256(''.join(f'{a}\n' for a in sorted(werte)).encode()).hexdigest()   # hash of the sorted values

def fund(m):   # fundamental solution (h0, k0) of x² − m·y² = 1 by the continued fraction of √m
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
def T(m, k):                            # x-part T_k of (T1 + U1·sqrt m)^k, by repeated squaring
    T1, U1 = fund(m); X, Y, bx, by = 1, 0, T1, U1
    while k:
        if k & 1: X, Y = X*bx + m*Y*by, X*by + Y*bx
        bx, by = bx*bx + m*by*by, 2*bx*by; k >>= 1
    return X

def werte(dateien, mmax3=None):   # `werte` = values: list of (n − 1, class, kernel m, digits) over the hits of the given files
    a = []                        # `dateien` = list of (class, file name); `mmax3` = optional kernel limit for class 3
    for cls, name in dateien:
        d = json.load(open(ERG/name, encoding='utf-8'))
        for m, mp, k, st, zg in d['hits']:   # hit = (kernel m, m' (`mp`), exponent k, digit count `st`, witnesses `zg`)
            if st > 70: continue
            if cls == 3 and mmax3 and m > mmax3: continue
            n = T(m, k); assert len(str(n)) == st, (cls, m, k, st)
            a.append((n - 1, cls, m, st))
    return a

K7 = (7, 'w9_v31_M100000000_H2000_result.json')   # `K7` = (class, file) for the class-7 kernels (m ≡ 7 mod 8) up to 10^8
# `p6`, `p8` = values with class 3 up to 10^6 and up to 10^8; `w6`, `w8` = the values up to 67 digits that enter the hash
p6 = werte([K7, (3, 'w9_v32_c3_M1000000_H2000_result.json')])
w6 = [x[0] for x in p6 if x[3] <= 67]
e1 = sha(w6) == SHA_ALEKSEYEV_33 and len(w6) == 33
sag(f'PK E1 {"✅" if e1 else "❌"}  Klasse 3 bis 10⁶ + Klasse 7 bis 10⁸: {len(w6)} Werte bis 67 Stellen, Hash wie w116: {sha(w6) == SHA_ALEKSEYEV_33}')
assert e1

p8 = werte([K7, (3, 'w9_v32_c3_M100000000_H2000_result.json')])
# `drei6` = class-3 hits (kernel, value) of the 10^6 run; `drei8_bis6` = the same from the 10^8 run, restricted to m ≤ 10^6
drei6 = sorted((x[2], x[0]) for x in p6 if x[1] == 3)
drei8_bis6 = sorted((x[2], x[0]) for x in p8 if x[1] == 3 and x[2] <= 10**6)
e2 = drei6 == drei8_bis6
sag(f'PK E2 {"✅" if e2 else "❌"}  Klasse-3-Liste bis 10⁸, eingeschraenkt auf m ≤ 10⁶, gleich der Liste bis 10⁶ ({len(drei6)} Treffer bis 70 Stellen)')
assert e2

w8 = [x[0] for x in p8 if x[3] <= 67]
# `neu` = additional values (digits, kernel) from class-3 kernels with 10^6 < m ≤ 10^8
neu = sorted((x[3], x[2]) for x in p8 if x[3] <= 67 and x[1] == 3 and x[2] > 10**6)
e3 = sha(w8) == SHA_ALEKSEYEV_33
sag(f'E3 {"✅" if e3 else "⚠️"}  Klasse 3 bis 10⁸: {len(w8)} Werte bis 67 Stellen, Hash wie Alekseyevs 33: {e3}; '
    f'zusaetzliche Werte aus Klasse-3-Kernen 10⁶ < m ≤ 10⁸ bis 67 Stellen: {len(neu)}' + (f' (Stellen, Kern: {neu})' if neu else ''))
# `kl8` = class-3 kernels with 10^6 < m ≤ 10^8 that have a pair of up to 70 digits
kl8 = sorted({x[2] for x in p8 if x[1] == 3 and x[2] > 10**6})
sag(f'          Klasse-3-Kerne 10⁶ < m ≤ 10⁸ mit einem Paar bis 70 Stellen: {len(kl8)}')

e4 = sha(w8[1:]) != SHA_ALEKSEYEV_33
sag(f'NK E4 {"✅" if e4 else "❌"}  ein Wert weniger gibt einen anderen Hash')
assert e4

res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), hash_wie_alekseyev_33_klasse3_1e6=e1,
           hash_wie_alekseyev_33_klasse3_1e8=e3, werte_bis_67_klasse3_1e8=len(w8), zusaetzliche_werte_1e6_1e8=len(neu),
           kern_klasse3_max=10**8 if e3 else 10**6, kern_klasse7_max=10**8,
           quelle='eigene W9-Listen (w9_v31, w9_v32_c3 bis 10^6 und bis 10^8); Alekseyevs Tabelle nur per Hash aus w116 verglichen')
(ERG/'w159_alekseyev_klasse3_1e8_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w159_alekseyev_klasse3_1e8_output.txt').write_text(f'w159 · {time.strftime("%Y-%m-%d %H:%M")}\n' + '\n'.join(aus) + '\n', encoding='utf-8')
sag('Ergebnis: w159_alekseyev_klasse3_1e8_result.json')
