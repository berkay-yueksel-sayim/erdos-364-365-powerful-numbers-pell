# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w16_a135735_klasse3_kern_2026-09-08.py
# Purpose: examine the two class-3 kernels (m = 3 mod 8) among the 21 known terms of OEIS A135735 (squarefree d with d | y,
#   y from the fundamental solution of x^2 - d*y^2 = 1; Reinhart 2024: list exhaustive up to 5.325e13): m = 1752299 and
#   m = 20256129307923. The fundamental solution is computed modulo m from the continued fraction (`fund_fast`):
#   m | U1, parity of T1, log10 T1 and the period length.
# Consequence for Prop. 6.4: a class-3 kernel with m | U1 below the Reinhart bound yields pairs from 10^99133 on at the
#   earliest (T1 even, 99134 digits), or none at all (T1 odd).
# Reads: no files. Writes: no files (stdout only). Usage: python w16_a135735_klasse3_kern_2026-09-08.py
# Controls: positive: the class-7 kernel 4099215 must give m | U1, T1 even, 228 digits (as in run w11); 209991 must give T1 odd.
#   Negative: m = 7 must NOT satisfy m | U1.
import math, time
from math import isqrt
# `fund_fast` = fast fundamental solution via continued fraction of sqrt(m); returns
#   (T1 mod m, U1 mod m, T1 mod 2, log10 T1, period length); the float convergent is rescaled by 1e100 to avoid overflow
def fund_fast(m):
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    hm1, hm0 = 1 % m, a0 % m; km1, km0 = 0, 1 % m
    h21, h20 = 1, a0 % 2
    hf1, hf0 = 1.0, float(a0); s = 0; i = 0
    while True:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q; i += 1
        if Q == 1 and i % 2 == 0: return hm0, km0, h20, math.log10(hf0) + s, i - 1
        hm1, hm0 = hm0, (a*hm0 + hm1) % m; km1, km0 = km0, (a*km0 + km1) % m
        h21, h20 = h20, (a*h20 + h21) % 2
        hf1, hf0 = hf0, a*hf0 + hf1
        if hf0 > 1e100: hf0 /= 1e100; hf1 /= 1e100; s += 100
t0 = time.time()
# Positive and negative controls.
# Result names: `tm` = T1 mod m, `um` = U1 mod m, `par` = parity of T1, `l10` = log10 T1, `per` = period length
tm, um, par, l10, per = fund_fast(4099215); assert um == 0 and par == 0 and int(l10) + 1 == 228, ('PK 4099215', um, par, l10)
tm, um, par, l10, per = fund_fast(209991);  assert um == 0 and par == 1, ('PK 209991', um, par)
tm, um, par, l10, per = fund_fast(7);       assert um != 0 and (tm, um, par) == (1, 3, 0), ('Negativ-PK m=7', tm, um, par)
print("PK ok: 4099215 (m|U1, T1 gerade, 228 Stellen), 209991 (m|U1, T1 ungerade), m=7 (7 teilt U1=3 nicht)")
A135735 = [46, 430, 1817, 58254, 209991, 1752299, 3124318, 4099215, 5374184665, 6459560882, 16466394154, 20565608894, 25666082990,
           117477414815, 125854178626, 1004569189366, 1188580642033, 15826129757609, 18803675974841, 20256129307923, 39028039587479]
print("A135735 (21 Glieder, OEIS 08.09.2026): ungerade Glieder nach Restklasse mod 8:",
      {r: [d for d in A135735 if d % 2 and d % 8 == r] for r in (1, 3, 5, 7)})
for m in [d for d in A135735 if d % 8 == 3]:
    tm, um, par, l10, per = fund_fast(m)
    print(f"Klasse 3: m = {m}: m | U1 = {um == 0}, T1 {'gerade' if par == 0 else 'ungerade'}, log10 T1 = {l10:.3f} ({int(l10)+1} Stellen), Periode {per}"
          + (f"  -> erste Stufe (Paar) bei 10^{l10:.1f}, weit ueber (5.325e13)^1.5 = 3.9e20 und ueber H = 2000" if um == 0 and par == 0 else "  -> keine Paare (T1 ungerade)"))
print(f"Folge: Klasse-3-Paarliste unbedingt vollstaendig bis min(3.88e20, 10^99133) = 3.88e20 (Reinhart-Schranke), fuer m nicht teilt U1 bis 2.7e37. Zeit {time.time()-t0:.1f}s")
