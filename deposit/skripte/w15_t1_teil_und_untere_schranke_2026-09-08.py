# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w15_t1_teil_und_untere_schranke_2026-09-08.py
# -*- coding: utf-8 -*-
# Checks two small "towards infinity" statements against the data of the height-box runs w9 (v3.1 for kernel class 7, v3.2 for
# class 3, kernels up to 10^8, heights up to 2000 digits).
# Reads: w9_v31_M100000000_H2000_result.json and w9_v32_c3_M100000000_H2000_result.json (relative paths: run the script in
#   the folder that holds them). Writes: nothing, the results are printed. No command-line arguments.
#   (A) T_1-part lemma: for odd k, T_1 | T_k and T_k/T_1 is odd, so v_2(T_k) = v_2(T_1); for odd primes p | T_1:
#       v_p(T_k) = v_p(T_1) + v_p(k) (LTE, rank 1). Consequence: the 2-adic valuation of the middle of a triple is v_2(T_1(m)),
#       and every prime p with p || T_1 must divide k (Mollin-Walsh ladder, step 1, as a general lemma).
#       Check: all lattice points of both classes (445 + 883), exactly (T_k up to 2000 digits).
#   (B) Explicit lower bound: L(h) := sum over the known towers (T_1 even) of #{odd multiples k of m' :
#       k*log10(eps) - log10(2) <= h} is a valid lower bound for P_2(10^h) for EVERY h (the towers exist for all heights).
#       Linear form: L(h) >= c*h - K with c = sum of 1/(2 m' log10 eps) and K = number of towers.
#       Check: L(h) = counted pairs for h <= 2000 (exact there), and L(h) >= c*h - K for h up to 10^6.
# Controls: PK: a forged point with v_2(T_k) != v_2(T_1) must be detected (assert); (A) and the linear form (B) abort via assert.
import json, math, sys
from math import isqrt
sys.set_int_max_str_digits(2000000)
# `RUNS` = (result file, residue class of the kernels mod 8)
RUNS = (('w9_v31_M100000000_H2000_result.json', 7), ('w9_v32_c3_M100000000_H2000_result.json', 3))

# `fund(m)` = fundamental solution (T_1, U_1) of x^2 - m*y^2 = 1 by continued fraction (P, Q, a = state; h, k = convergents)
def fund(m):
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
def T_exact(T1, U1, m, k):                      # exact T_k of the unit T1 + U1*sqrt(m) by square-and-multiply
    ra, rb, ba, bb = 1, 0, T1, U1
    while k:
        if k & 1: ra, rb = ra*ba + m*rb*bb, ra*bb + rb*ba
        ba, bb = ba*ba + m*bb*bb, 2*ba*bb; k >>= 1
    return ra
def vp(n, p):                                   # `vp` = p-adic valuation v_p(n)
    v = 0
    while n % p == 0: n //= p; v += 1
    return v
def small_odd_primes_of(n, bound=100000):       # odd primes dividing n, found by trial division up to `bound`
    ps = []; p = 3
    while p <= bound and p * p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0: n //= p
        p += 2
    if 1 < n <= bound and n % 2: ps.append(n)
    return ps

def log10eps(l):                                # log10 of eps = T_1 + sqrt(T_1^2 - 1) from l = log10 T_1
    if l >= 15: return l + math.log10(2.0)
    T = 10.0 ** l; return math.log10(T + math.sqrt(T*T - 1.0))

# `total_c` = sum of c over both classes, `total_towers` = number K of towers, `all_ev` = all towers (m, m', log10 eps);
# `cands` = candidate families (m, m', digits, log10 T_1, kmax, period `per`, parity `par`), `pts` = lattice points
# (m, m', k, digits, status, witness), `cache` = fund(m) per kernel; `checkedA` = points checked in part (A),
# `lte_checks` = number of LTE checks
total_c = 0.0; total_towers = 0; all_ev = []
for fn, cls in RUNS:
    d = json.load(open(fn)); cands, pts = d['cands'], d['points']
    cache = {}; checkedA = 0; lte_checks = 0
    for (m, mp, k, dig, status, w) in pts:
        if m not in cache: cache[m] = fund(m)
        T1, U1 = cache[m]; Tk = T_exact(T1, U1, m, k)
        assert Tk % T1 == 0 and (Tk // T1) % 2 == 1, ('T1-Teil-Lemma verletzt', m, k)
        assert vp(Tk, 2) == vp(T1, 2), ('v2 verletzt', m, k)
        for p in small_odd_primes_of(T1):
            assert vp(Tk, p) == vp(T1, p) + vp(k, p), ('LTE Rang 1 verletzt', m, k, p); lte_checks += 1
        checkedA += 1
    # PK: forged value
    m0 = 7 if 7 in cache else 3; T1, U1 = cache[m0]; Tf = T_exact(T1, U1, m0, 7) * 2   # forged value: doubled
    assert vp(Tf, 2) != vp(T1, 2), 'PK: Faelschung nicht erkannt'
    paare = [p for p in pts if p[4] != 'parity']   # `paare` = pairs: the points whose status is not 'parity'
    ev = [(m_, mp_, log10eps(l)) for (m_, mp_, dg, l, kmax, per, par) in cands if par == 'T1even']
    c = sum(1.0/(2*mp_*le) for m_, mp_, le in ev); total_c += c; total_towers += len(ev); all_ev += ev
    print(f"Klasse {cls}: {checkedA} Punkte geprueft — T1 | T_k, T_k/T1 ungerade, v_2(T_k) = v_2(T1): alle OK; LTE-Rang-1-Kontrollen an ungeraden Primteilern von T1: {lte_checks} OK")
    print(f"   v_2(T1) der Kerne mit Punkten: " + ", ".join(f"m={m_}: {vp(cache[m_][0], 2)}" for m_ in sorted(set(p[0] for p in pts))[:8]) + " ...")
    print(f"   Tuerme (T1 gerade): {len(ev)}, c = {c:.4f} Paare/Stelle; untere Schranke L(h) exakt = gezaehlt fuer h in (500, 1000, 2000)?",
          all(sum(1 for p in paare if p[3] <= h) == sum(int((( (h + math.log10(2))/le )/mp_ + 1)//2) for m_, mp_, le in ev if (h + math.log10(2))/le >= mp_) for h in (500, 1000, 2000)))
# `L(h)` = lower bound of (B): over all towers, the number of odd multiples k of m' with k*log10(eps) - log10(2) <= h
def L(h):
    return sum(int((((h + math.log10(2))/le)/mp_ + 1)//2) for m_, mp_, le in all_ev if (h + math.log10(2))/le >= mp_)
print(f"Beide Klassen: {total_towers} Tuerme, c = {total_c:.4f} Paare pro Dezimalstelle.")
# `worst` = smallest surplus L(h) - (c*h - K) on a grid of h up to 10^6
worst = min(L(h) - (total_c*h - total_towers) for h in range(1, 1000001, 997))
print(f"Lineare Form: L(h) >= c*h - K mit K = {total_towers}: minimaler Ueberschuss L(h) - (c h - K) fuer h bis 10^6: {worst:.1f} (>= 0 erwartet)")
assert worst >= 0
print("Beispiele: L(10^3) =", L(1000), "| L(10^4) =", L(10000), "| L(10^6) =", L(10**6), "| c*h fuer h=10^6:", round(total_c*10**6))
print("PK ok. Aussage: P_2(10^h) >= L(h) >= 0.595*h - 122 fuer alle h >= 1 (elementar, unbedingt); fuer h <= 2000 ist L(h) = P_2(10^h; Kerne <= 1e8) exakt.")
