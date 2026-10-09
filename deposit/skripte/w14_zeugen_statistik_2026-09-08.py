# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w14_zeugen_statistik_2026-09-08.py
# Purpose: witness statistics (`zeugen_statistik`) for the 371 class-7 lattice points that carry a witness (kernel <= 1e8,
#   n < 10^2000). For each point (m, m', k, digits, witness p) from the "hits" list of w9_v31_M100000000_H2000_result.json:
#     alpha(p) = smallest j >= 1 with p | T_j(m)   (rank of apparition in the T sequence; Bennett-Walsh, Lemma 3.3:
#                p | T_k iff k is an odd multiple of alpha)
#     primitive := (alpha(p) == k);  cofactor k/alpha(p) (odd);  v_p(T_alpha) exactly;
#     LTE check (lifting the exponent): v_p(T_k) = v_p(T_alpha) + v_p(k/alpha) = 1.
# Reads: w9_v31_M100000000_H2000_result.json (current directory); writes nothing. Usage: python w14_zeugen_statistik_2026-09-08.py
# Controls: positive: the Mollin-Walsh point T_7(7) with witness 29 (alpha(29) = 7, primitive); the LTE value must be 1 at ALL
#   points (abort otherwise). Negative: for p = 197 | T_7(7) and k = 7*197 = 1379, LTE must predict v_197(T_1379) = 2, and the
#   exact computation must agree.
import json, math, sys, time
from math import isqrt
from collections import Counter
sys.set_int_max_str_digits(2000000)
FN = 'w9_v31_M100000000_H2000_result.json'

# `fund` = fundamental solution (T1, U1) of x^2 - m*y^2 = 1, from the continued fraction of sqrt(m)
def fund(m):
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0

# exact T_k(m): k-th power of T1 + U1*sqrt(m) by binary exponentiation; (ra, rb) = result, (ba, bb) = current square
def T_exact(T1, U1, m, k):
    ra, rb, ba, bb = 1, 0, T1, U1
    while k:
        if k & 1: ra, rb = ra*ba + m*rb*bb, ra*bb + rb*ba
        ba, bb = ba*ba + m*bb*bb, 2*ba*bb; k >>= 1
    return ra

def alpha(T1, p, kmax):
    # `alpha` = rank of apparition: smallest j in 1..kmax with T_j = 0 (mod p); recursion T_{j+1} = 2 T1 T_j - T_{j-1}
    t0, t1 = 1 % p, T1 % p
    if t1 == 0: return 1
    for j in range(2, kmax + 1):
        t0, t1 = t1, (2*T1*t1 - t0) % p
        if t1 == 0: return j
    return None

# `vp` = p-adic valuation v_p(n)
def vp(n, p):
    v = 0
    while n % p == 0: n //= p; v += 1
    return v

t0 = time.time()
d = json.load(open(FN)); hits = d['hits']       # (m, m', k, digits, witness)
print(f"W14 Zeugen-Statistik: {len(hits)} tote Punkte (Zeuge) der Klasse 7, Kern <= {d['MMAX']}, n < 10^{d['H']}")
cache = {}   # kernel m -> (T1, U1)
# `rows`: one tuple per point: (m, m', k, digits, p, alpha, k/alpha, v_p(T_alpha), LTE value v_p(T_alpha) + v_p(k/alpha))
rows = []
for (m, mp, k, dig, p) in hits:
    if m not in cache: cache[m] = fund(m)
    T1, U1 = cache[m]
    a = alpha(T1, p, k); assert a is not None and k % a == 0 and (k // a) % 2 == 1, ('B-W Lemma 3.3 verletzt?', m, k, p, a)
    Ta = T_exact(T1, U1, m, a); va = vp(Ta, p)
    lte = va + vp(k // a, p)
    rows.append((m, mp, k, dig, p, a, k // a, va, lte))
# Positive control 1: the Mollin-Walsh point (`mw` = rows of the point m = 7, k = 7)
mw = [r for r in rows if r[0] == 7 and r[2] == 7]; assert mw and mw[0][4] == 29 and mw[0][5] == 7, ('PK M-W-Punkt', mw)
# Positive control 2: LTE value = 1 at all points (by definition the witness has exponent 1 in T_k)
bad = [r for r in rows if r[8] != 1]; assert not bad, ('LTE-Kontrolle fehlgeschlagen', bad[:3])
# Negative control: p = 197 | T_7(7); v_197(T_{7*197}) must be 2 (LTE), checked by exact computation
T1, U1 = cache[7]; a197 = alpha(T1, 197, 7*197); assert a197 == 7
v_pred = vp(T_exact(T1, U1, 7, 7), 197) + vp(197, 197); v_true = vp(T_exact(T1, U1, 7, 7*197), 197)
assert v_pred == 2 == v_true, ('Negativ-Kontrolle LTE', v_pred, v_true)
print("PK ok: M-W-Punkt (7,7) Zeuge 29 mit alpha(29) = 7 (primitiv); LTE = 1 an allen Punkten; Negativ-Kontrolle v_197(T_1379(7)) = 2 exakt bestaetigt")
prim = [r for r in rows if r[5] == r[2]]   # `prim` = points whose witness is primitive (alpha == k)
print(f"Zeuge primitiv (alpha(p) = k): {len(prim)} von {len(rows)} ({100*len(prim)/len(rows):.1f} %)")
print("Verteilung der kleinsten Zeugen-Primzahl p:", sorted(Counter(r[4] for r in rows).items())[:16], "...")
print(f"groesster kleinster Zeuge: p = {max(r[4] for r in rows)}; Median: {sorted(r[4] for r in rows)[len(rows)//2]}")
print("Kofaktor k/alpha(p) (ungerade) Verteilung (Top 8):", Counter(r[6] for r in rows).most_common(8))
print("v_p(T_alpha) = 1 in allen Faellen:", all(r[7] == 1 for r in rows), "| Punkte mit p | k (Index-Boost moeglich, aber Exponent trotzdem 1):", sum(1 for r in rows if r[2] % r[4] == 0))
print("Beispiele (m, m', k, Stellen, p, alpha, k/alpha):", [(r[0], r[1], r[2], r[3], r[4], r[5], r[6]) for r in rows[:6]])
print(f"Zeit: {time.time()-t0:.1f}s")
