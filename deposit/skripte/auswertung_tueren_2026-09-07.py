# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# auswertung_tueren_2026-09-07.py
# Purpose: three small calculations on existing numbers (`auswertung` = evaluation):
#   (1) 2-adic asymmetry: lattice points (m | U_k, T_k even) for kernels m = 7 mod 8 (middle = 0 mod 4, relevant for triples, i.e.
#       the Erdos-Mollin-Walsh conjecture) versus m = 3 mod 8 (middle = 2 mod 4, not relevant) at m <= 1e6, height 10^1000.
#       Same logic as w9.
#   (2) Exponent statistics: for living families m <= 1e4 (m = 7 mod 8) and primes 3 <= p <= 3000 with rank of apparition
#       alpha(p): counts of v_p(T_alpha) >= 2 (eps-Wieferich) and >= 3, against the heuristic sums of 1/p and 1/p^2.
#   (3) Tower structure of the 371 pairs in w9_v31_M100000000_H2000_result.json: pairs per kernel, k/m', heights.
# Reads: w9_v31_M100000000_H2000_result.json in the current directory. Writes: nothing (printed output only). No arguments.
# Controls: positive controls on the basic functions (fund(7) = (8, 3), T_7(7) = 2^3 * 29 * 197 * 2857, v_29 = 1)
#   abort the run on failure.

import sys, json, math, time
from math import isqrt
import numpy as np
sys.set_int_max_str_digits(2000000)

def sieve(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
# `PRIMES` = primes up to 200000 (`sieve` returns the list of primes up to n)
PRIMES = sieve(200000)
# `factor_small` = factorization {prime: exponent} by trial division with `PRIMES`
def factor_small(n):
    f = {}
    for p in PRIMES:
        if p * p > n: break
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
# `fund` = exact fundamental solution (T1, U1) of x^2 - m*y^2 = 1 by continued fraction (big integers);
# `fund_fast` = (T1 mod m, U1 mod m, log10 T1, period) without big integers
def fund(m):
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
def fund_fast(m):
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    hm1, hm0 = 1 % m, a0 % m; km1, km0 = 0, 1 % m; hf1, hf0 = 1.0, float(a0); s = 0; i = 0
    while True:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q; i += 1
        if Q == 1 and i % 2 == 0: return hm0, km0, math.log10(hf0) + s, i - 1
        hm1, hm0 = hm0, (a*hm0 + hm1) % m; km1, km0 = km0, (a*km0 + km1) % m; hf1, hf0 = hf0, a*hf0 + hf1
        if hf0 > 1e100: hf0 /= 1e100; hf1 /= 1e100; s += 100
# `eps_pow_mod` = T_k mod M;  `T_exact` = exact T_k (both by binary powering of T1 + U1*sqrt(m))
def eps_pow_mod(T1, U1, m, k, M):
    ra, rb = 1 % M, 0; ba, bb = T1 % M, U1 % M
    while k:
        if k & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        k >>= 1
    return ra
def T_exact(T1, U1, m, k):
    ra, rb, ba, bb = 1, 0, T1, U1
    while k:
        if k & 1: ra, rb = ra*ba + m*rb*bb, ra*bb + rb*ba
        ba, bb = ba*ba + m*bb*bb, 2*ba*bb
        k >>= 1
    return ra

# ---------- positive controls ----------
assert fund(7) == (8, 3) and T_exact(8, 3, 7, 7) == 2**3 * 29 * 197 * 2857
assert eps_pow_mod(8, 3, 7, 7, 29**2) % 29 == 0 and eps_pow_mod(8, 3, 7, 7, 29**2) != 0, "v_29(T_7(7)) = 1"
assert eps_pow_mod(8, 3, 7, 1, 8) == 0, "T_1(7) = 8 = 2^3"
print("PK ok", flush=True)

# ---------- (1) 2-adic asymmetry ----------
# `box_count` for the kernels m = residue (mod 8), m <= MMAX: `fam` = squarefree kernels, `cand` = candidate families
# (m' <= kmax), `par` = candidates with T1 odd, `pts` = lattice points with T_k even and at most H digits
def box_count(residue, MMAX, H):
    fam = par = cand = pts = 0
    for m in range(residue, MMAX + 1, 8):
        f = factor_small(m)
        if any(e > 1 for e in f.values()): continue
        fam += 1
        tm, um, l10, per = fund_fast(m)
        mp = 1
        for p in f:
            if um % p: mp *= p
        kmax = int((H + 0.31) / l10) + 1
        if mp > kmax: continue
        cand += 1
        T1, U1 = fund(m)
        if T1 % 2 == 1: par += 1; continue  # whole family without even middles
        k = mp
        while k <= kmax:
            if k % 2 == 1 and len(str(T_exact(T1, U1, m, k))) <= H: pts += 1
            k += mp
    return fam, cand, par, pts
t0 = time.time()
print("=== (1) 2-adische Asymmetrie, Kern <= 1e6, Hoehe 10^1000 ===")
for r, name in ((7, "m = 7 mod 8 (Mitte = 0 mod 4, EMW-relevant)"), (3, "m = 3 mod 8 (Mitte = 2 mod 4, irrelevant)")):
    fam, cand, par, pts = box_count(r, 10**6, 1000)
    print(f"  {name}: Familien {fam}, Kandidaten {cand} (davon T1 ungerade {par}), Paare (T_k gerade) {pts}   [{time.time()-t0:.0f}s]", flush=True)

# ---------- (2) exponent statistics ----------
print("=== (2) Exponent an der Apparition: lebende Familien m <= 1e4 (m = 7 mod 8), Primzahlen 3..3000 ===")
# `PR` = odd primes up to 3000; `fams` = living families: squarefree m = 7 (mod 8), m <= 10^4, with T1 even
PR = [p for p in PRIMES if 2 < p <= 3000]
Pn = np.array(PR, dtype=np.int64)
fams = []
for m in range(7, 10**4 + 1, 8):
    f = factor_small(m)
    if any(e > 1 for e in f.values()): continue
    T1, U1 = fund(m)
    if T1 % 2 == 0: fams.append((m, T1, U1))
# `n_alpha` = number of pairs (m, p) with a rank of apparition; `n_w2`, `n_w3` = number with v_p(T_alpha) >= 2 resp. >= 3;
# `exp2`, `exp3` = the corresponding lists of (m, p, alpha); `h1`, `h2` = heuristic sums of 1/p and 1/p^2
n_alpha = n_w2 = n_w3 = 0; exp2 = []; exp3 = []
h1 = h2 = 0.0
for m, T1, U1 in fams:
    # rank of apparition `a[j]` of each prime = first k <= 3002 with T_k = 0 (mod p) (0 if none), by the recurrence
    # T_{k+1} = 2*T1*T_k - T_{k-1} mod p; `tp`, `tc` = previous and current T_k mod p; `c` = 2*T1 mod p
    T1m = np.array([T1 % p for p in PR], dtype=np.int64); c = (2 * T1m) % Pn
    a = np.zeros(len(PR), dtype=np.int64); tp = np.ones(len(PR), dtype=np.int64); tc = T1m.copy()
    for k in range(1, 3003):
        z = (tc == 0) & (a == 0)
        if z.any(): a[z] = k
        tp, tc = tc, np.mod(c * tc - tp, Pn)
    for j, p in enumerate(PR):
        al = int(a[j])
        if al == 0 or m % p == 0: continue
        n_alpha += 1; h1 += 1.0 / p; h2 += 1.0 / (p * p)
        t2 = eps_pow_mod(T1, U1, m, al, p * p)
        if t2 == 0:
            n_w2 += 1; exp2.append((m, p, al))
            if eps_pow_mod(T1, U1, m, al, p**3) == 0: n_w3 += 1; exp3.append((m, p, al))
print(f"  Familien {len(fams)}, (m,p)-Paare mit Apparitionsrang: {n_alpha}")
print(f"  v_p(T_alpha) >= 2 (eps-Wieferich an der Apparition): {n_w2}   Heuristik sum 1/p: {h1:.0f}")
print(f"  v_p(T_alpha) >= 3: {n_w3}   Heuristik sum 1/p^2: {h2:.1f}")
print("  Beispiele Exponent >= 3 (m, p, alpha):", exp3[:12])
print("  Exponent >= 2, kleinste p:", sorted(exp2, key=lambda t: t[1])[:8], f"  [{time.time()-t0:.0f}s]", flush=True)

# ---------- (3) tower structure of the pairs ----------
print("=== (3) Turm-Struktur der 371 Paare (Kern <= 1e8, n < 10^2000) ===")
# reads the w9 result file; a point is (m, m', k, digits, status, witness); `pairs` = points with status 'witness'
d = json.load(open('w9_v31_M100000000_H2000_result.json'))
pairs = [p for p in d['points'] if p[4] == 'witness']
from collections import Counter
perk = Counter(p[0] for p in pairs)
print("  Paare pro Kern (Top 12):", perk.most_common(12))
print("  Kerne mit genau 1 Paar:", sum(1 for v in perk.values() if v == 1), "von", len(perk))
ratio = Counter(p[2] // p[1] for p in pairs)  # k / m' (odd multiples)
print("  k/m' Verteilung (Top 8):", ratio.most_common(8))
digs = sorted(p[3] for p in pairs)
print("  Hoehen (Stellen) Quantile 10/50/90 %:", digs[len(digs)//10], digs[len(digs)//2], digs[9*len(digs)//10])
# count N(X) = pairs with at most X digits against 0.16*ln(10^X) (half of the corrected heuristic, only middle = 0 mod 4)
for X in (30, 100, 300, 1000, 2000):
    N = sum(1 for p in pairs if p[3] <= X)
    print(f"  Paare mit <= {X} Stellen (Kern <= 1e8): {N}   Heuristik 0.16*ln(10^X) = {0.16*X*math.log(10):.0f}")
print(f"Zeit gesamt: {time.time()-t0:.0f}s")
