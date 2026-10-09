# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w8_pk_negativseite_2026-09-06.py
# Negative-side control for w8: the claim "alpha(p) does not exist" (3.7 million exclusions) and the values alpha(p) > 0
# are checked by an INDEPENDENT route.
#
# Here alpha(p) is the least k >= 1 with T_k = 0 (mod p), where T_k + U_k sqrt(m) = eps^k and eps = T_1 + U_1 sqrt(m).
# Route A (as in w8): the recurrence T_{k+1} = 2 T_1 T_k - T_{k-1} mod p, vectorized over all primes p.
# Route B (independent):
#   T_k = 0 (mod p)  <=>  eps^(2k) = -1 in (Z[sqrt m]/p)^*   (Norm(eps) = 1, hence eps_bar = eps^-1)
#   =>  such a k exists  <=>  4 | ord_p(eps), and then alpha(p) = ord_p(eps)/4.
# ord_p(eps) is determined from the factorization of p-1 (m a quadratic residue mod p) or p+1 (m a non-residue).
# Sample: 300 random squarefree kernels m < 1e5 with m = 7 (mod 8), each with all primes 3 <= p <= 2000 not dividing m. Own code.
#
# Reads:    nothing. Writes: nothing (printed output only; fixed random seed).
# Usage:    python w8_pk_negativseite_2026-09-06.py   (no arguments)
# Controls: positive control of the check itself (m = 7: ord_29(eps) = 28, ord_11(eps) = 12); any disagreement between the
#           two routes ("negative PK fired") aborts the run with an assertion error.
import sys, random, math, time
from math import isqrt
import numpy as np
random.seed(20260906)
MMAX, PMAX = 100000, 2000   # bound for the kernels m and for the primes p

def sieve(n):   # sieve of Eratosthenes: all primes <= n
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
PR = sieve(PMAX)   # `PR` = primes up to `PMAX`
def factor_small(n):   # factorization of n by trial division with the primes `PR`
    f = {}
    for p in PR:
        if p * p > n: break
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
def factor_full(n):          # n <= 2001, trial division suffices
    f = {}; d = 2
    while d * d <= n:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
def fund(m):   # `fund` = fundamental solution (T_1, U_1) of x^2 - m y^2 = 1 from the continued fraction of sqrt(m)
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1
        k1, k0 = k0, a*k0 + k1
    return h0, k0
def eps_pow(T1, U1, m, k, M):
    """(a, b) with eps^k = a + b sqrt(m) mod M."""
    ra, rb = 1 % M, 0; ba, bb = T1 % M, U1 % M
    while k:
        if k & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        k >>= 1
    return ra, rb
def order_eps_mod_p(T1, U1, m, p):
    """Order of eps in (Z[sqrt m]/p)^*; it divides p-1 (m a quadratic residue) or p+1 (m a non-residue)."""
    n = p - 1 if pow(m % p, (p - 1) // 2, p) == 1 else p + 1
    o = n
    for q, e in factor_full(n).items():
        for _ in range(e):
            if eps_pow(T1, U1, m, o // q, p) == (1, 0): o //= q
            else: break
    return o

# Route A (as in w8): recurrence, vectorized
PRn = np.array([p for p in PR if p > 2], dtype=np.int64)   # `PRn` = odd primes as a numpy array
def alphas_recurrence(T1):   # array of alpha(p) for all p in `PRn` (0 if no k was found within the loop range)
    T1m = np.array([T1 % int(p) for p in PRn], dtype=np.int64); c = (2 * T1m) % PRn
    a = np.zeros(len(PRn), dtype=np.int64); t_prev = np.ones(len(PRn), dtype=np.int64); t_cur = T1m.copy()
    for k in range(1, int(PRn.max()) + 3):
        z = (t_cur == 0) & (a == 0)
        if z.any(): a[z] = k
        t_prev, t_cur = t_cur, np.mod(c * t_cur - t_prev, PRn)
    return a

# Positive control of the check itself: m = 7, p = 29: ord = 28 -> alpha = 7 (Mollin-Walsh number); p = 7 | m: no order needed
T1, U1 = fund(7)
assert order_eps_mod_p(8, 3, 7, 29) == 28 and order_eps_mod_p(8, 3, 7, 11) == 12
print("PK ok: ord_29(eps_7) = 28 -> alpha = 7 (M-W); ord_11(eps_7) = 12 -> alpha = 3", flush=True)

fams = []   # `fams` = the sampled families (kernels m)
while len(fams) < 300:
    m = random.randrange(7, MMAX, 8)
    if all(e == 1 for e in factor_small(m).values()): fams.append(m)
t0 = time.time()
# counters: `n_pairs` = pairs (m, p) tested; `n_none_ok`/`n_none_bad` = "alpha does not exist" confirmed/refuted by route B;
#   `n_val_ok`/`n_val_bad` = values alpha(p) = ord/4 confirmed/refuted
n_pairs = n_none_ok = n_none_bad = n_val_ok = n_val_bad = 0
for m in fams:
    T1, U1 = fund(m)
    a = alphas_recurrence(T1)
    for i, p in enumerate(PRn):
        p = int(p)
        if m % p == 0: continue
        n_pairs += 1
        o = order_eps_mod_p(T1, U1, m, p)
        alpha_B = o // 4 if o % 4 == 0 else 0   # `alpha_B` = alpha(p) from route B (0 = does not exist)
        if a[i] == 0:
            if alpha_B == 0: n_none_ok += 1
            else: n_none_bad += 1; print("  WIDERSPRUCH none:", m, p, "ord", o)
        else:
            if alpha_B == a[i]: n_val_ok += 1
            else: n_val_bad += 1; print("  WIDERSPRUCH wert:", m, p, "rek", int(a[i]), "ord/4", alpha_B)
print(f"Stichprobe: {len(fams)} Familien m <= {MMAX}, Primzahlen 3..{PMAX}: {n_pairs} Paare (m, p)")
print(f"'alpha existiert nicht' (Rekursion) bestaetigt durch 4 !| ord: {n_none_ok}   widerlegt: {n_none_bad}")
print(f"alpha(p) = ord/4 bestaetigt: {n_val_ok}   widerlegt: {n_val_bad}")
print(f"Zeit: {time.time()-t0:.0f}s")
assert n_none_bad == 0 and n_val_bad == 0, "NEGATIV-PK GEFEUERT"
