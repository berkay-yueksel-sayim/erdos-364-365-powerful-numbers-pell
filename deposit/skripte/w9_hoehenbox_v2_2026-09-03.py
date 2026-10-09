# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w9_hoehenbox_v2_2026-09-03.py
# W9 v2, "height box" (kernel m <= MMAX, middle n = T_k(m) < 10^H), fast route: for the lattice points of the box it looks for a
#   witness (a prime p with v_p(T_k) = 1) or a parity obstruction (T_k odd) and lists the points that stay open.
#   m' = prod{p | m : p does not divide U_1} does not need the whole fundamental unit, only U_1 mod m.
#   The continued fraction of sqrt(m) runs with small numbers (P, Q, a); the convergents are carried only mod m (for m') and as
#   scaled floats (for log10 T_1). The exact unit is computed only for the few families with candidate points (m' <= kmax).
# Framework: Mollin-Walsh 1986 (criterion p. 110, lemma p. 111). Own code.
# Reads: nothing (kernels m = 7 mod 8, squarefree, are generated). Writes: nothing, the results are printed.
# Usage: python w9_hoehenbox_v2_2026-09-03.py [MMAX] [H]     (defaults 1000000 and 1000)
# Controls (abort on failure): positive controls on the exact and fast routes and on the test function (block below), among them
#   T_7(7) = 2^3 * 29 * 197 * 2857 (Mollin-Walsh 1986, p. 111); a cross-check is built in: the exact route must agree with the
#   fast route mod m and in log10 (abort otherwise).
import sys, time, math
from math import isqrt

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000   # `MMAX` = kernel limit
H    = int(sys.argv[2]) if len(sys.argv) > 2 else 1000      # `H` = height limit in digits of the middle
PMAX = 200000                                               # `PMAX` = prime limit of the witness search

def sieve(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
PRIMES = sieve(PMAX)

# `factor_small(n)` = factorization {prime: exponent} by trial division with PRIMES; a remaining cofactor counts as one factor
def factor_small(n):
    f = {}
    for p in PRIMES:
        if p * p > n: break
        while n % p == 0:
            f[p] = f.get(p, 0) + 1; n //= p
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def fund(m):
    """Exact fundamental solution (T1, U1) of x^2 - m y^2 = 1 (continued fraction)."""
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1
        k1, k0 = k0, a*k0 + k1
    return h0, k0

def fund_fast(m):
    """(T1 mod m, U1 mod m, log10 T1, period length) without large numbers.
    Identity: h_{i-1}^2 - m k_{i-1}^2 = (-1)^i Q_i  =>  Q_i = 1 at even i gives the +1 solution."""
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
        hf1, hf0 = hf0, a*hf0 + hf1
        if hf0 > 1e100:
            hf0 /= 1e100; hf1 /= 1e100; s += 100

def T_exact(T1, U1, m, k):               # exact T_k of the unit T1 + U1*sqrt(m) by square-and-multiply
    ra, rb, ba, bb = 1, 0, T1, U1
    while k:
        if k & 1: ra, rb = ra*ba + m*rb*bb, ra*bb + rb*ba
        ba, bb = ba*ba + m*bb*bb, 2*ba*bb
        k >>= 1
    return ra

def witness(T):                          # smallest prime p in PRIMES with v_p(T) = 1, else None
    for p in PRIMES:
        if T % p == 0:
            if (T // p) % p != 0: return p
    return None

def test_point(T):                       # 'parity' (T odd), ('witness', p), or 'open'
    if T % 2 == 1: return 'parity'
    w = witness(T)
    return ('witness', w) if w else 'open'

# ---------- positive controls (abort on failure) ----------
assert fund(7) == (8, 3)
tm, um, l10, per = fund_fast(7)
assert (tm, um) == (8 % 7, 3 % 7) and abs(l10 - math.log10(8)) < 1e-9, (tm, um, l10)
tm, um, l10, per = fund_fast(2)          # odd period: the +1 solution only appears in period 2 -> (3, 2)
assert (tm, um) == (3 % 2, 2 % 2) and abs(l10 - math.log10(3)) < 1e-9, (tm, um, l10)
for mm in (7, 15, 23, 87, 319, 4999):    # exact vs. fast route on samples
    T1, U1 = fund(mm); tm, um, l10, per = fund_fast(mm)
    assert (T1 % mm, U1 % mm) == (tm, um) and abs(math.log10(T1) - l10) < 1e-6, mm
assert test_point(2**3 * 3**2 * 5**2) == 'open'
assert test_point(2**3 * 3**2 * 5**2 * 7) == ('witness', 7)
assert test_point(3**2 * 5**3) == 'parity'
T7 = T_exact(8, 3, 7, 7)
assert T7 == 2**3 * 29 * 197 * 2857 and test_point(T7) == ('witness', 29)
print("PK ok: fund/fund_fast stimmen (7,2,15,23,87,319,4999); T_7(7)=2^3*29*197*2857 (M-W 1986 S.111); Testfunktion ok", flush=True)

# ---------- run ----------
t0 = time.time()
# `fam` = families (squarefree kernels), `pts` = lattice points tested, `dead_par`/`dead_wit` = points killed by parity / by
# a witness, `cand_fam` = families with candidate points, `open_pts` = points staying open, `hits` = points with a witness
fam = pts = dead_par = dead_wit = 0
cand_fam = 0; max_period = 0
open_pts, hits = [], []
for m in range(7, MMAX + 1, 8):
    if m % 500000 == 7 and m > 7:
        print(f"  ... m = {m}  Familien {fam}  Punkte {pts}  offen {len(open_pts)}  {time.time()-t0:.0f}s", flush=True)
    f = factor_small(m)
    if any(e > 1 for e in f.values()): continue
    fam += 1
    tm, um, l10, per = fund_fast(m)
    if per > max_period: max_period = per
    mp = 1
    for p in f:
        if um % p: mp *= p
    kmax = int((H + 0.31) / l10) + 1      # l10 = log10 T1 <= log10 eps  =>  upper bound for k
    if mp > kmax: continue
    cand_fam += 1
    T1, U1 = fund(m)                       # exact only here
    assert (T1 % m, U1 % m) == (tm, um) and abs(math.log10(T1) - l10) < 1e-6, ("Cross-Check", m)
    k = mp
    while k <= kmax:
        if k % 2 == 1:
            T = T_exact(T1, U1, m, k)
            if len(str(T)) <= H:
                pts += 1
                r = test_point(T)
                if r == 'parity': dead_par += 1
                elif r == 'open': open_pts.append((m, k, len(str(T))))
                else:
                    dead_wit += 1; hits.append((m, mp, k, len(str(T)), r[1]))
        k += mp

print(f"Box: m <= {MMAX}, n < 10^{H}   Familien(squarefree, 7 mod 8): {fam}   Kandidaten-Familien (exakt gerechnet): {cand_fam}   max Periode: {max_period}")
print(f"Gitterpunkte in der Box: {pts}   tot(Paritaet): {dead_par}   tot(Zeuge): {dead_wit}   OFFEN: {len(open_pts)}")
print(f"Zeit: {time.time()-t0:.1f}s")
print("Familien mit Punkten in der Box (m, m', k, Stellen, Zeuge):")
for h in hits: print("  ", h)
print("OFFEN:", open_pts)
