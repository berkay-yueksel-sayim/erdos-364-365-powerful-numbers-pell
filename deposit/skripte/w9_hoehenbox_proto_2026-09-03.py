# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w9_hoehenbox_proto_2026-09-03.py
# Prototype of the "height box" search for triples of consecutive powerful numbers: kernels m <= MMAX (squarefree, m = 7 mod
#   8),
# middle n = T_k(m) < 10^H. Every lattice point (m, k) is killed by parity or by a witness, or stays open.
# Framework (Mollin and Walsh 1986): middle n = T_k, m | U_k, k odd, T_k even and powerful.
# Lemma L (= Mollin and Walsh, p. 111): m | U_k  <=>  m' | k, where m' = prod{p | m : p does not divide U_1}.
# Height: T_k ~ eps^k / 2  =>  k <= (H + 0.31) / log10(eps).
# Test per point: parity (T_k must be even) + a witness p with v_p(T_k) = 1.
# Usage: python w9_hoehenbox_proto_2026-09-03.py [MMAX] [H]   (defaults 50000 and 400). Reads no file; prints to stdout only
# (ASCII output, for a cp1252 console).
# Controls: positive controls at the output-determining step (asserted, abort on failure): fund(7) = (8, 3), the
#   classification of
# test numbers (powerful -> open, exponent-1 prime found, odd -> parity), and T_7(7) = 2^3*29*197*2857 (Mollin and Walsh 1986,
# p. 111) with witness 29. There is no negative control. Own code, no external source.
import sys, time, math
from math import isqrt

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 50000
H    = int(sys.argv[2]) if len(sys.argv) > 2 else 400
PMAX = 200000   # witness primes up to this bound

# Sieve of Eratosthenes: all primes <= n.
def sieve(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

PRIMES = sieve(PMAX)

# Trial division by the primes up to PMAX; returns {prime: exponent} (a remaining cofactor > 1 is stored as one factor).
def factor_small(n):
    f = {}
    for p in PRIMES:
        if p * p > n: break
        while n % p == 0:
            f[p] = f.get(p, 0) + 1; n //= p
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def fund(m):
    """Fundamental solution (T1, U1) of x^2 - m y^2 = 1 by continued fraction (exact)."""
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1
        k1, k0 = k0, a*k0 + k1
    return h0, k0

def T_exact(T1, U1, m, k):
    """T_k exactly: (T1 + U1 sqrt m)^k = T_k + U_k sqrt m."""
    ra, rb, ba, bb = 1, 0, T1, U1
    while k:
        if k & 1: ra, rb = ra*ba + m*rb*bb, ra*bb + rb*ba
        ba, bb = ba*ba + m*bb*bb, 2*ba*bb
        k >>= 1
    return ra

def witness(T):
    """First prime p < PMAX with v_p(T) = 1, else None."""
    for p in PRIMES:
        if T % p == 0:
            if (T // p) % p != 0: return p
    return None

def test_point(T):
    """'parity' | ('witness', p) | 'open'"""
    if T % 2 == 1: return 'parity'
    w = witness(T)
    return ('witness', w) if w else 'open'

# ---------- Positive controls at the output-determining step (abort on failure) ----------
assert fund(7) == (8, 3), fund(7)
assert test_point(2**3 * 3**2 * 5**2) == 'open', "PK: powerful muss 'open' geben"
assert test_point(2**3 * 3**2 * 5**2 * 7) == ('witness', 7), "PK: Exponent-1-Prim muss gefunden werden"
assert test_point(3**2 * 5**3) == 'parity', "PK: ungerade muss 'parity' geben"
T7 = T_exact(8, 3, 7, 7)
assert T7 == 2**3 * 29 * 197 * 2857, "PK: T_7(7) = 2^3*29*197*2857 (Mollin-Walsh 1986, S.111)"
assert test_point(T7) == ('witness', 29)
print("PK ok: fund(7)=(8,3); T_7(7)=2^3*29*197*2857; Testfunktion trennt powerful/Zeuge/Paritaet")

# ---------- Run ----------
t0 = time.time()
# `fam` = number of kernels (families), `pts` = lattice points in the box, `dead_par` = killed by parity,
# `dead_wit` = killed by a witness, `open_pts` = open points (m, k, digits), `hits` = points killed by a witness.
fam = pts = dead_par = dead_wit = 0
open_pts, hits = [], []
for m in range(7, MMAX + 1, 8):
    f = factor_small(m)
    if any(e > 1 for e in f.values()): continue     # squarefree
    fam += 1
    T1, U1 = fund(m)
    # `mp` = m', the product of the primes p | m with p not dividing U_1; k runs over the odd multiples of `mp`.
    mp = 1
    for p in f:
        if U1 % p: mp *= p
    l10 = math.log10(T1) + math.log10(1 + math.sqrt(1 - 1/(T1*T1)))  # log10(eps)
    # `kmax` = largest k allowed by the height bound.
    kmax = int((H + 0.31) / l10) + 1
    if mp > kmax: continue
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
                    dead_wit += 1
                    hits.append((m, mp, k, len(str(T)), r[1]))
        k += mp
        if mp % 2 == 0: break   # mp even => k never odd ... (mp is a product of odd primes; safeguard only)

print(f"Box: m <= {MMAX}, n < 10^{H}   Familien(squarefree, 7 mod 8): {fam}")
print(f"Gitterpunkte in der Box: {pts}   tot(Paritaet): {dead_par}   tot(Zeuge): {dead_wit}   OFFEN: {len(open_pts)}")
print(f"Zeit: {time.time()-t0:.1f}s")
print("Familien mit Punkten in der Box (m, m', k, Stellen, Zeuge):")
for h in hits: print("  ", h)
if len(hits) > 60: print("   ... insgesamt", len(hits))
print("OFFEN:", open_pts)
