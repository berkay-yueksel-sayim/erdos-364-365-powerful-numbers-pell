# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w18_rang_kongruenz_2026-09-09.py
# Purpose: check the rank congruence and evaluate the "escape route" of the witness criterion (Corollary 4.8(iii)).
# Claim 1 (classical, easy to verify): if p has rank d in the T-sequence (smallest j with p | T_j), then its rank in the
#   associated Lucas sequence is 2d and divides p − (m|p) (Legendre symbol). So 2d | p ∓ 1, i.e. p ≡ ±1 (mod 2d); in
#   particular p ≥ 2d − 1. The script checks the sharper form p ≡ ±1 (mod 4d): the order of ε = T1 + U1√m modulo p is 4d
#   and divides p − (m|p).
# Claim 2 (consequence for Corollary 4.8(iii)): the escape route "p | k/d" needs p ≤ k/d, hence k/d ≥ 4d − 1, hence
#   k ≥ d(4d − 1). For every divisor d of k with d(4d − 1) > k the escape route is closed, and the condition becomes the
#   pure Wieferich condition p² | T_d.
# The script verifies Claim 1 on the prime/rank pairs (p, d) of the class-7 kernels with the most lattice points, and counts,
#   over all class-7 lattice points, how many divisors d ≥ 5 of k have no escape route.
# Reads:  w9_v31_M100000000_H2000_result.json (from the current working directory). Writes: nothing (standard output only).
# Usage:  python w18_rang_kongruenz_2026-09-09.py [P=1000000] [NKERNE=15]   (P = prime bound, NKERNE = number of kernels)
# Controls (asserted): positive: for m = 7 the primes 29, 197, 2857 have rank 7 and satisfy p ≡ ±1 (mod 28 = 4d);
#   negative: the prime 3 has no rank for m = 7.
import json, math, sys, time
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
from math import isqrt
from collections import Counter
P = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000          # `P` = bound for the primes p
NKERNE = int(sys.argv[2]) if len(sys.argv) > 2 else 15          # `NKERNE` = number of kernels (those with most lattice points)

def sieve_spf(n):   # smallest-prime-factor sieve up to n
    spf = list(range(n + 1))
    for i in range(2, isqrt(n) + 1):
        if spf[i] == i:
            for j in range(i*i, n + 1, i):
                if spf[j] == j: spf[j] = i
    return spf
def factor(n, spf):   # prime factorization of n as {prime: exponent}, using the sieve `spf`
    f = {}
    while n > 1:
        p = spf[n]; f[p] = f.get(p, 0) + 1; n //= p
    return f
def divisors(n, spf):   # sorted list of all divisors of n
    ds = [1]
    for q, e in factor(n, spf).items():
        ds = [d * q**i for d in ds for i in range(e + 1)]
    return sorted(ds)
def fund(m):   # fundamental solution (h0, k0) = (T1, U1) of x² − m·y² = 1, by the continued fraction of √m
    a0 = isqrt(m); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (m - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
def eps_pow(T1, U1, m, e, M):   # (T_e, U_e) mod M: the e-th power of ε = T1 + U1·√m, by repeated squaring
    ra, rb, ba, bb = 1 % M, 0, T1 % M, U1 % M
    while e:
        if e & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        e >>= 1
    return ra, rb
# multiplicative order of ε modulo p: starts from p − (m|p), then removes prime factors while possible
def order(T1, U1, m, p, spf):
    N = p - 1 if pow(m % p, (p - 1)//2, p) == 1 else p + 1
    o = N
    for q in factor(N, spf):
        while o % q == 0 and eps_pow(T1, U1, m, o//q, p) == (1, 0): o //= q
    return o

t0 = time.time(); spf = sieve_spf(P + 2)
primes = [i for i in range(3, P + 1) if spf[i] == i]   # odd primes up to P
T1, U1 = fund(7)
for p, d in ((29, 7), (197, 7), (2857, 7)):
    o = order(T1, U1, 7, p, spf); assert o % 4 == 0 and o//4 == d and p % (4*d) in (1, 4*d - 1), ('PK Kongruenz', p, o, p % (4*d))
assert order(T1, U1, 7, 3, spf) % 4 != 0, 'Negativ-PK: 3 hat keinen Rang fuer m = 7'
print(f'W18: P = {P} ({len(primes)} Primzahlen). PK ok: 29, 197, 2857 haben Rang 7 und sind ≡ ±1 (mod 28 = 4d).', flush=True)

d7 = json.load(open('w9_v31_M100000000_H2000_result.json'))   # `d7` = class-7 lattice points: (m, m', k, digits, status, witness)
# `KERNE` = the NKERNE kernels with the most lattice points
cnt = Counter(pt[0] for pt in d7['points']); KERNE = [m for m, _ in cnt.most_common(NKERNE)]
# `paare` = pairs (p, d) examined; `ok_kong` = those satisfying the congruence; `verstoss` = violations (must stay empty)
paare = 0; ok_kong = 0; verstoss = []
for m in KERNE:
    T1, U1 = fund(m)
    for p in primes:
        if m % p == 0: continue
        o = order(T1, U1, m, p, spf)
        if o % 4: continue
        d = o // 4; paare += 1
        if p % (4*d) in (1, 4*d - 1): ok_kong += 1
        else: verstoss.append((m, p, d, p % (4*d)))
print(f'Behauptung 1 VERSCHAERFT (p ≡ ±1 mod 4d, da ord_p(eps) = 4d | p - (m|p)) an {paare} (p, d)-Paaren ueber {len(KERNE)} Kerne: erfuellt {ok_kong}, Verstoesse {len(verstoss)} {verstoss[:5]}')
assert not verstoss, 'Rang-Kongruenz verletzt'

# Claim 2: escape-route statistics over all class-7 lattice points
# `tot` = divisors d ≥ 5 in total; `zu` = escape route closed; `offen` = escape route possible;
# `klein` = number of excluded divisors (d = 1, 3)
tot = zu = offen = klein = 0; per_point = []
for (m, mp, k, dig, status, w) in d7['points']:
    ds = [d for d in divisors(k, spf) if d >= 5]
    z = sum(1 for d in ds if d*(4*d - 1) > k)
    tot += len(ds); zu += z; offen += len(ds) - z; per_point.append((len(ds), z))
    klein += len(divisors(k, spf)) - len(ds)
print(f'Behauptung 2 ueber alle {len(d7["points"])} Gitterpunkte der Klasse 7: Teiler d >= 5 insgesamt {tot} '
      f'(ausgeschlossen d in 1,3: {klein}); davon **ausweglos** (k < d(4d-1), Fluchtweg p | k/d unmoeglich): {zu} ({100*zu/tot:.1f} %), '
      f'mit moeglichem Fluchtweg: {offen}')
q = sorted(z for _, z in per_point)
print(f'ausweglose Bedingungen je Punkt: min {q[0]}, Median {q[len(q)//2]}, Mittel {sum(q)/len(q):.1f}, max {q[-1]}')
q2 = sorted(a for a, _ in per_point)
print(f'Bedingungen je Punkt insgesamt (Teiler d >= 5): min {q2[0]}, Median {q2[len(q2)//2]}, max {q2[-1]}')
print(f'Zeit {time.time()-t0:.0f}s')
