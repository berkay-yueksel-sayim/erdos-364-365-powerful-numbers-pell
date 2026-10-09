# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w17c_rang_abdeckung_tief_2026-09-08.py
# Purpose: deep rank coverage for the witness search at the class-7 lattice points (m, k). For each of the NKERNE kernels with
#   the most lattice points, the primes p <= P are sorted by their rank d (the smallest j with p | T_j, computed as ord_p(ε)/4
#   with ε = T1 + U1·√m). For each rank the first NP primes are kept, and for every divisor d of k the question is: DOES rank d
#   have a witness, i.e. a kept prime p with v_p(T_d) = 1 and p not dividing k/d? (not merely: is the smallest prime of that
#   rank one?)
# Observation that simplifies the evaluation: k is always odd, so every divisor d of k is odd and k/d is odd. The set of ranks
#   that can contribute to T_k is therefore exactly the set of divisors of k, and its size is tau(k). A counterexample with
#   index k would have to fail at every occupied rank d | k simultaneously.
# Output per kernel: sum of tau(k), occupied ranks (a prime <= P was found), ranks with at least one witness, ranks at which
#   ALL examined primes fail (p² | T_d or p | k/d; these structural exceptions would be the interesting cases), and finally
#   the distribution of tau(k) over the points.
# Reads:  w9_v31_M100000000_H2000_result.json (from the current working directory). Writes: nothing (standard output only).
# Usage:  python w17c_rang_abdeckung_tief_2026-09-08.py [P=2000000] [NKERNE=25] [NP=3]
# Controls (asserted): positive: for m = 7 the primes 29, 197, 2857 have rank 7, the prime 3 has no rank, and v_29(T_7) = 1;
#   negative: a value artificially set to a multiple of p² must be recognized as a non-witness.
import json, math, sys, time
from math import isqrt
from collections import Counter
sys.set_int_max_str_digits(2000000)
P = int(sys.argv[1]) if len(sys.argv) > 1 else 2000000      # `P` = bound for the primes p
NKERNE = int(sys.argv[2]) if len(sys.argv) > 2 else 25      # `NKERNE` = number of kernels (those with most lattice points)
NP = int(sys.argv[3]) if len(sys.argv) > 3 else 3           # `NP` = number of primes kept per rank

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

t0 = time.time()
spf = sieve_spf(P + 2)
primes = [i for i in range(3, P + 1) if spf[i] == i]   # odd primes up to P
print(f'W17c: P = {P} ({len(primes)} ungerade Primzahlen), NKERNE = {NKERNE}, NP = {NP} Primzahlen je Rang; Sieb {time.time()-t0:.1f}s', flush=True)
d7 = json.load(open('w9_v31_M100000000_H2000_result.json'))   # `d7` = class-7 lattice points: (m, m', k, digits, status, witness)
cnt = Counter(pt[0] for pt in d7['points'])
KERNE = [m for m, _ in cnt.most_common(NKERNE)]   # `KERNE` = kernels with the most lattice points
print('Kerne:', KERNE, flush=True)
T1, U1 = fund(7)
for p, exp in ((29, 7), (197, 7), (2857, 7)):
    o = order(T1, U1, 7, p, spf); assert o % 4 == 0 and o//4 == exp, ('PK Rang', p, o)
assert order(T1, U1, 7, 3, spf) % 4 != 0, 'Negativ-PK Rang 3'
Td = eps_pow(T1, U1, 7, 7, 29*29)[0]; assert Td % 29 == 0 and Td % (29*29) != 0, 'PK v_29(T_7) = 1'   # `Td` = T_d mod p²
assert not ((29*29) % (29*29) != 0), 'Negativ-PK: ein durch p^2 teilbarer Wert ist kein Zeuge'
print('PK ok (Raenge, Exponent, Negativ-PK)', flush=True)

# `TOT` = running totals ('tau' = sum of tau(k), 'occ' = occupied ranks, 'wit' = ranks with a witness, 'onlyw' = occupied ranks
# without a witness, 'pts' = points); `taus` = tau(k) of every point
TOT = Counter(); taus = []
print(f"{'m':>8} {'Pkt':>4} {'tau-Summe':>10} {'besetzt':>8} {'mit Zeuge':>10} {'nur Wieferich':>14} {'unbesetzt':>10}  Zeit", flush=True)
for m in KERNE:
    tm = time.time(); T1, U1 = fund(m)
    rank_ps = {}   # `rank_ps` = rank d -> the first NP primes of that rank
    for p in primes:
        if m % p == 0: continue
        o = order(T1, U1, m, p, spf)
        if o % 4: continue
        d = o // 4
        lst = rank_ps.setdefault(d, [])
        if len(lst) < NP: lst.append(p)
    pts = [pt for pt in d7['points'] if pt[0] == m]
    # per kernel: sum of tau(k); occupied ranks; ranks with a witness; occupied ranks without one
    tau_sum = occ = wit = onlyw = 0
    for (mm, mp, k, dig, status, w) in pts:
        ds = divisors(k, spf); tau_sum += len(ds); taus.append(len(ds))
        for d in ds:
            ps = rank_ps.get(d)
            if not ps: continue
            occ += 1
            has = False   # `has` = a witness exists at rank d: v_p(T_d) = 1 and p does not divide k/d, so v_p(T_k) = 1
            for p in ps:
                Td = eps_pow(T1, U1, m, d, p*p)[0]
                if Td % (p*p) != 0 and (k // d) % p != 0: has = True; break
            if has: wit += 1
            else: onlyw += 1
    TOT['tau'] += tau_sum; TOT['occ'] += occ; TOT['wit'] += wit; TOT['onlyw'] += onlyw; TOT['pts'] += len(pts)
    print(f'{m:>8} {len(pts):>4} {tau_sum:>10} {occ:>8} {wit:>10} {onlyw:>14} {tau_sum-occ:>10}  {time.time()-tm:.0f}s', flush=True)
print(f"\nGESAMT: {TOT['pts']} Punkte, {TOT['tau']} Raenge (= Summe tau(k)), besetzt {TOT['occ']} ({100*TOT['occ']/TOT['tau']:.1f} %), "
      f"mit Zeuge {TOT['wit']} ({100*TOT['wit']/max(TOT['occ'],1):.2f} % der besetzten), nur Wieferich {TOT['onlyw']}, unbesetzt (Primzahl > P) {TOT['tau']-TOT['occ']}")
ts = sorted(taus)
print(f"tau(k) ueber die Punkte: min {ts[0]}, Median {ts[len(ts)//2]}, Mittel {sum(ts)/len(ts):.1f}, max {ts[-1]}  (= Zahl der Bedingungen je Punkt)")
print(f"Verteilung tau(k):", sorted(Counter(taus).items())[:12], '...')
print(f"Lesart: ein Gegenbeispiel mit Index k muesste an allen tau(k) Raengen gleichzeitig scheitern; besetzte Raenge ohne Zeugen: {TOT['onlyw']}.")
print(f"Verzerrung: unbesetzte Raenge haben ihre kleinste Primzahl oberhalb von {P}; sie sind unbekannt, nicht zeugenfrei. Zeit gesamt {time.time()-t0:.0f}s")
