# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w19_teilfaelle_woerterbuch_2026-09-09.py
# Purpose: dictionary of sub-cases (`woerterbuch` = dictionary, `teilfaelle` = sub-cases), two parts.
#   Part A: index types of the 445 lattice points. Part B: middle = power of two.
# Reads: w9_v31_M100000000_H2000_result.json in the current directory (part A). Writes: nothing (printed output only).
# Usage: python w19_teilfaelle_woerterbuch_2026-09-09.py [AMAX=100000] [PMAX=1000000]
#   AMAX = largest exponent a, PMAX = bound for the primes p in the witness search of part B.
# Controls: PK and negative PK in part B (listed below); the run aborts if one fails.
#
# Part A. If k = q is prime (q >= 5), then q is the only divisor >= 5 and k/q = 1, so the escape route p | k/d in
#   Corollary 4.8(iii) is impossible: the condition becomes the pure Wieferich condition p_q^2 | T_q(m).
#
# Part B: middle = power of two (the simplest conceivable sub-case). A triple (2^a - 1, 2^a, 2^a + 1) requires 2^a - 1 AND 2^a + 1
#   to be powerful. Same structure as in our setting, only with the sequence 2^n instead of a Pell unit:
#     p | 2^a - 1  <=>  e := ord_p(2) divides a;  and then  v_p(2^a - 1) = v_p(2^e - 1) + v_p(a/e)  (LTE, p odd).
#     p | 2^a + 1  <=>  e is even and a = e/2 (mod e); and then v_p(2^a + 1) = v_p(2^{e/2} + 1) + v_p(2a/e).
#   So p is a witness (exponent exactly 1) if and only if p^2 does not divide the first term (that is the Wieferich condition
#   to base 2) and p does not divide the cofactor. The loop therefore runs over the PRIMES and marks the a, not the other
#   way round.
#   Controls: 2^3 + 1 = 9 = 3^2 is powerful and must have no witness; 2^3 - 1 = 7 and 2^6 - 1 = 63 = 3^2*7 have the witness 7;
#   the only Wieferich primes below 10^6 are 1093 and 3511 (check of the order computation).

import json, sys, time
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
from math import isqrt
from collections import Counter
AMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 1000000

# `sieve_spf` = smallest-prime-factor sieve; `factor` = factorization {prime: exponent}; `divisors_of` = sorted divisors;
# `is_prime` = Miller-Rabin test with the first twelve primes as bases
def sieve_spf(n):
    spf = list(range(n + 1))
    for i in range(2, isqrt(n) + 1):
        if spf[i] == i:
            for j in range(i*i, n + 1, i):
                if spf[j] == j: spf[j] = i
    return spf
def factor(n, spf):
    f = {}
    while n > 1:
        p = spf[n]; f[p] = f.get(p, 0) + 1; n //= p
    return f
def divisors_of(n):
    ds = []; i = 1
    while i * i <= n:
        if n % i == 0:
            ds.append(i)
            if i != n // i: ds.append(n // i)
        i += 1
    return sorted(ds)
def is_prime(n):
    if n < 2: return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0: return n == p
    d = n - 1; r = 0
    while d % 2 == 0: d //= 2; r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

t0 = time.time()
# ---------------- part A ----------------
d7 = json.load(open('w9_v31_M100000000_H2000_result.json'))
# `kinds` = counts of the index types of k; `prim_idx` = number of points with k prime, k >= 5
kinds = Counter(); prim_idx = 0
for (m, mp, k, dig, status, w) in d7['points']:
    if k == 1: kinds['k = 1 (erste Stufe, AAC-Kern)'] += 1
    elif is_prime(k):
        kinds['k = q prim'] += 1
        if k >= 5: prim_idx += 1
    else: kinds[f'zusammengesetzt, tau(k) = {len(divisors_of(k))}'] += 1
print('Teil A - Index-Typen der 445 Gitterpunkte (Klasse 7, Kern <= 1e8, n < 10^2000):')
for kk, vv in sorted(kinds.items(), key=lambda t: -t[1])[:8]: print(f'   {kk}: {vv}')
print(f'   PRIMZAHL-Index k = q >= 5: {prim_idx} Punkte. Dort ist k/q = 1, der Fluchtweg p | k/d also unmoeglich:')
print('   Corollary 4.8(iii) wird zur reinen Wieferich-Bedingung p_q^2 | T_q(m). Das ist der schaerfste Teilfall.', flush=True)

# ---------------- part B ----------------
spf = sieve_spf(PMAX)
primes = [i for i in range(3, PMAX + 1) if spf[i] == i]
hasw_minus = bytearray(AMAX + 1)  # 1 = witness found for 2^a - 1  (`hasw_minus`; `hasw_plus` is the same for 2^a + 1)
hasw_plus = bytearray(AMAX + 1)
# `wieferich` = Wieferich primes to base 2 found in the loop
wieferich = []
for p in primes:
    e = p - 1
    for q in factor(p - 1, spf):
        while e % q == 0 and pow(2, e // q, p) == 1: e //= q
    if pow(2, e, p * p) == 1:  # p^2 | 2^e - 1: Wieferich-like, no witness via this rank
        wieferich.append(p)
    else:
        a = e
        while a <= AMAX:
            if (a // e) % p: hasw_minus[a] = 1
            a += e
    if e % 2 == 0:  # 2^{e/2} = -1 (mod p)
        h = e // 2
        first_ok = pow(2, h, p * p) != (p * p - 1) % (p * p)  # p^2 does not divide 2^h + 1
        if first_ok:
            a = h
            while a <= AMAX:
                if ((a // h) % 2 == 1) and ((a // h) % p): hasw_plus[a] = 1
                a += 2 * h
print(f'\nTeil B - Mitte = Zweierpotenz. Primzahlen <= {PMAX} ({len(primes)}), a <= {AMAX}. Wieferich-Primzahlen zur Basis 2 gefunden: {wieferich}')
assert wieferich == [1093, 3511], ('PK Wieferich-Liste', wieferich)
assert hasw_minus[3] and hasw_minus[6], 'PK: 2^3-1 und 2^6-1 muessen Zeugen haben'
assert not hasw_plus[3], 'Negativ-PK: 2^3+1 = 9 ist powerful und darf keinen Zeugen haben'
# `om`, `op` = exponents a without a witness for 2^a - 1 resp. 2^a + 1; `beide` = without witness on both sides
om = [a for a in range(2, AMAX + 1) if not hasw_minus[a]]
op = [a for a in range(2, AMAX + 1) if not hasw_plus[a]]
beide = sorted(set(om) & set(op))
print(f'   2^a - 1 ohne Zeugen <= {PMAX}: {len(om)} von {AMAX-1} Exponenten ({100*len(om)/(AMAX-1):.2f} %); kleinste: {om[:15]}')
print(f'   2^a + 1 ohne Zeugen <= {PMAX}: {len(op)} von {AMAX-1} ({100*len(op)/(AMAX-1):.2f} %); kleinste: {op[:15]}')
print(f'   BEIDE ohne Zeugen (einzige Tripel-Kandidaten dieser Form bei dieser Suchtiefe): {len(beide)} -> {beide[:25]}')
prim_unter_ihnen = [a for a in beide[:200] if is_prime(2**a - 1) or is_prime(2**a + 1)]
print(f'   davon mit 2^a - 1 oder 2^a + 1 PRIM (dann trivial nicht powerful, also kein Tripel): {len(prim_unter_ihnen)} -> {prim_unter_ihnen[:15]}')
print('   Lesart: jeder Zeuge schliesst ein Tripel mit dieser Mitte aus. Die verbleibenden a sind bei dieser Suchtiefe unentschieden,')
print('   nicht Kandidaten: ihre Primfaktoren liegen oberhalb der Schranke. Struktur = dieselbe Wieferich-Wand wie in Cor. 4.6/4.8,')
print('   hier zur Basis 2: ein Tripel mit Mitte 2^a verlangt, dass JEDER Primteiler von 2^a - 1 und 2^a + 1 Wieferich-artig ist')
print('   oder den Kofaktor teilt. Unter 10^6 gibt es genau zwei Wieferich-Primzahlen (1093, 3511).')
print(f'Zeit {time.time()-t0:.0f}s')
