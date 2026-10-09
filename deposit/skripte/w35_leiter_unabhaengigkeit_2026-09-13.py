# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w35_leiter_unabhaengigkeit_2026-09-13.py
# Purpose: the "third axis", independence along the LADDER. Tests whether the Wieferich status of rank d is correlated with the
#   status of rank q*d, and measures how strong a coupling would have to be for this test to detect it.
# Context: the earlier measurements looked across primes or across kernels, none along the divisor structure:
#   (1) clustering of Wieferich events WITHIN a kernel (across the primes), z = -0.87;
#   (2) kernel split m = m'*d: do the prime factors of a kernel divide U1 independently? z = +0.05 / -0.10.
# The question: our ladder is k -> k*q. The index gets new divisors, and every divisor d brings its own rank and thereby its own
#   Wieferich condition. Is the status of rank d CORRELATED with the status of rank q*d?
#     Status(d) := 1 if the smallest prime p of rank d satisfies v_p(T_d) >= 2 (i.e. is "eps-Wieferich"), else 0.
#   A coupling would be a new constraint on the ladder.
# Power of the test: a null result without power is silence that looks like knowledge. Therefore, after the test, the script
#   simulates how strong a coupling must be for this test to SEE it (marginal-preserving: with probability rho the status of
#   rank q*d is copied from rank d, otherwise it is drawn independently).
# Reads:  ergebnisse/w9_v31_M100000000_H2000_result.json (class-7 lattice points).
# Writes: ergebnisse/w35_leiter_P<P>_N<NKERNE>_result.json.
# Usage:  python w35_leiter_unabhaengigkeit_2026-09-13.py [P=2000000] [NKERNE=25] [SIM=2000]
#   (P = bound for the primes, NKERNE = number of kernels, SIM = number of simulations per coupling strength)
# Controls (asserted): PK+: (a) m = 7: the primes 29, 197, 2857 have rank 7 and v_p = 1 (not Wieferich);
#   (b) Cor. 4.8(i): every prime found of rank d satisfies p = +-1 (mod 4d), checked at run time.
#   PK-: a value artificially set to p^2 MUST be counted as Wieferich (status 1).

import sys, json, math, random, pathlib, time
from math import isqrt
from collections import Counter, defaultdict

P       = int(sys.argv[1]) if len(sys.argv) > 1 else 2000000   # `P` = bound for the primes p
NKERNE  = int(sys.argv[2]) if len(sys.argv) > 2 else 25        # `NKERNE` = number of kernels (those with most lattice points)
SIM     = int(sys.argv[3]) if len(sys.argv) > 3 else 2000      # `SIM` = simulations per coupling strength

HIER = pathlib.Path(__file__).resolve().parent
ERG  = HIER.parent / 'ergebnisse'
PUNKTE = ERG / 'w9_v31_M100000000_H2000_result.json'   # `PUNKTE` = file with the class-7 lattice points

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
print(f'W35: P = {P} ({len(primes)} Primzahlen), NKERNE = {NKERNE}, SIM = {SIM}; Sieb {time.time()-t0:.0f}s', flush=True)

# ---------------------------------------------------------------- PK
T1_7, U1_7 = fund(7)
for p in (29, 197, 2857):
    o = order(T1_7, U1_7, 7, p, spf)
    assert o % 4 == 0 and o//4 == 7, ('PK Rang', p)
    Td = eps_pow(T1_7, U1_7, 7, 7, p*p)[0]
    assert Td % (p*p) != 0, ('PK: v_p = 1, also KEIN Wieferich', p)
assert (29*29) % (29*29) == 0, 'PK-: ein durch p^2 teilbarer Wert MUSS als Wieferich zaehlen'
print('PK ok: m=7 Raenge und Nicht-Wieferich-Status; Negativseite.\n', flush=True)

d7 = json.load(open(PUNKTE, encoding='utf-8'))   # `d7` = class-7 lattice points
cnt = Counter(pt[0] for pt in d7['points'])
KERNE = [m for m, _ in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[:NKERNE]]   # `KERNE` = kernels with the most points
print('Kerne:', KERNE, flush=True)

# ---------------------------------------------------------------- status for each (kernel, rank)
status = {}          # `status`: (m, d) -> 0/1, the Wieferich status of rank d for kernel m
for m in KERNE:
    T1, U1 = fund(m)
    kleinste = {}    # `kleinste` = rank d -> the smallest prime of that rank
    for p in primes:
        if m % p == 0: continue
        o = order(T1, U1, m, p, spf)
        if o % 4: continue
        d = o // 4
        if d not in kleinste: kleinste[d] = p
    for d, p in kleinste.items():
        assert p % (4*d) in (1, 4*d - 1), ('Cor. 4.8(i) verletzt', m, d, p)
        Td = eps_pow(T1, U1, m, d, p*p)[0]
        status[(m, d)] = 1 if Td % (p*p) == 0 else 0
    print(f'  m = {m:>9}: {len(kleinste)} besetzte Raenge, davon Wieferich {sum(status[(m,d)] for d in kleinste)}', flush=True)

# ---------------------------------------------------------------- pairs (d, q*d) along the ladder
paare = []   # `paare` = pairs (Status(d), Status(q*d)) for all (m, d) where rank q*d is occupied too
for (m, d), s in status.items():
    for q in (3, 5, 7, 11, 13):
        if (m, q*d) in status:
            paare.append((s, status[(m, q*d)]))
# 2x2 table of the pairs: `n11` = both statuses 1, `n10` = first only, `n01` = second only, `n00` = neither; `N` = number of pairs
n11 = sum(1 for a, b in paare if a and b)
n10 = sum(1 for a, b in paare if a and not b)
n01 = sum(1 for a, b in paare if not a and b)
n00 = sum(1 for a, b in paare if not a and not b)
N = len(paare)
print(f'\n=== LEITER-PAARE (d, q*d), beide Raenge besetzt: {N} ===')
print(f'   Status(d)=1 & Status(qd)=1 : {n11}')
print(f'   Status(d)=1 & Status(qd)=0 : {n10}')
print(f'   Status(d)=0 & Status(qd)=1 : {n01}')
print(f'   Status(d)=0 & Status(qd)=0 : {n00}')

p_d  = (n11 + n10) / N if N else 0     # `p_d` = marginal frequency of Status(d) = 1
p_qd = (n11 + n01) / N if N else 0     # `p_qd` = marginal frequency of Status(q*d) = 1
erw11 = N * p_d * p_qd                 # `erw11` = expected n11 under independence
print(f'\n  Randraten: P(Status(d)=1) = {p_d:.5f}   P(Status(qd)=1) = {p_qd:.5f}')
print(f'  n11 beobachtet {n11}, unter UNABHAENGIGKEIT erwartet {erw11:.3f}')

def phi_koeff(a, b, c, d_):   # phi coefficient (correlation of two binary variables) of the 2x2 table (a, b; c, d_)
    nen = math.sqrt((a+b)*(c+d_)*(a+c)*(b+d_))
    return (a*d_ - b*c)/nen if nen > 0 else 0.0
phi = phi_koeff(n11, n10, n01, n00)
chi2 = (N*(n11*n00 - n10*n01)**2) / max(1, (n11+n10)*(n01+n00)*(n11+n01)*(n10+n00))
print(f'  phi-Koeffizient: {phi:+.4f}   chi^2 = {chi2:.3f}   (chi^2 > 3.84 waere p < 0.05)')

# ---------------------------------------------------------------- POWER of the test (detection strength)
print(f'\n=== TRENNSCHAERFE: wie stark muesste eine Kopplung sein, damit DIESER Test sie sieht? ===')
rng = random.Random(20260913)
def sim(rho):
    """Marginal-preserving: with probability rho, Status(qd) is set to Status(d), otherwise drawn independently.
    Returns the fraction of `SIM` simulations in which the chi-square test (> 3.84) detects the coupling."""
    tr = 0
    for _ in range(SIM):
        a = b = c = e = 0
        for (s1, _) in paare:
            s2 = s1 if rng.random() < rho else (1 if rng.random() < p_qd else 0)
            if s1 and s2: a += 1
            elif s1: b += 1
            elif s2: c += 1
            else: e += 1
        n = a+b+c+e
        x2 = (n*(a*e - b*c)**2) / max(1, (a+b)*(c+e)*(a+c)*(b+e))
        if x2 > 3.84: tr += 1
    return tr / SIM
for rho in (0.0, 0.05, 0.1, 0.2, 0.5):
    print(f'   rho = {rho:<5} -> in {100*sim(rho):5.1f} % der Simulationen erkannt', flush=True)

print(f'\n{"="*70}')
print('LESART: ein kleiner phi-Wert bedeutet NUR dann "keine Kopplung", wenn die Trennschaerfe zeigt,')
print('dass eine Kopplung ueberhaupt sichtbar waere.')

out = ERG / f'w35_leiter_P{P}_N{NKERNE}_result.json'
out.write_text(json.dumps(dict(
    skript='w35_leiter_unabhaengigkeit_2026-09-13.py', P=P, NKERNE=NKERNE, SIM=SIM, kerne=KERNE,
    n_paare=N, tafel=dict(n11=n11, n10=n10, n01=n01, n00=n00),
    p_d=p_d, p_qd=p_qd, erwartet_n11=erw11, phi=phi, chi2=chi2,
    trennschaerfe={str(r): sim(r) for r in (0.0, 0.1, 0.2)},
    laufzeit_s=time.time()-t0, python=sys.version.split()[0]), indent=1), encoding='utf-8')
print(f'Ergebnis: {out.name}   Zeit {time.time()-t0:.0f}s')
