# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w22_yokoi_darstellung_2026-09-09.py
# Purpose: test the Yokoi representation of the fundamental unit on our kernels.
# Source: Mollin-Walsh 1986 (C. R. Math. Rep. Acad. Sci. Canada 8(2), p. 112): if a fundamental unit satisfies condition (3) of
#   their theorem, then U1 = 0 (mod m); since m = 7 (mod 8), some prime p | m has p = 3 (mod 4), and Yokoi 1970 (Prop. 1, p. 107)
#   gives a positive integer l with T1 = p^2*l +/- 2 and (U1/p)^2 * m = p^2*l^2 +/- 4l, where T1^2 - U1^2*m = 4.
#
# Two preconditions that the theorem states only implicitly:
#   (V1) NORMALIZATION: Yokoi works with t^2 - u^2*m = 4. For m = 3 (mod 4) the ring of integers is Z[sqrt(m)], hence
#        t = 2*T1, u = 2*U1 with our norm-1 solution T1^2 - m*U1^2 = 1.
#   (V2) p MUST divide u, otherwise (u/p) is not an integer. The theorem is stated in the context U1 = 0 (mod m), i.e. in the AAC
#        (Ankeny-Artin-Chowla) case, or at least in the fast-AAC case. m = 7 has U1 = 3 and 7 does not divide 3, so it is
#        NOT a test point.
#   Therefore the test runs only over pairs (m, p) with p | m, p = 3 (mod 4) and p | u. These pairs are exactly the fast-AAC
#   structure.
#
# Reads: nothing (the kernels m = 7 mod 8 are generated up to MMAX).
# Writes: w22_yokoi_M<MMAX>_result.json in the current directory.
# Usage: python w22_yokoi_darstellung_2026-09-09.py [MMAX]     (default 20000)
# Controls (both abort the run on failure):
#   PK+ : the smallest admissible pair (m, p) found is printed explicitly and the representation must hold for it.
#   PK- : the same pair with a falsified t (t+1) must NOT satisfy the representation.
#
# Own code. Yokoi and Mollin-Walsh are compared only as STATEMENTS; no foreign code was consulted (clean-room rule).
import sys, json
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
from math import isqrt

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20000

def fund(m):
    """Fundamental solution of x^2 - m y^2 = 1 via the continued fraction expansion of sqrt(m)."""
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0

def primfaktoren(n):                            # `primfaktoren` = distinct prime factors of n by trial division
    f = []; d = 2
    while d*d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0: n //= d
        d += 1 if d == 2 else 2
    if n > 1: f.append(n)
    return f

def quadratfrei(n):                             # `quadratfrei` = squarefree test
    d = 2
    while d*d <= n:
        if n % (d*d) == 0: return False
        if n % d == 0: n //= d
        d += 1 if d == 2 else 2
    return True

def yokoi(t, u, m, p):
    """Signs s and l for which BOTH Yokoi equations hold exactly (norm-4 quantities t, u)."""
    out = []
    if u % p: return out
    for s in (1, -1):
        rest = t - s*2
        if rest % (p*p): continue
        l = rest // (p*p)
        if (u//p)**2 * m == p*p*l*l + s*4*l:
            out.append((s, l))
    return out

# ---------------- collect admissible pairs ----------------
paare = []          # `paare` = pairs: (m, p, T1, U1)
for m in range(7, MMAX+1, 8):
    if not quadratfrei(m): continue
    T1, U1 = fund(m)
    if T1 % 2: continue                       # only living families (T1 even)
    t, u = 2*T1, 2*U1                         # Yokoi normalization (t, u) = (2*T1, 2*U1)
    for p in primfaktoren(m):
        if p % 4 == 3 and u % p == 0:
            paare.append((m, p, T1, U1))

print(f'Kerne m = 7 (mod 8), quadratfrei, T1 gerade, m <= {MMAX}')
print(f'Zulaessige Paare (m, p) mit p | m, p = 3 (mod 4) UND p | u  (= Fast-AAC-Struktur): {len(paare)}')
if not paare:
    sys.exit('ABBRUCH: kein zulaessiges Paar im Suchbereich - MMAX erhoehen, bevor irgendetwas behauptet wird.')

# ---------------- PK+ / PK- on the smallest pair ----------------
m0, p0, T0, U0 = paare[0]
r0 = yokoi(2*T0, 2*U0, m0, p0)
print(f'\nPK+ am kleinsten zulaessigen Paar: m={m0}, p={p0}, T1={T0}, U1={U0}  ->  {r0 or "leer"}')
if not r0:
    sys.exit('ABBRUCH PK+: die Darstellung geht am kleinsten zulaessigen Paar nicht auf. '
             'Zitat oder Normierung falsch -> Yokoi 1970 im Volltext, bevor etwas darauf baut.')
if yokoi(2*T0 + 1, 2*U0, m0, p0):
    sys.exit('ABBRUCH PK-: verfaelschtes t geht trotzdem auf - der Test ist entartet.')
print('PK- ok: verfaelschtes t geht nicht auf.\n')

# ---------------- run ----------------
treffer = 0; fehl = []   # `treffer` = hits (representation holds), `fehl` = failures
for (m, p, T1, U1) in paare:
    if yokoi(2*T1, 2*U1, m, p): treffer += 1
    else: fehl.append((m, p, T1, U1))
print(f'Yokoi-Darstellung geht auf: {treffer} von {len(paare)} Paaren.')
if fehl:
    print(f'  Ausnahmen ({len(fehl)}), erste fuenf: {fehl[:5]}')
else:
    print('  Keine Ausnahme.')
print('  Stichproben (m, p, T1, U1, (s, l)):')
for (m, p, T1, U1) in paare[:5]:
    print(f'    m={m:6d} p={p:5d} T1={T1} U1={U1} -> {yokoi(2*T1, 2*U1, m, p)}')

json.dump({'MMAX': MMAX, 'paare': len(paare), 'treffer': treffer, 'fehl': len(fehl)},
          open(f'w22_yokoi_M{MMAX}_result.json', 'w'), indent=1)
print(f'\nErgebnis in w22_yokoi_M{MMAX}_result.json')
