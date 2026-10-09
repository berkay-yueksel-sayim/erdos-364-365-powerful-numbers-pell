# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w26_gap1_hoehe_und_dritte_zahl_2026-09-10.py
# Purpose: pairs (a, a+1) of consecutive powerful numbers up to height 10^H via the Pell tower of the kernel D, and then the
#   THIRD number: if two consecutive powerful numbers are found, the neighbour (a-1 or a+2) is tested for a witness; if it were
#   powerful too, there would be a triple.
# Reads: nothing (the kernels D are generated). Writes: w26_gap1_D<DMAX>_H<H>_result.json in the current directory.
# Usage: python w26_gap1_hoehe_und_dritte_zahl_2026-09-10.py [DMAX] [H] [P]   (defaults 100000 2000 200000)
#   D <= DMAX, height 10^H, witness search among the primes below P (`PMAX`).
# Controls (all abort on failure): PK-1 to PK-4, described below.
#
# Structure: (a, a+1) both powerful  <=>  T := 2a+1, T^2 - D U^2 = 1 with D = kernel(T^2-1) and D | U.
#   a = (T-1)/2, a+1 = (T+1)/2. This is the odd branch of the same Pell family as the pairs at distance two.
#
# Three filters, from cheap to expensive:
#   (F1) INDEX:   {k : D | U_k} = multiples of D' := product of the p | D with p not dividing U_1  (analogue of Lemma L).
#   (F2) 2-ADIC:  T_k = +-1 (mod 8). Exactly one of the two numbers is even; its 2-adic valuation is v2(T_k -+ 1) - 1
#         and must be >= 2. T = 3 (pair (1,2)) correctly drops out, all known pairs pass.
#   (F3) WITNESS: an odd prime p with T_k = 1 (mod p) but not (mod p^2)  =>  v_p(a) = 1  =>  a is not powerful.
#         Analogously with -1 for a+1. Factoring is impossible at 10^2000; one witness suffices.
#
# Third number: if a pair (a, a+1) survives all filters, it would extend to a triple (a-1, a, a+1) or (a, a+1, a+2),
#   with a-1 = (T_k-3)/2 and a+2 = (T_k+3)/2. For odd p, v_p(a-1) = v_p(T_k-3) and v_p(a+2) = v_p(T_k+3), so it is the
#   same witness search, against 3 instead of 1.
#
# Controls:
#   PK-1: analogue of Lemma L at the known pairs (k must be a multiple of D').
#   PK-2: the known pairs must PASS F1 and F2 and SURVIVE the pair witness search.
#   PK-3 (negative): D = 2, k = 1 (T = 3, pair (1,2)) must fail at F2.
#   PK-4: for EVERY known pair the THIRD number must have a witness, since there is no triple there. If the script finds none,
#         the third-number machinery is broken.

import sys, math, json, time
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
from math import isqrt

# `DMAX` = largest kernel D, `H` = height bound (log10), `PMAX` = bound for the witness primes
DMAX = int(float(sys.argv[1])) if len(sys.argv) > 1 else 100000
H    = float(sys.argv[2]) if len(sys.argv) > 2 else 2000.0
PMAX = int(float(sys.argv[3])) if len(sys.argv) > 3 else 200000

# `sieb` = sieve: the odd primes up to n; `PRIMES` = odd primes up to `PMAX`
def sieb(n):
    s = bytearray(b'\x01') * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(3, n + 1, 2) if s[i]]
PRIMES = sieb(PMAX)

# `kern` = kernel: the squarefree part of n
def kern(n):
    m = 1; d = 2
    while d * d <= n:
        e = 0
        while n % d == 0: n //= d; e += 1
        if e % 2: m *= d
        d += 1 if d == 2 else 2
    return m * (n if n > 1 else 1)

# `primfaktoren` = the distinct prime factors of n
def primfaktoren(n):
    f = []; d = 2
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0: n //= d
        d += 1 if d == 2 else 2
    if n > 1: f.append(n)
    return f

def fund_mod(D, M):
    """(T1 mod M, U1 mod M, log10 T1) without big integers; None if D is a square."""
    a0 = isqrt(D)
    if a0 * a0 == D: return None
    P, Q, a = 0, 1, a0
    hm1, hm0 = 1 % M, a0 % M; km1, km0 = 0, 1 % M
    hf1, hf0 = 1.0, float(a0); s = 0; i = 0
    while True:
        P = a * Q - P; Q = (D - P * P) // Q; a = (a0 + P) // Q; i += 1
        if Q == 1 and i % 2 == 0:
            return hm0, km0, math.log10(hf0) + s
        hm1, hm0 = hm0, (a * hm0 + hm1) % M; km1, km0 = km0, (a * km0 + km1) % M
        hf1, hf0 = hf0, a * hf0 + hf1
        if hf0 > 1e100: hf0 /= 1e100; hf1 /= 1e100; s += 100

def pot_mod(T1, U1, D, k, M):
    """(T_k mod M, U_k mod M) by fast exponentiation of T1 + U1*sqrt(D)."""
    ra, rb = 1 % M, 0
    ba, bb = T1 % M, U1 % M
    while k:
        if k & 1: ra, rb = (ra*ba + D*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + D*bb*bb) % M, (2*ba*bb) % M
        k >>= 1
    return ra, rb

def zeuge(T1, U1, D, k, ziel, primes):
    """Search for p with T_k = ziel (mod p) but not (mod p^2)  ->  v_p((T_k - ziel)/2) = 1.
    Returns p or None.  (`zeuge` = witness, `ziel` = target value)"""
    for p in primes:
        if D % p == 0: continue
        pp = p * p
        t, _ = pot_mod(T1, U1, D, k, pp)
        if (t - ziel) % p == 0 and (t - ziel) % pp != 0:
            return p
    return None

# `fund_exakt` = exact fundamental solution (T1, U1) of x^2 - D*y^2 = 1 by continued fraction (big integers)
def fund_exakt(D):
    a0 = isqrt(D); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - D*k0*k0 != 1:
        P = a*Q - P; Q = (D - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0

t0 = time.time()
print(f'W26 -- Abstand-1 in die Hoehe + dritte Zahl.  D <= {DMAX}, Hoehe 10^{H:.0f}, Zeugen unter {PMAX}')

# ================= CONTROLS at the known pairs =================
# `BEKANNT` = known pairs: (a, kernel D of (2a+1)^2 - 1)
BEKANNT = [(8, 2), (288, 2), (675, 3), (9800, 2), (12167, 46), (235224, 6), (332928, 2), (465124, 5)]
print('\n--- PK-1 bis PK-4 an den bekannten Abstand-1-Paaren ---')
for a, D_soll in BEKANNT:
    T = 2*a + 1
    D = kern(T*T - 1)
    assert D == D_soll, f'PK-1: Kern von {T}^2-1 ist {D}, erwartet {D_soll}'
    T1, U1 = fund_exakt(D)
    Dp = 1
    for p in primfaktoren(D):
        if U1 % p: Dp *= p
    # find the index k with T_k = T
    tt, uu, k = T1, U1, 1
    while tt < T: tt, uu = T1*tt + D*U1*uu, T1*uu + U1*tt; k += 1
    assert tt == T, f'PK-1: {T} liegt nicht im Turm von D={D}'
    assert k % Dp == 0, f'PK-1: k={k} ist kein Vielfaches von D\'={Dp} (D={D})'
    assert T % 8 in (1, 7), f'PK-2: T={T} faellt am 2-adischen Filter'
    zp  = zeuge(T1, U1, D, k, 1, PRIMES)  # a not powerful?
    zp1 = zeuge(T1, U1, D, k, -1, PRIMES)  # a+1 not powerful?
    assert zp is None and zp1 is None, f'PK-2: fuer das echte Paar a={a} wurde faelschlich ein Zeuge gefunden ({zp}, {zp1})'
    # PK-4: the THIRD number must have a witness (`z_links` / `z_rechts` = left / right)
    z_links  = zeuge(T1, U1, D, k, 3, PRIMES)    # a-1 = (T-3)/2
    z_rechts = zeuge(T1, U1, D, k, -3, PRIMES)   # a+2 = (T+3)/2
    assert z_links is not None, f'PK-4: kein Zeuge fuer a-1 bei a={a}'
    assert z_rechts is not None, f'PK-4: kein Zeuge fuer a+2 bei a={a}'
    print(f'  a={a:>8} D={D:>4} D\'={Dp:>4} k={k:>3} T mod 8={T%8}  |  Paar ueberlebt ✓  |  a-1 tot durch p={z_links}, a+2 tot durch p={z_rechts}')
# PK-3 (negative)
T1, U1 = fund_exakt(2)
t3, _ = pot_mod(T1, U1, 2, 1, 8)
assert t3 % 8 not in (1, 7), 'PK-3: T=3 muesste am 2-adischen Filter scheitern'
print(f'  PK-3 ok (negativ): D=2, k=1 gibt T=3 = {t3 % 8} (mod 8) -> korrekt abgewiesen (Paar (1,2), 2 ist nicht powerful).')
print(f'  Alle Kontrollen bestanden ({time.time()-t0:.0f}s).')

# ================= SCAN =================
print(f'\n--- Scan: D <= {DMAX}, Hoehe 10^{H:.0f} ---')
LOG2 = math.log10(2.0)
# `fam` = squarefree kernels, `kand` = candidate pairs, `paar_tot` = pairs with one side dead through 2,
# `ueberlebt` = pairs whose remaining side got a witness, `offen` = remaining side without a witness (`offen` = open)
fam = 0; kand = 0; paar_tot = 0; ueberlebt = []; offen = []
for D in range(2, DMAX + 1):
    if kern(D) != D: continue
    M = 8 * D
    fm = fund_mod(D, M)
    if fm is None: continue
    fam += 1
    tM, uM, l10 = fm
    # `leps` = log10(eps), where eps = T1 + U1*sqrt(D)
    leps = l10 + LOG2 if l10 >= 15 else math.log10(round(10.0**l10) + math.sqrt(round(10.0**l10)**2 - 1))
    Dp = 1
    for p in primfaktoren(D):
        if uM % p: Dp *= p
    # `kmax` = largest index k with height at most 10^H
    kmax = int(H / leps) + 1
    if Dp > kmax: continue
    T1m, U1m = tM, uM
    k = Dp
    while k <= kmax:
        t8, _ = pot_mod(T1m, U1m, D, k, 8)
        if t8 % 8 in (1, 7):  # F2 -> (a, a+1) is a genuine pair (F1+F2 are sufficient)
            kand += 1
            # F2b (free): one side is always already dead through the prime 2.
            #   T = 1 (mod 8): T+3 = 4 (mod 8) => v2((T+3)/2) = 1 => a+2 not powerful.
            #   T = 7 (mod 8): T-3 = 4 (mod 8) => v2((T-3)/2) = 1 => a-1 not powerful.
            # Only the other side needs a witness search. Its middle is always = 0 (mod 4) and would thus be
            # a point of our distance-2 map.
            paar_tot += 1
            ziel   = 3 if t8 % 8 == 1 else -3  # the LIVING side: a-1 resp. a+2
            seite  = 'a-1' if t8 % 8 == 1 else 'a+2'  # `seite` = side, `tot_gratis` = the side that is dead for free
            tot_gratis = 'a+2' if t8 % 8 == 1 else 'a-1'
            T1e, U1e = fund_exakt(D)
            z = zeuge(T1e, U1e, D, k, ziel, PRIMES)
            stellen = k * leps  # `stellen` = digits (log10 of T_k)
            if z is not None:
                ueberlebt.append((D, k, round(stellen, 1), seite, z, tot_gratis))
            else:
                offen.append((D, k, round(stellen, 1), seite, None, tot_gratis))
        k += Dp

print(f'  quadratfreie D: {fam} | **Abstand-1-Paare** (F1+F2, hinreichend): {kand}')
print(f'  davon eine Tripel-Seite gratis tot durch p = 2: {paar_tot} (= alle)')
print(f'  ✅ die verbleibende Seite durch einen Zeugen erledigt: {len(ueberlebt)}')
for r in ueberlebt[:12]:
    print(f'     D={r[0]:>7} k={r[1]:>5} ~10^{r[2]:<9} {r[3]} tot durch p={r[4]:<9} ({r[5]} war schon durch 2 tot)')
print(f'  🔴 OFFEN (verbleibende Seite ohne Zeugen unter {PMAX}): {len(offen)}')
for r in offen[:12]:
    print(f'     D={r[0]:>7} k={r[1]:>5} ~10^{r[2]:<9} {r[3]}: kein Zeuge  ({r[5]} war schon durch 2 tot)')
print(f'  Lesart: "offen" heisst NICHT "Tripel", sondern "Zeuge unter {PMAX} nicht gefunden".')
print(f'  Erwartung aus unserer Abstand-2-Karte: eine ueberlebende Seite waere ein Punkt mit Mitte = 0 (mod 4),')
print(f'  also genau ein Punkt jener Karte -- und die ist bis 10^2000 fuer Kerne <= 10^8 leer.')

res = {'DMAX': DMAX, 'H': H, 'PMAX': PMAX, 'familien': fam, 'kandidaten': kand,
       'paar_tot': paar_tot, 'ueberlebt': ueberlebt, 'offen': offen, 'sekunden': round(time.time()-t0, 1)}
fn = f'w26_gap1_D{DMAX}_H{int(H)}_result.json'
json.dump(res, open(fn, 'w'), indent=1)
print(f'\nErgebnis in {fn}  ({time.time()-t0:.0f}s)')
