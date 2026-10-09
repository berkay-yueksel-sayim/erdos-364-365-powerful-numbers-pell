# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# -*- coding: utf-8 -*-
# w163_pocklington_beweise_2026-10-02.py
# Primality proofs (Pocklington) for the probable primes of up to 80 digits in the inventory of w162.
# Motivation: the inventory w162 lists 24 numbers that are only probably prime; 18 of them have at most 75 digits.
# Method (own code): Pocklington's theorem (1914) in the form: N - 1 = F*R with F completely factored, F^2 > N, and for every
#   prime q | F an a with a^(N-1) = 1 (mod N) and gcd(a^((N-1)/q) - 1, N) = 1 imply that N is prime. (Every prime divisor p of N
#   then satisfies p = 1 mod F, hence p > sqrt(N).) The prime factors q of F are handled below psi_12 by the deterministic
#   12-base test and above it RECURSIVELY by the same method. The factorization of N - 1 is delivered by sympy.factorint
#   (library call); this is uncritical, because the certificate checks every factor and the product itself. A number for which
#   no certificate is found stays "probable" and is reported.
# SECOND CHECKER (`pruefe`, written separately, factors nothing): takes a certificate as stored in the JSON and checks only
#   the product, the powers, the gcd, F^2 > N and the recursion.
# Reads: ergebnisse/w162_inventar_wahrscheinliche_primzahlen_result.json.
# Writes: ergebnisse/w163_pocklington_beweise_result.json and the matching _output.txt. No command-line arguments.
# Controls (expectations fixed in advance):
#   E1 (PK)  The next prime after 10^40 gets a certificate that the second checker accepts.
#   E2 (NK)  The product of two 20-digit primes gets NO certificate; a certificate with one exponent changed and a certificate
#            with a wrong N (N + 2) are rejected by the second checker.
#   E3       Every number of up to 80 digits with status "only probable" in w162 gets a certificate, or is reported as not
#            achieved (printed, not asserted).
import sys, json, time, pathlib, random
from math import gcd, isqrt
import sympy
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'
# `aus` = collected output lines; `sag` (= say) prints a line and records it for the output file.
aus = []
def sag(s=''):
    print(s, flush=True); aus.append(s)
# `PSI12` = bound below which the Miller-Rabin test with the first 12 prime bases (`BASEN`) is deterministic;
# `mr12` = that test (with trial division by the bases first).
PSI12 = 318665857834031151167461
BASEN = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
def mr12(n):
    if n < 2: return False
    for p in BASEN:
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in BASEN:
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

# `zertifikat(N, tiefe)` = certificate (`tiefe` = recursion depth). Keys: `art` = kind ('mr12' or 'pocklington'), `F` = the
#   prime
# powers [q, e] used for F, `a` = the base chosen for each q, `unter` = the certificates of the primes q.
def zertifikat(N, tiefe=0):
    """Pocklington certificate as a dict, or None. Small N (< psi_12): kind 'mr12'."""
    if N < PSI12:
        return {'N': str(N), 'art': 'mr12'} if mr12(N) else None
    if not mr12(N): return None
    # `fs` = factorization of N - 1 (prime -> exponent); `teile` = its prime powers, smallest first.
    fs = sympy.factorint(N - 1)
    F, teile = 1, []
    for q in sorted(fs, reverse=False):
        F *= q ** fs[q]; teile.append((q, fs[q]))
    assert F == N - 1
    # Use only as many prime powers as needed (F^2 > N), the smallest first: less recursion. `benutzt` = the prime powers
    #   used.
    F, benutzt = 1, []
    for q, e in teile:
        if F * F > N: break
        F *= q ** e; benutzt.append((q, e))
    if F * F <= N: return None
    zeugen = []
    for q, e in benutzt:
        for a in range(2, 1000):
            if pow(a, N - 1, N) != 1: return None              # Fermat test violated => composite
            if gcd(pow(a, (N - 1) // q, N) - 1, N) == 1:
                zeugen.append(a); break
        else: return None
    # `unter` = certificates of the prime factors q of F.
    unter = []
    for q, e in benutzt:
        z = zertifikat(q, tiefe + 1)
        if z is None: return None
        unter.append(z)
    return {'N': str(N), 'art': 'pocklington', 'F': [[str(q), e] for q, e in benutzt], 'a': zeugen, 'unter': unter}

# ---------------------------------------------- second checker (factors nothing)
def pruefe(z):
    N = int(z['N'])
    if z['art'] == 'mr12':
        return N < PSI12 and mr12(N)
    if z['art'] != 'pocklington': return False
    F = 1
    for q, e in z['F']:
        F *= int(q) ** int(e)
    if (N - 1) % F != 0 or F * F <= N: return False
    if len(z['a']) != len(z['F']) or len(z['unter']) != len(z['F']): return False
    for (q, e), a, u in zip(z['F'], z['a'], z['unter']):
        q = int(q)
        if int(u['N']) != q or not pruefe(u): return False
        if pow(a, N - 1, N) != 1 or gcd(pow(a, (N - 1) // q, N) - 1, N) != 1: return False
    return True

# E1 (PK)
N1 = int(sympy.nextprime(10**40)); z1 = zertifikat(N1)
e1 = z1 is not None and pruefe(z1)
sag(f'PK E1 {"✅" if e1 else "❌"}  naechste Primzahl nach 10^40: Zertifikat erzeugt und vom zweiten Pruefer bestaetigt'); assert e1
# E2 (NK): `falle` = trap (product of two primes); `z_falsch` / `z_falsch2` = tampered certificates
p20, q20 = int(sympy.nextprime(10**19)), int(sympy.nextprime(3 * 10**19))
falle = p20 * q20
e2a = zertifikat(falle) is None
z_falsch = json.loads(json.dumps(z1)); z_falsch['F'][0][1] = int(z_falsch['F'][0][1]) + 1
e2b = not pruefe(z_falsch)
z_falsch2 = json.loads(json.dumps(z1)); z_falsch2['N'] = str(N1 + 2)
e2c = not pruefe(z_falsch2)
e2 = e2a and e2b and e2c
sag(f'NK E2 {"✅" if e2 else "❌"}  Produkt zweier 20-stelliger Primzahlen: kein Zertifikat ({e2a}); verfaelschter Exponent abgelehnt ({e2b}); falsches N abgelehnt ({e2c})'); assert e2

# `inv` = factor inventory of w162; `ziel` = its only-probable numbers of at most 80 digits, smallest first; `ergebnisse` =
#   results.
inv = json.load(open(ERG/'w162_inventar_wahrscheinliche_primzahlen_result.json', encoding='utf-8'))['inventar']
ziel = sorted({x['zahl']: x for x in inv if x['status'] == 'nur_wahrscheinlich' and x['stellen'] <= 80}.values(), key=lambda x: x['stellen'])
ergebnisse = []
for x in ziel:
    N = int(x['zahl']); t0 = time.time()
    z = zertifikat(N); dt = time.time() - t0
    ok = z is not None and pruefe(z)
    ergebnisse.append(dict(m=x['m'], d=x['d'], stellen=x['stellen'], zahl=x['zahl'], bewiesen=ok, sek=round(dt, 1), zertifikat=z))
    sag(f'  ({x["m"]}, {x["d"]}) {x["stellen"]:>3} St.: ' + ('✅ bewiesen prim (Pocklington), zweiter Pruefer bestaetigt' if ok else '🔴 KEIN Beweis') + f'  {dt:.1f} s')
n_ok = sum(1 for r in ergebnisse if r['bewiesen'])
sag(f'E3 {"✅" if n_ok == len(ergebnisse) else "⚠️"}  {n_ok} von {len(ergebnisse)} Zahlen bis 80 Stellen bewiesen')
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), psi12=str(PSI12), n=len(ergebnisse), bewiesen=n_ok,
           ergebnisse=ergebnisse)
(ERG/'w163_pocklington_beweise_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w163_pocklington_beweise_output.txt').write_text(f'w163 · {time.strftime("%Y-%m-%d %H:%M")}\n' + '\n'.join(aus) + '\n', encoding='utf-8')
sag('Ergebnis: w163_pocklington_beweise_result.json')
