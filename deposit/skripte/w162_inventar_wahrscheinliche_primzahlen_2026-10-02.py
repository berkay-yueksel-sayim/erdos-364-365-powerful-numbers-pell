# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# -*- coding: utf-8 -*-
# w162_inventar_wahrscheinliche_primzahlen_2026-10-02.py
# Inventory: which numbers in the 41 building blocks of the classes "voll" / "prpzeuge" are ONLY probably prime?
# Purpose: before any proof is attempted, a complete list, independent of which older file stored which factors.
# Method (own code): for each block (m, d) of w138 with class voll or prpzeuge, M_d is built exactly (T_d by recursion,
#   Moebius product over the divisors of d). Then (1) primes below 10^6 are divided out; (2) all known factors from the sources
#   are divided out (w30c, w83, w89, w136, and the w158 parts via gcd); (3) the remainder must be 1 or probably prime (BPSW of
#   sympy AND 12 fixed Miller-Rabin bases), otherwise an assertion fails (a small composite remainder of at most 70 digits is
#   first split with the library function sympy.factorint). Status per factor: "proven" below psi_12 = 318665857834031151167461
#   (deterministic 12-base test, Sorenson and Webster 2017) or with a certificate from w95; otherwise "only probably prime". For
#   the block (7, 1449) of class prpzeuge only the part B (the witness) counts; A is composite.
# Reads (ergebnisse/): w138_karte_961_daten.json, w30c_ueberschneidung_result.json, w83_primtest_result.json,
#   w89_gmp_ecm_result.json, w136_stufe2_zweite_route_result.json, w95_primzertifikate_result.json.
# Writes: ergebnisse/w162_inventar_wahrscheinliche_primzahlen_result.json and the matching _output.txt. No command-line
#   arguments.
# Controls (expectations fixed in advance):
#   E1       41 blocks (w138: voll 40 + prpzeuge 1); per block the product of the found factors (with exponents) equals M_d
#            (voll), resp. B | M_d with exponent 1 (prpzeuge).
#   E2 (PK)  For (255, 75) the method yields the two large factors known from w83 (26 and 35 digits). (255, 75) is one of the
#            41 blocks; a control on an object outside the 41 would check nothing.
#   E3       No remainder stays composite (otherwise "completely factored" would be wrong).
#   E4       Output: list of the only probably prime numbers with digit count and block; count per size class.
import sys, json, time, math, pathlib
from math import isqrt, gcd
import sympy
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'
# `J` = load a JSON file from the results folder; `aus` / `sag` (= say) collect and print the output lines.
def J(n): return json.load(open(ERG/n, encoding='utf-8'))
aus = []
def sag(s=''):
    print(s, flush=True); aus.append(s)
# `PSI12` = bound psi_12 below which the Miller-Rabin test with the first 12 prime bases (`BASEN`) is deterministic.
PSI12 = 318665857834031151167461
BASEN = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
# `mr(n, a)` = Miller-Rabin round for base a; `wahrscheinlich_prim` = probable-prime test (trial division by the bases,
# then all 12 Miller-Rabin rounds and sympy.isprime).
def mr(n, a):
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    x = pow(a, d, n)
    if x in (1, n - 1): return True
    for _ in range(s - 1):
        x = x * x % n
        if x == n - 1: return True
    return False
def wahrscheinlich_prim(n):
    if n < 2: return False
    for p in BASEN:
        if n % p == 0: return n == p
    return all(mr(n, a) for a in BASEN) and sympy.isprime(n)
# `fund` = fundamental solution (T_1, U_1) of x^2 - m*y^2 = 1 by continued fractions; `mu` = Moebius function.
def fund(m):
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
def mu(n):
    r, x, f = 1, n, 2
    while f*f <= x:
        if x % f == 0:
            x //= f
            if x % f == 0: return 0
            r = -r
        f += 1
    return -r if x > 1 else r
# `baustein(m, d)` = building block: returns (M_d, T_d), where M_d is the primitive part of T_d, computed as the Moebius
#   product of
# T_e over the divisors e of d, and T_n follows the recursion T_(n+1) = 2*T_1*T_n - T_(n-1).
def baustein(m, d):
    T1, _ = fund(m); t = [1, T1]
    for _ in range(d - 1): t.append(2*T1*t[-1] - t[-2])
    num = den = 1
    for e in range(1, d + 1):
        if d % e: continue
        s = mu(d // e)
        if s == 1: num *= t[e]
        elif s == -1: den *= t[e]
    assert num % den == 0
    return num // den, t[d]
# `sieb` = sieve of Eratosthenes; `KLEIN` = primes up to 10^6.
def sieb(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
KLEIN = sieb(10**6)

# Known factors from the sources, per (m, d): `bekannt` collects them, `dazu` adds one.
bekannt = {}
def dazu(m, d, f):
    f = int(f)
    if f > 1: bekannt.setdefault((int(m), int(d)), set()).add(f)
for u in J('w30c_ueberschneidung_result.json')['unsichtbare']: dazu(u['m'], u['d'], u['kleinster_primfaktor'])
w83 = J('w83_primtest_result.json')
for x in w83['voll']:
    for f, e in x[2]: dazu(x[0], x[1], f)
for x in w83['rest_teilfaktorisiert']:
    for f in x[2]: dazu(x[0], x[1], f)
for k, v in J('w89_gmp_ecm_result.json')['je_zahl'].items():
    for f in v['funde']:
        dazu(v['m'], v['d'], f['faktor'])
for b in J('w136_stufe2_zweite_route_result.json')['probable']: dazu(b[0], b[1], b[2])
w95 = J('w95_primzertifikate_result.json')
zertifiziert = set()
for k, v in w95.get('zertifikate', {}).items() if isinstance(w95.get('zertifikate'), dict) else []:
    pass
def sammle_zahlen(o):
    if isinstance(o, dict):
        for x in o.values(): yield from sammle_zahlen(x)
    elif isinstance(o, list):
        for x in o: yield from sammle_zahlen(x)
    elif isinstance(o, str) and o.isdigit() and len(o) >= 20: yield int(o)
# all numbers with at least 20 digits occurring in w95 (witnesses with certificate AND the three cofactors)
w95_zahlen = set(sammle_zahlen(w95))
w95_prob = {int(x) for x in sammle_zahlen(w95.get('probable', []))}

# `ziel` = the target blocks (m, d, class, route) of the classes voll and prpzeuge; `inventar` = list of all factors found.
w138 = J('w138_karte_961_daten.json')
ziel = [(int(r['m']), int(r['d']), r['klasse'], r['route']) for r in w138['raenge'] if r['klasse'] in ('voll', 'prpzeuge')]
e1a = len(ziel) == 41 and sum(1 for z in ziel if z[2] == 'prpzeuge') == 1
sag(f'PK E1a {"✅" if e1a else "❌"}  {len(ziel)} Bausteine (voll {sum(1 for z in ziel if z[2] == "voll")} + prpzeuge {sum(1 for z in ziel if z[2] == "prpzeuge")})'); assert e1a

inventar, t0 = [], time.time()
for m, d, kl, route in sorted(ziel, key=lambda z: (z[1], z[0])):
    M, Td = baustein(m, d)
    if kl == 'prpzeuge':
        s = isqrt(Td + 1); assert s*s == Td + 1
        B = max(gcd(M, s - 1), gcd(M, s + 1), key=lambda x: wahrscheinlich_prim(x))
        assert wahrscheinlich_prim(B) and M % B == 0 and (M // B) % B != 0, ('prpzeuge', m, d)
        faktoren = [(B, 1)]
    else:
        R, faktoren = M, []
        for p in KLEIN:
            if R % p == 0:
                e = 0
                while R % p == 0: R //= p; e += 1
                faktoren.append((p, e))
            if p * p > R: break
        for f in sorted(bekannt.get((m, d), ())):
            if R % f == 0:
                e = 0
                while R % f == 0: R //= f; e += 1
                faktoren.append((f, e))
        # Parts from w158 for (23, 161) and other kernels with T_1 + 1 = square.
        s = isqrt(Td + 1)
        if s*s == Td + 1 and R > 1:
            for g in (gcd(R, s - 1), gcd(R, s + 1)):
                if 1 < g < R and wahrscheinlich_prim(g):
                    e = 0
                    while R % g == 0: R //= g; e += 1
                    faktoren.append((g, e))
        if R > 1 and not wahrscheinlich_prim(R):
            # Remainder composite but small: the older runs found these factors, but the sources do not store them. Split with
            #   sympy.factorint (library call); above 70 digits "completely factored" would not be established, so assert.
            assert len(str(R)) <= 70, ('E3 verletzt: Rest zusammengesetzt und gross', m, d, len(str(R)))
            for f, e in sympy.factorint(R).items():
                faktoren.append((int(f), e)); R //= int(f)**e
            assert R == 1
        if R > 1:
            faktoren.append((R, 1))
        prod = 1
        for f, e in faktoren: prod *= f**e
        assert prod == M, ('E1 verletzt: Produkt ≠ M_d', m, d)
    for f, e in faktoren:
        assert wahrscheinlich_prim(f), ('Faktor nicht wahrscheinlich prim', m, d, len(str(f)))
        # Status per factor: `bewiesen_psi12` (proven by the deterministic test), `zertifikat_w95` (certificate),
        # `nur_wahrscheinlich` (only probably prime).
        status = 'bewiesen_psi12' if f < PSI12 else ('zertifikat_w95' if (f in w95_zahlen and f not in w95_prob) else 'nur_wahrscheinlich')
        inventar.append(dict(m=m, d=d, klasse=kl, route=route, stellen=len(str(f)), exponent=e, status=status, zahl=str(f)))
# control on the right object: (255, 75) is among the 41, its factors are in w83
pk = [x for x in inventar if (x['m'], x['d']) == (255, 75)]
e2 = sorted(x['stellen'] for x in pk if x['stellen'] > 6) == [26, 35]
sag(f'PK E2 {"✅" if e2 else "❌"}  (255, 75): grosse Faktoren mit {sorted(x["stellen"] for x in pk if x["stellen"] > 6)} Stellen (w83: 26 und 35)'); assert e2
sag(f'E3 ✅  kein zusammengesetzter Rest; alle 41 Produkte = M_d (bzw. B | M_d genau einmal) ({time.time() - t0:.0f} s)')
# `nur` = the only probably prime factors; `klassen` = their count per size class (number of digits).
nur = [x for x in inventar if x['status'] == 'nur_wahrscheinlich']
from collections import Counter
klassen = Counter(('≤ 40' if x['stellen'] <= 40 else '41–80' if x['stellen'] <= 80 else '81–120' if x['stellen'] <= 120 else '121–250' if x['stellen'] <= 250 else '> 250') for x in nur)
sag(f'E4     Faktoren insgesamt {len(inventar)}; bewiesen unter psi_12: {sum(1 for x in inventar if x["status"] == "bewiesen_psi12")}; '
    f'mit w95-Zertifikat: {sum(1 for x in inventar if x["status"] == "zertifikat_w95")}; NUR WAHRSCHEINLICH PRIM: {len(nur)} in {len({(x["m"], x["d"]) for x in nur})} Bausteinen')
sag(f'       nach Stellen: {dict(sorted(klassen.items()))}')
for x in sorted(nur, key=lambda x: x['stellen']):
    sag(f'       ({x["m"]}, {x["d"]}) {x["klasse"]}/{x["route"]}: {x["stellen"]} Stellen')
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), bausteine=len(ziel), faktoren=len(inventar),
           nur_wahrscheinlich=len(nur), bausteine_mit_unbewiesenem=len({(x['m'], x['d']) for x in nur}), inventar=inventar)
(ERG/'w162_inventar_wahrscheinliche_primzahlen_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w162_inventar_wahrscheinliche_primzahlen_output.txt').write_text(f'w162 · {time.strftime("%Y-%m-%d %H:%M")}\n' + '\n'.join(aus) + '\n', encoding='utf-8')
sag('Ergebnis: w162_inventar_wahrscheinliche_primzahlen_result.json')
