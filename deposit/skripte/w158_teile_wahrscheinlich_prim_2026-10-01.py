# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w158_teile_wahrscheinlich_prim_2026-10-01.py
# -*- coding: utf-8 -*-
# Closes two open blocks: the algebraic parts of the blocks M_d are (probably) prime.
# Background: w155 split the blocks on kernels with T_1 + 1 = x^2 into M_d = A * B with A = gcd(M_d, s - 1), B = gcd(M_d, s + 1),
#   gcd(A, B) = 1, but did not test whether A or B is itself prime. If a part is prime, it divides M_d exactly once (the parts are
#   coprime), and M_d is not powerful. ECM can never find this (a prime has no small factor).
# What it does, for EACH of the ten open blocks of the stored file w147 (second route, checked against the stored data rather than
#   against w155):
#   (1) unit from our own continued fraction, T_d exact; if T_1 + 1 is a square: s = sqrt(T_d + 1), A, B; A * B == stored number;
#       gcd(A, B) = 1;
#   (2) rank: A | T_d and gcd(A, T_e) = 1 for all divisors e < d of d (likewise B), so every prime divisor has rank exactly d;
#   (3) primality by two routes: SymPy isprime (BPSW) and gmpy2 is_prime with 50 rounds; exponent in the block: v_P(M_d) = 1 for
#       every probable-prime part P;
#   (4) without an algebraic decomposition: M_d itself by both routes (expectation: composite, otherwise the chain would have
#       treated it as prime).
#   Classes: "both parts probably prime" = complete factorization into probable primes; "one part probably prime" = a witness that
#   is a probable prime, with exponent 1; otherwise open. NO primality proof (not feasible with w95 for 111 to 477 digits):
#   the blocks close under the same condition as the blocks of the class "complete factorization into probable primes",
#   namely if the probable primes are prime.
# Reads: ergebnisse/w147_offene_bloecke_dezimal.txt. Writes: ergebnisse/w158_teile_wahrscheinlich_prim_result.json and
#   ergebnisse/w158_teile_wahrscheinlich_prim_output.txt. Requires SymPy and gmpy2. No command-line arguments.
# Expectations: (23, 161): both parts probably prime (111 and 112 digits); (7, 1449): B probably prime (477 digits), A composite;
#   the other eight: no decomposition, M_d composite.
# Controls: PK: (7, 5): M = 64 261 = 179 * 359, both parts prime (provable here by trial division), hence class "both";
#   PK primality test: 2^127 - 1 prime, (2^127 - 1)(2^61 - 1) not. NK: m = 15 (T_1 + 1 not a square) is not decomposed.
# Own code.
import sys, math, json, time, pathlib, re
from sympy import isprime
import gmpy2
sys.set_int_max_str_digits(0)
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent / 'ergebnisse'   # `ERG` = `ergebnisse` (results)
aus = []   # `aus` = output lines
# `sag` = say: print a line and record it
def sag(s=''): print(s, flush=True); aus.append(s)

# `fund(m)` = fundamental solution (T_1, U_1) of x^2 - m*y^2 = 1 by continued fraction; `T(m, k)` = T_k exactly
def fund(m):
    a0 = math.isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0 * h0 - m * k0 * k0 != 1:
        P = a * Q - P; Q = (m - P * P) // Q; a = (a0 + P) // Q
        h1, h0 = h0, a * h0 + h1; k1, k0 = k0, a * k0 + k1
    return h0, k0
def T(m, k):
    T1, U1 = fund(m); a, b = 1, 0; x, y = T1, U1
    while k:
        if k & 1: a, b = a * x + m * b * y, a * y + b * x
        x, y = x * x + m * y * y, 2 * x * y; k >>= 1
    return a
# `prp(n)` = probable-prime test by two routes: (SymPy isprime, gmpy2 is_prime with 50 rounds)
def prp(n): return isprime(n), bool(gmpy2.is_prime(n, 50))
# `rang_genau` = rank exactly d: P | T_d and gcd(P, T_e) = 1 for every divisor e < d of d
def rang_genau(P, m, d, Td):
    return Td % P == 0 and all(math.gcd(P, T(m, e)) == 1 for e in range(1, d) if d % e == 0)

# `pruefe(m, d, M)` = check one block M = M_d of kernel m: decomposition A * B (if T_1 + 1 is a square), primality of the
# parts, rank, exponent. `klasse`: `beide_prp` = both parts probably prime, `zeuge_prp` = one part probably prime, `offen` =
# open, `M_prim` = no decomposition but M itself probably prime. Result fields: `zerlegt` = decomposed (T_1 + 1 is a square),
# `stellen` / `stellen_A` / `stellen_B` = digits of M / A / B, `A_prp` / `B_prp` = parts probably prime, `M_prp` = the two
# tests on M, `rang_genau_A` / `rang_genau_B` = rank exactly d, `exponent_eins` = every probable-prime part divides M exactly
# once.
def pruefe(m, d, M):
    T1, _ = fund(m); x = math.isqrt(T1 + 1); Td = T(m, d)
    assert Td % M == 0, f'({m}, {d}): hinterlegter Block teilt T_d nicht'
    if x * x != T1 + 1:
        r = prp(M); return dict(m=m, d=d, stellen=len(str(M)), zerlegt=False, M_prp=r, klasse='offen' if not any(r) else 'M_prim')
    s = math.isqrt(Td + 1); assert s * s == Td + 1
    A, B = math.gcd(M, s - 1), math.gcd(M, s + 1)
    assert A * B == M and math.gcd(A, B) == 1, f'({m}, {d}): A·B ≠ M oder nicht teilerfremd'
    rA, rB = prp(A), prp(B)
    assert rA[0] == rA[1] and rB[0] == rB[1], f'({m}, {d}): die zwei Primtest-Routen widersprechen sich {rA} {rB}'
    rangA, rangB = rang_genau(A, m, d, Td), rang_genau(B, m, d, Td)
    # `expo` = for each probable-prime part: it divides M exactly once
    expo = {n: (M % n == 0 and (M // n) % n != 0) for n, r in (('A', rA), ('B', rB)) if r[0] for n in [A if n == 'A' else B]}
    klasse = 'beide_prp' if rA[0] and rB[0] else ('zeuge_prp' if rA[0] or rB[0] else 'offen')
    return dict(m=m, d=d, stellen=len(str(M)), zerlegt=True, stellen_A=len(str(A)), stellen_B=len(str(B)), A_prp=rA[0], B_prp=rB[0],
                rang_genau_A=rangA, rang_genau_B=rangB, exponent_eins=all(expo.values()), klasse=klasse)

# controls: `P_PK` = 2^127 - 1 (prime), `Z_PK` = product of two Mersenne primes (composite)
sag('=' * 100); sag(f'w158 — algebraische Teile der offenen Bloecke: wahrscheinlich prim? · {time.strftime("%Y-%m-%d %H:%M")}'); sag('=' * 100)
P_PK, Z_PK = 2**127 - 1, (2**127 - 1) * (2**61 - 1)
assert prp(P_PK) == (True, True) and prp(Z_PK) == (False, False); sag('PK ✅ Primtest: 2^127 − 1 prim, (2^127 − 1)(2^61 − 1) nicht (beide Routen)')
M = T(7, 5)
# rank-5 part: divide out everything shared with T_1 (1 is the only proper divisor of 5)
while math.gcd(M, T(7, 1)) > 1: M //= math.gcd(M, T(7, 1))
pk = pruefe(7, 5, M); assert M == 64261 == 179 * 359 and pk['klasse'] == 'beide_prp', pk
assert all(179 % q and 359 % q for q in range(2, 19)), 'Probedivision'
sag(f'PK ✅ (7, 5): M = 64 261 = 179 · 359, Klasse „{pk["klasse"]}" (179 und 359 durch Probedivision prim)')
nk = pruefe(15, 7, T(15, 7))   # M = T_7 suffices: this only checks that nothing is decomposed when T_1 + 1 is not a square
assert nk['zerlegt'] is False; sag('NK ✅ m = 15: T_1 + 1 kein Quadrat, keine Zerlegung')
# the ten open blocks from the stored file: a line 'm = .. d = .. digits = ..' followed by the decimal number
# (`bloecke` = blocks)
txt = (ERG / 'w147_offene_bloecke_dezimal.txt').read_text(encoding='utf-8').split('\n')
bloecke = []
for i, z in enumerate(txt):
    mm = re.match(r'^m = (\d+)\s+d = (\d+)\s+digits = (\d+)', z)
    if mm:
        n = int(txt[i + 1].strip()); assert len(str(n)) == int(mm.group(3)); bloecke.append((int(mm.group(1)), int(mm.group(2)), n))
assert len(bloecke) == 10, len(bloecke)
sag(f'\nOffene Bloecke aus w147 (hinterlegt): {len(bloecke)}')
zeilen = []   # `zeilen` = result rows, one per block
for m, d, M in bloecke:
    z = pruefe(m, d, M); zeilen.append(z)
    if z['zerlegt']:
        sag(f'  ({m}, {d}) {z["stellen"]} St. = A ({z["stellen_A"]} St., {"wahrsch. prim" if z["A_prp"] else "zusammengesetzt"}) · '
            f'B ({z["stellen_B"]} St., {"wahrsch. prim" if z["B_prp"] else "zusammengesetzt"}); Rang genau d: A {z["rang_genau_A"]}, B {z["rang_genau_B"]}; '
            f'Exponent 1: {z["exponent_eins"]} ⇒ {z["klasse"]}')
    else:
        sag(f'  ({m}, {d}) {z["stellen"]} St.: keine algebraische Zerlegung; M_d wahrsch. prim: {z["M_prp"]} ⇒ {z["klasse"]}')
# `klassen` = the blocks (m, d) per class
klassen = {k: [(z['m'], z['d']) for z in zeilen if z['klasse'] == k] for k in ('beide_prp', 'zeuge_prp', 'offen', 'M_prim')}
assert klassen['beide_prp'] == [(23, 161)] and klassen['zeuge_prp'] == [(7, 1449)] and len(klassen['offen']) == 8 and not klassen['M_prim'], klassen
assert all(z['rang_genau_A'] and z['rang_genau_B'] and z['exponent_eins'] for z in zeilen if z['zerlegt'])
sag(f'\nGESAMT ✅ {len(klassen["beide_prp"])} Block vollstaendig in wahrscheinliche Primzahlen zerlegt {klassen["beide_prp"]}, '
    f'{len(klassen["zeuge_prp"])} Block mit wahrscheinlich primem Zeugen {klassen["zeuge_prp"]}; offen bleiben {len(klassen["offen"])}. Kein Primbeweis.')
# JSON keys: `skript` = script, `datum` = date, `quelle` = source, `bloecke` = blocks, `klassen` = classes,
# `offen_danach` = number of blocks still open afterwards
(ERG / 'w158_teile_wahrscheinlich_prim_result.json').write_text(json.dumps(dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'),
    quelle='w147_offene_bloecke_dezimal.txt', bloecke=zeilen, klassen={k: [list(x) for x in v] for k, v in klassen.items()},
    offen_danach=len(klassen['offen'])), indent=1, ensure_ascii=False), encoding='utf-8')
(ERG / 'w158_teile_wahrscheinlich_prim_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
sag('Ergebnis: w158_teile_wahrscheinlich_prim_result.json')
