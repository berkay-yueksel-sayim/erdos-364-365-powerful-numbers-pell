# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w153_t1plus1_quadrat_bloecke_2026-09-30.py
# Purpose: kernels with T_1 + 1 = square (`quadrat`): on them the building block (`baustein`) M_d, the primitive part of T_d, has
#   at least two distinct prime factors.
# Reads: ergebnisse/w138_karte_961_daten.json, ergebnisse/w78_omega_heuristik_result.json. Requires gmpy2.
# Writes: ergebnisse/w153_t1plus1_quadrat_bloecke_result.json and ergebnisse/w153_t1plus1_quadrat_bloecke_output.txt.
# Usage: python w153_t1plus1_quadrat_bloecke_2026-09-30.py (no arguments; the results folder is found relative to the script).
# Controls: negative-type control (PK-): on kernels WITHOUT T_1 + 1 = square there are prime building blocks (w78: 18); the test
#   must recognize them as "prime" (`prim`), i.e. it discriminates. NK: no prime block of w78 lies on any of the 26 square
#   kernels.

# Background: Remark I.6.11 (a second primitive prime factor is not available, so the theorem does not apply to any block) was
#   shown for ONE Lehmer pair only. If T_1 + 1 = x^2 is a square, then eps = gamma^2 with gamma + 1/gamma = x*sqrt(2), and
#   (gamma, 1/gamma) is a DIFFERENT real Lehmer pair with L' = 2x^2, M' = 1. For this pair kappa = k(M*max{L' - 4M, L'}) =
#   k(2x^2) = 2 and eta = 2 (Def. 1.29 in Juricevic 2007), so eta*kappa = 4. For n = 4d with d odd, n/(eta*kappa) = d is odd,
#   and Thm 1.1 (real pairs: n > 4, n != 6, n/(eta*kappa) odd, triple not in Table 1.2) gives at least TWO primitive divisors
#   of u_{4d}(gamma, 1/gamma). The primitive divisors are exactly the primes of order 4d of eps, i.e. the primes of rank d, so
#   M_d has at least two distinct prime factors there (for every odd d >= 5; (n, L', M') = (20, ., .) occurs in Table 1.2 only
#   for (1, -2) and (9, 2)). Elementary: T_d + 1 = (T_1 + 1)*w_d^2 = (x*w_d)^2 (d odd), hence T_d = (s - 1)(s + 1) with
#   s = x*w_d.
# What: (1) among the 79 kernels of the 961 building blocks (w138), find those with T_1 + 1 = square; (2) for them: T_d + 1 is a
#   square for d = 3, 5, 7, 9; (3) for EVERY block (m, d) on these kernels (at most 3000 digits) form M_d exactly (T_d without
#   all prime factors of the product of T_e over the divisors e < d of d) and check that M_d is composite (Miller-Rabin, gmpy2)
#   and not a perfect power, which gives at least two distinct prime factors (strictly: a composite number with only one prime
#   factor would be a proper prime power); (4) every prime kernel m = 7 (mod 8) up to 5000 has T_1 + 1 = square (odd class
#   number, so x^2 - m*y^2 = 2 is solvable).
# Expectations (stated before the run): 26 kernels; 361 building blocks on them; NO violation.
# No external code is used; Juricevic (2007) is used for its theorem only.
import sys, json, math, time, pathlib
sys.set_int_max_str_digits(0)
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
import gmpy2
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent / 'ergebnisse'   # script folder; results folder
aus = []   # `aus` = output lines, written to the _output.txt file
def sag(s=''):   # `sag` = say: print a line and record it in `aus`
    print(s, flush=True); aus.append(s)

# fundamental solution (T_1, U_1) of x^2 - m*y^2 = 1 (continued fraction of sqrt(m); an odd period length gives a solution of -1,
#   which is squared)
def fund(m):
    a0 = math.isqrt(m); P, Q, a = 0, 1, a0; p1, p = 1, a0; q1, q = 0, 1; l = 0
    while True:
        P = a * Q - P; Q = (m - P * P) // Q; a = (a0 + P) // Q; l += 1
        p1, p = p, a * p + p1; q1, q = q, a * q + q1
        if Q == 1: break
    T, U = p1, q1
    if l % 2 == 1: T, U = T * T + m * U * U, 2 * T * U
    assert T * T - m * U * U == 1
    return T, U

# `tfolge` = T sequence: list [T_0, ..., T_dmax] with T_0 = 1 and T_j = 2 T1 T_{j-1} - T_{j-2}
def tfolge(T1, dmax):
    T = [1, T1]
    for j in range(2, dmax + 1): T.append(2 * T1 * T[-1] - T[-2])
    return T

sag('=' * 100); sag(f'w153 — Kerne mit T_1 + 1 = Quadrat: der Baustein M_d hat mindestens zwei Primteiler · {time.strftime("%Y-%m-%d %H:%M")}'); sag('=' * 100)
w138 = json.load(open(ERG / 'w138_karte_961_daten.json', encoding='utf-8'))
w78 = json.load(open(ERG / 'w78_omega_heuristik_result.json', encoding='utf-8'))
# `bloecke` = building blocks as (m, d, lg); `lg` = log10 of the block
bloecke = [(r['m'], r['d'], r['lg']) for r in w138['raenge']]
kerne = sorted({m for m, d, lg in bloecke}); assert len(kerne) == 79 and len(bloecke) == 961   # `kerne` = the kernels m
T1 = {m: fund(m)[0] for m in kerne}
quad = [m for m in kerne if math.isqrt(T1[m] + 1) ** 2 == T1[m] + 1]   # `quad` = square kernels (T_1 + 1 a square)
sag(f'(1) Kerne mit T_1 + 1 = Quadrat: {len(quad)} von {len(kerne)}: {quad}')
assert len(quad) == 26, len(quad)
# (2) identity
for m in quad:
    T = tfolge(T1[m], 9)
    for d in (3, 5, 7, 9):
        assert math.isqrt(T[d] + 1) ** 2 == T[d] + 1, (m, d)
sag('(2) ✅ T_d + 1 ist Quadrat fuer d = 3, 5, 7, 9 auf allen 26 Kernen (T_d = (s − 1)(s + 1))')
# (3) building blocks
# `teiler_ungerade` = divisors (of an odd number): the divisors e < d of d
def teiler_ungerade(d): return [e for e in range(1, d) if d % e == 0]
# `baustein` = building block M_d: T_d with all common factors with `Lw` (the product of T_e over the divisors e < d) divided out
def baustein(m, d, Tl):
    Lw = 1
    for e in teiler_ungerade(d): Lw *= Tl[e]
    M = Tl[d]
    while True:
        g = math.gcd(M, Lw)
        if g == 1: break
        M //= g
    return M
# `zweiprimteiler` = two-prime-factor test; returns 'eins' (M = 1), 'prim' (prime), 'potenz' (perfect power),
#   or 'mind. zwei' (at least two distinct prime factors)
def zweiprimteiler(M):
    # strict: composite and not a perfect power  =>  at least two distinct prime factors
    if M == 1: return 'eins'
    if gmpy2.is_prime(M, 25): return 'prim'
    if gmpy2.is_power(M): return 'potenz'
    return 'mind. zwei'
auf_quad = [(m, d, lg) for m, d, lg in bloecke if m in quad]   # `auf_quad` = blocks lying on the square kernels
sag(f'(3) Bausteine auf diesen Kernen: {len(auf_quad)} von {len(bloecke)}')
assert len(auf_quad) == 361, len(auf_quad)
# `ergebnis` = counts per verdict of `zweiprimteiler`; `uebersprungen` = blocks skipped (more than 3000 digits);
#   `lgabw` = largest deviation of log10(M_d) from the column `lg` (checked only for M_d < 10^300)
ergebnis = {}; uebersprungen = 0; lgabw = 0.0; t0 = time.time(); cache = {}
for m, d, lg in auf_quad:
    if m not in cache: cache[m] = tfolge(T1[m], max(dd for mm, dd, _ in auf_quad if mm == m))
    M = baustein(m, d, cache[m]); ziffern = len(str(M)) if M.bit_length() < 12000 else None   # `ziffern` = digits of M_d
    lgabw = max(lgabw, abs(math.log10(M) - lg) if M.bit_length() < 1000 else 0.0) if lg else lgabw
    if ziffern is None or ziffern > 3000: uebersprungen += 1; continue
    r = zweiprimteiler(M); ergebnis[r] = ergebnis.get(r, 0) + 1
    assert r == 'mind. zwei', f'VERSTOSS: m = {m}, d = {d}: {r}'
sag(f'    geprueft {sum(ergebnis.values())} Bausteine (Stellenzahl ≤ 3000; uebersprungen {uebersprungen}): {ergebnis} — {time.time() - t0:.0f} s')
sag(f'    Vergleich log10(M_d) mit der Spalte lg aus w138 (nur M_d < 10^300): groesste Abweichung {lgabw:.3f}')
assert sum(ergebnis.values()) + uebersprungen == 361
# PK-: the prime blocks of w78 on kernels without T_1 + 1 = square are recognized as prime
pk_prim = 0   # `pk_prim` = number of prime blocks recognized
for m, d in w78['phi_prim']:
    assert m not in quad, f'NK: prime Baustein auf einem Quadrat-Kern: {m}'
    if m not in T1: T1[m] = fund(m)[0]
    Tl = tfolge(T1[m], d); M = baustein(m, d, Tl)
    assert zweiprimteiler(M) == 'prim', f'PK−: ({m}, {d}) nicht als prim erkannt'
    pk_prim += 1
sag(f'PK− ✅ die {pk_prim} primen Bausteine aus w78 (Kerne {sorted({m for m, d in w78["phi_prim"]})}) werden als „prim" erkannt — der Test unterscheidet')
sag('NK ✅ kein prime Baustein aus w78 liegt auf einem der 26 Kerne')
# (4) prime kernels = 7 (mod 8)
PRIM_GRENZE = 5000   # `PRIM_GRENZE` = prime limit
ps = [p for p in range(7, PRIM_GRENZE) if p % 8 == 7 and gmpy2.is_prime(p)]   # `ps` = primes p = 7 (mod 8) below the limit
nq = [p for p in ps if math.isqrt(fund(p)[0] + 1) ** 2 != fund(p)[0] + 1]   # `nq` = exceptions (T_1 + 1 not a square)
sag(f'(4) prime Kerne m ≡ 7 (mod 8), m < 5000: {len(ps)}; T_1 + 1 ist Quadrat bei {len(ps) - len(nq)}; Ausnahmen: {nq}')
assert not nq
sag('GESAMT ✅')
# result keys: `skript` = script name, `datum` = date, `kerne_quadrat` = the square kernels, `n_kerne` = number of kernels,
#   `n_bloecke_quadrat` = number of blocks on them, `geprueft` = counts per verdict of `zweiprimteiler`, `uebersprungen` =
#   skipped, `pk_prim` = number of prime blocks recognized, `prim_grenze` = prime limit, `prim_kerne` = number of prime kernels,
#   `prim_kerne_ausnahmen` = number of exceptions
(ERG / 'w153_t1plus1_quadrat_bloecke_result.json').write_text(json.dumps(dict(
    skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), kerne_quadrat=quad, n_kerne=len(kerne), n_bloecke_quadrat=len(auf_quad),
    geprueft=ergebnis, uebersprungen=uebersprungen, pk_prim=pk_prim, prim_grenze=5000, prim_kerne=len(ps), prim_kerne_ausnahmen=len(nq)), indent=1, ensure_ascii=False), encoding='utf-8')
(ERG / 'w153_t1plus1_quadrat_bloecke_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
