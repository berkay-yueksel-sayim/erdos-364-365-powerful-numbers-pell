# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w136_stufe2_zweite_route_zertifikate_2026-09-29.py
# Purpose: stage 2, second route and primality proof for ALL new witnesses found by w133, w134 and w135. Every witness p of a pair
#   (m, d) is re-verified by two independent routes (A: powers of eps; L: Lucas sequence) and proved prime by a Pocklington
#   certificate; a second, separately written checker re-verifies the certificates from the file.
# Reads (ergebnisse/): w9_v31_M100000000_H2000_result.json, w133_probelauf_79_kerne_result.json, w134_stufe2_rangsieb_result.json,
#   w135_stufe2_ecm_checkpoint.json, and (unless in probe mode) w135_stufe2_ecm_result.json. Requires SymPy.
# Writes (ergebnisse/): w136_zertifikate.json, w136_stufe2_zweite_route_result.json, w136_stufe2_zweite_route_output.txt.
# Usage: python w136_stufe2_zweite_route_zertifikate_2026-09-29.py [probe]   (`probe`: read the w135 checkpoint, write nothing)
# Controls: E1 to E5 below (positive and negative controls of both routes, of the prime proof and of its checker, and the
#   partition of the 150 new pairs); the run aborts if a control fails or if the two routes disagree.

# Note: without `probe`, the script refuses to run while the lock file ECM_LAEUFT.lock exists two folders above the script folder
#   (a w135 run is still in progress).
# Background: the 150 new pairs (m, d) of the 54 kernels outside the 25 got their witnesses by three methods, one per script:
#   w133 (sieve <= 3·10⁶)  — Lucas sequence V_n(2T₁, 1) by doubling             → independent of it is ROUTE A
#   w134 (sieve <= 10⁹)    — binary exponentiation of eps (eps_pow)              → independent of it is ROUTE L
#   w135 (GMP-ECM <= t35)  — linear recurrence `T_rek` for the rank               → independent of it are A AND L
#   Here EVERY witness gets BOTH routes; so each one has at least one route that is not the one it was found by.
# ROUTE A (as in w90): eps^(4d) = 1 mod p, eps^(4d/q) ≠ 1 for every prime q | 4d (order exactly 4d), p | T_d, v_p(T_d) = 1 via
#   eps^d mod p³, p ≡ ±1 (mod 4d), negative control p ∤ T_{d+2}.
# ROUTE L: V_d ≡ 0 (mod p), V_{d/r} ≢ 0 for every prime r | d (rank exactly d), V_d ≢ 0 (mod p²), p ≡ ±1 (mod 4d). (V_n = 2 T_n.)
# PRIME PROOF per witness: below ψ₁₂ = 318665857834031151167461 a deterministic Miller–Rabin test on the first twelve prime
#   bases (Sorenson–Webster 2017); above it POCKLINGTON (1914) with a PARTIALLY factored N − 1 = F·R: F completely factored,
#   F² > N, for every prime q | F a base a with a^(N−1) ≡ 1 and gcd(a^((N−1)/q) − 1, N) = 1 ⇒ every prime divisor of N is
#   ≡ 1 (mod F), hence > √N, hence N is prime. (w95 used the special case F = N − 1; Brillhart–Lehmer–Selfridge 1975,
#   Math. Comp. 29(130), 620–647.) The primes q | F are certified recursively by the same procedure.
#   The FACTORIZATION of N − 1 comes from sympy.factorint (library used as a black box; with a time limit, in a separate
#   process). It is not critical: the SECOND, separately written checker reads the certificate from the file and only
#   recomputes (product, F² > N, powers, gcd), so a wrong factorization fails there.
# COUNTING RULE (fixed before the real run; the same TWO classes as Prop. 6.10 (a) for the 811 pairs — 766 proved by certified
#   witnesses, 38 by complete factorization into probable primes):
#   (a) CERTIFIED — at least ONE candidate passes both routes AND is proved prime.
#   (b) COMPLETELY FACTORED, PROBABLE — no proof succeeds, but the block is COMPLETELY factored into probable primes (route `prim`
#       or `voll` in w135) and a candidate passes both routes; three tests (12 bases, 12 other bases, BPSW) are recorded.
#   Otherwise (only one probable factor, the rest not factored): OPEN.
#
# Expectations (stated before the run):
#   E1  PK routes: (7, 5, 179) passes A and L.
#   E2 NK routes: (7, 7, 179) [rank 5, not 7], (7, 5, 29) [rank 7], (15, 45, 181) and (143, 13, 311) [Wieferich, v = 2, from w77]
#       fail in A AND L.
#   E3  PK proof: the next prime after 10³⁰ and after 10⁴⁰ proved, confirmed by the second checker. NK: a 30-digit semiprime and
#       the Carmichael number 5394826801 rejected. Second checker: a partial certificate built by hand (R ≠ 1) is accepted; the
#       same with too small an F (F² < N) and with a wrong base is rejected.
#   E4  Partition: w133 (new, with witness) + w134 settled + w135 settled + w135 open = the 150 new pairs, disjoint.
#   E5 Every witness reported by the finding route passes A and L (otherwise: do NOT book the finding). A and L agree on EVERY
#       candidate.
#   Measured: number of pairs with proof / probable; size of the witnesses.
import json, sys, time, pathlib, multiprocessing as mp
from math import isqrt, gcd
from collections import Counter, defaultdict

PSI12 = 318665857834031151167461   # `PSI12` = ψ₁₂: Miller–Rabin with the first 12 primes is deterministic below this bound
BASEN12 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]        # `BASEN12` = the first twelve prime bases
BASEN_B = (41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89)    # `BASEN_B` = twelve other bases (second probable-prime test)

# fundamental solution (T1, U1) of x^2 - m*y^2 = 1 from the continued fraction of sqrt(m)
def fund(m):
    a0 = isqrt(m); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (m - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
# `primfaktoren_klein` = set of the prime factors of a small n by trial division
def primfaktoren_klein(n):
    f = set(); x = n; q = 2
    while q*q <= x:
        while x % q == 0: f.add(q); x //= q
        q += 1 if q == 2 else 2
    if x > 1: f.add(x)
    return f

# ---- ROUTE A: powers of eps
# (T1 + U1*sqrt(m))^e modulo M by binary exponentiation; returns (rational part, sqrt(m) part)
def eps_pow(T1, U1, m, e, M):
    ra, rb, ba, bb = 1 % M, 0, T1 % M, U1 % M
    while e:
        if e & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        e >>= 1
    return ra, rb
# route A for the witness p of the pair (m, d); `eins` = the identity (1, 0)
def route_A(T1, U1, m, d, p):
    M = 4*d; eins = (1 % p, 0)
    if eps_pow(T1, U1, m, M, p) != eins: return False
    if any(eps_pow(T1, U1, m, M//q, p) == eins for q in primfaktoren_klein(M)): return False
    if eps_pow(T1, U1, m, d, p)[0] % p: return False
    Td3 = eps_pow(T1, U1, m, d, p**3)[0]
    if Td3 % (p*p) == 0: return False                                   # v_p >= 2
    if p % M not in (1, M - 1): return False
    return eps_pow(T1, U1, m, d + 2, p)[0] % p != 0                    # negative control: p must not divide T_{d+2}

# ---- ROUTE L: Lucas V_n(P, 1), P = 2 T₁, doubling from the left (over the binary digits of n)
def V(n, P, mod):
    a, b = 2 % mod, P % mod
    for bit in bin(n)[2:]:
        if bit == '1': a, b = (a*b - P) % mod, (b*b - 2) % mod
        else: a, b = (a*a - 2) % mod, (a*b - P) % mod
    return a
def route_L(T1, d, p):
    if V(d, 2*(T1 % p), p) != 0: return False
    if any(V(d//r, 2*(T1 % p), p) == 0 for r in primfaktoren_klein(d)): return False
    p2 = p*p
    if V(d, 2*(T1 % p2), p2) == 0: return False
    return p % (4*d) in (1, 4*d - 1)

# ---- prime proof
# Miller–Rabin with the twelve bases `BASEN12` (deterministic below PSI12)
def mr12(n):
    assert n < PSI12
    if n < 2: return False
    for q in BASEN12:
        if n % q == 0: return n == q
    dd, s = n - 1, 0
    while dd % 2 == 0: dd //= 2; s += 1
    for a in BASEN12:
        x = pow(a, dd, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x*x % n
            if x == n - 1: break
        else: return False
    return True
# factorization in a subprocess: puts {prime (as a string): exponent} into the queue `q`; `limit` = trial-division limit or None
def _zerlege(n, limit, q):
    import sympy
    q.put({str(k): v for k, v in (sympy.factorint(n, limit=limit) if limit else sympy.factorint(n)).items()})
# `zerlege_mit_zeit` = factor n in a separate process with a time limit of `sek` seconds; returns {prime: exponent} or None on
#   timeout
def zerlege_mit_zeit(n, limit, sek):
    q = mp.Queue(); pr = mp.Process(target=_zerlege, args=(n, limit, q)); pr.start()
    try: erg = q.get(timeout=sek)
    except Exception: erg = None
    pr.join(1)
    if pr.is_alive(): pr.terminate(); pr.join()
    return None if erg is None else {int(k): v for k, v in erg.items()}
# `zertifikat` = certificate for the primality of N, or None: for N < PSI12 a Miller–Rabin certificate ('mr12'); otherwise a
#   Pocklington certificate ('pocklington'): F = product of the smallest prime powers of N − 1 with F² > N, R = (N − 1)/F, a base
#   per prime q | F, and (recursively) the certificates of the primes q | F
def zertifikat(N, sek=120):
    import sympy
    if N < PSI12:
        return {'N': str(N), 'art': 'mr12'} if mr12(N) else None
    for limit in (10**6, None):
        fk = zerlege_mit_zeit(N - 1, limit, sek)
        if fk is None: continue
        prim = sorted(q for q in fk if sympy.isprime(q))
        F, gewaehlt = 1, []   # `gewaehlt` = chosen: the primes q taken into F
        for q in prim:                                                   # smallest first, until F² > N
            if F*F > N: break
            F *= q**fk[q]; gewaehlt.append(q)
        if F*F <= N: continue
        R = (N - 1) // F
        assert F*R == N - 1
        basen, unter = {}, []   # `basen` = base per prime q | F; `unter` = sub-certificates of these primes
        for q in gewaehlt:
            for a in range(2, 200):
                if pow(a, N - 1, N) != 1: return None                    # Fermat test fails ⇒ N is composite
                if gcd(pow(a, (N - 1)//q, N) - 1, N) == 1: basen[str(q)] = a; break
            else: break
            z = zertifikat(q, sek)
            if z is None: break
            unter.append(z)
        else:
            return {'N': str(N), 'art': 'pocklington', 'F': [[str(q), fk[q]] for q in gewaehlt], 'R': str(R),
                    'basen': basen, 'faktoren': unter}
    return None

# ---- second checker, written separately: factors nothing, only recomputes
# Miller–Rabin test of n for the single base a
def _mr_basis(n, a):
    if n % a == 0: return n == a
    t = n - 1; r = 0
    while not t & 1: t >>= 1; r += 1
    y = pow(a, t, n)
    if y == 1 or y == n - 1: return True
    for _ in range(r - 1):
        y = pow(y, 2, n)
        if y == n - 1: return True
    return False
# `pruefe` = check a certificate z (recursively): True if it proves N prime
def pruefe(z):
    N = int(z['N'])
    if z['art'] == 'mr12':
        if N < 10**12:
            if N < 2: return False
            i = 2
            while i*i <= N:
                if N % i == 0: return False
                i += 1
            return True
        return N < PSI12 and all(_mr_basis(N, a) for a in BASEN12)
    if z['art'] != 'pocklington': return False
    F = 1
    for q, e in z['F']: F *= int(q)**e
    if F*int(z['R']) != N - 1 or F*F <= N: return False
    kinder = {k['N']: k for k in z['faktoren']}   # `kinder` = sub-certificates by N
    for q, _ in z['F']:
        a = z['basen'].get(q)
        if a is None or q not in kinder: return False
        if pow(a, N - 1, N) != 1 or gcd(pow(a, (N - 1)//int(q), N) - 1, N) != 1: return False
        if not pruefe(kinder[q]): return False
    return True

def main():
    import sympy
    sys.stdout.reconfigure(encoding='utf-8')
    # script folder; results folder; folder above it
    HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'; BOX = HIER.parent.parent
    PROBE = len(sys.argv) > 1 and sys.argv[1] == 'probe'
    aus = []   # `aus` = output lines, written to the _output.txt file
    def sag(s=''):   # `sag` = say: print a line and record it in `aus`
        print(s, flush=True); aus.append(s)
    t0 = time.time()
    sag('='*100); sag(f'w136 — ZWEITE ROUTE (A + L) UND PRIMBEWEIS FUER ALLE NEUEN ZEUGEN   {time.strftime("%Y-%m-%d %H:%M")}'
                      + ('   [PROBE: w135-Checkpoint, schreibt nichts]' if PROBE else '')); sag('='*100)
    # E1 / E2
    fc = {}   # cache of fundamental solutions
    def F_(m):
        if m not in fc: fc[m] = fund(m)
        return fc[m]
    def beide(m, d, p):   # `beide` = both: the results of route A and route L for the witness p of (m, d)
        T1, U1 = F_(m); return route_A(T1, U1, m, d, p), route_L(T1, d, p)
    e1 = beide(7, 5, 179) == (True, True)
    nk = {(7, 7, 179): beide(7, 7, 179), (7, 5, 29): beide(7, 5, 29), (15, 45, 181): beide(15, 45, 181), (143, 13, 311): beide(143, 13, 311)}
    e2 = all(v == (False, False) for v in nk.values())
    sag(f'PK E1 {"✅" if e1 else "❌"}  (7, 5, 179) besteht A und L.')
    sag(f'NK E2 {"✅" if e2 else "❌"}  abgelehnt in A und L: ' + ', '.join(f'{k} → {v}' for k, v in nk.items()))
    assert e1 and e2, 'E1/E2 VERLETZT'
    # E3
    N30, N40 = sympy.nextprime(10**30), sympy.nextprime(10**40)
    z30, z40 = zertifikat(N30), zertifikat(N40)
    zus = sympy.nextprime(10**14) * sympy.nextprime(3 * 10**15)   # `zus` = composite (semiprime) for the negative control
    carm = 7 * 13 * 17 * 23 * 31 * 67 * 73                         # `carm` = Carmichael number
    assert carm == 5394826801 and all((carm - 1) % (q - 1) == 0 for q in (7, 13, 17, 23, 31, 67, 73)) and pow(2, carm - 1, carm) == 1
    # partial certificate by hand: N = 2·A·C + 1 with A = product of small primes, C = q1·q2 (composite, stays as R)
    A = 1; q_ = 3
    while (2*A)**2 < 10**90: A *= q_; q_ = sympy.nextprime(q_)
    q1 = sympy.nextprime(10**14); q2 = sympy.nextprime(2*10**14)
    while not sympy.isprime(2*A*q1*q2 + 1): q2 = sympy.nextprime(q2)
    Nt = 2*A*q1*q2 + 1; Ft = sympy.factorint(2*A)
    def hand(N, fk):   # `hand` = build a Pocklington certificate by hand from the given factorization `fk` of the factored part
        basen, unter = {}, []
        for q in fk:
            a = next(a for a in range(2, 200) if gcd(pow(a, (N - 1)//q, N) - 1, N) == 1)
            basen[str(q)] = a; unter.append(zertifikat(q))
        Fh = 1
        for q, e in fk.items(): Fh *= q**e
        return {'N': str(N), 'art': 'pocklington', 'F': [[str(q), e] for q, e in fk.items()], 'R': str((N - 1)//Fh), 'basen': basen, 'faktoren': unter}
    zh = hand(Nt, Ft)
    klein = dict(list(Ft.items())[:3]); zk = hand(Nt, klein)                     # too small an F: F² < N
    zb = json.loads(json.dumps(zh)); qb = zb['F'][-1][0]; zb['basen'][qb] = 1    # wrong base: gcd(1 − 1, N) = N
    e3 = (z30 is not None and pruefe(z30) and z40 is not None and pruefe(z40) and zertifikat(zus) is None and not mr12(carm)
          and zertifikat(carm) is None and int(zh['R']) > 1 and pruefe(zh) and not pruefe(zk) and not pruefe(zb))
    sag(f'PK E3 {"✅" if e3 else "❌"}  10³⁰+{N30 - 10**30} und 10⁴⁰+{N40 - 10**40} bewiesen und bestaetigt; Semiprim + Carmichael abgelehnt; '
        f'Teil-Zertifikat ({len(str(Nt))} St., R mit {len(zh["R"])} St.) angenommen, mit F² < N und mit falscher Basis abgelehnt.')
    assert e3, 'E3 VERLETZT'

    # E4: partition of the 150 new pairs
    # `pts` = lattice points of the w9 run; `cnt` = points per kernel; `rang` = kernels by decreasing count; `K25` = the 25
    #   kernels with the most points; `paare` = the new pairs (m, d): odd d >= 5 dividing k, for lattice points (m, m', k, ...) of
    #   kernels outside K25
    pts = json.loads((ERG/'w9_v31_M100000000_H2000_result.json').read_text(encoding='utf-8'))['points']
    cnt = Counter(p[0] for p in pts); rang = sorted(cnt, key=lambda m: (-cnt[m], m)); K25 = set(rang[:25])
    paare = set()
    for m, mp_, k, st, status, z in pts:
        for d in range(5, k + 1, 2):
            if k % d == 0 and m not in K25: paare.add((m, d))
    w133 = json.loads((ERG/'w133_probelauf_79_kerne_result.json').read_text(encoding='utf-8'))
    w134 = json.loads((ERG/'w134_stufe2_rangsieb_result.json').read_text(encoding='utf-8'))
    ck = json.loads((ERG/'w135_stufe2_ecm_checkpoint.json').read_text(encoding='utf-8'))
    if PROBE:
        w135_erl = [[e['m'], e['d'], e['zeuge'], e['route']] for e in ck.values() if e['fertig']]
        w135_off = [[e['m'], e['d']] for e in ck.values() if not e['fertig']]
    else:
        assert not (BOX/'ECM_LAEUFT.lock').exists(), 'w135 laeuft noch (Sperrdatei) — erst nach dem Lauf, oder „probe"'
        w135 = json.loads((ERG/'w135_stufe2_ecm_result.json').read_text(encoding='utf-8'))
        w135_erl = w135['erledigt']; w135_off = [[o[0], o[1]] for o in w135['offen']]
    # `quelle` = source: (m, d) -> (finding route, candidate witnesses[, reported witness]); `offen` = pairs still open in w135
    quelle = {}
    for k, p in w133['zeugen'].items():
        m, d = map(int, k.split('_'))
        if m not in K25: quelle[(m, d)] = ('w133', [int(p)])
    for m, d, p, r in w134['erledigt']:
        assert (m, d) not in quelle; quelle[(m, d)] = ('w134', [int(p)])
    for m, d, z, r in w135_erl:
        assert (m, d) not in quelle
        e = ck[f'{m}_{d}']
        quelle[(m, d)] = ('w135/' + r, sorted({int(x) for x in e['faktoren']} | {int(z)}), int(z))
    offen = {(m, d) for m, d in w135_off}
    e4 = (set(quelle) | offen == paare and not (set(quelle) & offen) and len(paare) == 150)
    n_w = Counter(v[0].split('/')[0] for v in quelle.values())
    sag(f'PK E4 {"✅" if e4 else "❌"}  Partition: w133 {n_w["w133"]} + w134 {n_w["w134"]} + w135 {n_w["w135"]} + offen {len(offen)} = '
        f'{len(quelle) + len(offen)} von {len(paare)} neuen Paaren, disjunkt.')
    assert e4, ('E4 VERLETZT', len(paare), len(quelle), len(offen))

    # main run
    # lists: `bewiesen` = proved witnesses, `probable` = class (b), `nur_probable` = only one probable factor, `fehler` = errors,
    #   `ablehn` = other factors that are not witnesses, `zert` = certificates by pair
    bewiesen, probable, nur_probable, fehler, ablehn, zert = [], [], [], [], [], {}
    for (m, d) in sorted(quelle):
        # `weg` = finding route, `kand` = candidates, `gemeldet` = reported witness
        q = quelle[(m, d)]; weg, kand = q[0], q[1]; gemeldet = q[2] if len(q) > 2 else kand[0]
        ok_liste = []   # candidates that pass both routes
        for p in kand:
            a, l = beide(m, d, p)
            if a != l: fehler.append((m, d, p, 'A ≠ L')); continue
            if not a:
                if p == gemeldet: fehler.append((m, d, p, 'gemeldeter Zeuge faellt in A und L durch'))
                else: ablehn.append((m, d, len(str(p))))
                continue
            ok_liste.append(p)
        beweis = None
        for p in ok_liste:
            ta = time.time(); z = zertifikat(p)
            if z is not None and pruefe(z):
                beweis = (p, z['art'], round(time.time() - ta, 1)); zert[f'{m}_{d}'] = {str(p): z}; break
        if beweis:
            bewiesen.append([m, d, str(beweis[0]), weg, beweis[1]])
            if len(str(beweis[0])) >= 20 or weg.startswith('w135'):
                sag(f'  ({m}, {d}) {weg:<11} Zeuge {len(str(beweis[0])):>3} St.: ✅ A + L, bewiesen prim ({beweis[1]}, {beweis[2]}s)')
        elif ok_liste and weg in ('w135/prim', 'w135/voll'):
            p = ok_liste[0]
            # `drei` = three: the three probable-prime tests
            drei = (all(_mr_basis(p, a) for a in BASEN12), all(_mr_basis(p, a) for a in BASEN_B), bool(sympy.isprime(p)))
            assert all(drei), ('probable-Zeuge besteht nicht alle drei Tests', m, d)
            probable.append([m, d, str(p), weg, list(drei)])
            sag(f'  ({m}, {d}) {weg:<11} Zeuge {len(str(p)):>3} St.: ✅ A + L · ⚠️ KEIN BEWEIS — Klasse (b) vollzerlegt, probable '
                f'(12 Basen {drei[0]}, 12 andere {drei[1]}, BPSW {drei[2]})')
        elif ok_liste:
            nur_probable.append([m, d, str(ok_liste[0]), weg])
            sag(f'  ({m}, {d}) {weg:<11}: ⚠️ nur ein wahrscheinlicher Faktor, Block nicht vollzerlegt — OFFEN')
        else:
            sag(f'  ({m}, {d}) {weg:<11}: 🔴 kein Kandidat besteht A und L')
    sag()
    # `groesse` = size classes of the proved witnesses
    groesse = Counter(('≤10⁹' if int(b[2]) <= 10**9 else ('<ψ₁₂' if int(b[2]) < PSI12 else '>ψ₁₂')) for b in bewiesen)
    sag(f'ERGEBNIS: {len(quelle)} Paare mit Fund · (a) zertifiziert {len(bewiesen)} (Zeugen {dict(groesse)}) · (b) vollzerlegt, probable {len(probable)} · '
        f'nur wahrscheinlicher Faktor {len(nur_probable)} · Fehler {len(fehler)} · weitere Faktoren, die keine Zeugen sind: {len(ablehn)} · '
        f'offen (w135) {len(offen)} ⇒ erledigt {len(bewiesen) + len(probable)} von {len(paare)}, offen {len(offen) + len(nur_probable)}')
    for f in fehler: sag(f'  🔴 {f}')
    sag(f'Laufzeit {time.time() - t0:.0f}s.')
    if PROBE:
        sag('[PROBE — nichts geschrieben]'); return
    # certificate file: `psi12` = ψ₁₂, `verfahren` = method, `zertifikate` = certificates by pair "m_d"
    (ERG/'w136_zertifikate.json').write_text(json.dumps(dict(skript=pathlib.Path(__file__).name, psi12=str(PSI12),
        verfahren='Pocklington mit F·R = N−1, F voll zerlegt, F² > N, rekursiv; unter psi12 deterministischer 12-Basen-MR', zertifikate=zert),
        ensure_ascii=False), encoding='utf-8')
    # `geladen` = certificates read back from the file
    geladen = json.loads((ERG/'w136_zertifikate.json').read_text(encoding='utf-8'))['zertifikate']
    n_ok = sum(1 for k in geladen for p in geladen[k] if pruefe(geladen[k][p]))
    assert n_ok == len(bewiesen), 'Zertifikate aus der Datei nicht alle bestaetigt'
    sag(f'Zweiter Pruefer: {n_ok} Zertifikate aus der Datei gelesen und bestaetigt.')
    # `res` = result dict: `skript` = script name, `datum` = date, `neue_paare` = number of new pairs, `mit_fund` = with a
    #   finding, `zertifiziert` = certified (a), `vollzerlegt_probable` = completely factored, probable (b), `erledigt` = settled,
    #   `bewiesen` = proved witnesses, `nur_probable` = only a probable factor, `offen` = open pairs, `fehler` = errors,
    #   `weitere_nichtzeugen` = other factors that are not witnesses, `laufzeit_s` = run time in seconds
    res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), neue_paare=len(paare),
               mit_fund=len(quelle), zertifiziert=len(bewiesen), vollzerlegt_probable=len(probable),
               erledigt=len(bewiesen) + len(probable), bewiesen=bewiesen, probable=probable, nur_probable=nur_probable,
               offen=sorted([list(x) for x in offen] + [[b[0], b[1]] for b in nur_probable]),
               fehler=[list(map(str, f)) for f in fehler], weitere_nichtzeugen=len(ablehn), laufzeit_s=round(time.time() - t0, 1))
    (ERG/'w136_stufe2_zweite_route_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
    (ERG/'w136_stufe2_zweite_route_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
    assert not fehler, 'w136: Gegenrechnung faellt durch — Funde NICHT verbuchen'

if __name__ == '__main__':
    main()
