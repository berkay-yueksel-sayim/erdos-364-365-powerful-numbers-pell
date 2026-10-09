# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w95_primzertifikate_2026-09-24.py
# -*- coding: utf-8 -*-
# Primality proof for the witnesses from w89 that lie above psi_12.
# Background: w74 certifies every witness with the Miller-Rabin test on the first twelve prime bases, which is deterministic
#   ONLY below psi_12 = 318665857834031151167461 (24 digits; Sorenson-Webster 2017). Of the 23 ranks that w89 closed, 8 have even
#   their SMALLEST witness above it (25-35 digits). Without a proof they would only be "probable", and the statement
#   "every witness is certified" would be false.
# Method (own code): Pocklington (1914) with FULLY factored N - 1:
#   If N - 1 = prod q_i^e_i and for each q_i there is an a with a^(N-1) = 1 (mod N) and gcd(a^((N-1)/q_i) - 1, N) = 1,
#   then every prime divisor of N is = 1 (mod N - 1), hence N is prime. (This is the special case F = N - 1 of Pocklington's
#   theorem; general case: Brillhart-Lehmer-Selfridge 1975, Math. Comp. 29(130), 620-647.) The q_i themselves are proved below
#   psi_12 by the deterministic 12-base test, above it RECURSIVELY by the same method.
#   The FACTORIZATION of N - 1 is delivered by sympy.factorint (library, used but not read). It is uncritical: the certificate
#   contains the product and every factor is checked itself; a wrong factorization fails there.
# Reads: ergebnisse/w89_gmp_ecm_result.json. Writes: ergebnisse/w95_primzertifikate.json (the certificates),
#   ergebnisse/w95_primzertifikate_result.json and ergebnisse/w95_primzertifikate_output.txt. Requires SymPy. No arguments.
# Certificates are produced for the smallest witness of every rank and for every further witness above psi_12 up to 40 digits.
#
# Controls (expectations, each an assert; the printed lines are named E1 ... E5):
#   E1  PK positive: the next prime after 10^30 (sympy.nextprime, BPSW) gets a certificate that the second checker confirms.
#   E2  PK negative: a 30-digit product of two primes and a Carmichael number are REJECTED.
#   E3  Each of the 8 smallest witnesses above psi_12 gets a certificate, and so does every further witness up to 40 digits.
#   E4  A SECOND, separately written checker reads the certificates from the JSON and confirms them without factoring
#       anything (only product, powers, gcd).
#   E5  The cofactors with 66, 70 and 342 digits (for "completely factored") stay probable primes, like the factors of the
#       complete factorizations in the text. They are not needed for the count: every rank has a smaller witness. They get the
#       same three tests as the factors in Prop. 6.10: the first twelve bases, twelve OTHER bases (41 to 89), BPSW; all must say
#       "probable prime" (w89 itself only tests with the first twelve).
import json, sys, pathlib, time
from math import gcd
import sympy
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent
ERG  = HIER.parent/'ergebnisse'   # `ERG` = `ergebnisse` (results)
aus = []   # `aus` = output lines
# `sag` = say: print a line and record it
def sag(s=''):
    print(s, flush=True); aus.append(s)
PSI12 = 318665857834031151167461   # psi_12: below it the strong test to the first twelve prime bases (`BASEN12`) is a proof
BASEN12 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

def mr12(n):
    """Miller-Rabin on the first twelve prime bases -- deterministic for n < psi_12."""
    assert n < PSI12
    if n < 2: return False
    for p in BASEN12:
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in BASEN12:
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

# certificate = dict with `art` (type: 'mr12' or 'pocklington'), `n_minus_1` (factorization of N - 1), `basen` (bases a per q),
# `faktoren` (certificates of the q); `tiefe` = recursion depth
def zertifikat(N, tiefe=0):
    """Returns a certificate or None (then N is NOT proved prime)."""
    if N < PSI12:
        return {'N': str(N), 'art': 'mr12'} if mr12(N) else None
    fk = sympy.factorint(N - 1)
    prod = 1
    for q, e in fk.items(): prod *= q**e
    if prod != N - 1: return None
    # `basen` = bases a per q, `unter` = sub-certificates (lower levels)
    basen, unter = {}, []
    for q in sorted(fk):
        for a in range(2, 200):
            if pow(a, N - 1, N) != 1: return None           # Fermat fails => N composite
            if gcd(pow(a, (N - 1)//q, N) - 1, N) == 1: basen[str(q)] = a; break
        else: return None
        z = zertifikat(q, tiefe + 1)
        if z is None: return None
        unter.append(z)
    return {'N': str(N), 'art': 'pocklington', 'n_minus_1': [[str(q), e] for q, e in sorted(fk.items())],
            'basen': basen, 'faktoren': unter}

# ---- second checker, separate: factors nothing, only recomputes (`pruefe` = check)
def pruefe(z):
    N = int(z['N'])
    if z['art'] == 'mr12':
        # independent form: trial division up to sqrt(N) is cheap for N < 10^12, otherwise Miller-Rabin with the same set of bases
        # but as a separate implementation (witness loop written differently)
        if N < 10**12:
            if N < 2: return False
            i = 2
            while i * i <= N:
                if N % i == 0: return False
                i += 1
            return True
        return N < PSI12 and all(_mr_basis(N, a) for a in BASEN12)
    prod = 1
    for q, e in z['n_minus_1']: prod *= int(q)**e
    if prod != N - 1: return False
    # `kinder` = children: the sub-certificates by their N
    kinder = {k['N']: k for k in z['faktoren']}
    for q, _ in z['n_minus_1']:
        a = z['basen'][q]
        if pow(a, N - 1, N) != 1 or gcd(pow(a, (N - 1)//int(q), N) - 1, N) != 1: return False
        if not pruefe(kinder[q]): return False
    return True
def _mr_basis(n, a):                       # strong probable-prime test of n to the single base a
    if n % a == 0: return n == a
    t = n - 1; r = 0
    while not t & 1: t >>= 1; r += 1
    y = pow(a, t, n)
    if y == 1 or y == n - 1: return True
    for _ in range(r - 1):
        y = pow(y, 2, n)
        if y == n - 1: return True
    return False

sag('='*96); sag('w95 — PRIMALITAETS-BEWEIS (Pocklington, rekursiv) fuer die w89-Zeugen ueber ψ₁₂'); sag('='*96)
t0 = time.time()
# E1 (positive control)
N1 = sympy.nextprime(10**30)
z = zertifikat(N1); assert z and pruefe(z), 'E1 VERLETZT'
sag(f'PK E1 ✅  naechste Primzahl nach 10³⁰ (sympy) = 10³⁰ + {N1 - 10**30}: bewiesen und vom zweiten Pruefer bestaetigt.')
# E2 (negative controls: `zus` = product of two primes, `carm` = Carmichael number)
zus = sympy.nextprime(10**14) * sympy.nextprime(3 * 10**15)
carm = 7 * 13 * 17 * 23 * 31 * 67 * 73                                      # = 5394826801, Carmichael (Korselt: p-1 | n-1)
assert carm == 5394826801 and all((carm - 1) % (q - 1) == 0 for q in (7, 13, 17, 23, 31, 67, 73))
assert pow(2, carm - 1, carm) == 1                                          # it passes the Fermat test to base 2
assert zertifikat(zus) is None and not mr12(carm) and zertifikat(carm) is None, 'E2 VERLETZT'
sag(f'PK E2 ✅  30-stelliges Produkt abgelehnt; Carmichael-Zahl 5394826801 (besteht Fermat zur Basis 2) abgelehnt.')

w89 = json.load(open(ERG/'w89_gmp_ecm_result.json', encoding='utf-8'))
# `je_zahl` = per number (rank); `zeugen` = witnesses; `kleinster` = smallest; `rolle` = role; `ergeb` = certificates per rank;
# `fehl` = ranks whose smallest witness could not be proved
ergeb, fehl = {}, []
for k, v in w89['je_zahl'].items():
    zeugen = sorted({int(x) for x in v['zeugen']})
    if not zeugen: continue
    kleinster = zeugen[0]
    for p in zeugen:
        if p < PSI12 and p != kleinster: continue
        if len(str(p)) > 40 and p != kleinster: continue
        ta = time.time()
        z = zertifikat(p)
        ok = z is not None and pruefe(z)
        rolle = 'kleinster' if p == kleinster else 'weiterer'
        sag(f'  {k:<10} {rolle:<9} {len(str(p)):>3} St.: ' + ('✅ bewiesen prim (' + z['art'] + ')' if ok else '🔴 KEIN BEWEIS')
            + f'  {time.time()-ta:.1f}s')
        if ok: ergeb.setdefault(k, {})[str(p)] = z
        elif p == kleinster: fehl.append((k, p))
# `ueber` = ranks whose smallest witness lies above psi_12
ueber = [k for k, v in w89['je_zahl'].items() if v['zeugen'] and min(int(x) for x in v['zeugen']) >= PSI12]
sag()
sag(f'E3      Raenge, deren kleinster Zeuge ueber ψ₁₂ liegt: {len(ueber)} — ' +
    ('✅ alle mit Beweis.' if not fehl else f'🔴 ohne Beweis: {fehl}'))
# JSON keys: `skript` = script, `datum` = date, `psi12`, `verfahren` = method, `zertifikate` = certificates (by rank and witness);
# `json_pfad` = path of the JSON file
json_pfad = ERG/'w95_primzertifikate.json'
json_pfad.write_text(json.dumps(dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d'), psi12=str(PSI12),
                                     verfahren='Pocklington mit voll zerlegtem N-1, rekursiv; unter psi12 deterministischer 12-Basen-MR',
                                     zertifikate=ergeb), ensure_ascii=False), encoding='utf-8')
# E4: read afresh from the file (`geladen` = loaded certificates) and confirm with the second checker ONLY
geladen = json.load(open(json_pfad, encoding='utf-8'))['zertifikate']
n4 = sum(1 for k in geladen for p in geladen[k] if pruefe(geladen[k][p]))
n_ges = sum(len(v) for v in geladen.values())
assert n4 == n_ges, 'E4 VERLETZT'
sag(f'PK E4 ✅  {n4} Zertifikate aus der Datei gelesen und vom zweiten Pruefer bestaetigt (ohne Zerlegung).')
# E5: the factors of the complete factorizations WITHOUT certificate get the same three tests as the factors in the text of
#   Prop. 6.10 -- twelve bases, twelve OTHER bases (41..89), BPSW. w89 itself only tests with the first twelve.
#   Expectation: all three tests say "probable prime" for each. (`voll` = completely factored, `drei` = the three test results)
BASEN_B = (41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89)
probable = []
for k, v in w89['je_zahl'].items():
    if not v['voll']: continue
    for f in dict.fromkeys(v['zeugen']):
        if f in geladen.get(k, {}) or int(f) < PSI12: continue
        n = int(f)
        drei = (all(_mr_basis(n, a) for a in BASEN12), all(_mr_basis(n, a) for a in BASEN_B), bool(sympy.isprime(n)))
        sag(f'E5      {k:<10} Kofaktor {len(f):>3} St.: 12 Basen {drei[0]}, 12 andere Basen {drei[1]}, BPSW {drei[2]}')
        assert all(drei), ('E5 VERLETZT', k, len(f))
        probable.append([k, len(f)])
sag(f'        Laufzeit {time.time()-t0:.0f}s.')
# result keys: `raenge_ueber_psi12` = ranks above psi_12, `ohne_beweis` = ranks without proof, `n_zertifikate` = number of
# certificates, `probable_drei_tests` = [rank, digits] of the probable primes tested three ways
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d'), raenge_ueber_psi12=ueber, ohne_beweis=[list(map(str, f)) for f in fehl],
           n_zertifikate=n_ges, probable_drei_tests=probable)
(ERG/'w95_primzertifikate_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w95_primzertifikate_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
assert not fehl, 'w95: nicht alle kleinsten Zeugen bewiesen'
