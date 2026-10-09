# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w144_abstand1_stufe2_2026-09-29.py
# -*- coding: utf-8 -*-
# Part II, Sec. 7, stage 2 of the distance-1 route: the 918 open cases from w27.
# Starting point: w26 builds all pairs (a, a + 1) of powerful numbers with D = sqfree(a(a + 1)) <= 20 000 and T = 2a + 1 < 10^2000
#   (6063 pairs) and searches a witness p < 5*10^5 for the third number (a - 1 or a + 2; the other side is already settled by 2);
#   w27 raises the bound to 5*10^6. 918 cases remain "without witness below 5*10^6", which does NOT mean triples. w27 itself names
#   one of them, (D, k) = (3, 24), whose third number is a prime: there only the SEARCH was open, not the question.
# Stage 2 (as in Part I, Prop. 6.8), for each third number N:
#   (1) N exactly from the Pell solution: T = T_k(D), N = (T - 3)/2 resp. (T + 3)/2.
#   (2) divide out all primes p < P = 5*10^6 (gcd with the primorial, then the multiplicity of each). Each must have
#       exponent >= 2, otherwise w27 would be wrong (a prime with exponent 1 is reported as class `kleiner_zeuge`, small witness).
#   (3) Cofactor C:  C = 1 => N is powerful => TRIPLE CANDIDATE => loud abort.
#                    C (probably) prime => C divides N exactly once => witness C; certified if C < psi_12 (a strong test to the
#                      first twelve prime bases is a proof there, Sorenson-Webster 2017), otherwise "probably prime" (BPSW,
#                      gmpy2 + sympy, two routes).
#                    C a square or higher power => open (noted).
#                    otherwise, for C < 10^90: ECM (sympy; B1 = 50000, B2 = 5000000, at most 400 curves); a prime factor q found
#                      with v_q(N) = 1 is a witness.
#                    otherwise open.
# Reads: ergebnisse/w26_gap1_D20000_H2000_result.json and ergebnisse/w27_tiefpass_P5000000_result.json.
# Writes: ergebnisse/w144_abstand1_stufe2_result.json. Requires gmpy2 and SymPy. Usage: python w144_abstand1_stufe2_2026-09-29.py
# Controls (a failure stops the run before the main loop):
#   PK: (D, k) = (3, 24), side a - 1 => N itself is prime and certified (the head of w27 names exactly this case).
#   NK: for five cases already settled by w26 (witness p < 5*10^5), step (2) must report the witness p as "exponent 1";
#       otherwise step (2) would be useless.
# Own code.
import json, math, sys, time, pathlib
import gmpy2
from sympy import factorint, isprime as sympy_isprime
from sympy.ntheory import ecm
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'   # `ERG` = `ergebnisse` (results)
P = 5_000_000                                  # bound for the small primes of step (2)
PSI12 = 318665857834031151167461               # psi_12: below it the strong test to the first twelve prime bases is a proof
BASEN12 = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)

# `fund(D)` = fundamental solution (T_1, U_1) of x^2 - D*y^2 = 1 by continued fraction; `T_exakt(D, k)` = exact (T_k, U_k)
def fund(D):
    a0 = math.isqrt(D); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - D*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (D - Pp*Pp)//Q; a = (a0 + Pp)//Q; h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
def T_exakt(D, k):
    t1, u1 = fund(D); ra, rb, ba, bb, n = 1, 0, t1, u1, k
    while n:
        if n & 1: ra, rb = ra*ba + D*rb*bb, ra*bb + rb*ba
        ba, bb = ba*ba + D*bb*bb, 2*ba*bb; n >>= 1
    return ra, rb
# `dritte` = third number N for the pair (a, a + 1) given by (D, k); `seite` = side: 'a-1' gives N = a - 1, otherwise N = a + 2
def dritte(D, k, seite):
    T, U = T_exakt(D, k); assert T*T - D*U*U == 1 and T % 2 == 1
    a = (T - 1)//2
    # plausibility: (a, a + 1) is a pair, a(a + 1) = D*V^2 with D | V (Part II, distance-one proposition)
    assert U % 2 == 0 and (U//2) % D == 0, (D, k)
    return (T - 3)//2 if seite == 'a-1' else (T + 3)//2
# `mr12` = strong probable-prime test (Miller-Rabin) to the first twelve prime bases `BASEN12`
def mr12(n):
    d, r = n - 1, 0
    while d % 2 == 0: d //= 2; r += 1
    for b in BASEN12:
        if n == b: return True
        x = pow(b, d, n)
        if x in (1, n - 1): continue
        for _ in range(r - 1):
            x = x*x % n
            if x == n - 1: break
        else: return False
    return True

PRIMORIAL = gmpy2.primorial(P - 1)
# `klein_abteilen` = divide out the small primes: returns (C, exp1), C = cofactor, exp1 = primes dividing N exactly once
def klein_abteilen(N):
    g = int(gmpy2.gcd(PRIMORIAL, N)); C = N; exp1 = []
    for p in (factorint(g) if g > 1 else {}):
        v = 0
        while C % p == 0: C //= p; v += 1
        if v == 1: exp1.append(p)
    return C, exp1

# `pruefe` = check one case. Result class (`klasse`): `kleiner_zeuge` = small witness, `TRIPEL_KANDIDAT` = triple candidate,
# `prim_zertifiziert` / `prim_wahrscheinlich` = cofactor is a certified / probable prime witness, `offen_potenz` = open (power),
# `ecm_zeuge_*` = witness found by ECM, `offen` = open. Result fields: `stellen` = digits of N, `kofaktor_stellen` = digits of C,
# `zeuge` = witness, `zeuge_ist_N` = the witness is N itself, `ecm_s` = ECM seconds, `ecm_versucht` = ECM was tried
def pruefe(D, k, seite, ecm_s=20.0):
    N = dritte(D, k, seite); C, exp1 = klein_abteilen(N)
    if exp1: return dict(klasse='kleiner_zeuge', zeuge=str(min(exp1)), stellen=len(str(N)))
    r = dict(stellen=len(str(N)), kofaktor_stellen=len(str(C)))
    if C == 1:
        return dict(r, klasse='TRIPEL_KANDIDAT')
    if gmpy2.is_prime(C, 30) and sympy_isprime(C):
        zert = C < PSI12 and mr12(C)
        return dict(r, klasse='prim_zertifiziert' if zert else 'prim_wahrscheinlich', zeuge=str(C) if C < 10**40 else f'{len(str(C))}-stellig',
                    zeuge_ist_N=(C == N))
    if gmpy2.is_power(C):
        return dict(r, klasse='offen_potenz')
    if C < 10**90:
        t0 = time.time()
        try:
            fs = ecm(C, B1=50000, B2=5000000, max_curve=400, seed=1)
        except Exception:
            fs = set()
        for q in sorted(fs):
            v, x = 0, N
            while x % q == 0: x //= q; v += 1
            if v == 1 and gmpy2.is_prime(q, 30) and sympy_isprime(q):
                zert = q < PSI12 and mr12(q)
                return dict(r, klasse='ecm_zeuge' + ('_zertifiziert' if zert else '_wahrscheinlich'), zeuge=str(q) if q < 10**40 else f'{len(str(q))}-stellig',
                            ecm_s=round(time.time() - t0, 1))
        return dict(r, klasse='offen', ecm_versucht=True, ecm_s=round(time.time() - t0, 1))
    return dict(r, klasse='offen')

if __name__ == '__main__':
    t0 = time.time()
    w26 = json.loads((ERG/'w26_gap1_D20000_H2000_result.json').read_text(encoding='utf-8'))
    w27 = json.loads((ERG/'w27_tiefpass_P5000000_result.json').read_text(encoding='utf-8'))
    assert w27['P'] == P
    print('=' * 100); print(f'w144 — Abstand-1-Route, Stufe 2 auf den {len(w27["bleibt_offen"])} offenen Faellen aus w27   {time.strftime("%Y-%m-%d %H:%M")}')
    print('=' * 100)
    pk = pruefe(3, 24, 'a-1'); pk_ok = pk['klasse'] == 'prim_zertifiziert' and pk.get('zeuge_ist_N')
    print(f'PK (3, 24, a − 1): {pk} → {"✅" if pk_ok else "❌"}')
    # `nk_faelle` = NK cases: the first five entries of w26 with a witness x[4] > 1000 (`faelle` = cases)
    nk_faelle = [x for x in w26['ueberlebt'] if x[4] and x[4] > 1000][:5]
    nk = [(x, pruefe(x[0], x[1], x[3])) for x in nk_faelle]
    nk_ok = all(r['klasse'] == 'kleiner_zeuge' and int(r['zeuge']) <= x[4] for x, r in nk)
    for x, r in nk: print(f'NK {x[:4]} Zeuge w26 = {x[4]}: {r}')
    print(f'NK: {"✅" if nk_ok else "❌"}')
    if not (pk_ok and nk_ok): print('❌ Kontrollen gescheitert — kein Lauf.'); sys.exit(1)
    erg = []   # `erg` = results per open case
    for i, (D, k, st, seite, _) in enumerate(w27['bleibt_offen']):
        r = pruefe(D, k, seite); r.update(D=D, k=k, seite=seite); erg.append(r)
        if r['klasse'] == 'TRIPEL_KANDIDAT':
            print(f'🔴 ({D}, {k}, {seite}): dritte Zahl POWERFUL — TRIPEL-KANDIDAT. ABBRUCH.'); break
        if (i + 1) % 50 == 0: print(f'  {i + 1} / {len(w27["bleibt_offen"])} · {time.time() - t0:.0f} s', flush=True)
    from collections import Counter
    z = Counter(r['klasse'] for r in erg)
    print(f'\nERGEBNIS: {dict(z)} · {time.time() - t0:.0f} s')
    # JSON keys: `skript` = script, `datum` = date, `quelle` = source, `faelle` = number of cases, `klassen` = counts per class,
    # `ergebnisse` = results per case, `sekunden` = seconds
    (ERG/'w144_abstand1_stufe2_result.json').write_text(json.dumps(dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'),
        P=P, quelle='w27_tiefpass_P5000000_result.json', pk=pk, pk_ok=pk_ok, nk=[[x, r] for x, r in nk], nk_ok=nk_ok, faelle=len(erg),
        klassen=dict(z), ergebnisse=erg, sekunden=round(time.time() - t0, 1)), ensure_ascii=False, indent=1), encoding='utf-8')
