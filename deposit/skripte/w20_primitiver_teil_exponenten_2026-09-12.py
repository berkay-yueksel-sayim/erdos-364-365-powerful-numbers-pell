# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w20_primitiver_teil_exponenten_2026-09-12.py
# Purpose: can the primitive part Phi_d of T_d be powerful? For the primes p of rank d (alpha_m(p) = d) it measures the exact
#   exponent v_p(Phi_d) = v_p(T_d) (capped at 3 via T_d mod p^3) and compares the counts with the heuristic sum of 1/p.
# Reads: ergebnisse/w9_v31_M100000000_H2000_result.json (lattice points (m, k)).
# Writes: ergebnisse/w20_primteil_P<P>_N<NKERNE>_NP<NP>_result.json; prints the tables (A), (A'), (B) and (C).
# Usage: python w20_primitiver_teil_exponenten_2026-09-12.py [P=3000000] [NKERNE=25] [NP=5]
#   P = sieve bound for the primes p, NKERNE = number of kernels m (those with the most lattice points),
#   NP = number of primes kept per rank d.
# Controls: positive controls PK1-PK6 and negative controls PK-a, PK-b (listed below) run first; any failure aborts the run.
#
# Background. Murty-Seguin 2019, Lemma 2.3, bounds the exponent of an INDEX prime (p | n with p | Phi_n(a)) by 1. The
#   primitive part considered here has no such primes: by Cor. 4.8(i), p = +-1 (mod 4d), hence p >= 4d - 1 > d. The
#   literature says nothing about the exponent of a PRIMITIVE prime, which is the quantity measured here.
#   Unlike script w17c, which computes T_d mod p^2 and keeps only a Boolean ("does the rank have a witness"),
#   this script keeps the number (T_d mod p^3): v = 1, v = 2 or v >= 3.
#
# Two populations, kept strictly separate:
#   P1  question population: d odd, d | k for at least one lattice point (m, k). Reported: v_p(Phi_d), size log10 Phi_d,
#       known fraction, ranks without witness.
#       The values d = 1 and d = 3 divide k but lie outside the question (Phi_1 is the part of T_1, and Cor. 4.8 requires
#       d >= 5, Carmichael); they are tallied separately (`P1k_V`, `P1k_je_p`), neither averaged in nor discarded.
#   P2  broad population (in the style of Prop. 6.1): every rank hit by some p <= P. Only the v distribution and the 1/p
#       comparison. It is not the object of the question (d need not divide any k, and d may be even), but it is a much
#       larger sample for the Wieferich rate than Prop. 6.1 (m <= 1e4, p <= 3000).
#
# Controls (the run aborts if one fails):
#   PK1  ranks at m = 7: alpha(29) = alpha(197) = alpha(2857) = 7; 3 has no rank (negative side, PK-b).
#   PK2  Mollin-Walsh: T_7(7) = 2^3 * 29 * 197 * 2857, and v_29 = v_197 = v_2857 = 1.
#   PK3  congruence (Cor. 4.8(i)): 29, 197, 2857 are all +-1 (mod 28); checked at run time for EVERY (d, p).
#   PK4  size formula on the exact number: Phi_7(7) = 29*197*2857 = 16322041, log10 = 7.2128. Prediction
#        phi(2d) * log10(eps_7) = 6 * 1.20255 = 7.2145, deviation < 0.01.
#        PK4b: the variant "phi(2d)/2 * log eps" (it gives 3.6) MUST fail. PK4c: sum_{e | d} phi(2e) = d for odd d.
#   PK5  the totals of script w17c at P = 3e6, NKERNE = 25, NP = 5:
#        378 points, 2536 ranks, 2025 occupied, 1988 with witness, 37 Wieferich only, 511 unoccupied.
#   PK6  P1 contains only odd d, and the number of occupied P1 ranks (m, d) must not exceed the number of occupied
#        rank occurrences counted per point (PK5).
#   PK-a  an artificial value 29^2 must give v = 2 and 29^3 must give v = 3 (meaning >= 3), so the counter can see a defect.
#   PK-b  rank 3 at m = 7 does not exist, so it must not be counted as occupied.

import json, math, sys, time, pathlib
from math import isqrt
from collections import Counter

sys.set_int_max_str_digits(2000000)
P      = int(sys.argv[1]) if len(sys.argv) > 1 else 3000000
NKERNE = int(sys.argv[2]) if len(sys.argv) > 2 else 25
NP     = int(sys.argv[3]) if len(sys.argv) > 3 else 5

HIER = pathlib.Path(__file__).resolve().parent
ERG  = HIER.parent / 'ergebnisse'
PUNKTE_JSON = ERG / 'w9_v31_M100000000_H2000_result.json'

# ---------------------------------------------------------------- building blocks (as in w17c)
def sieve_spf(n):
    spf = list(range(n + 1))  # `spf[n]` = smallest prime factor of n (sieve)
    for i in range(2, isqrt(n) + 1):
        if spf[i] == i:
            for j in range(i*i, n + 1, i):
                if spf[j] == j: spf[j] = i
    return spf

def factor(n, spf):
    f = {}
    while n > 1:
        p = spf[n]; f[p] = f.get(p, 0) + 1; n //= p
    return f

def divisors(n, spf):
    ds = [1]
    for q, e in factor(n, spf).items():
        ds = [d * q**i for d in ds for i in range(e + 1)]
    return sorted(ds)

def phi_von(n, spf):  # Euler's totient phi(n)
    r = n
    for q in factor(n, spf): r = r // q * (q - 1)
    return r

def fund(m):  # fundamental solution (T1, U1) of x^2 - m*y^2 = 1: convergents (h0, k0) of sqrt(m) until h0^2 - m*k0^2 = 1
    a0 = isqrt(m); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (m - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0

def eps_pow(T1, U1, m, e, M):
    """Return (ra, rb) with eps^e = ra + rb*sqrt(m) (mod M); ra = T_e (mod M)."""
    ra, rb, ba, bb = 1 % M, 0, T1 % M, U1 % M
    while e:
        if e & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        e >>= 1
    return ra, rb

def order(T1, U1, m, p, spf):
    N = p - 1 if pow(m % p, (p - 1)//2, p) == 1 else p + 1
    o = N
    for q in factor(N, spf):
        while o % q == 0 and eps_pow(T1, U1, m, o//q, p) == (1, 0): o //= q
    return o

def log10_int(n):
    if n <= 0: return float('-inf')
    b = n.bit_length()
    if b <= 900: return math.log10(n)
    s = b - 900
    return math.log10(n >> s) + s * math.log10(2)

def log10_eps(T1, U1, m):
    """log10(eps), eps = T1 + U1*sqrt(m). Exact for small values, asymptotic (eps ~ 2*U1*sqrt(m)) for large ones."""
    if T1.bit_length() < 900:
        return math.log10(T1 + U1 * math.sqrt(m))
    return log10_int(U1) + 0.5*math.log10(m) + math.log10(2)

def v_p_von_T(T1, U1, m, d, p):
    """v_p(T_d), capped at 3 (return value 3 means '>= 3'). Assumes p | T_d."""
    Td = eps_pow(T1, U1, m, d, p**3)[0]
    if Td % p != 0: return 0
    if Td % (p*p) != 0: return 1
    if Td % (p**3) != 0: return 2
    return 3

# ---------------------------------------------------------------- positive/negative controls before the run
t0 = time.time()
spf = sieve_spf(P + 2)
primes = [i for i in range(3, P + 1) if spf[i] == i]
print(f'W20: P = {P} ({len(primes)} ungerade Primzahlen), NKERNE = {NKERNE}, NP = {NP}; Sieb {time.time()-t0:.1f}s', flush=True)

T1_7, U1_7 = fund(7)  # (T1_7, U1_7) = (8, 3) is the fundamental solution of x^2 - 7*y^2 = 1
assert (T1_7, U1_7) == (8, 3), ('PK fund(7)', T1_7, U1_7)
for p, soll in ((29, 7), (197, 7), (2857, 7)):
    o = order(T1_7, U1_7, 7, p, spf)
    assert o % 4 == 0 and o//4 == soll, ('PK1 Rang', p, o)
    assert p % 28 in (1, 27), ('PK3 Kongruenz mod 4d', p, p % 28)
    assert v_p_von_T(T1_7, U1_7, 7, 7, p) == 1, ('PK2 Exponent 1', p)
assert order(T1_7, U1_7, 7, 3, spf) % 4 != 0, 'PK-b: Rang 3 existiert bei m = 7 nicht'
T7_exakt = eps_pow(T1_7, U1_7, 7, 7, 10**30)[0]
assert T7_exakt == 8 * 29 * 197 * 2857, ('PK2 Mollin-Walsh T_7(7)', T7_exakt)
phi14  = phi_von(14, spf)
vorher = phi14 * log10_eps(T1_7, U1_7, 7)  # `vorher` = predicted log10 Phi_7(7) = phi(2d) * log10(eps_7)
echt   = log10_int(29 * 197 * 2857)  # `echt` = actual log10 of 29*197*2857
assert abs(vorher - echt) < 0.01, ('PK4 Groessenformel phi(2d)*log eps', vorher, echt)
assert abs(phi14/2 * log10_eps(T1_7, U1_7, 7) - echt) > 3, 'PK4b: die halbierte Fassung MUSS durchfallen'
assert sum(phi_von(2*e, spf) for e in divisors(7, spf)) == 7, 'PK4c: sum_{e|d} phi(2e) = d'
# PK-a: an artificial value 29^2 (resp. 29^3) must give v = 2 (resp. 3)
for kunst, soll in ((29**2, 2), (29**3, 3)):
    Td = kunst % (29**3)
    v = 1 if Td % (29*29) != 0 else (2 if Td % (29**3) != 0 else 3)
    assert v == soll, ('PK-a', kunst, v, soll)
print(f'PK ok: Raenge, Exponent 1, Kongruenz mod 4d, Mollin-Walsh T_7(7), '
      f'Groessenformel ({vorher:.4f} gegen {echt:.4f}), halbierte Fassung faellt durch, PK- Zaehler', flush=True)

# ---------------------------------------------------------------- data
# each point is a list (m, m-prime, k, digits, status, witness), unpacked below
d7 = json.load(open(PUNKTE_JSON, encoding='utf-8'))
cnt = Counter(pt[0] for pt in d7['points'])
KERNE = [m for m, _ in cnt.most_common(NKERNE)]  # the NKERNE kernels with the most lattice points
print('Kerne:', KERNE, flush=True)

# P1 = question population, P2 = broad population
P1_V, P2_V = Counter(), Counter()
P1_je_p, P2_je_p = {}, {}
# d = 1 and d = 3 lie OUTSIDE the question: Phi_1 is the part of T_1, and Cor. 4.8 requires d >= 5
# (Carmichael). They are tallied separately, neither averaged in nor discarded.
P1k_V, P1k_je_p = Counter(), {}
P1_ohne_zeuge = []  # `P1_ohne_zeuge` = ranks without witness: (m, d, [p...], [v...], klein); all found primes have v >= 2
P1_groessen   = []  # `P1_groessen` = sizes: (m, d, log10 Phi_d, log10 of the known part, number of primes)
TOT = Counter()  # totals for PK5 / PK6

# `kopf` = table header line
kopf = (f"{'m':>9} {'Pkt':>4} {'P1-d':>6} {'P1(m,d,p)':>10} {'v=1':>7} {'v=2':>4} {'v>=3':>5} "
        f"{'P2(m,d,p)':>10} {'v>=2':>5} {'bes.':>6} {'Zeuge':>6} {'nurW':>5}  Zeit")
print(kopf, flush=True)

for m in KERNE:
    tm = time.time()
    T1, U1 = fund(m)
    lg_eps = log10_eps(T1, U1, m)

    # collect the ranks once (serves both populations): `rank_ps` maps a rank d to up to NP primes p of that rank
    rank_ps = {}
    for p in primes:
        if m % p == 0: continue
        o = order(T1, U1, m, p, spf)
        if o % 4: continue
        d = o // 4
        lst = rank_ps.setdefault(d, [])
        if len(lst) < NP: lst.append(p)

    pts = [pt for pt in d7['points'] if pt[0] == m]

    # ---- P1: the divisors of the indices k of genuine lattice points (the object of the question)
    d_frage = set()
    for (mm, mp, k, dig, status, w) in pts:
        assert k % 2 == 1, ('k muss ungerade sein', m, k)
        d_frage.update(divisors(k, spf))
    assert all(d % 2 == 1 for d in d_frage), 'PK6: P1 enthaelt nur ungerade d'

    m1V = Counter(); besetzt_p1 = 0  # `m1V` = v histogram of this kernel in P1; `besetzt_p1` = occupied P1 ranks of this kernel
    for d in sorted(d_frage):
        ps = rank_ps.get(d)
        if not ps: continue
        klein = (d < 5)  # `klein` = small d, outside the question: Phi_1 belongs to T_1, Cor. 4.8 requires d >= 5
        if not klein: besetzt_p1 += 1
        vs = []
        for p in ps:
            v = v_p_von_T(T1, U1, m, d, p)
            assert v >= 1, ('Rang-Widerspruch: p teilt T_d nicht', m, d, p)
            assert p % (4*d) in (1, 4*d - 1), ('Cor. 4.8(i) verletzt', m, d, p, p % (4*d))
            vs.append(v)
            if klein:
                P1k_V[v] += 1; P1k_je_p[p] = P1k_je_p.get(p, 0) + 1
            else:
                P1_V[v] += 1; m1V[v] += 1; P1_je_p[p] = P1_je_p.get(p, 0) + 1
        if min(vs) >= 2:
            P1_ohne_zeuge.append((m, d, ps[:], vs[:], klein))
        if not klein:
            lgPhi  = phi_von(2*d, spf) * lg_eps  # predicted log10 Phi_d = phi(2d) * log10(eps), see PK4
            # log10 of the known part (product of p^v over the found primes)
            lg_bek = sum(v * math.log10(p) for p, v in zip(ps, vs))
            P1_groessen.append((m, d, lgPhi, lg_bek, len(ps)))

    # ---- P2: all hit ranks (a DIFFERENT question, only the v distribution)
    m2V = Counter()  # v histogram of this kernel in P2
    for d, ps in rank_ps.items():
        for p in ps:
            v = v_p_von_T(T1, U1, m, d, p)
            P2_V[v] += 1; m2V[v] += 1
            P2_je_p[p] = P2_je_p.get(p, 0) + 1

    # ---- the w17c totals as PK5 (counted per point, not per rank); `tau_sum` = total number of divisors of the k,
    # `occ` = occupied (point, d) pairs, `wit` = of these with a witness, `onlyw` = Wieferich only;
    # a witness is a prime p with p^2 not dividing T_d and p not dividing k/d
    tau_sum = occ = wit = onlyw = 0
    for (mm, mp, k, dig, status, w) in pts:
        ds = divisors(k, spf); tau_sum += len(ds)
        for d in ds:
            ps = rank_ps.get(d)
            if not ps: continue
            occ += 1
            has = False
            for p in ps:
                Td = eps_pow(T1, U1, m, d, p*p)[0]
                if Td % (p*p) != 0 and (k // d) % p != 0: has = True; break
            if has: wit += 1
            else: onlyw += 1
    TOT['tau'] += tau_sum; TOT['occ'] += occ; TOT['wit'] += wit
    TOT['onlyw'] += onlyw; TOT['pts'] += len(pts); TOT['p1d'] += besetzt_p1
    print(f'{m:>9} {len(pts):>4} {besetzt_p1:>6} {sum(m1V.values()):>10} '
          f'{m1V[1]:>7} {m1V[2]:>4} {m1V[3]:>5} {sum(m2V.values()):>10} {m2V[2]+m2V[3]:>5} '
          f'{occ:>6} {wit:>6} {onlyw:>5}  {time.time()-tm:.0f}s', flush=True)

# ---------------------------------------------------------------- PK5 / PK6
# `SOLL` = expected w17c totals (reference values); `ist` = totals of this run
SOLL = dict(pts=378, tau=2536, occ=2025, wit=1988, onlyw=37, unbes=511)
ist  = dict(pts=TOT['pts'], tau=TOT['tau'], occ=TOT['occ'], wit=TOT['wit'],
            onlyw=TOT['onlyw'], unbes=TOT['tau'] - TOT['occ'])
if (P, NKERNE, NP) == (3000000, 25, 5):
    assert ist == SOLL, ('PK5 W17c nicht reproduziert', ist, SOLL)
    print(f'\nPK5 ok: W17c Ziffer fuer Ziffer reproduziert {ist}', flush=True)
else:
    print(f'\nPK5 uebersprungen (andere Argumente); W17c-Zaehlung dieses Laufs: {ist}', flush=True)
print(f"PK6: P1 = {TOT['p1d']} besetzte (m, d)-Raenge; W17c zaehlt {TOT['occ']} besetzte Rang-VORKOMMEN "
      f"ueber die Punkte (ein d zaehlt bei jedem Punkt neu, daher >= P1).", flush=True)
assert TOT['p1d'] <= TOT['occ'], 'PK6: P1 darf nicht mehr Raenge haben als Vorkommen'

# ---------------------------------------------------------------- evaluation
N1 = sum(P1_V.values()); N2 = sum(P2_V.values())
print(f"\n=== (A) P1 — FRAGE-POPULATION: v_p(Phi_d) fuer d | k, ueber {N1} Paare (m, d, p) ===")
for v in (1, 2, 3):
    lab = 'v >= 3' if v == 3 else f'v = {v}'
    print(f"  {lab:>7}: {P1_V[v]:>7}  ({100*P1_V[v]/max(N1,1):.4f} %)")
# heuristic: the expected number of pairs with v >= 2 is the sum of c_p / p, with v >= 3 the sum of c_p / p^2
# (c_p = number of pairs (m, d, p) with this prime p)
e2_1 = sum(c/p    for p, c in P1_je_p.items())
e3_1 = sum(c/p**2 for p, c in P1_je_p.items())
print(f"  Heuristik ueber DIESE Paare: v >= 2 erwartet {e2_1:.2f}, beobachtet {P1_V[2]+P1_V[3]}"
      f"  |  v >= 3 erwartet {e3_1:.4f}, beobachtet {P1_V[3]}")
Nk   = sum(P1k_V.values())
e2_k = sum(c/p for p, c in P1k_je_p.items())
print(f"  Getrennt gehalten (AUSSERHALB der Frage): d = 1 und d = 3, {Nk} Paare "
      f"[v=1: {P1k_V[1]}, v=2: {P1k_V[2]}, v>=3: {P1k_V[3]}], Erwartung v>=2 dort {e2_k:.2f}.")
print(f"    -> Diese {Nk} Paare tragen {e2_k:.2f} der Erwartung, die Frage-Population {e2_1:.2f}. "
      f"Zusammen {e2_k + e2_1:.2f}.")
print(f"    DESHALB faellt die Erwartung beim Schnitt so stark: sum 1/p wird von KLEINEN p beherrscht,")
print(f"    und die kleinsten Primzahlen sitzen genau bei d = 1 und d = 3 (bei m = 39 ist es p = 5).")
print(f"    Grund: Phi_1 ist der Teil von T_1, und Cor. 4.8 setzt d >= 5 voraus (Carmichael).")
print(f"    Zaehlte man sie mit, stuende der Eintrag (m = 39, d = 1, p = 5, v = 2)")
print(f"    dadurch faelschlich in der Liste (C) der zeugenlosen Raenge.")

print(f"\n=== (A') P2 — BREITE POPULATION (andere Frage): {N2} Paare (m, d, p) ===")
for v in (1, 2, 3):
    lab = 'v >= 3' if v == 3 else f'v = {v}'
    print(f"  {lab:>7}: {P2_V[v]:>7}  ({100*P2_V[v]/max(N2,1):.5f} %)")
e2_2 = sum(c/p    for p, c in P2_je_p.items())
e3_2 = sum(c/p**2 for p, c in P2_je_p.items())
print(f"  Heuristik ueber DIESE Paare: v >= 2 erwartet {e2_2:.2f}, beobachtet {P2_V[2]+P2_V[3]}"
      f"  |  v >= 3 erwartet {e3_2:.4f}, beobachtet {P2_V[3]}")
print(f"  Reichweite: P2 ist NICHT Wand 04s Objekt (d muss kein Teiler eines k sein, d auch gerade).")
print(f"  Prop. 6.1 misst am Apparitionsrang ueber m <= 1e4, p <= 3000: 470 gegen 486 erwartet.")
print(f"  P2 ist dieselbe Art Groesse in viel groesserer Stichprobe -> auf KONSISTENZ lesen, nicht auf Gleichheit.")

print(f"\n=== (B) GROESSE UND BEKANNTER ANTEIL (nur P1) ===")
gl = sorted(P1_groessen, key=lambda x: x[2])
ges_lg  = sum(g[2] for g in P1_groessen)
ges_bek = sum(g[3] for g in P1_groessen)
print(f"  {len(P1_groessen)} Raenge (m, d) mit mindestens einer gefundenen Primzahl.")
print(f"  log10 Phi_d: min {gl[0][2]:.2f}, Median {gl[len(gl)//2][2]:.2f}, max {gl[-1][2]:.2f}")
print(f"  Bekannter Anteil: {100*ges_bek/ges_lg:.4f} % der Dezimalstellen von Phi_d.")
print(f"  -> Der unbekannte Rest ist die Schranke jeder Aussage 'Phi_d ist nicht powerful'.")

print(f"\n=== (C) P1-RAENGE (d >= 5), DEREN GEFUNDENE PRIMZAHLEN ALLE v >= 2 HABEN ===")
echte = [x for x in P1_ohne_zeuge if not x[4]]
ausser = [x for x in P1_ohne_zeuge if x[4]]
print(f"  Anzahl: {len(echte)} von {len(P1_groessen)} Raengen")
for (m, d, ps, vs, _) in echte[:15]:
    print(f"    m = {m:>9}, d = {d:>7}, Primzahlen {ps}, v = {vs}")
if ausser:
    print(f"  Zusaetzlich {len(ausser)} mit d < 5, AUSSERHALB der Frage (nur zur Vollstaendigkeit):")
    for (m, d, ps, vs, _) in ausser[:5]:
        print(f"    m = {m:>9}, d = {d:>7}, Primzahlen {ps}, v = {vs}")
print(f"\n  Lesart: das sind KEINE powerful Phi_d — nur Raenge, bei denen die ersten {NP} Primzahlen")
print(f"  keinen Zeugen lieferten. Der unbekannte Rest aus (B) entscheidet.")

# result record; keys: `erwartet_v_ge_2` = expected number with v >= 2, `beobachtet_v_ge_2` = observed number,
# `n_paare` = number of pairs, `n_raenge` = number of ranks, `anteil_bekannt_prozent` = known fraction in percent,
# `raenge_ohne_zeugen` = ranks without witness, `ausserhalb_d_klein` = tallies for d = 1 and d = 3
erg = dict(
    skript='w20_primitiver_teil_exponenten_2026-09-12.py', fassung=2, P=P, NKERNE=NKERNE, NP=NP,
    kerne=KERNE,
    P1=dict(n_paare=N1, v={str(k): v for k, v in P1_V.items()},
            erwartet_v_ge_2=e2_1, beobachtet_v_ge_2=P1_V[2]+P1_V[3],
            erwartet_v_ge_3=e3_1, beobachtet_v_ge_3=P1_V[3],
            n_raenge=len(P1_groessen), log10_summe=ges_lg, log10_bekannt=ges_bek,
            anteil_bekannt_prozent=100*ges_bek/ges_lg,
            raenge_ohne_zeugen=[[m, d, ps, vs] for (m, d, ps, vs, kl) in P1_ohne_zeuge if not kl],
            ausserhalb_d_klein=dict(n_paare=sum(P1k_V.values()),
                                    v={str(k): v for k, v in P1k_V.items()},
                                    raenge_ohne_zeugen=[[m, d, ps, vs] for (m, d, ps, vs, kl) in P1_ohne_zeuge if kl])),
    P2=dict(n_paare=N2, v={str(k): v for k, v in P2_V.items()},
            erwartet_v_ge_2=e2_2, beobachtet_v_ge_2=P2_V[2]+P2_V[3],
            erwartet_v_ge_3=e3_2, beobachtet_v_ge_3=P2_V[3]),
    w17c_reproduktion=ist, w17c_soll=SOLL,
    laufzeit_s=time.time()-t0, python=sys.version.split()[0],
)
out = ERG / f'w20_primteil_P{P}_N{NKERNE}_NP{NP}_result.json'
out.write_text(json.dumps(erg, indent=1), encoding='utf-8')
print(f"\nErgebnis: {out.name}   Zeit gesamt {time.time()-t0:.0f}s")
