# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w137_wieferich_suche_79kerne_2026-09-29.py
# -*- coding: utf-8 -*-
# Purpose: Wieferich-type search (as in w77) on all 79 kernels m (`kerne` = kernels). For each kernel and each odd prime
#   p <= 3*10^6 with p not dividing m: order of eps mod p, rank d = order/4 (if 4 | order), v_p(T_d) from eps^d mod p^3.
# Reads: ergebnisse/w9_v31_M100000000_H2000_result.json (lattice points), ergebnisse/w77_wieferich_ereignisse_result.json.
# Writes: ergebnisse/w137_wieferich_79kerne_result.json and ergebnisse/w137_wieferich_79kerne_output.txt. No arguments.
# Controls: E1 (reproduces w77 on the 25 old kernels), E2 (congruence), E3 (size test), E4 (negative control); described below.
#
# Background. Remark 6.11 and question Q3 depend on eps_m Wieferich events: a triple (m, d, p) with rank alpha_m(p) = d and
#   v_p(T_d) >= 2. If p were the ONLY prime of rank d, then Phi_d = p^2 would be powerful.
#   The earlier run w77 covered only 25 kernels.
# Counted TWICE (`gedeckelt` = capped, `ungedeckelt` = uncapped):
#   capped   - only the first NP = 5 primes per rank (the counting of w20/w77, so the numbers are comparable),
#   uncapped - every prime (the full population that Remark 6.11 needs for "below 3*10^6").
#   All 79 kernels, not only the 54 new ones: one common base population for both counts; the 25 old kernels serve as
#   positive control.
#
# Expectations and controls (checked at run time):
#   E1  PK: capped on the 25 old kernels exactly w77: v=1: 1 932 631, v=2: 16, v>=3: 0, and the same 16 triples.
#   E2  Every event (capped as well as uncapped) satisfies p = +-1 (mod 4d)  (Cor. 4.8(i)).
#   E3  For every event with d >= 5: log10 Phi_d > 2*log10 p + 1, hence Phi_d != p^2. If E3 fails, that is a COUNTEREXAMPLE
#       CANDIDATE to Q3: abort loudly, record nothing, check by hand.
#   E4  NK: the step function `stufe` detects artificially set values: p*c => 1, p^2*c => 2, 0 => >= 3, c => 0 (c coprime to p).
# Measured, not expected: the number of events on the 54 new kernels (capped, uncapped), of these those in the question population
#   (d >= 5, d | k of a lattice point), and per rank the number of witnesses v = 1 below 3*10^6; plus the heuristic sum of 1/p.
# Building blocks: our own code, taken from w77/w20.


import json, math, sys, time, pathlib
from math import isqrt
from collections import Counter, defaultdict
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'
# `P_MAX` = bound for the primes p, `NP` = number of primes per rank counted in the capped count
P_MAX, NP = 3_000_000, 5

# `sieve_spf` = smallest-prime-factor sieve;  `factor` = prime factorization {prime: exponent}
def sieve_spf(n):
    spf = list(range(n + 1))
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
# `phi_von` = Euler's totient of n
def phi_von(n, spf):
    r = n
    for q in factor(n, spf): r = r // q * (q - 1)
    return r
# `fund` = fundamental solution (T1, U1) of x^2 - m*y^2 = 1 by continued fraction (big integers)
def fund(m):
    a0 = isqrt(m); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (m - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
# `eps_pow` = (ra, rb) with eps^e = ra + rb*sqrt(m) (mod M); ra = T_e (mod M)
def eps_pow(T1, U1, m, e, M):
    ra, rb, ba, bb = 1 % M, 0, T1 % M, U1 % M
    while e:
        if e & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        e >>= 1
    return ra, rb
# `order` = multiplicative order of eps modulo the prime p (a divisor of p-1 or p+1)
def order(T1, U1, m, p, spf):
    N = p - 1 if pow(m % p, (p - 1)//2, p) == 1 else p + 1
    o = N
    for q in factor(N, spf):
        while o % q == 0 and eps_pow(T1, U1, m, o//q, p) == (1, 0): o //= q
    return o
def log10_int(n):
    b = n.bit_length()
    if b <= 900: return math.log10(n)
    s = b - 900; return math.log10(n >> s) + s*math.log10(2)
def log10_eps(T1, U1, m):
    if T1.bit_length() < 900: return math.log10(T1 + U1*math.sqrt(m))
    return log10_int(U1) + 0.5*math.log10(m) + math.log10(2)
# `stufe` = step (level) function: the exponent v_p of Td, capped at 3
def stufe(Td, p):
    # v_p from T_d mod p^3: 0, 1, 2 or 3 (= at least 3)
    if Td % p != 0: return 0
    if Td % (p*p) != 0: return 1
    if Td % (p**3) != 0: return 2
    return 3
def v_p_von_T(T1, U1, m, d, p):
    return stufe(eps_pow(T1, U1, m, d, p**3)[0], p)
# `teiler` = the set of divisors of n
def teiler(n):
    ds = [i for i in range(1, isqrt(n)+1) if n % i == 0]; return set(ds + [n//i for i in ds])

# `aus` = output lines (also written to the output file); `sag` = say: print a line and keep it
t0 = time.time(); aus = []
def sag(s=''):
    print(s, flush=True); aus.append(s)
sag('='*100); sag('w137 — WIEFERICH-SUCHE (w77-Art) AUF ALLEN 79 KERNEN   ' + time.strftime('%Y-%m-%d %H:%M')); sag('='*100)

# E4 first: does the step function detect a defect?
p_t = 2857
e4 = (stufe(p_t*7, p_t) == 1 and stufe(p_t*p_t*7 % p_t**3, p_t) == 2 and stufe(0, p_t) == 3 and stufe(7, p_t) == 0)
assert e4, 'E4 VERLETZT — Stufenfunktion'
sag('NK E4 ✅  Stufenfunktion: p·7 → 1, p²·7 → 2, 0 → ≥3, 7 → 0.')

spf = sieve_spf(P_MAX + 2)
primes = [i for i in range(3, P_MAX + 1) if spf[i] == i]
d7 = json.load(open(ERG/'w9_v31_M100000000_H2000_result.json', encoding='utf-8'))
cnt = Counter(pt[0] for pt in d7['points'])
# `RANG` = ranking: the kernels sorted by number of lattice points (descending), then by m
RANG = [m for m, _ in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))]
assert len(RANG) == 79, ('Kernzahl', len(RANG))
# `K25` = the 25 old kernels (those of w77), `NEU` = the 54 new ones
K25, NEU = RANG[:25], RANG[25:]
# `frage` = question population: pairs (m, d) with d >= 5 dividing the index k of a lattice point of m
frage = set()
for pt in d7['points']:
    for d in teiler(pt[2]):
        if d >= 5: frage.add((pt[0], d))

# `V_ged` = per kernel, histogram of v in the capped count; `ev_ged`, `ev_alle` = events (m, d, p, v), capped / all
V_ged = {m: Counter() for m in RANG}; ev_ged, ev_alle = [], []
# `zeugen_je_rang` = number of witnesses (v = 1) per rank (m, d); `primes_je_rang` = number of primes per rank;
# `heur` = heuristic sum of 1/p per kernel
zeugen_je_rang = defaultdict(int); primes_je_rang = defaultdict(int); heur = defaultdict(float)
for m in RANG:
    T1, U1 = fund(m)
    # `rank_all` = maps a rank d to ALL primes p of that rank
    rank_all = defaultdict(list)
    for p in primes:
        if m % p == 0: continue
        o = order(T1, U1, m, p, spf)
        if o % 4: continue
        rank_all[o // 4].append(p)
    for d, ps in rank_all.items():
        for j, p in enumerate(ps):
            v = v_p_von_T(T1, U1, m, d, p)
            primes_je_rang[(m, d)] += 1; heur[m] += 1/p
            if v == 1: zeugen_je_rang[(m, d)] += 1
            if v >= 2: ev_alle.append((m, d, p, v))
            # capped count: the first NP primes of the rank
            if j < NP:
                V_ged[m][v] += 1
                if v >= 2: ev_ged.append((m, d, p, v))
    sag(f'   m = {m:>9} ({"alt" if m in K25 else "neu"}): {sum(len(v) for v in rank_all.values()):>8} Rang-Primzahlen, '
        f'Ereignisse gedeckelt {sum(1 for e in ev_ged if e[0] == m):>2}, ungedeckelt {sum(1 for e in ev_alle if e[0] == m):>2}   ({time.time()-t0:.0f}s)')

sag()
# E1: reproduction of w77 on the 25 (`soll_trip` = expected triples, `ist_trip` = triples found here)
w77 = json.load(open(ERG/'w77_wieferich_ereignisse_result.json', encoding='utf-8'))
V25 = Counter()
for m in K25: V25.update(V_ged[m])
soll_trip = {(e['m'], e['d'], e['p'], e['v']) for e in w77['ereignisse']}
ist_trip = {e for e in ev_ged if e[0] in set(K25)}
e1 = (V25[1] == 1_932_631 and V25[2] == 16 and V25[3] == 0 and ist_trip == soll_trip)
sag(f'PK E1 {"✅" if e1 else "❌"}  gedeckelt auf den 25: v=1 {V25[1]:,} · v=2 {V25[2]} · v>=3 {V25[3]} · Tripel wie w77: {ist_trip == soll_trip}')
assert e1, ('E1 VERLETZT — w77 nicht reproduziert', dict(V25), sorted(ist_trip ^ soll_trip))
e2 = all(p % (4*d) in (1, 4*d - 1) for m, d, p, v in ev_alle)
sag(f'PK E2 {"✅" if e2 else "❌"}  alle {len(ev_alle)} Ereignisse (ungedeckelt) erfuellen p ≡ ±1 (mod 4d).')
assert e2, 'E2 VERLETZT'

# E3 + tally per event (uncapped; the capped ones are a subset); `gegen` = counterexample candidates
# (`gegen` = against), `bilanz` = tally, `fund_c` = cache of fundamental solutions
gegen = []; bilanz = []
fund_c = {}
for m, d, p, v in sorted(ev_alle):
    if m not in fund_c: fund_c[m] = fund(m)
    T1, U1 = fund_c[m]
    # `lgPhi` = log10 Phi_d by the size formula phi(2d) * log10(eps); `lg2p` = log10 of p^2
    lgPhi = phi_von(2*d, spf) * log10_eps(T1, U1, m); lg2p = 2*math.log10(p)
    ist_p2 = lgPhi <= lg2p + 1.0
    if ist_p2 and d >= 5: gegen.append((m, d, p))
    bilanz.append(dict(m=m, d=d, p=p, v=v, kern='alt' if m in K25 else 'neu', gedeckelt=(m, d, p, v) in set(ev_ged),
                       log10_phi=round(lgPhi, 2), log10_p2=round(lg2p, 2), phi_ist_p2=ist_p2,
                       zeugen_unter_3e6=zeugen_je_rang[(m, d)], rangprimzahlen_unter_3e6=primes_je_rang[(m, d)],
                       in_frage_population=(m, d) in frage))
sag(f"{'m':>9} {'d':>7} {'p':>9} {'v':>2} {'Kern':>4} {'ged.':>4} {'log10 Phi_d':>12} {'2·log10 p':>10} {'Zeugen<3e6':>11} {'Frage-Pop.':>10}")
for b in bilanz:
    sag(f"{b['m']:>9} {b['d']:>7} {b['p']:>9} {b['v']:>2} {b['kern']:>4} {'ja' if b['gedeckelt'] else '':>4} {b['log10_phi']:>12.1f} "
        f"{b['log10_p2']:>10.1f} {b['zeugen_unter_3e6']:>11} {'ja' if b['in_frage_population'] else '':>10}")
sag()
if gegen:
    sag('🔴🔴🔴 E3 VERLETZT — Phi_d = p² moeglich bei: ' + str(gegen) + '   GEGENBEISPIEL-KANDIDAT ZU Q3. NICHT VERBUCHEN, VON HAND PRUEFEN.')
    raise SystemExit('GEGENBEISPIEL-KANDIDAT')
sag('PK E3 ✅  Kein Ereignis mit d >= 5 hat Phi_d = p² (Groessenvergleich).')

neu_ged = [b for b in bilanz if b['kern'] == 'neu' and b['gedeckelt']]
neu_alle = [b for b in bilanz if b['kern'] == 'neu']
alt_alle = [b for b in bilanz if b['kern'] == 'alt']
neu_F = [b for b in neu_alle if b['in_frage_population']]
heur_neu = sum(heur[m] for m in NEU)
V54 = Counter()
for m in NEU: V54.update(V_ged[m])
sag('='*100)
sag(f'  54 NEUE KERNE: Ereignisse gedeckelt (NP = {NP}) {len(neu_ged)} · ungedeckelt {len(neu_alle)} · davon in der Frage-Population {len(neu_F)}')
sag(f'                 gedeckelte Zaehlung v=1 {V54[1]:,} · v=2 {V54[2]} · v>=3 {V54[3]}')
sag(f'                 Heuristik (ungedeckelt): Summe 1/p ueber alle Rang-Primzahlen = {heur_neu:.2f} erwartete Ereignisse')
sag(f'  25 ALTE KERNE: ungedeckelt {len(alt_alle)} (gedeckelt 16 = w77)')
sag(f'  Frage-Population, neue Kerne: ' + (', '.join(f"({b['m']}, {b['d']}, {b['p']}) Zeugen<3e6: {b['zeugen_unter_3e6']}" for b in neu_F) or 'keine'))
sag('='*100)
# result record: `gedeckelt_25` / `gedeckelt_54` = capped v counts of the 25 old / 54 new kernels,
# `ereignisse` = events, `gegenbeispiele` = counterexample candidates, `neu_frage_population` = events of the new kernels
# in the question population, `heuristik_summe_1_durch_p_neu` = heuristic sum of 1/p for the new kernels
res = dict(skript=pathlib.Path(__file__).name, datum='2026-09-29', P=P_MAX, NP=NP, kerne=len(RANG),
           gedeckelt_25=dict(V25), gedeckelt_54=dict(V54), ereignisse=bilanz, gegenbeispiele=gegen,
           neu_gedeckelt=len(neu_ged), neu_ungedeckelt=len(neu_alle), neu_frage_population=len(neu_F),
           alt_ungedeckelt=len(alt_alle), heuristik_summe_1_durch_p_neu=round(heur_neu, 2), laufzeit_s=round(time.time()-t0, 1))
(ERG/'w137_wieferich_79kerne_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w137_wieferich_79kerne_output.txt').write_text('\n'.join(aus)+'\n', encoding='utf-8')
print(f'\nErgebnis: w137_wieferich_79kerne_result.json   Zeit {time.time()-t0:.0f}s')
