# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# -*- coding: utf-8 -*-
# w138_karte_961_daten_2026-09-29.py
# The 961 pairs (m, d): map of the 811 pairs of w92 plus the 150 new pairs of the 54 further kernels (79 kernels in total).
# Reads only result files (no new search). The 811 pairs come UNCHANGED from w92; the 150 new ones get their route from the
#   path by which the witness was found (w133 "erste", w134 "sieb8"/"sieb9", w135 "prim"/"ecm135"/"voll") - but ONLY if w136
#   confirmed the witness by two routes. CLASS per pair as in Prop. 6.10 (a): "zert" (witness proven prime), "voll" (completely
#   factored into probable primes), "offen" (open). For the 811 the class follows from the route (erste/sieb8/sieb9/ecm81/gmpecm
#   -> zert, voll/prim -> voll); for the 150 it comes from w136 (a pair with route "prim" can be CERTIFIED there if the
#   primality proof succeeded). Two later updates: blocks that w158 closes via their algebraic parts (class "voll" or the new
#   class "prpzeuge" = witness that is a probable prime), and blocks whose factors are all proven prime by w162-w164 (class
#   "zert").
# Unequal effort: the 811 pairs additionally received GMP-ECM up to t40 (w89); the 150 only up to t35 (w135).
# Reads (ergebnisse/): w9_v31_M100000000_H2000_result.json, w92_karte_811_daten.json, w133_probelauf_79_kerne_result.json,
#   w134_stufe2_rangsieb_result.json, w135_stufe2_ecm_result.json, w135_stufe2_ecm_checkpoint.json,
#   w136_stufe2_zweite_route_result.json, w137_wieferich_79kerne_result.json, w158_teile_wahrscheinlich_prim_result.json,
#   w162_inventar_wahrscheinliche_primzahlen_result.json, w163_pocklington_beweise_result.json, w164_ecpp_pruefer_result.json.
# Writes: ergebnisse/w138_karte_961_daten.json. No command-line arguments.
# Controls (expectations fixed in advance):
#   E1 (PK)  961 pairs and 79 kernels, rebuilt from w9 (d >= 5, d | k of a lattice point) = 811 (25 kernels) + 150 (54
#            kernels).
#   E2 (PK)  The 811 are exactly as in w92 (set of pairs and count per route).
#   E3       The 150: zert + voll + offen = 150 and the numbers agree with w136; the 811: zert / voll / offen as in Prop. 6.10
#            (a), checked deliberately at the state before the w158 update.
#   E4 (PK)  Size: log10 D_d from the Moebius formula gives the EXACT digit count (floor + 1) for the blocks of w135.
#   E5       Every Wieferich event of w137 with `in_frage_population` = True lies in the set of 961 pairs.
#   E6       The blocks closed by w158 were open before. E7: count of blocks upgraded by primality proofs (printed, not
#            asserted).
import json, math, pathlib, sys
from math import isqrt
from collections import Counter
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'
# `J` = load a JSON file from the results folder.
def J(n): return json.load(open(ERG/n, encoding='utf-8'))
# `fund` = fundamental solution (T_1, U_1) of x^2 - m*y^2 = 1 by continued fractions; `teiler` = divisors of n;
# `mu` = Moebius function.
def fund(m):
    a0 = isqrt(m); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (m - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
def teiler(n):
    ds = [i for i in range(1, isqrt(n)+1) if n % i == 0]; return sorted(set(ds + [n//i for i in ds]))
def mu(n):
    r, x, p = 1, n, 2
    while p*p <= x:
        if x % p == 0:
            x //= p
            if x % p == 0: return 0
            r = -r
        p += 1
    return -r if x > 1 else r
# `lg_T(L, e)` = log10 T_e for L = log10 of the fundamental unit; `lg_Phi(L, d)` = log10 of the primitive part D_d by Moebius
# inversion over the divisors of d; `lg_eps(m)` = log10 of the fundamental unit T_1 + U_1*sqrt(m).
def lg_T(L, e): return e*L - math.log10(2) + math.log10(1 + 10**(-2*e*L))
def lg_Phi(L, d): return sum(mu(d//e) * lg_T(L, e) for e in teiler(d))
def lg_eps(m):
    T1, U1 = fund(m)
    return math.log10(T1 + U1*math.sqrt(m)) if T1 < 10**15 else math.log10(2*T1)

# `cnt` = number of lattice points per kernel; `RANG` = kernels ordered by decreasing count; `K25` = the 25 kernels with the
#   most points;
# `menge` = all pairs (m, d) with d | k (d >= 5) for a lattice point (m, k); `alt` = pairs of the 25 kernels, `neu` = pairs of
#   the other 54.
w9 = J('w9_v31_M100000000_H2000_result.json')
cnt = Counter(p[0] for p in w9['points']); RANG = sorted(cnt, key=lambda m: (-cnt[m], m)); K25 = set(RANG[:25])
menge = set()
for m, _, k, *_ in w9['points']: menge |= {(m, d) for d in teiler(k) if d >= 5}
alt = {md for md in menge if md[0] in K25}; neu = menge - alt
e1 = len(menge) == 961 and len(RANG) == 79 and len(alt) == 811 and len(neu) == 150
print(f'PK E1 {"✅" if e1 else "❌"}  {len(menge)} Paare, {len(RANG)} Kerne = {len(alt)} (25 Kerne) + {len(neu)} (54 Kerne).')
assert e1

w92 = J('w92_karte_811_daten.json'); assert not w92['w89_vorlaeufig']
e2 = {(r['m'], r['d']) for r in w92['raenge']} == alt and dict(Counter(r['route'] for r in w92['raenge'])) == w92['routen']
print(f'PK E2 {"✅" if e2 else "❌"}  die 811 wie w92: {w92["routen"]}')
assert e2

# w133 ... w136: the earlier runs; `bew` / `prob` = witnesses of w136 proven prime / only probably prime, keyed by (m, d);
# `route` = how the witness was found, `klasse` = class (zert / voll / offen, later also prpzeuge), `KL_ALT` = route -> class
#   for the 811.
w133, w134, w135, w136 = J('w133_probelauf_79_kerne_result.json'), J('w134_stufe2_rangsieb_result.json'), \
                         J('w135_stufe2_ecm_result.json'), J('w136_stufe2_zweite_route_result.json')
assert not w136['fehler']
bew = {(b[0], b[1]): b for b in w136['bewiesen']}; prob = {(b[0], b[1]): b for b in w136['probable']}
route = {}
for k in w133['zeugen']:
    m, d = map(int, k.split('_'))
    if (m, d) in neu: route[(m, d)] = 'erste'
for m, d, p, r in w134['erledigt']: route[(m, d)] = r
for m, d, z, r in w135['erledigt']: route[(m, d)] = r
klasse = {}
KL_ALT = dict(erste='zert', sieb8='zert', sieb9='zert', ecm81='zert', gmpecm='zert', voll='voll', prim='voll', offen='offen')
for r in w92['raenge']: klasse[(r['m'], r['d'])] = KL_ALT[r['route']]
for md in neu:
    if md in bew: klasse[md] = 'zert'
    elif md in prob: klasse[md] = 'voll'
    else: klasse[md] = 'offen'; route[md] = 'offen'
    assert md in route, ('erledigt, aber ohne Fundweg', md)
k_neu = Counter(klasse[md] for md in neu); k_alt = Counter(klasse[md] for md in alt)
e3 = (sum(k_neu.values()) == 150 and k_neu['zert'] == w136['zertifiziert'] and k_neu['voll'] == w136['vollzerlegt_probable']
      and k_neu['offen'] == len(w136['offen']) and dict(k_alt) == dict(zert=767, voll=38, offen=6))
# E3 expects 767 / 38 / 6 for the 811: the pair (95, 475) was settled by a 38-digit witness found in a later ECM block, which
#    w92 took over as route "gmpecm" (second route w90, proof w95).
print(f'PK E3 {"✅" if e3 else "❌"}  die 150: zert {k_neu["zert"]} + voll {k_neu["voll"]} + offen {k_neu["offen"]} — wie w136 · '
      f'die 811: {dict(k_alt)} — wie Prop. 6.10 (a).')
assert e3
# Update from w158: blocks closed via the algebraic parts M_d = A*B (kernels with T_1 + 1 = square). Both parts probable prime
#   = complete factorization into probable primes => class "voll", route "teile_voll"; one probable-prime part with exponent 1
#   => NEW class "prpzeuge" (witness that is a probable prime), route "teile_zeuge". E3 above deliberately checks the state
#   BEFORE w158 (Prop. 6.10 (a) as published).
w158 = J('w158_teile_wahrscheinlich_prim_result.json')
ext = {tuple(x): ('voll', 'teile_voll') for x in w158['klassen']['beide_prp']}
ext.update({tuple(x): ('prpzeuge', 'teile_zeuge') for x in w158['klassen']['zeuge_prp']})
for md, (kl_, ro_) in ext.items():
    assert md in menge and klasse[md] == 'offen', ('w158 schliesst einen Block, der nicht offen war', md, klasse.get(md))
    klasse[md] = kl_; route[md] = ro_
print(f'PK E6 ✅  w158 schliesst {sorted(ext)} als {[ext[x][0] for x in sorted(ext)]} (vorher offen)')
# Primality proofs: w162 lists every factor of the building blocks of the classes voll / prpzeuge; w163 (Pocklington, own
#   checker) and w164 (ECPP from PARI, own checker) prove those that are only probable primes. A block becomes "zert" when EVERY
#   one of its factors is proven (deterministically below psi_12, otherwise by certificate).
# `bewiesen` = numbers proven by w163 or w164; `je_block` = per block the list of flags "factor proven"; `hoch` = blocks
#   upgraded to zert.
w162 = J('w162_inventar_wahrscheinliche_primzahlen_result.json')
bewiesen = {e['zahl'] for e in J('w163_pocklington_beweise_result.json')['ergebnisse'] if e['bewiesen']} \
         | {e['zahl'] for e in J('w164_ecpp_pruefer_result.json')['ergebnisse'] if e['bewiesen']}
je_block = {}
for f in w162['inventar']:
    ok = f['status'] == 'bewiesen_psi12' or f['zahl'] in bewiesen
    je_block.setdefault((f['m'], f['d']), []).append(ok)
hoch = sorted(md for md, oks in je_block.items() if all(oks))
for md in hoch:
    assert klasse[md] in ('voll', 'prpzeuge'), ('w162-Baustein nicht in voll/prpzeuge', md, klasse[md])
    klasse[md] = 'zert'
print(f'PK E7 {"✅" if len(hoch) == len(je_block) == 41 else "⚠️"}  Primbeweise: {len(hoch)} von {len(je_block)} Bausteinen der Klassen voll/prpzeuge '
      f'haben jetzt nur bewiesene Faktoren ⇒ zert (Faktoren: {len(w162["inventar"])}, davon per Zertifikat {len(bewiesen)})')
# `z_neu` = count of routes among the 150 new pairs.
z_neu = Counter(route[md] for md in neu)

# `L` = log10 of the fundamental unit per kernel; `lg` = log10 D_d per pair; `ck` = checkpoint of w135 with the digit counts
#   (`stellen`)
# of its blocks, used for E4.
L = {m: lg_eps(m) for m in RANG}
lg = {md: lg_Phi(L[md[0]], md[1]) for md in menge}
ck = J('w135_stufe2_ecm_checkpoint.json')
falsch = [(e['m'], e['d'], e['stellen'], round(lg[(e['m'], e['d'])], 3)) for e in ck.values()
          if math.floor(lg[(e['m'], e['d'])]) + 1 != e['stellen']]
print(f'PK E4 {"✅" if not falsch else "❌"}  log10 D_d stellengenau fuer {len(ck)} Bloecke aus w135' + (f' — FALSCH: {falsch}' if falsch else '.'))
assert not falsch

w137 = J('w137_wieferich_79kerne_result.json')
# `wief` = Wieferich events (m, d, p) of the question population.
wief = [(e['m'], e['d'], e['p']) for e in w137['ereignisse'] if e['in_frage_population']]
e5 = all((m, d) in menge for m, d, _ in wief)
print(f'PK E5 {"✅" if e5 else "❌"}  {len(wief)} Wieferich-Ereignisse in der Frage-Population, alle in der 961-Menge: {wief}')
assert e5

# `alt_r` = rows of w92 by (m, d); `raenge` = final rows (rank d of each kernel m), sorted by (d, m).
alt_r = {(r['m'], r['d']): r for r in w92['raenge']}
raenge = []
for m, d in sorted(menge, key=lambda x: (x[1], x[0])):
    if (m, d) in alt_r:
        r = dict(alt_r[(m, d)]); r['kern'] = 'alt'
        if (m, d) in ext: r['route'] = route[(m, d)]   # the route from w158 instead of "offen" from w92
    else:
        r = dict(m=m, d=d, lg=round(lg[(m, d)], 2), route=route[(m, d)], kern='neu')
        b = bew.get((m, d)) or prob.get((m, d))
        if b: r['zeuge_stellen'] = len(b[2])
        if (m, d) in bew: r['beweis'] = bew[(m, d)][4]
    r['klasse'] = klasse[(m, d)]
    w = [p for mm, dd, p in wief if (mm, dd) == (m, d)]
    if w: r['wieferich_p'] = w[0]
    raenge.append(r)
zahl = Counter(r['route'] for r in raenge); kl = Counter(r['klasse'] for r in raenge)
assert sum(kl.values()) == 961
erl = kl['zert'] + kl['voll'] + kl['prpzeuge']
print(f'          961: zert {kl["zert"]} + voll {kl["voll"]} + prpzeuge {kl["prpzeuge"]} = erledigt {erl} · offen {kl["offen"]} · Routen {dict(zahl)}')
out = dict(skript=pathlib.Path(__file__).name, datum='2026-10-02', n=len(raenge), kerne=len(RANG), routen=dict(zahl),
           zert=kl['zert'], voll=kl['voll'], prpzeuge=kl['prpzeuge'], erledigt=erl, offen=kl['offen'],
           # pairs upgraded to zert by primality proofs
           zert_durch_primbeweis=len(hoch), primbeweis_bausteine=[list(md) for md in hoch],
           # class counts of the 811 / 150 after all updates
           alt=dict(n=811, **dict(Counter(klasse[md] for md in alt))), neu=dict(n=150, **dict(Counter(klasse[md] for md in neu))),
           aufwand='die 811 zusaetzlich GMP-ECM bis t40 (w89); die 150 bis t35 (w135)', raenge=raenge)
(ERG/'w138_karte_961_daten.json').write_text(json.dumps(out, ensure_ascii=False), encoding='utf-8')
print('Ergebnis: w138_karte_961_daten.json')
