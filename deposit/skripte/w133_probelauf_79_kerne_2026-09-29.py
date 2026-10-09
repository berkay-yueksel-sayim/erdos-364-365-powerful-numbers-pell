# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w133_probelauf_79_kerne_2026-09-29.py
# Trial run: the FIRST SEARCH (witness <= 3·10⁶) on all 961 pairs (m, d) of the 79 kernels.
#
# Reads:    in ../ergebnisse: w9_v31_M100000000_H2000_result.json (points of the box run), w92_karte_811_daten.json (routes of
#           the 811 pairs of the 25 kernels with the most points), and w131_wieferich_zweite_route_result.json.
# Writes:   ../ergebnisse/w133_probelauf_79_kerne_output.txt and ../ergebnisse/w133_probelauf_79_kerne_result.json.
# Usage:    python w133_probelauf_79_kerne_2026-09-29.py   (no arguments); exit code 0 if all controls pass, otherwise 1.
# Controls: PK (doubling = recurrence; all pairs of route `erste` of w92 refound), NK (no witness for the open pairs); see below.
#
# Decision rule fixed before the run: extend the search to the 54 further kernels if there is no structural blocker and
#   the share of the 150 new pairs with a witness <= 3·10⁶ is at least 37 % (half of 606/811 = 74.7 %). This script
#   delivers NUMBERS, no decision.
# WITNESS for (m, d): a prime q <= 3·10⁶ of rank exactly d in T_k(m) with v_q(T_d) = 1 (then v_q(Φ_d) = 1
#   and Φ_d is not powerful).
# Own code: T₁ from the continued fraction of √m; V_n = 2 T_n as a Lucas sequence (P, Q) = (2 T₁, 1), evaluated by
#   doubling modulo q; q | T_d ⇔ q | V_d (q odd); rank exactly d ⇔ in addition q ∤ V_{d/r} for every prime r | d (d
#   odd); v_q = 1 ⇔ q² ∤ V_d.
# The search runs over the classes q ≡ ±1 (mod 4d) (Corollary 4.8(i)); a witness that is found is checked directly anyway. If the
#   congruence were wrong, witnesses would be missing; the PK catches this: all 606 pairs of route `erste` (w92) must be refound.
# NK: the 7 open pairs (searched far beyond 3·10⁶) must not get a witness.
# NOT done here: the heavier routes (sieve up to 10⁸/10⁹, ECM, factorization) for the new pairs; that
#   would be the extension itself.
import sys, json, math, pathlib, time
from collections import Counter, defaultdict
sys.stdout.reconfigure(encoding='utf-8')
ERG = pathlib.Path(__file__).resolve().parents[1] / 'ergebnisse'   # `ERG` = results directory
t0 = time.time(); aus = []   # `aus` = lines of the output file
def sag(s=''): print(s, flush=True); aus.append(s)   # `sag` = print a line and record it for the output file

# returns T_1 of the fundamental solution of x^2 - m y^2 = 1 (continued fraction of sqrt(m)); m must not be a square
def pell(m):
    a0 = math.isqrt(m); assert a0 * a0 != m
    mm, dd, a = 0, 1, a0; h1, h = 1, a0; k1, k = 0, 1
    while h * h - m * k * k != 1:
        mm = dd * a - mm; dd = (m - mm * mm) // dd; a = (a0 + mm) // dd; h1, h = h, a * h + h1; k1, k = k, a * k + k1
    return h

def V(n, P, mod):
    # Lucas V_n(P, 1) modulo `mod`, doubling from the left over the binary digits of n; the pair (V_j, V_{j+1}) is kept
    a, b = 2 % mod, P % mod
    for bit in bin(n)[2:]:
        if bit == '1': a, b = (a * b - P) % mod, (b * b - 2) % mod
        else: a, b = (a * a - 2) % mod, (a * b - P) % mod
    return a

def primteiler(n):   # `primteiler` = set of the prime divisors of n (trial division)
    r, q = set(), 2
    while q * q <= n:
        while n % q == 0: r.add(q); n //= q
        q += 1
    if n > 1: r.add(n)
    return r

N = 3_000_000   # search bound for the witness q
sieb = bytearray([1]) * (N + 1); sieb[0] = sieb[1] = 0   # `sieb` = sieve of Eratosthenes up to N (1 = prime)
for i in range(2, math.isqrt(N) + 1):
    if sieb[i]: sieb[i * i::i] = bytearray(len(range(i * i, N + 1, i)))

def zeuge(T1, d):
    # `zeuge` = witness search for a pair (m, d) whose fundamental solution has T_1 = `T1`: returns the smallest prime q <= N
    #   of rank exactly d with q^2 not dividing V_d (q in the classes +-1 mod 4d), or None.
    #   `pr` = prime divisors of d, `M` = 4d, `kand` = candidates q.
    pr = primteiler(d); M = 4 * d
    kand = sorted(q for r in (1, M - 1) for q in range(r, N + 1, M) if q > 2 and sieb[q])
    for q in kand:
        P = 2 * (T1 % q)
        if V(d, P, q) != 0: continue
        if any(V(d // r, P, q) == 0 for r in pr): continue          # rank smaller than d
        q2 = q * q
        if V(d, 2 * (T1 % q2), q2) == 0: continue                    # q² | V_d (Wieferich-type case): not a witness
        return q
    return None

# Check of the doubling against the recurrence (small): T_k(7) exact, V_k = 2 T_k
T1_7 = pell(7); t = [1, T1_7]
for _ in range(20): t.append(2 * T1_7 * t[-1] - t[-2])
pk_v = all(V(k, 2 * T1_7, 10 ** 30) == (2 * t[k]) % 10 ** 30 for k in range(22))

pts = json.loads((ERG / 'w9_v31_M100000000_H2000_result.json').read_text(encoding='utf-8'))['points']   # `pts` = points (m, m', k, digits, status, witness)
# `cnt` = number of points per kernel; `rang` = kernels sorted by that number (descending), then by m; `K25` = the
#   25 kernels with most points
cnt = Counter(p[0] for p in pts); rang = sorted(cnt, key=lambda m: (-cnt[m], m)); K25 = set(rang[:25])
paare = defaultdict(set)   # `paare` = pairs: kernel m -> set of the odd divisors d >= 5 of the k values of its points
for m, mp, k, st, status, z in pts:
    for d in range(5, k + 1, 2):
        if k % d == 0: paare[m].add(d)
w92 = json.loads((ERG / 'w92_karte_811_daten.json').read_text(encoding='utf-8'))
# `route` = route by which each of the 811 pairs was settled in w92
route = {(r['m'], r['d']): r['route'] for r in w92['raenge']}
sag('=' * 100); sag(f'w133 — PROBELAUF ERSTE SUCHE (Zeuge ≤ 3·10⁶) AUF ALLEN 79 KERNEN   {time.strftime("%Y-%m-%d %H:%M")}'); sag('=' * 100)
sag(f'PK Verdopplung = Rekursion (m = 7, k ≤ 21): {"✅" if pk_v else "❌"}')

T1 = {}; blocker = []   # `T1` = fundamental T_1 per kernel; `blocker` = kernels for which it could not be computed
for m in rang:
    try: T1[m] = pell(m)
    except Exception as e: blocker.append((m, str(e)))
sag(f'T₁ berechnet fuer {len(T1)} von {len(rang)} Kernen; groesstes T₁: {max(len(str(v)) for v in T1.values())} Stellen; Blocker: {blocker or "keine"}')
ungerade = sorted(m for m in T1 if T1[m] % 2 == 1)   # `ungerade` = kernels with T_1 odd
sag(f'Kerne mit T₁ ungerade: in den 25 {sum(1 for m in ungerade if m in K25)}, in den uebrigen 54 {sum(1 for m in ungerade if m not in K25)}')

erg = {}   # `erg` = (m, d) -> witness q or None
for m in rang:
    for d in sorted(paare[m]):
        erg[(m, d)] = zeuge(T1[m], d)
# `alt` = pairs of the 25 kernels in `K25` (old), `neu` = pairs of the other 54 kernels (new)
alt = [md for md in erg if md[0] in K25]; neu = [md for md in erg if md[0] not in K25]
assert len(alt) == 811 and set(alt) == set(route), ('Paarmenge der 25 ≠ w92', len(alt))
# `erste` = old pairs on route `erste` in w92, `offen` = old pairs that are open in w92 (route `offen`)
erste = [md for md in alt if route[md] == 'erste']; offen = [md for md in alt if route[md] == 'offen']
pk_erste = sum(1 for md in erste if erg[md]); nk_offen = [md for md in offen if erg[md]]   # refound pairs; open pairs that wrongly got a witness
andere = Counter(route[md] for md in alt if erg[md] and route[md] != 'erste')   # witnesses on pairs of other routes
# Route `erste` in w92 also covers small blocks that were COMPLETELY FACTORED, so a pair can be on this route without having a
#   witness <= 3·10⁶: (319, 5) has none, its witnesses are 1375502561 and 158585313121 (w131). A missing pair therefore counts
#   as explained only if the INDEPENDENT route w131 documents "no witness below 3·10⁶" for exactly this pair; no exception is
#   hard-coded.
w131 = json.loads((ERG / 'w131_wieferich_zweite_route_result.json').read_text(encoding='utf-8'))
belegt_ohne = {(b['m'], b['d']) for b in w131['bloecke_reichweite'] if not b['zeugen_unter_3e6']}   # pairs w131 documents without witness
fehlt = [md for md in erste if not erg[md]]; unerklaert = [md for md in fehlt if md not in belegt_ohne]   # missing / unexplained pairs
sag(f'PK: Route „erste" (w92) wiedergefunden: {pk_erste} von {len(erste)} · fehlend {fehlt}, davon durch w131 als „kein Zeuge ≤ 3·10⁶" belegt: '
    f'{[md for md in fehlt if md in belegt_ohne]} · unerklaert: {unerklaert or "keine"} {"✅" if not unerklaert else "❌"}')
pk_erste_ok = not unerklaert
sag(f'NK: Zeugen fuer die {len(offen)} offenen: {nk_offen or "keine"} {"✅" if not nk_offen else "❌ — PRUEFEN"}')
sag(f'Kalibrierung: auf den 811 hier mit Zeugen ≤ 3·10⁶ insgesamt {sum(1 for md in alt if erg[md])} '
    f'(zusaetzlich zu „erste": {dict(andere) or "keine"} — dort fand die erste Suche keinen, diese Suche schon)')
n_neu = sum(1 for md in neu if erg[md]); anteil = n_neu / len(neu)   # `anteil` = share of the new pairs with a witness
sag(f'\nNEUE PAARE (54 Kerne): {len(neu)} · mit Zeugen ≤ 3·10⁶: {n_neu} = {100 * anteil:.1f} % · ohne: {len(neu) - n_neu}')
sag(f'Vergleich: die 811 in der ersten Suche 606 = 74,7 %; Schwelle aus dem Kriterium (vorher festgelegt): 37 %')
ohne = sorted(md for md in neu if not erg[md])   # `ohne` = new pairs without a witness
nach_kern = Counter(m for m, d in ohne)   # `nach_kern` = number of pairs without a witness per kernel
sag(f'ohne Zeugen, je Kern (m: Anzahl): {dict(sorted(nach_kern.items()))}')
sag(f'davon auf Kernen mit T₁ ungerade: {sum(1 for m, d in ohne if T1[m] % 2 == 1)}')
sag(f'Laufzeit {time.time() - t0:.0f} s')
sag('Hinweis: diese Suche ist STRENGER als „erste" (keine Zerlegung kleiner Bloecke) — der Anteil der neuen Paare ist eher unter- als ueberschaetzt.')
ok = pk_v and pk_erste_ok and not nk_offen and not blocker   # `ok` = all controls passed
(ERG / 'w133_probelauf_79_kerne_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
(ERG / 'w133_probelauf_79_kerne_result.json').write_text(json.dumps(dict(
    skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d'), schwelle_prozent=37, vergleich_811=[606, 811],
    neu_paare=len(neu), neu_mit_zeuge=n_neu, neu_anteil=round(anteil, 4), neu_ohne=[list(x) for x in ohne],
    pk_erste=[pk_erste, len(erste)], nk_offen=[list(x) for x in nk_offen], zusaetzlich_auf_811=dict(andere),
    kerne_T1_ungerade=ungerade, blocker=blocker, zeugen={f'{m}_{d}': q for (m, d), q in erg.items() if q}, alle_kontrollen_ok=ok),
    ensure_ascii=False, indent=1), encoding='utf-8')
sys.exit(0 if ok else 1)
