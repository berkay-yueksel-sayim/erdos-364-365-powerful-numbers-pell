# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w91_wieferich_landschaft_2026-09-24.py
# Purpose: data for a coordinate plot of Wieferich quotients (the "Wieferich landscape") for the 25 kernels of w67, and a
#   genuine test of uniformity: is the Wieferich quotient equidistributed in [0, 1)?
# The quantity: for every prime p not dividing 2·m·U₁, and the unit ε = T₁ + U₁√m (norm 1), ε^(p−e) ≡ 1 (mod p) with e = (m/p).
#   Write ε^(p−e) = A + B√m (mod p²). Then p | B, and the WIEFERICH QUOTIENT is
#       q_p = (B/p) mod p ∈ {0, …, p−1},    y_p = q_p / p ∈ [0, 1).
#   p is Lucas–Wieferich (McIntosh–Roettger 2007, p. 2088; Prop. 6.11) exactly when q_p = 0. The usual heuristic treats y_p as
#   uniform in [0, 1), which gives the expectation ~1/p per prime. So hidden structure would have to show up as non-uniformity
#   of y_p; this is what the script measures.
# Reads:  ergebnisse/w67_beweismenge_result.json (key `kerne`, the 25 kernels) and
#         ergebnisse/w77_wieferich_ereignisse_result.json (key `ereignisse`, the 16 known Wieferich events).
# Writes: ergebnisse/w91_wieferich_landschaft_result.json, ..._daten.json (plot data and histograms) and ..._output.txt.
# Usage:  python w91_wieferich_landschaft_2026-09-24.py [P_MAX=3000000]   (P_MAX = bound for the primes p)
# Controls and expectations (stated before the run):
#   E1  structure control: for EVERY prime, B ≡ 0 and A ≡ 1 (mod p); otherwise the script computes wrongly (asserted).
#   E2  control against the record: all 16 Wieferich events of w77 (p < 3·10⁶, 25 kernels) appear as q_p = 0 (asserted).
#   E3  there are MORE Lucas–Wieferich primes than 16: the events of w77 count only primes dividing some T_k (a rank exists),
#       Lucas–Wieferich counts all primes. The difference is counted, not conjectured.
#   E4  uniformity: chi-square test with 20 classes, per kernel and overall, separately for p split (e = +1) and inert (e = −1).
#       Expectation: no significant deviation (the p-values over the 25 tests roughly uniform). A systematic deviation would
#       be a real finding and would be checked twice before any interpretation.
#   E5  number of hits against the heuristic Σ 1/p over the primes tested.

import sys, json, time, pathlib, math
from math import isqrt
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent
ERG  = HIER.parent/'ergebnisse'
P_MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3_000_000   # `P_MAX` = bound for the primes p
# `DUENNUNG` = thinning for the plot: every 20th prime per kernel is kept (hits: ALL); histograms use ALL primes
DUENNUNG = 20
aus = []   # `aus` = output lines, written to the output file at the end
def sag(s=''):   # `sag` = "say": print a line and keep it for the output file
    print(s, flush=True); aus.append(s)
def fund(m):   # fundamental solution (h0, k0) = (T1, U1) of x² − m·y² = 1, by the continued fraction of √m
    a0 = isqrt(m); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (m - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
def eps_pow(T1, U1, m, e, M):   # (A, B) with ε^e = A + B√m (mod M), ε = T1 + U1·√m, by repeated squaring
    ra, rb, ba, bb = 1 % M, 0, T1 % M, U1 % M
    while e:
        if e & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        e >>= 1
    return ra, rb
def primzahlen_bis(N):   # `primzahlen_bis` = primes up to N: the odd primes <= N, by a sieve
    s = bytearray([1])*(N+1); s[0] = s[1] = 0
    for i in range(2, isqrt(N)+1):
        if s[i]: s[i*i::i] = bytearray(len(range(i*i, N+1, i)))
    return [i for i in range(3, N+1) if s[i]]
def chi2_p(zaehl, k):
    """p-value of the chi-square test for uniformity (k classes), Wilson–Hilferty approximation; returns (chi-square, p-value)."""
    n = sum(zaehl); erw = n / k
    x2 = sum((z - erw)**2 / erw for z in zaehl); df = k - 1
    z = ((x2/df)**(1/3) - (1 - 2/(9*df))) / math.sqrt(2/(9*df))
    return x2, 0.5 * math.erfc(z / math.sqrt(2))

t0 = time.time()
sag('='*100)
sag(f'w91 — WIEFERICH-LANDSCHAFT der 25 Kerne, Primzahlen bis {P_MAX:,}   {time.strftime("%Y-%m-%d %H:%M")}')
sag('='*100)
KERNE = json.load(open(ERG/'w67_beweismenge_result.json', encoding='utf-8'))['kerne']   # `KERNE` = the 25 kernels
assert len(KERNE) == 25
PR = primzahlen_bis(P_MAX)   # `PR` = odd primes up to P_MAX
sag(f'{len(PR):,} ungerade Primzahlen, 25 Kerne ⇒ ~{len(PR)*25/1e6:.1f} Mio. Quotienten.')

w77 = json.load(open(ERG/'w77_wieferich_ereignisse_result.json', encoding='utf-8'))['ereignisse']
ereig = {(x['m'], x['p']) for x in w77}   # `ereig` = the known Wieferich events as pairs (kernel m, prime p)
assert len(ereig) == 16

K = 20   # `K` = number of classes (histogram bins for y_p)
# `daten` = per-kernel data; `treffer` = hits (q_p = 0); `gesamt_split`, `gesamt_traege` = overall histograms for split
# (e = +1) and inert (e = −1) primes; `erw_summe` = heuristic expectation Σ 1/p; `e1_verletzt` = number of violations of E1
daten, treffer, gesamt_split, gesamt_traege, erw_summe, e1_verletzt = {}, [], [0]*K, [0]*K, 0.0, 0
for m in KERNE:
    T1, U1 = fund(m); tm = time.time()
    # `hist` = histograms of y_p by e; `punkte` = plot points [p, y, e]; `n_m` = number of primes tested for this kernel
    hist = {+1: [0]*K, -1: [0]*K}; punkte = []; n_m = 0
    for i, p in enumerate(PR):
        if m % p == 0 or U1 % p == 0: continue
        e = 1 if pow(m % p, (p-1)//2, p) == 1 else -1
        A, B = eps_pow(T1, U1, m, p - e, p*p)
        if A % p != 1 or B % p != 0: e1_verletzt += 1; continue
        q = (B // p) % p; y = q / p
        hist[e][min(K-1, int(y*K))] += 1; n_m += 1; erw_summe += 1/p
        if q == 0:
            treffer.append(dict(m=m, p=p, e=e, im_record_w77=(m, p) in ereig))
        if i % DUENNUNG == 0 or q == 0:
            punkte.append([p, round(y, 4), e])
    for k_ in range(K): gesamt_split[k_] += hist[1][k_]; gesamt_traege[k_] += hist[-1][k_]
    x2s, ps = chi2_p(hist[1], K); x2t, pt = chi2_p(hist[-1], K)
    daten[m] = dict(n=n_m, hist_zerfallend=hist[1], hist_traege=hist[-1], chi2_p_zerfallend=round(ps, 4),
                    chi2_p_traege=round(pt, 4), punkte=punkte)
    sag(f'  m = {m:>8}: {n_m:>7,} Primzahlen · Treffer {sum(1 for t in treffer if t["m"] == m):>2} · '
        f'χ²-p zerfallend {ps:.3f}, traege {pt:.3f} · {time.time()-tm:.0f}s')

sag()
assert e1_verletzt == 0, ('E1 VERLETZT: A ≢ 1 oder B ≢ 0 (mod p) in', e1_verletzt, 'Faellen')
sag(f'PK E1 ✅  fuer alle {sum(d["n"] for d in daten.values()):,} Primzahlen: ε^(p−e) ≡ 1 (mod p). Die Arithmetik stimmt.')
im_record = [t for t in treffer if t['im_record_w77']]   # `im_record` = hits that are known events of w77
fehlend = ereig - {(t['m'], t['p']) for t in treffer}   # `fehlend` = known events that were not found (must be empty)
assert not fehlend, ('E2 VERLETZT: Record-Ereignisse fehlen', fehlend)
sag(f'PK E2 ✅  alle 16 Wieferich-Ereignisse aus w77 erscheinen als q_p = 0.')
sag(f'E3      Lucas-Wieferich-Primzahlen insgesamt: {len(treffer)} — davon {len(im_record)} im Record w77, '
    f'{len(treffer)-len(im_record)} ZUSAETZLICH (Primzahlen, die kein T_k teilen).')
for t in treffer:
    if not t['im_record_w77']: sag(f'          zusaetzlich: m = {t["m"]}, p = {t["p"]}, e = {t["e"]:+d}')
x2s, ps = chi2_p(gesamt_split, K); x2t, pt = chi2_p(gesamt_traege, K)
ps_liste = [d['chi2_p_zerfallend'] for d in daten.values()] + [d['chi2_p_traege'] for d in daten.values()]
unter5 = sum(1 for x in ps_liste if x < 0.05)   # `unter5` = number of single tests with p-value below 0.05
sag(f'E4      Gleichverteilung gesamt: zerfallend χ² = {x2s:.1f} (p = {ps:.3f}), traege χ² = {x2t:.1f} (p = {pt:.3f}).')
sag(f'        Einzeltests: {unter5} von {len(ps_liste)} unter 0,05 — erwartet bei reiner Zufaelligkeit: ~{0.05*len(ps_liste):.1f}.')
sag(f'E5      Treffer {len(treffer)} gegen die Heuristik Σ 1/p = {erw_summe:.1f}.')
sag()
sag('='*100)
# `urteil` = verdict text (stored in the result file)
urteil = ('KEINE Abweichung von der Gleichverteilung sichtbar — die Wolke verhaelt sich, wie die Heuristik es sagt.'
          if ps > 0.01 and pt > 0.01 and unter5 <= max(6, 0.05*len(ps_liste)*3) else
          'AUFFAELLIG — vor jeder Deutung zweimal pruefen (andere Klassenzahl, anderer Test, anderer Primzahlbereich).')
sag(f'  {urteil}')
sag('='*100)
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), p_max=P_MAX, klassen=K,
           duennung=DUENNUNG, kerne=KERNE, treffer=treffer, n_treffer=len(treffer), n_im_record=len(im_record),
           heuristik_summe=round(erw_summe, 3), chi2_gesamt=dict(zerfallend=[round(x2s, 2), round(ps, 4)],
           traege=[round(x2t, 2), round(pt, 4)]), einzeltests_unter_5prozent=unter5, n_einzeltests=len(ps_liste),
           urteil=urteil, laufzeit_s=round(time.time()-t0))
(ERG/'w91_wieferich_landschaft_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w91_wieferich_landschaft_daten.json').write_text(json.dumps({str(k): v for k, v in daten.items()}, ensure_ascii=False), encoding='utf-8')
(ERG/'w91_wieferich_landschaft_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
print('\nErgebnis: w91_wieferich_landschaft_result.json + _daten.json')
