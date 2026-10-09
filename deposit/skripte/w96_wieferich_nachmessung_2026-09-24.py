# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w96_wieferich_nachmessung_2026-09-24.py
# -*- coding: utf-8 -*-
# Three follow-up measurements to w91 (Wieferich landscape), in one run:
#   (a) fine structure near q = 0 ("near hits"), which the 20 classes of w91 cannot see;
#   (b) dependence on the order of ε modulo p;
#   (c) cross-correlation between kernels (as in w94) at ALL primes instead of every 20th prime.
#
# Reads:    ../ergebnisse/w67_beweismenge_result.json (key `kerne`, the kernels), w91_wieferich_landschaft_daten.json and
#           w91_wieferich_landschaft_result.json (same directory).
# Writes:   ../ergebnisse/w96_wieferich_nachmessung_result.json and ../ergebnisse/w96_wieferich_nachmessung_output.txt
# Usage:    python w96_wieferich_nachmessung_2026-09-24.py   (no arguments; uses numpy)
# Controls: PK1 to PK4 (below); each failed control aborts the run with an assertion.
#
# Quantity as in w91: ε^(p−e) = A + B√m (mod p²), q_p = (B/p) mod p, y_p = q_p/p. In addition, for each prime the 2-adic valuation
# s = v₂(ord(ε mod p)): ε^(odd part of n) is squared until it equals 1 (n = p − e). Classes: s = 0 (odd order),
# s = 1 (order ≡ 2 mod 4), s ≥ 2 (order ≡ 0 mod 4 ⟺ p divides some T_k).
#
# EXPECTATIONS (set before the run):
#   PK1  A ≡ 1, B ≡ 0 (mod p) for every prime (as in w91).
#   PK2  REPRODUCTION: per kernel the 20-class histograms of w91 (split/inert) EXACTLY, and the 56 hits as a set.
#   PK3  the 16 record hits lie in s ≥ 2, the 40 additional hits in s ≤ 1 (measured one by one beforehand).
#   PK4  (c) positive control: an artificially planted correlation (2 % of the values of kernel 2 replaced by those of kernel 1,
#        r ≈ 0.02) is detected with |z| > 5; without planting the same pair is unremarkable.
#   (a)  For thresholds t = 10⁻², …, 10⁻⁷: observed number of y < t against the EXACT expectation Σ ⌈t·p⌉/p
#        (q uniform on 0..p−1). Expectation: all |z| < 3. For t < 1/(3·10⁶) the count equals the number of hits (56).
#        In addition the lowest percent in 10 subclasses, χ² against the exact discrete expectation.
#   (b)  Per class s (0 / 1 / ≥2) and per e (±1): χ² uniformity test (20 classes) and hits against Σ 1/p.
#        Expectation: no deviation (6 tests; remarkable only for p < 0.01/6 after Bonferroni).
#   (c)  300 kernel pairs at all common primes: |z| > 2.576 for about 3 pairs, maximum around 3, Bonferroni p > 0.01.
#        Detection limit at n ≈ 216 000: |r| ≈ 4/√n ≈ 0.009.
# A deviation would be a finding, and would be checked a second time (on a different range) before any interpretation.

import sys, json, time, pathlib, math, itertools
from math import isqrt
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent   # `HIER` = this directory
ERG  = HIER.parent/'ergebnisse'                  # `ERG` = results directory
P_MAX = 3_000_000                                # primes up to P_MAX are used
aus = []   # `aus` = lines of the output file
def sag(s=''):   # `sag` = print a line and record it for the output file
    print(s, flush=True); aus.append(s)
def fund(m):   # `fund` = fundamental solution (T_1, U_1) of x^2 - m y^2 = 1 from the continued fraction of sqrt(m)
    a0 = isqrt(m); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (m - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
def eps_pow(T1, U1, m, e, M):   # (a, b) with ε^e = a + b√m modulo M, by square-and-multiply
    ra, rb, ba, bb = 1 % M, 0, T1 % M, U1 % M
    while e:
        if e & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        e >>= 1
    return ra, rb
def primzahlen_bis(N):   # `primzahlen_bis` = all odd primes <= N (sieve of Eratosthenes)
    s = bytearray([1])*(N+1); s[0] = s[1] = 0
    for i in range(2, isqrt(N)+1):
        if s[i]: s[i*i::i] = bytearray(len(range(i*i, N+1, i)))
    return [i for i in range(3, N+1) if s[i]]
# chi-square statistic of observed vs expected counts, and its upper-tail p-value (Wilson-Hilferty approximation)
def chi2_p(beob, erw):
    x2 = float(sum((b - e)**2 / e for b, e in zip(beob, erw))); df = len(beob) - 1
    z = ((x2/df)**(1/3) - (1 - 2/(9*df))) / math.sqrt(2/(9*df))
    return x2, 0.5 * math.erfc(z / math.sqrt(2))

t0 = time.time()
sag('='*100); sag(f'w96 — WIEFERICH-NACHMESSUNGEN (a) Beinahe-Treffer · (b) Ordnung · (c) Kreuzkorrelation, p ≤ {P_MAX:,}   {time.strftime("%Y-%m-%d %H:%M")}'); sag('='*100)
KERNE = json.load(open(ERG/'w67_beweismenge_result.json', encoding='utf-8'))['kerne']   # `KERNE` = the kernels m
d91 = json.load(open(ERG/'w91_wieferich_landschaft_daten.json', encoding='utf-8'))   # `d91` = w91 histograms per kernel
r91 = json.load(open(ERG/'w91_wieferich_landschaft_result.json', encoding='utf-8'))   # `r91` = w91 result (list of hits)
PR = primzahlen_bis(P_MAX); NP = len(PR); idx = {p: i for i, p in enumerate(PR)}   # `PR` = odd primes, `idx` = index of each prime
# arrays over (kernel, prime): `Y` = y_p = q_p/p (NaN where p | m or p | U_1), `Qm` = q_p, `E` = e = (m/p) as
#   +-1, `S` = s = v_2(order)
Y = np.full((len(KERNE), NP), np.nan); Qm = np.zeros((len(KERNE), NP), dtype=np.int64)
E = np.zeros((len(KERNE), NP), dtype=np.int8); S = np.full((len(KERNE), NP), -1, dtype=np.int8)
pk1, treffer = 0, set()   # `pk1` = number of violations of PK1; `treffer` = hits: set of (m, p) with q = 0
for j, m in enumerate(KERNE):
    T1, U1 = fund(m); tm = time.time()
    for i, p in enumerate(PR):
        if m % p == 0 or U1 % p == 0: continue
        e = 1 if pow(m % p, (p-1)//2, p) == 1 else -1
        A, B = eps_pow(T1, U1, m, p - e, p*p)
        if A % p != 1 or B % p != 0: pk1 += 1; continue
        q = (B // p) % p
        n = p - e; v = (n & -n).bit_length() - 1
        a, b = eps_pow(T1, U1, m, n >> v, p); s = 0
        while (a, b) != (1 % p, 0):
            a, b = (a*a + m*b*b) % p, (2*a*b) % p; s += 1
            assert s <= v, ('Ordnung passt nicht zu p − e', m, p)
        Y[j, i] = q / p; Qm[j, i] = q; E[j, i] = e; S[j, i] = s
        if q == 0: treffer.add((m, p))
    sag(f'  m = {m:>8}: {int(np.sum(~np.isnan(Y[j]))):>7,} Primzahlen · {time.time()-tm:.0f}s')
sag()
assert pk1 == 0, ('PK1 VERLETZT', pk1)
sag(f'PK1 ✅  ε^(p−e) ≡ 1 (mod p) fuer alle {int(np.sum(~np.isnan(Y))):,} Paare (m, p).')
# PK2: reproduce w91 exactly
for j, m in enumerate(KERNE):
    ok = ~np.isnan(Y[j])
    for e, key in ((1, 'hist_zerfallend'), (-1, 'hist_traege')):   # keys: split (e = +1) and inert (e = -1) histograms
        h = np.bincount(np.minimum(19, (Y[j][ok & (E[j] == e)] * 20).astype(int)), minlength=20).tolist()
        assert h == d91[str(m)][key], ('PK2 VERLETZT: Histogramm', m, key)
t91 = {(t['m'], t['p']) for t in r91['treffer']}   # `t91` = the hits of w91 as a set
assert treffer == t91 and len(treffer) == 56, ('PK2 VERLETZT: Treffer', len(treffer))
sag('PK2 ✅  alle 50 Histogramme von w91 Klasse fuer Klasse gleich; dieselben 56 Treffer.')
rec = {(t['m'], t['p']) for t in r91['treffer'] if t['im_record_w77']}   # `rec` = the record hits (flag `im_record_w77` in w91)
for (m, p) in treffer:
    s = int(S[KERNE.index(m), idx[p]])
    assert (s >= 2) == ((m, p) in rec), ('PK3 VERLETZT', m, p, s)
sag('PK3 ✅  die 16 Record-Treffer haben Ordnung ≡ 0 (mod 4), die 40 zusaetzlichen nicht.')

# flattened views over all valid (kernel, prime) pairs: `Pv` = p, `Yv` = y, `Qv` = q, `Ev` = e, `Sv` = s
OK = ~np.isnan(Y)
Pgrid = np.broadcast_to(np.array(PR, dtype=np.float64), Y.shape)
Pv, Yv, Qv, Ev, Sv = Pgrid[OK], Y[OK], Qm[OK], E[OK], S[OK]
N_ALL = Yv.size   # `N_ALL` = total number of valid (m, p) pairs

# ---------------- (a) near hits
sag(); sag('(a) BEINAHE-TREFFER — Anzahl y = q/p < t gegen die exakte Erwartung Σ ⌈t·p⌉/p')
ergebnis_a = []   # `ergebnis_a` = rows (threshold, observed, expected, z)
for t in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
    pi = np.ceil(t * Pv) / Pv                       # P(q < t·p) = #{0 ≤ q < t·p}/p = ⌈t·p⌉/p
    erw = float(pi.sum()); var = float((pi * (1 - pi)).sum())   # `erw` = expected count, `var` = its variance
    beob = int(np.sum(Qv < t * Pv))                 # `beob` = observed count
    z = (beob - erw) / math.sqrt(var)
    ergebnis_a.append(dict(t=t, beobachtet=beob, erwartet=round(erw, 2), z=round(z, 2)))
    sag(f'    t = {t:.0e}: beobachtet {beob:>8,}  erwartet {erw:>12,.2f}  z = {z:+.2f}')
assert ergebnis_a[-1]['beobachtet'] == 56, 'bei t = 1e-7 muss die Anzahl die Treffer sein'
# the lowest percent in 10 subclasses, exactly discrete (`kanten` = class edges, `beob10`/`erw10` = observed / expected counts)
kanten = np.linspace(0, 0.01, 11)
beob10 = [int(np.sum((Qv >= np.ceil(kanten[i]*Pv)) & (Qv < np.ceil(kanten[i+1]*Pv)))) for i in range(10)]
erw10 = [float(((np.ceil(kanten[i+1]*Pv) - np.ceil(kanten[i]*Pv)) / Pv).sum()) for i in range(10)]
x2a, pa = chi2_p(beob10, erw10)   # chi-square statistic and p-value of the 10-class test
sag(f'    unterstes Prozent in 10 Teilklassen: χ² = {x2a:.1f}, p = {pa:.3f}   (beobachtet {beob10})')
za_max = max(abs(x['z']) for x in ergebnis_a)   # largest |z| over the thresholds

# ---------------- (b) by order and splitting behavior
sag(); sag('(b) NACH ORDNUNG — Klasse s = v₂(ord ε mod p) und e = (m/p)')
NAME = {0: 'ungerade Ordnung', 1: 'Ordnung ≡ 2 (mod 4)', 2: 'Ordnung ≡ 0 (mod 4)'}   # class labels (printed)
ergebnis_b = []   # `ergebnis_b` = one row per (class s, e)
for sk in (0, 1, 2):
    for e in (1, -1):
        sel = (np.minimum(Sv, 2) == sk) & (Ev == e)
        n = int(sel.sum())
        if n == 0:
            sag(f'    {NAME[sk]:<22} e = {e:+d}: keine Primzahlen'); continue
        h = np.bincount(np.minimum(19, (Yv[sel]*20).astype(int)), minlength=20)
        x2, pw = chi2_p(h.tolist(), [n/20]*20)
        tr = int(np.sum(Qv[sel] == 0)); erw_tr = float((1/Pv[sel]).sum())   # hits and expected hits Σ 1/p in this class
        ergebnis_b.append(dict(klasse=NAME[sk], e=e, n=n, chi2=round(x2, 2), p=round(pw, 4), treffer=tr, erwartet=round(erw_tr, 2)))
        sag(f'    {NAME[sk]:<22} e = {e:+d}: {n:>9,} Primzahlen · χ² = {x2:5.1f}, p = {pw:.3f} · Treffer {tr:>2} gegen {erw_tr:5.1f}')
pb_min = min(x['p'] for x in ergebnis_b)   # smallest p-value over the classes

# ---------------- (c) cross-correlation at ALL primes
sag(); sag('(c) KREUZKORRELATION zwischen den Kernen, an ALLEN gemeinsamen Primzahlen')
def korr(a, b):   # `korr` = Pearson correlation of two y-vectors on the primes where both are defined; returns (r, n)
    ok = ~np.isnan(a) & ~np.isnan(b); x, y = a[ok], b[ok]
    r = float(np.corrcoef(x, y)[0, 1]); return r, int(ok.sum())
# self-correlation (must be 1) and a one-prime shift (must be ~0)
r_s, n_s = korr(Y[0], Y[0]); r_v, n_v = korr(Y[0][:-1], Y[0][1:])
assert abs(r_s - 1) < 1e-12 and abs(r_v*math.sqrt(n_v)) < 4
rng = np.random.default_rng(20260924)
maske = rng.random(NP) < 0.02   # `maske` = the 2 % of primes at which kernel 1 is planted into kernel 2
y2p = np.where(maske & ~np.isnan(Y[0]), Y[0], Y[1])   # `y2p` = kernel 2 with the planted correlation
r_o, n_o = korr(Y[0], Y[1]); r_p, n_p = korr(Y[0], y2p)   # unplanted and planted correlation of the same pair
assert abs(r_p*math.sqrt(n_p)) > 5 and abs(r_o*math.sqrt(n_o)) < 4, ('PK4 VERLETZT', r_p, r_o)
sag(f'PK4 ✅  gepflanzt (2 %): z = {r_p*math.sqrt(n_p):+.1f} erkannt; dasselbe Paar ohne Pflanzung z = {r_o*math.sqrt(n_o):+.2f}; '
    f'Selbst r = 1, um eine Primzahl verschoben z = {r_v*math.sqrt(n_v):+.2f}.')
zs = []   # `zs` = (z, m1, m2, n) for every pair of kernels, z = r*sqrt(n)
for a_, b_ in itertools.combinations(range(len(KERNE)), 2):
    r, n = korr(Y[a_], Y[b_]); zs.append((r*math.sqrt(n), KERNE[a_], KERNE[b_], n))
# `N` = number of kernel pairs; `u1` / `u5` = pairs with |z| > 2.576 / > 1.96; `zmax` = the pair with the largest |z|;
#   `pbon` = its Bonferroni p-value
N = len(zs); u1 = sum(1 for z in zs if abs(z[0]) > 2.576); u5 = sum(1 for z in zs if abs(z[0]) > 1.96)
zmax = max(zs, key=lambda t: abs(t[0])); pbon = min(1.0, math.erfc(abs(zmax[0])/math.sqrt(2)) * N)
zm = float(np.mean([z[0] for z in zs])); zsd = float(np.std([z[0] for z in zs], ddof=1))   # mean and standard deviation of the z values
sag(f'    {N} Paare, je ~{zs[0][3]:,} Primzahlen: |z| > 2,576 bei {u1} (erwartet ~3), |z| > 1,96 bei {u5} (~15); '
    f'groesstes |z| = {abs(zmax[0]):.2f} (m = {zmax[1]} / {zmax[2]}), Bonferroni p = {pbon:.3f}; Mittel {zm:+.3f}, Streuung {zsd:.3f}.')

sag(); sag('='*100)
auff = []   # `auff` = the parts (a), (b), (c) that look remarkable
if za_max >= 3 or pa < 0.01: auff.append('(a)')
if pb_min < 0.01/6: auff.append('(b)')
if pbon < 0.01 or u1 > 9: auff.append('(c)')
urteil = ('KEINE der drei Nachmessungen weicht ab — auch nahe q = 0, auch nach Ordnung getrennt, auch an allen Primzahlen.'
          if not auff else f'AUFFAELLIG in {", ".join(auff)} — vor jeder Deutung an einem anderen Bereich wiederholen.')
sag('  ' + urteil); sag('='*100)
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), p_max=P_MAX, n_paare=N_ALL,
           a_schwellen=ergebnis_a, a_unterstes_prozent=dict(beobachtet=beob10, erwartet=[round(x, 2) for x in erw10], chi2=round(x2a, 2), p=round(pa, 4)),
           b_klassen=ergebnis_b, c=dict(paare=N, ueber_2_576=u1, ueber_1_96=u5, z_max=[round(zmax[0], 3), zmax[1], zmax[2]],
           p_bonferroni=round(pbon, 4), z_mittel=round(zm, 4), z_streuung=round(zsd, 4), pk_gepflanzt_z=round(r_p*math.sqrt(n_p), 2)),
           urteil=urteil, laufzeit_s=round(time.time()-t0))
(ERG/'w96_wieferich_nachmessung_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w96_wieferich_nachmessung_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
