# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w108_gedaechtnis_entlang_p_2026-09-25.py
# Purpose: does the Wieferich quotient y_p = ((B // p) % p) / p (with eps^(p-e) = A + B*sqrt(m) mod p^2, e = (m/p)) have a memory
#   along the primes? For each of the 25 kernels: autocorrelation r(k), k = 1 ... 20, of the sequence y_p over all primes p <=
#   3e6, z = r*sqrt(n), and the Ljung-Box statistic Q = n(n+2) * sum r(k)^2/(n-k) ~ chi^2(20).
# Reads: ergebnisse/w67_beweismenge_result.json (key `kerne`: the 25 kernels).
# Writes: ergebnisse/w108_gedaechtnis_entlang_p_result.json and ergebnisse/w108_gedaechtnis_entlang_p_output.txt.
# Usage: python w108_gedaechtnis_entlang_p_2026-09-25.py (no arguments; the results folder is found relative to the script).
# Controls: PK: an artificially built-in memory (y'_i = y_i, replaced by y_(i-1) with probability 5 %) must be detected at k = 1
#   with |z| > 5; the same sequence shuffled (order destroyed) must show nothing (|z| < 4). Expectation E1 below.

# Background: earlier runs measured the distribution (w91), fine structure / order / cross-kernel behavior (w96) and distances
#   between values (w105). Not yet measured: whether y_p at the NEXT prime (or the k-th next) is related to y_p, i.e. an order in
#   the SEQUENCE rather than in the distribution.
# Expectations (stated before the run):
#   E1  no memory: of 500 z-values about 5 exceed |z| > 2.576 (1 %), the largest |z| is around 3.5; Bonferroni p > 0.01;
#       the 25 Ljung-Box p-values are roughly uniformly distributed (about 1-2 below 0.05).
#   PK  an artificially built-in memory (y'_i = y_i, replaced by y_(i-1) with probability 5 %) is detected at k = 1 with
#       |z| > 5; the same sequence shuffled (order destroyed) shows nothing.

import sys, json, math, time, pathlib
from math import isqrt
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
ERG = pathlib.Path(__file__).resolve().parent.parent/'ergebnisse'
aus = []   # `aus` = output lines, written to the _output.txt file
def sag(s=''):   # `sag` = say: print a line and record it in `aus`
    print(s, flush=True); aus.append(s)
t0 = time.time()
sag('='*100); sag('w108 — GEDAECHTNIS ENTLANG p: Autokorrelation der Wieferich-Quotienten, 25 Kerne, Abstand 1 … 20   ' + time.strftime('%Y-%m-%d %H:%M')); sag('='*100)
# fundamental solution (T1, U1) of x^2 - m*y^2 = 1 from the continued fraction of sqrt(m)
def fund(m):
    a0 = isqrt(m); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (m - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
# (T1 + U1*sqrt(m))^e modulo M by binary exponentiation; returns (rational part, sqrt(m) part)
def eps_pow(T1, U1, m, e, M):
    ra, rb, ba, bb = 1 % M, 0, T1 % M, U1 % M
    while e:
        if e & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        e >>= 1
    return ra, rb
N3 = 3_000_000
sieb = bytearray([1]) * (N3 + 1); sieb[0] = sieb[1] = 0   # `sieb` = sieve of Eratosthenes up to N3
for p in range(2, isqrt(N3) + 1):
    if sieb[p]: sieb[p*p::p] = bytearray(len(range(p*p, N3 + 1, p)))
PR = [p for p in range(3, N3 + 1) if sieb[p]]   # odd primes up to N3
KERNE = json.load(open(ERG/'w67_beweismenge_result.json', encoding='utf-8'))['kerne']   # `KERNE` = the kernels m
LAGS = 20
# `autokorr` = autocorrelation: r(k), k = 1 ... lags, of the mean-centered sequence y; returns (list of r, n = length)
def autokorr(y, lags=LAGS):
    y = np.asarray(y, dtype=float); y = y - y.mean(); v = float(np.dot(y, y)); n = y.size
    return [float(np.dot(y[:-k], y[k:]) / v) for k in range(1, lags + 1)], n
# Ljung-Box statistic Q and its upper-tail p-value (Wilson-Hilferty normal approximation of chi^2 with df = number of lags)
def ljung_box(r, n):
    Q = n*(n + 2)*sum(rk*rk/(n - k) for k, rk in enumerate(r, 1))
    df = len(r); z = ((Q/df)**(1/3) - (1 - 2/(9*df))) / math.sqrt(2/(9*df))
    return Q, 0.5*math.erfc(z/math.sqrt(2))
# `res` = result dict: `skript` = script name, `datum` = date, `lags` = number of lags, `kerne` = per-kernel results (`z` =
#   z-values, `ljung_box_Q`, `ljung_box_p`), `pk` = positive control (`z_eingebaut` = z with the built-in memory, `z_gemischt` = z
#   of the shuffled sequence), `n_tests`, `ueber_2_576` = number of |z| > 2.576, `z_max` = largest |z|, `p_bonferroni`,
#   `ljung_box_unter_5` = number of Ljung-Box p < 0.05, `urteil` = verdict, `laufzeit_s` = run time in seconds
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), lags=LAGS, kerne={})
# `alle_z` = list of (|z|, m, lag) over all kernels and lags; `lb_p` = Ljung-Box p-values; `pk_gemacht` = flag: positive control
#   done
alle_z, lb_p, pk_gemacht = [], [], False
rng = np.random.default_rng(20260925)
for m in KERNE:
    T1, U1 = fund(m); ys = []
    for p in PR:
        if m % p == 0 or U1 % p == 0: continue
        e = 1 if pow(m % p, (p-1)//2, p) == 1 else -1   # e = Legendre symbol (m/p)
        A, B = eps_pow(T1, U1, m, p - e, p*p); ys.append(((B // p) % p) / p)
    r, n = autokorr(ys); z = [rk*math.sqrt(n) for rk in r]; Q, pQ = ljung_box(r, n)
    alle_z.extend((abs(zz), m, k) for k, zz in enumerate(z, 1)); lb_p.append(pQ)
    res['kerne'][str(m)] = dict(n=n, z=[round(x, 3) for x in z], ljung_box_Q=round(Q, 2), ljung_box_p=round(pQ, 4))
    if not pk_gemacht:                          # positive control, run once at the first kernel
        # built-in memory: each y_i is replaced by y_(i-1) with probability 5 % (`maske` = mask of replaced positions)
        y = np.array(ys); y2 = y.copy(); maske = rng.random(y.size) < 0.05; y2[1:][maske[1:]] = y[:-1][maske[1:]]
        r2, n2 = autokorr(y2, 1); r3, _ = autokorr(rng.permutation(y2), 1)
        assert abs(r2[0]*math.sqrt(n2)) > 5 and abs(r3[0]*math.sqrt(n2)) < 4, ('PK VERLETZT', r2, r3)
        sag(f'PK ✅  eingebaute 5-%-Erinnerung: z = {r2[0]*math.sqrt(n2):+.1f} erkannt; dieselbe Folge gemischt: z = {r3[0]*math.sqrt(n2):+.2f}.')
        res['pk'] = dict(z_eingebaut=round(r2[0]*math.sqrt(n2), 2), z_gemischt=round(r3[0]*math.sqrt(n2), 2)); pk_gemacht = True
    sag(f'  m = {m:>8}: n = {n:,} · groesstes |z| = {max(abs(x) for x in z):.2f} (Abstand {1 + int(np.argmax(np.abs(z)))}) · Ljung–Box p = {pQ:.3f}')
# summary over all kernels: N tests, `u1` = number with |z| > 2.576, `zmax` = largest (|z|, m, lag), `pbon` = Bonferroni p of the
#   maximum
N = len(alle_z); u1 = sum(1 for a, _, _ in alle_z if a > 2.576); zmax = max(alle_z)
pbon = min(1.0, math.erfc(zmax[0]/math.sqrt(2)) * N)
lb_u5 = sum(1 for p in lb_p if p < 0.05)   # number of kernels with Ljung-Box p < 0.05
sag(); sag(f'E1  {N} Tests: |z| > 2,576 bei {u1} (erwartet ~{0.01*N:.0f}); groesstes |z| = {zmax[0]:.2f} (m = {zmax[1]}, Abstand {zmax[2]}), '
           f'Bonferroni p = {pbon:.3f}; Ljung–Box unter 0,05 bei {lb_u5} von {len(lb_p)} Kernen (erwartet ~1,25).')
# `urteil` = verdict: no memory if the Bonferroni p exceeds 0.01, at most 12 exceedances of 2.576 and at most 4 Ljung-Box p below
#   0.05
urteil = ('KEIN Gedaechtnis entlang der Primzahlen — die Quotienten folgen einander wie unabhaengige Wuerfe.' if pbon > 0.01 and u1 <= 12 and lb_u5 <= 4
          else 'AUFFAELLIG — vor jeder Deutung an einem zweiten Bereich wiederholen.')
sag('    ' + urteil)
res.update(n_tests=N, ueber_2_576=u1, z_max=[round(zmax[0], 3), zmax[1], zmax[2]], p_bonferroni=round(pbon, 4), ljung_box_unter_5=lb_u5,
           urteil=urteil, laufzeit_s=round(time.time()-t0))
(ERG/'w108_gedaechtnis_entlang_p_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w108_gedaechtnis_entlang_p_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
