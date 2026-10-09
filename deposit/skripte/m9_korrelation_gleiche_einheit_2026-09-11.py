# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# m9_korrelation_gleiche_einheit_2026-09-11.py
# Purpose: are the Wieferich conditions of different primes independent, or does the shared unit eps_m couple them? Tests whether
#   the hits (events with p^2 | T_alpha) cluster within ONE kernel more strongly than independence allows (statistic S,
#   simulated null).
# Reads: nothing (kernels m = 7 mod 8, m <= 10^4, with T1 even, and odd primes p <= 3000 are generated).
# Writes: m9_korrelation_result.json in the current directory.
# Usage: python m9_korrelation_gleiche_einheit_2026-09-11.py   (no arguments)
# Controls (both abort the run): PK+ must reproduce the published counts exactly; PK- (synthetic independent data) must not
#   show up as overdispersed.
#
# Background. All conditions depend on the SAME unit, yet Cor. 4.8 counts tau(k)-1 conditions and treats them as independent, and
#   Prop. 6.1 measures only the MARGINAL rate (470 observed against 486 expected), never the JOINT distribution.
#   Statistic S = sum over kernels of (k_m - lambda_m)^2, where k_m is the number of hits in kernel m and lambda_m the expected
#   number (sum of 1/p); null distribution by simulation (REPS = 4000 independent Bernoulli(1/p) draws per pair).
# Method: the apparition loop is taken verbatim from auswertung_tueren_2026-09-07.py; only the bookkeeping changes: instead of
#   counters, every pair (m, p, alpha, w2, w3) is kept (`w2`, `w3`: p^2 resp. p^3 divides T_alpha).
# PK+ (at the comparison step, not at the extraction): the rebuild MUST reproduce the published numbers exactly: 750 families,
#   120373 pairs, 470 (w2), 52 (w3), h1 ~ 486, h2 ~ 49.1 (h1 = sum of 1/p, h2 = sum of 1/p^2). If anything deviates, the
#   computation was altered, and the data have not said anything new.
# PK- : a SYNTHETICALLY independent data set (Bernoulli(1/p) per pair) must pass through exactly the same evaluation path and must
#   NOT show up as overdispersed (|z| <= 4). If it does, the statistic is broken.

import sys, json, math, time
from math import isqrt
import numpy as np
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.set_int_max_str_digits(2000000)
t0 = time.time()

# ---------- helper functions: verbatim from auswertung_tueren_2026-09-07.py ----------
def sieve(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
# `PRIMES` = primes up to 200000 (`sieve` returns the list of primes up to n)
PRIMES = sieve(200000)
# `factor_small` = factorization {prime: exponent} by trial division with `PRIMES`
def factor_small(n):
    f = {}
    for p in PRIMES:
        if p * p > n: break
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
# `fund` = fundamental solution (T1, U1) of x^2 - m*y^2 = 1 by continued fraction (big integers)
def fund(m):
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
# `eps_pow_mod` = T_k mod M;  `T_exact` = exact T_k (both by binary powering of T1 + U1*sqrt(m))
def eps_pow_mod(T1, U1, m, k, M):
    ra, rb = 1 % M, 0; ba, bb = T1 % M, U1 % M
    while k:
        if k & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        k >>= 1
    return ra
def T_exact(T1, U1, m, k):
    ra, rb, ba, bb = 1, 0, T1, U1
    while k:
        if k & 1: ra, rb = ra*ba + m*rb*bb, ra*bb + rb*ba
        ba, bb = ba*ba + m*bb*bb, 2*ba*bb
        k >>= 1
    return ra

# ---------- positive controls of the basic functions (taken over verbatim) ----------
# fund(7) = (8, 3); T_7(7) = 2^3 * 29 * 197 * 2857 (Mollin-Walsh); v_29(T_7(7)) = 1
assert fund(7) == (8, 3) and T_exact(8, 3, 7, 7) == 2**3 * 29 * 197 * 2857
assert eps_pow_mod(8, 3, 7, 7, 29**2) % 29 == 0 and eps_pow_mod(8, 3, 7, 7, 29**2) != 0, "v_29(T_7(7)) = 1"
print("PK Grundfunktionen ok", flush=True)

# ---------- apparition loop: mathematics unchanged, bookkeeping extended ----------
# `PR` = odd primes up to 3000; `fams` = families: squarefree kernels m = 7 (mod 8), m <= 10^4, with T1 even
PR = [p for p in PRIMES if 2 < p <= 3000]
Pn = np.array(PR, dtype=np.int64)
fams = []
for m in range(7, 10**4 + 1, 8):
    f = factor_small(m)
    if any(e > 1 for e in f.values()): continue
    T1, U1 = fund(m)
    if T1 % 2 == 0: fams.append((m, T1, U1))

# per pair (m, p): kernel index `kern_id`, prime `prim`, flags `w2_flag`, `w3_flag` (p^2 resp. p^3 divides T_alpha);
# `n_alpha` = number of pairs, `n_w2`, `n_w3` = number of hits, `h1`, `h2` = sums of 1/p and 1/p^2 (expectations)
kern_id, prim, w2_flag, w3_flag = [], [], [], []
n_alpha = n_w2 = n_w3 = 0
h1 = h2 = 0.0
for idx, (m, T1, U1) in enumerate(fams):
    # rank of apparition `a[j]` of each prime = first k <= 3002 with T_k = 0 (mod p) (0 if none), by the recurrence
    # T_{k+1} = 2*T1*T_k - T_{k-1} mod p; `tp`, `tc` = previous and current T_k mod p; `c` = 2*T1 mod p
    T1m = np.array([T1 % p for p in PR], dtype=np.int64); c = (2 * T1m) % Pn
    a = np.zeros(len(PR), dtype=np.int64); tp = np.ones(len(PR), dtype=np.int64); tc = T1m.copy()
    for k in range(1, 3003):
        z = (tc == 0) & (a == 0)
        if z.any(): a[z] = k
        tp, tc = tc, np.mod(c * tc - tp, Pn)
    for j, p in enumerate(PR):
        al = int(a[j])
        if al == 0 or m % p == 0: continue
        n_alpha += 1; h1 += 1.0 / p; h2 += 1.0 / (p * p)
        t2 = eps_pow_mod(T1, U1, m, al, p * p)
        w2 = (t2 == 0)
        w3 = bool(w2 and eps_pow_mod(T1, U1, m, al, p**3) == 0)
        if w2: n_w2 += 1
        if w3: n_w3 += 1
        kern_id.append(idx); prim.append(p); w2_flag.append(1 if w2 else 0); w3_flag.append(1 if w3 else 0)
print(f"Schleife durch [{time.time()-t0:.0f}s]: Familien {len(fams)}, Paare {n_alpha}, w2 {n_w2}, w3 {n_w3}, h1 {h1:.0f}, h2 {h2:.1f}", flush=True)

# ---------- PK+ : the published numbers MUST come out exactly (`soll` = expected, `ist` = obtained) ----------
soll = {"fams": 750, "paare": 120373, "w2": 470, "w3": 52, "h1": 486, "h2": 49.1}
ist  = {"fams": len(fams), "paare": n_alpha, "w2": n_w2, "w3": n_w3, "h1": round(h1), "h2": round(h2, 1)}
if ist != soll:
    sys.exit(f"PK+ GEFEUERT: Umbau hat die Rechnung veraendert.\n  soll {soll}\n  ist  {ist}\nABBRUCH.")
print("PK+ ok: alle fuenf veroeffentlichten Zahlen exakt reproduziert.\n", flush=True)

kern_id = np.array(kern_id); prim = np.array(prim, dtype=np.float64)
w2 = np.array(w2_flag, dtype=np.float64)
pi = 1.0 / prim  # single probability per pair (the heuristic)
K = len(fams)

# ---------- statistic S: clustering WITHIN a kernel ----------
def S_von(labels):
    """S = sum_m (k_m - lambda_m)^2 over the kernels; k_m = hits in kernel m, lambda_m = expected hits (`labels` = 0/1 flags)."""
    k_m   = np.bincount(kern_id, weights=labels, minlength=K)
    lam_m = np.bincount(kern_id, weights=pi,     minlength=K)
    return float(np.sum((k_m - lam_m) ** 2)), k_m, lam_m

S_obs, k_obs, lam = S_von(w2)

# ---------- null distribution: simulated instead of guessing a formula ----------
REPS = 4000
rng = np.random.default_rng(20260911)
S_null = np.empty(REPS)
for r in range(REPS):
    S_null[r] = S_von((rng.random(len(pi)) < pi).astype(np.float64))[0]
mu, sd = S_null.mean(), S_null.std(ddof=1)
# z-score of the observed S; `p_emp` = empirical one-sided p-value (overdispersion)
z = (S_obs - mu) / sd
p_emp = float((np.sum(S_null >= S_obs) + 1) / (REPS + 1))

# ---------- PK- : a synthetically INDEPENDENT data set must not stand out ----------
syn = (rng.random(len(pi)) < pi).astype(np.float64)
S_syn = S_von(syn)[0]
z_syn = (S_syn - mu) / sd
if abs(z_syn) > 4:
    sys.exit(f"PK- GEFEUERT: synthetisch unabhaengige Daten geben z = {z_syn:.2f}. Die Statistik ist kaputt. ABBRUCH.")
print(f"PK- ok: synthetisch unabhaengiger Datensatz gibt z = {z_syn:+.2f} (liegt im Rauschen).\n", flush=True)

# ---------- result ----------
print("=== (A) Klumpung INNERHALB der Kerne ===")
print(f"  Kerne {K} | Paare {len(pi)} | lambda_m: min {lam.min():.3f}  median {np.median(lam):.3f}  max {lam.max():.3f}")
print(f"  Treffer je Kern: {np.bincount(k_obs.astype(int))}   (Index = Anzahl Treffer)")
print(f"  S_beobachtet = {S_obs:.1f}")
print(f"  S_null       = {mu:.1f} +/- {sd:.1f}   (4000 Simulationen des Unabhaengigkeits-Modells)")
print(f"  z = {z:+.2f}   empirisches p (einseitig, Ueberstreuung) = {p_emp:.4f}")

print("\n=== (B) Trennung: koppelt die EINHEIT, oder ist die PRIMZAHL oft Wieferich? ===")
# `tre` = per prime p: [observed hits, expected hits]; `auff` = conspicuous primes (expected >= 1), sorted by deviation
tre = {}
for p_, wv in zip(prim, w2):
    d = tre.setdefault(int(p_), [0, 0.0])
    d[0] += int(wv); d[1] += 1.0 / p_
auff = sorted(((p_, o, e) for p_, (o, e) in tre.items() if e >= 1.0),
              key=lambda t: -(t[1] - t[2]) / math.sqrt(t[2]))
print("  groesste POSITIVE Abweichungen (p, beobachtet, erwartet, z):")
for p_, o, e in auff[:8]:
    print(f"     p = {p_:5d}   beob {o:4d}   erw {e:7.2f}   z = {(o-e)/math.sqrt(e):+5.2f}")
print("  groesste NEGATIVE Abweichungen:")
for p_, o, e in auff[-5:]:
    print(f"     p = {p_:5d}   beob {o:4d}   erw {e:7.2f}   z = {(o-e)/math.sqrt(e):+5.2f}")

# result record; `treffer_je_kern` = number of kernels with 0, 1, 2, ... hits (index = number of hits)
res = {"fams": K, "paare": int(len(pi)), "w2": int(n_w2), "w3": int(n_w3),
       "S_obs": S_obs, "S_null_mu": float(mu), "S_null_sd": float(sd), "z": float(z),
       "p_emp": p_emp, "z_synthetisch": float(z_syn), "reps": REPS, "seed": 20260911,
       "treffer_je_kern": np.bincount(k_obs.astype(int)).tolist()}
json.dump(res, open("m9_korrelation_result.json", "w"), indent=1)
print(f"\nJSON: m9_korrelation_result.json   [{time.time()-t0:.0f}s]")
