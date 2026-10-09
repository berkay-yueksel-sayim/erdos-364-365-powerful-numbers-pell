# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# m9_trennschaerfe_2026-09-11.py
# Power check of the overdispersion statistic: a control of the control.
#
# Reads:    nothing (the scaffold of the main run is rebuilt: kernels m = 7 (mod 8), m <= 10^4, and primes 3 <= p <= 3000).
# Writes:   nothing (printed output only; fixed random seed).
# Usage:    python m9_trennschaerfe_2026-09-11.py   (no arguments; uses numpy)
# Controls: PK (power control): at rho = 0 the statistic z must lie within the noise, at rho = 1 it must be clearly positive;
#           if z(rho=1) is not clearly larger than z(rho=0), the run aborts ("PK GEFEUERT").
#           The scaffold must reproduce the main run exactly (750 kernels, 120373 pairs), otherwise an assertion fails.
#
# BACKGROUND: the main run gives z = -0.87, i.e. no overdispersion. A null result is only worth as much as the ability of the test
#   to SEE the opposite. So far it has only been shown that the statistic stays quiet for INDEPENDENT data (negative control).
#   Open question: does it respond to COUPLED data?
#
# COUPLING MODEL (marginal-preserving, which is the point):
#   For each kernel, with probability rho a COMMON shock u_m ~ U(0,1) is drawn and event_i := (u_m < p_i) is set, so that all
#   events of a kernel are then comonotone. With probability 1-rho the events are drawn independently.
#   In BOTH cases the marginal probability stays exactly p_i = 1/p. The coupling therefore changes ONLY the joint distribution,
#   which is exactly the quantity under test. A test that cannot see this could not support the main run's null result.
#
# PK: at rho = 0, z must be within the noise; at rho = 1, z must be clearly positive.
#   If z(rho=1) is not clearly larger than z(rho=0), the statistic is BLIND and the null result is worthless.
import sys, time
from math import isqrt
import numpy as np
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.set_int_max_str_digits(2000000)
t0 = time.time()

def sieve(n):   # sieve of Eratosthenes: all primes <= n
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
PRIMES = sieve(200000)
def factor_small(n):   # factorization of n by trial division with `PRIMES`
    f = {}
    for p in PRIMES:
        if p * p > n: break
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
def fund(m):   # `fund` = fundamental solution (T_1, U_1) of x^2 - m y^2 = 1 from the continued fraction of sqrt(m)
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0

# ---------- scaffold identical to the main run: only (kernel, p) is needed, no Wieferich tests ----------
PR = [p for p in PRIMES if 2 < p <= 3000]   # `PR` = odd primes up to 3000
Pn = np.array(PR, dtype=np.int64)
fams = []   # `fams` = families (m, T_1, U_1): squarefree kernels m = 7 (mod 8), m <= 10^4, with T_1 even
for m in range(7, 10**4 + 1, 8):
    f = factor_small(m)
    if any(e > 1 for e in f.values()): continue
    T1, U1 = fund(m)
    if T1 % 2 == 0: fams.append((m, T1, U1))
# `kern_id` = kernel index of each pair, `prim` = its prime p (a pair is kept if alpha(p) exists and p does not divide m)
kern_id, prim = [], []
for idx, (m, T1, U1) in enumerate(fams):
    T1m = np.array([T1 % p for p in PR], dtype=np.int64); c = (2 * T1m) % Pn
    a = np.zeros(len(PR), dtype=np.int64); tp = np.ones(len(PR), dtype=np.int64); tc = T1m.copy()
    for k in range(1, 3003):   # recurrence T_{k+1} = 2 T_1 T_k - T_{k-1} mod p; `a[j]` = first k with T_k = 0 (mod p), else 0
        z = (tc == 0) & (a == 0)
        if z.any(): a[z] = k
        tp, tc = tc, np.mod(c * tc - tp, Pn)
    for j, p in enumerate(PR):
        al = int(a[j])
        if al == 0 or m % p == 0: continue
        kern_id.append(idx); prim.append(p)
kern_id = np.array(kern_id); pi = 1.0 / np.array(prim, dtype=np.float64)   # `pi` = marginal event probability 1/p of each pair
K = len(fams)   # `K` = number of kernels
assert len(pi) == 120373 and K == 750, f"Geruest weicht ab: {K} Kerne, {len(pi)} Paare"
print(f"Geruest ok [{time.time()-t0:.0f}s]: {K} Kerne, {len(pi)} Paare (identisch zum Hauptlauf)\n", flush=True)

lam_m = np.bincount(kern_id, weights=pi, minlength=K)   # `lam_m` = expected number of events per kernel
def S_von(lab):   # `S_von` = statistic S: sum over kernels of (event count - expected count)^2, for the 0/1 labeling `lab`
    return float(np.sum((np.bincount(kern_id, weights=lab, minlength=K) - lam_m) ** 2))

rng = np.random.default_rng(20260911)
def ziehe(rho):   # `ziehe` = draw one 0/1 labeling of all pairs under the coupling model with parameter rho
    """marginal-preserving: per kernel with probability rho comonotone (common shock), otherwise independent."""
    u_gem = rng.random(K)[kern_id]                 # common shock per kernel
    u_ind = rng.random(len(pi))                    # independent per pair
    gek   = (rng.random(K) < rho)[kern_id]         # which kernels are coupled
    u = np.where(gek, u_gem, u_ind)
    return (u < pi).astype(np.float64)

# null distribution (rho = 0), as in the main run
REPS = 4000   # `REPS` = number of null repetitions
S0 = np.array([S_von((rng.random(len(pi)) < pi).astype(np.float64)) for _ in range(REPS)])   # `S0` = null sample of S
mu, sd = S0.mean(), S0.std(ddof=1)   # mean and standard deviation of S under the null
print(f"Null (rho=0): S = {mu:.1f} +/- {sd:.1f}\n")

print("=== Trennschaerfe: bewegt Kopplung die Statistik ueberhaupt? ===")
print("  rho   mittleres z   Anteil Laeufe mit z > 2   (40 Laeufe je rho)")
erg = {}   # `erg` = for each rho: (mean z, fraction of runs with z > 2)
for rho in (0.0, 0.05, 0.1, 0.2, 0.5, 1.0):
    zs = np.array([(S_von(ziehe(rho)) - mu) / sd for _ in range(40)])   # `zs` = z values of 40 runs
    erg[rho] = (zs.mean(), float(np.mean(zs > 2)))
    print(f"  {rho:4.2f}   {zs.mean():+9.2f}   {np.mean(zs > 2)*100:20.0f} %")

# ---------- PK: blind or not? ----------
if erg[1.0][0] <= erg[0.0][0] + 3:
    sys.exit(f"PK GEFEUERT: volle Kopplung gibt z = {erg[1.0][0]:+.2f}, unabhaengig gibt {erg[0.0][0]:+.2f}. "
             "Die Statistik ist BLIND — der Nullbefund des Hauptlaufs traegt nicht. ABBRUCH.")
print(f"\nPK ok: volle Kopplung hebt z von {erg[0.0][0]:+.2f} auf {erg[1.0][0]:+.2f}. Die Statistik SIEHT Kopplung.")

# smallest rho whose mean z exceeds 2
nach = [r for r in sorted(erg) if erg[r][0] > 2]   # `nach` = tested rho values with mean z > 2
print(f"\nKLEINSTE erkannte Kopplung: rho = {nach[0] if nach else 'keine der getesteten'}")
print(f"Beobachtet im Hauptlauf: z = -0.87  =>  eine Kopplung dieser Staerke oder groesser ist ausgeschlossen;")
print(f"schwaechere Kopplung bleibt moeglich und ist mit diesen Daten NICHT entscheidbar.")
print(f"\n[{time.time()-t0:.0f}s]")
