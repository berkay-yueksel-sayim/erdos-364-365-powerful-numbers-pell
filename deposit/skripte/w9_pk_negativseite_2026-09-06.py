# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w9_pk_negativseite_2026-09-06.py
# Negative-side control for w9: the rejection of more than 10^7 families by the w9 box run is itself checked here.
# Three independent checks on a sample of REJECTED families:
#   (1) m' from the fast route (U_1 mod m) against m' from the EXACT fundamental unit,
#   (2) log10 T_1 from the scaled floating-point value against the exact number of digits,
#   (3) RIGOROUS rejection without floating point: T_k >= T_1^k / 2 (k >= 1) and k >= m'
#       =>  digits(T_k) > m' * (digits(T_1) - 1) - 1.  If this exceeds H, the family is provably free of points.
# Addition: the A135735 kernels = 7 mod 8 (parity and digit count of T_1), for the open question on m = 209991.
# Own code, no outside source.
#
# Reads:    nothing. Writes: nothing (printed output only; fixed random seed).
# Usage:    python w9_pk_negativseite_2026-09-06.py   (no arguments)
# Controls: positive control of the check itself: the known candidate families m = 7 (m' = 7) and m = 4099215 (m' = 1) must NOT be
#           rejected. Any mismatch between the fast and the exact route is counted and printed; a family that is rejected but not
#           confirmed rigorously is printed as "RIGOROS NICHT BESTAETIGT".
import sys, random, math, time
from math import isqrt
sys.set_int_max_str_digits(2000000)

H = 2000; MMAX = 10**8; PMAX = 20000   # digit bound H for T_k, kernel bound MMAX, prime bound PMAX for the sieve
random.seed(20260906)

def sieve(n):   # sieve of Eratosthenes: all primes <= n
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
PRIMES = sieve(PMAX)

def factor_small(n):   # factorization of n by trial division with `PRIMES`
    f = {}
    for p in PRIMES:
        if p * p > n: break
        while n % p == 0:
            f[p] = f.get(p, 0) + 1; n //= p
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def fund(m):   # `fund` = exact fundamental solution (T_1, U_1) of x^2 - m y^2 = 1 from the continued fraction of sqrt(m)
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1
        k1, k0 = k0, a*k0 + k1
    return h0, k0

def fund_fast(m):
    # `fund_fast` = fast route: runs the continued fraction up to the end of an even period and returns
    # (T_1 mod m, U_1 mod m, log10 T_1), carrying T_1 only as a rescaled float (rescaled by 10^100 steps, offset kept in `s`)
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    hm1, hm0 = 1 % m, a0 % m; km1, km0 = 0, 1 % m
    hf1, hf0 = 1.0, float(a0); s = 0; i = 0
    while True:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q; i += 1
        if Q == 1 and i % 2 == 0: return hm0, km0, math.log10(hf0) + s
        hm1, hm0 = hm0, (a*hm0 + hm1) % m
        km1, km0 = km0, (a*km0 + km1) % m
        hf1, hf0 = hf0, a*hf0 + hf1
        if hf0 > 1e100: hf0 /= 1e100; hf1 /= 1e100; s += 100

# --- Positive control of the check: known candidate families must NOT be rejected ---
for m_known, mp_known in ((7, 7), (4099215, 1)):
    f = factor_small(m_known); tm, um, l10 = fund_fast(m_known)
    mp = 1   # `mp` = m' = product of the primes p | m with p not dividing U_1
    for p in f:
        if um % p: mp *= p
    assert mp == mp_known, (m_known, mp, mp_known)
    assert mp <= int((H + 0.31) / l10) + 1, ("PK: bekannte Kandidaten-Familie faellt durch", m_known)
print("PK ok: m=7 (m'=7) und m=4099215 (m'=1) werden als Kandidaten erkannt, nicht verworfen", flush=True)

def sample_sqfree(lo, hi, n):   # n distinct random squarefree m = 7 (mod 8) in [lo, hi), sorted
    out = set()
    while len(out) < n:
        m = random.randrange(lo, hi); m -= (m - 7) % 8
        if m < lo: continue
        if all(e == 1 for e in factor_small(m).values()): out.add(m)
    return sorted(out)

S = sample_sqfree(10**6, 10**7, 250) + sample_sqfree(10**7, 10**8, 250)   # `S` = the sample of kernels
top = []; m = MMAX - ((MMAX - 7) % 8)   # `top` = the 30 largest squarefree m = 7 (mod 8) below 10^8
while len(top) < 30:
    if all(e == 1 for e in factor_small(m).values()): top.append(m)
    m -= 8
S += sorted(top)

t0 = time.time()
# counters: `n_rej` = rejected, `n_cand` = candidates (not rejected), `n_rig_ok`/`n_rig_fail` = rigorous rejection
#   confirmed/not confirmed, `n_mp_mismatch`/`n_log_mismatch` = disagreements between the fast and the exact route;
#   `worst` = smallest margin (margin, m, m', digits)
n_rej = n_cand = n_rig_ok = n_rig_fail = n_mp_mismatch = n_log_mismatch = 0
maxdig = 0; worst = None
for i, m in enumerate(S):
    f = factor_small(m)
    tm, um, l10 = fund_fast(m)
    mp_fast = 1   # `mp_fast` = m' from the fast route
    for p in f:
        if um % p: mp_fast *= p
    kmax = int((H + 0.31) / l10) + 1   # `kmax` = largest k for which T_k can still have at most H digits (T_k >= T_1^k / 2)
    rejected = mp_fast > kmax
    T1, U1 = fund(m)
    mp_ex = 1   # `mp_ex` = m' from the exact fundamental unit
    for p in f:
        if U1 % p: mp_ex *= p
    dig = len(str(T1)); maxdig = max(maxdig, dig)
    if (T1 % m, U1 % m) != (tm, um): n_mp_mismatch += 1; print("  MOD-MISMATCH:", m)
    if abs(math.log10(T1) - l10) > 1e-6: n_log_mismatch += 1; print("  LOG-MISMATCH:", m, dig, l10)
    if mp_ex != mp_fast: n_mp_mismatch += 1; print("  M-STRICH-MISMATCH:", m, mp_fast, mp_ex)
    if rejected:
        n_rej += 1
        margin = mp_ex * (dig - 1) - 1 - H          # > 0  =>  provably free of points, without floating point
        if margin > 0: n_rig_ok += 1
        else: n_rig_fail += 1; print("  RIGOROS NICHT BESTAETIGT:", m, "m'=", mp_ex, "digits T1=", dig, "margin=", margin)
        if worst is None or margin < worst[0]: worst = (margin, m, mp_ex, dig)
    else:
        n_cand += 1; print("  Kandidat in der Stichprobe (kein Widerspruch):", m, "m'=", mp_ex, "digits T1=", dig, "kmax=", kmax)
    if (i + 1) % 100 == 0: print(f"  ... {i+1}/{len(S)}  {time.time()-t0:.0f}s", flush=True)

print()
print(f"Stichprobe: {len(S)} quadratfreie m = 7 mod 8 (250 aus [1e6,1e7], 250 aus [1e7,1e8], 30 groesste < 1e8)")
print(f"verworfen: {n_rej}   Kandidaten: {n_cand}")
print(f"Mismatch schnell/exakt (m' oder mod oder log): {n_mp_mismatch + n_log_mismatch}")
print(f"rigoros bestaetigt (ohne Float): {n_rig_ok}   NICHT bestaetigt: {n_rig_fail}")
print(f"knappste Reserve in log10-Einheiten ueber H=2000: {worst}   max Stellen T1: {maxdig}")
print(f"Zeit: {time.time()-t0:.0f}s")

print()
print("A135735-Kerne = 7 mod 8 (Reinhart 2024): Paritaet und Groesse von T_1")
for m in (209991, 4099215):
    T1, U1 = fund(m)
    print(f"  m = {m}: digits(T1) = {len(str(T1))}, T1 {'ungerade (paritaetstot)' if T1 % 2 else 'gerade'}, m | U1: {U1 % m == 0}")
print("  (117477414815 und 39028039587479 liegen ueber 1e8, also ausserhalb der W9-Box)")
