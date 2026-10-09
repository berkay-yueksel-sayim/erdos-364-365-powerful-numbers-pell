# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w8_bsignatur_box_2026-09-06.py
# "b-signature box" (ALL heights): no triple of consecutive powerful numbers with kernel m <= MMAX and
# squarefree part b <= BMAX of the middle.
#
# Foundations (credits):
#   Mollin-Walsh 1986 (C.R. Math. Rep. 8, pp. 110/111): triple <=> middle n = T_k(m), m squarefree = 7 mod 8,
#     k odd, T_k even and powerful, m | U_k;  lemma (p. 111): m | U_k <=> m | k*U_1  <=> m' | k.
#   Bennett-Walsh 1999 (Proc. AMS 127), Lemma 3.3: for a prime p (p not dividing m), {k : p | T_k} = {t*alpha(p): t odd};
#     Thm 1.2: for squarefree b > 1 there is at most ONE k with T_k = b*x^2, namely k = alpha(b).
#   Case b = 1 (square middle): Ljunggren 1942 / Cohn 1997: T_k is a square only for k in {1, 2}  =>  k = 1.
# Form: n powerful with squarefree part b  <=>  n = b^3 z^2.
# For each (m, b) there is exactly one candidate k = alpha_m(b) = lcm(alpha(p) : p | b) (only if all alpha(p) exist and have the
#   same 2-adic valuation; otherwise b never divides a T_k). Tests: k odd, m' | k, b^3 | T_k,
#   T_k / b^3 a square (quadratic-residue test modulo primes), parity (T_1 even). Own code.
#
# Reads:    nothing. Writes: w8_bsignatur_M<MMAX>_B<BMAX>_result.json in the working directory, and printed output.
# Usage:    python w8_bsignatur_box_2026-09-06.py [MMAX=10000] [BMAX=10000]   (uses numpy)
# Controls: positive controls at the start (abort on failure): fund(7) = (8, 3); T_7(7) = 2^3 * 29 * 197 * 2857;
#           alpha(p) for the primes of T_3(7) and T_7(7); Bennett-Walsh Lemma 3.3 at p = 11; the form test accepts 8 =
#           2^3 and rejects 2024 and 127.
import sys, math, time, json
from math import isqrt, gcd
import numpy as np
sys.set_int_max_str_digits(2000000)

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10000   # bound for the kernel m
BMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 10000   # bound for the squarefree part b of the middle
NQR  = 40   # square test modulo NQR primes; only if all pass: recompute exactly (if feasible) or 'open'

def sieve(n):   # sieve of Eratosthenes: all primes <= n
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

def spf_table(n):   # `spf_table` = table of the smallest prime factor of every number <= n
    spf = list(range(n + 1))
    for i in range(2, isqrt(n) + 1):
        if spf[i] == i:
            for j in range(i*i, n + 1, i):
                if spf[j] == j: spf[j] = i
    return spf

def factor_small(n, primes):   # factorization of n by trial division with the given primes
    f = {}
    for p in primes:
        if p * p > n: break
        while n % p == 0:
            f[p] = f.get(p, 0) + 1; n //= p
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def fund(m):   # `fund` = fundamental solution (T_1, U_1) of x^2 - m y^2 = 1 from the continued fraction of sqrt(m)
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1
        k1, k0 = k0, a*k0 + k1
    return h0, k0

def eps_pow_mod(T1, U1, m, k, M):
    """T_k mod M via square-and-multiply in Z[sqrt m]/(M)."""
    ra, rb = 1 % M, 0; ba, bb = T1 % M, U1 % M
    while k:
        if k & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        k >>= 1
    return ra

def T_exact(T1, U1, m, k):   # exact T_k by square-and-multiply
    ra, rb, ba, bb = 1, 0, T1, U1
    while k:
        if k & 1: ra, rb = ra*ba + m*rb*bb, ra*bb + rb*ba
        ba, bb = ba*ba + m*bb*bb, 2*ba*bb
        k >>= 1
    return ra

def v2(n):   # 2-adic valuation of n
    return (n & -n).bit_length() - 1

def is_square(n):
    if n < 0: return False
    r = isqrt(n); return r * r == n

PRIMES_ALL = sieve(max(BMAX, 200000))   # `PRIMES_ALL` = primes for factoring m and for the square test
PRIMES_B = [p for p in PRIMES_ALL if p <= BMAX]   # `PRIMES_B` = primes up to BMAX (the possible prime factors of b)
Pn = np.array(PRIMES_B, dtype=np.int64)
pidx = {p: i for i, p in enumerate(PRIMES_B)}   # `pidx` = index of each prime in `PRIMES_B`
QR = [p for p in PRIMES_ALL if 1000 < p][:NQR]          # fixed test primes for the square test
SPF = spf_table(BMAX)
def sqfree_list(B):   # `sqfree_list` = list of (b, [prime factors of b]) for all squarefree b <= B
    out = []
    for b in range(1, B + 1):
        n, ps, ok = b, [], True
        while n > 1:
            p = SPF[n]; n //= p
            if n % p == 0: ok = False; break
            ps.append(p)
        if ok: out.append((b, ps))
    return out
SQF = sqfree_list(BMAX)   # `SQF` = the squarefree b up to BMAX with their prime factors

def alphas_for_family(T1, U1, m):
    """alpha(p) for all p in PRIMES_B: smallest k >= 1 with T_k = 0 (mod p); 0 = does not exist.
    Recurrence T_{k+1} = 2 T_1 T_k - T_{k-1} (mod p), vectorized over all p."""
    T1m = np.array([T1 % p for p in PRIMES_B], dtype=np.int64)
    c = (2 * T1m) % Pn
    a = np.zeros(len(PRIMES_B), dtype=np.int64)
    t_prev = np.ones(len(PRIMES_B), dtype=np.int64)
    t_cur = T1m.copy()
    kmax_iter = int(Pn.max()) + 2         # alpha(p) <= p + 1 (the order of eps divides p -/+ 1)
    for k in range(1, kmax_iter + 1):
        z = (t_cur == 0) & (a == 0)
        if z.any(): a[z] = k
        t_prev, t_cur = t_cur, np.mod(c * t_cur - t_prev, Pn)
    return a

def legendre(a, p):
    return pow(a % p, (p - 1) // 2, p)     # 1, p-1 (= -1), or 0

def form_test(T1, U1, m, k, b):
    """Is T_k = b^3 z^2 ?  -> 'b3_fail' | 'nonsquare' | 'form_ok_exact' | 'form_ok_modular_open'"""
    b3 = b ** 3
    if eps_pow_mod(T1, U1, m, k, b3) != 0: return 'b3_fail'
    for q in QR:
        if b % q == 0 or m % q == 0: continue
        t = eps_pow_mod(T1, U1, m, k, q)
        if t == 0: continue                                  # q | T_k: no verdict
        r = (t * pow(b3 % q, q - 2, q)) % q                  # T_k / b^3 mod q
        if legendre(r, q) == q - 1: return 'nonsquare'
    # all QR tests passed: exact check if feasible
    if k * math.log10(T1 + 1) < 200000:
        T = T_exact(T1, U1, m, k)
        return 'form_ok_exact' if (T % b3 == 0 and is_square(T // b3)) else 'nonsquare_exact'
    return 'form_ok_modular_open'

# ---------- Positive controls (abort on failure) ----------
T1, U1 = fund(7); assert (T1, U1) == (8, 3)
assert T_exact(8, 3, 7, 7) == 2**3 * 29 * 197 * 2857                     # Mollin-Walsh 1986, p. 111
assert eps_pow_mod(8, 3, 7, 7, 10**6) == (2**3 * 29 * 197 * 2857) % 10**6
a7 = alphas_for_family(8, 3, 7)
assert a7[pidx[2]] == 1 and a7[pidx[11]] == 3 and a7[pidx[23]] == 3, "T_3(7) = 2^3*11*23"
assert a7[pidx[29]] == 7 and a7[pidx[197]] == 7 and a7[pidx[2857]] == 7, "T_7(7) = 2^3*29*197*2857"
assert a7[pidx[7]] == 0, "p | m teilt nie T_k"
assert eps_pow_mod(8, 3, 7, 9, 11) == 0 and eps_pow_mod(8, 3, 7, 6, 11) != 0 and eps_pow_mod(8, 3, 7, 15, 11) == 0, "B-W Lemma 3.3: 11 | T_k <=> k ungerades Vielfaches von 3"
assert form_test(8, 3, 7, 1, 2) == 'form_ok_exact', "PK Form: T_1(7) = 8 = 2^3 * 1^2"
assert form_test(8, 3, 7, 7, 2) in ('nonsquare', 'nonsquare_exact'), "T_7(7)/8 = 29*197*2857 ist kein Quadrat"
assert form_test(8, 3, 7, 3, 2) in ('b3_fail', 'nonsquare', 'nonsquare_exact'), "T_3(7) = 2024 = 2^3*11*23 ist nicht 8*z^2"
assert form_test(8, 3, 7, 2, 127) == 'b3_fail', "T_2(7) = 127 ist nicht 127^3 z^2"
print("PK ok: fund(7)=(8,3); T_7(7) = 2^3*29*197*2857 (M-W 1986 S.111); alpha(11)=alpha(23)=3, alpha(29)=alpha(197)=alpha(2857)=7;"
      " B-W Lemma 3.3 an p=11; Form-Test erkennt 8 = 2^3 und verwirft 2024, 127", flush=True)

# ---------- Run ----------
t0 = time.time()
# counters: `fam` = squarefree kernels m = 7 (mod 8) seen; `fam_par` = of these, dead by parity (T_1 odd);
#   `pairs` = pairs (m, b) considered
fam = fam_par = 0
pairs = 0
# `cnt` = number of pairs per exclusion reason: `no_alpha` = some alpha(p) does not exist, `twoadic` = different 2-adic
#   valuations, `k_even` = candidate k even, `lemmaL` = m' does not divide k, the rest = outcomes of `form_test`
cnt = dict(no_alpha=0, twoadic=0, k_even=0, lemmaL=0, b3_fail=0, nonsquare=0, nonsquare_exact=0, form_ok_exact=0, form_ok_modular_open=0)
survivors = []       # (m, b, k, result) for all pairs that survive Lemma L
open_pts = []        # the survivors whose form test did not exclude them
for m in range(7, MMAX + 1, 8):
    f = factor_small(m, PRIMES_ALL)
    if any(e > 1 for e in f.values()): continue
    fam += 1
    T1, U1 = fund(m)
    if T1 % 2 == 1:
        fam_par += 1; continue                       # whole family dead by parity (M-W p. 111)
    mp = 1   # `mp` = m' = product of the primes p | m with p not dividing U_1
    for p in f:
        if U1 % p: mp *= p
    a = alphas_for_family(T1, U1, m)
    for b, ps in SQF:
        pairs += 1
        if b == 1:
            k = 1                                        # Ljunggren/Cohn: square middle only for k = 1 (odd)
        else:
            al = [int(a[pidx[p]]) for p in ps]
            if any(x == 0 for x in al): cnt['no_alpha'] += 1; continue
            if len({v2(x) for x in al}) > 1: cnt['twoadic'] += 1; continue
            k = 1
            for x in al: k = k * x // gcd(k, x)
            if k % 2 == 0: cnt['k_even'] += 1; continue
        if k % mp != 0: cnt['lemmaL'] += 1; continue
        r = form_test(T1, U1, m, k, b)
        cnt[r] += 1
        survivors.append((m, b, k, r))
        if r.startswith('form_ok'): open_pts.append((m, b, k, r))
    if fam % 200 == 0:
        print(f"  ... m = {m}  Familien {fam}  Paare {pairs}  Lemma-L-Ueberlebende {len(survivors)}  offen {len(open_pts)}  {time.time()-t0:.0f}s", flush=True)

print(f"b-Box: m <= {MMAX}, b <= {BMAX} (squarefree: {len(SQF)}), ALLE Hoehen.  Familien: {fam} (davon paritaetstot: {fam_par})  Paare (lebende Familien): {pairs}")
print("Ausschluss-Gruende:", cnt)
print(f"Lemma-L-Ueberlebende (Kandidat k = alpha_m(b) ungerade, m' | k): {len(survivors)}   OFFEN: {len(open_pts)}   Zeit: {time.time()-t0:.0f}s")
print("Ueberlebende (m, b, k, Ergebnis):")
for s in survivors: print("  S", s)
print("OFFEN:", open_pts)
with open(f"w8_bsignatur_M{MMAX}_B{BMAX}_result.json", "w") as fjs:
    json.dump(dict(MMAX=MMAX, H='all', BMAX=BMAX, n_sqfree_b=len(SQF), fam=fam, fam_parity_dead=fam_par, pairs=pairs,
                   counts=cnt, survivors=survivors, open=open_pts, elapsed=time.time()-t0, NQR=NQR), fjs, indent=1)
