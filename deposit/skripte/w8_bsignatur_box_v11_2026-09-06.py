# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w8_bsignatur_box_v11_2026-09-06.py
# W8 v1.1, "b-signature box": for the squarefree kernels m <= MMAX (m = 7 mod 8, living families, i.e. T_1 even) and every
#   squarefree b <= BMAX, for ALL heights, test whether the middle n = T_k(m) can have squarefree part b (n powerful with
#   squarefree part b <=> n = b^3 z^2). Per pair (m, b) the filters below exclude the case or leave a survivor; survivors that
#   the form test does not exclude are listed as open.
# This version adds checkpoint/resume (every CKPT_SEC seconds, atomic; the same call continues) and a result JSON to the logic of
#   w8_bsignatur_box_2026-09-06.py (same logic and filters; the result for MMAX = BMAX = 10^4 must be identical).
# Foundations (credits): Mollin-Walsh 1986 (criterion p. 110, lemma p. 111: m | U_k <=> m' | k);
#   Bennett-Walsh 1999 Lemma 3.3 ({k : p | T_k} = odd multiples of alpha(p)) and Thm 1.2 (T_k = b x^2 at most for k = alpha(b),
#   b > 1 squarefree); b = 1: Ljunggren 1942 / Cohn 1997 (k in {1,2}). Own code.
# Reads: nothing, except the checkpoint file if it exists. Writes (current directory):
#   w8_v11_ckpt_M<MMAX>_B<BMAX>.json (checkpoint) and w8_v11_M<MMAX>_B<BMAX>_result.json.
# Usage: python w8_bsignatur_box_v11_2026-09-06.py MMAX BMAX [CKPT_SEC=900] [STOP_AT=0]     (defaults 10000 10000 900 0)
#   STOP_AT: if nonzero, save a checkpoint and stop once m >= STOP_AT (used to test the resume).
# Controls (abort on failure, block "positive controls" below): T_7(7) = 2^3 * 29 * 197 * 2857 (Mollin-Walsh 1986, p. 111);
#   alpha(11) = alpha(23) = 3 and alpha(29) = alpha(197) = alpha(2857) = 7 for kernel 7; Bennett-Walsh Lemma 3.3 at p = 11;
#   three examples for the form test.
import sys, math, time, json, os
from math import isqrt, gcd
import numpy as np
sys.set_int_max_str_digits(2000000)

MMAX     = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
BMAX     = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
CKPT_SEC = int(sys.argv[3]) if len(sys.argv) > 3 else 900
STOP_AT  = int(sys.argv[4]) if len(sys.argv) > 4 else 0
NQR  = 40   # `NQR` = number of primes q > 1000 used in the modular square test (`QR`)
# `CKPT` = checkpoint file, `RESULT` = result file
CKPT   = f"w8_v11_ckpt_M{MMAX}_B{BMAX}.json"
RESULT = f"w8_v11_M{MMAX}_B{BMAX}_result.json"

def sieve(n):                      # list of the primes up to n
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

def spf_table(n):                  # table of the smallest prime factor of every number up to n
    spf = list(range(n + 1))
    for i in range(2, isqrt(n) + 1):
        if spf[i] == i:
            for j in range(i*i, n + 1, i):
                if spf[j] == j: spf[j] = i
    return spf

def factor_small(n, primes):        # factorization {prime: exponent} by trial division with `primes`
    f = {}
    for p in primes:
        if p * p > n: break
        while n % p == 0:
            f[p] = f.get(p, 0) + 1; n //= p
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def fund(m):                        # fundamental solution (T1, U1) of x^2 - m*y^2 = 1 by continued fraction
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1
        k1, k0 = k0, a*k0 + k1
    return h0, k0

def eps_pow_mod(T1, U1, m, k, M):  # T_k modulo M of the unit T1 + U1*sqrt(m) (square-and-multiply)
    ra, rb = 1 % M, 0; ba, bb = T1 % M, U1 % M
    while k:
        if k & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        k >>= 1
    return ra

def T_exact(T1, U1, m, k):         # exact T_k
    ra, rb, ba, bb = 1, 0, T1, U1
    while k:
        if k & 1: ra, rb = ra*ba + m*rb*bb, ra*bb + rb*ba
        ba, bb = ba*ba + m*bb*bb, 2*ba*bb
        k >>= 1
    return ra

def v2(n):                         # 2-adic valuation
    return (n & -n).bit_length() - 1

def is_square(n):
    if n < 0: return False
    r = isqrt(n); return r * r == n

# `PRIMES_ALL`/`PRIMES_B` = primes up to max(BMAX, 200000) / up to BMAX, `Pn` = PRIMES_B as array, `pidx` = index of each prime,
# `QR` = the primes > 1000 for the square test, `SPF` = smallest-prime-factor table
PRIMES_ALL = sieve(max(BMAX, 200000))
PRIMES_B = [p for p in PRIMES_ALL if p <= BMAX]
Pn = np.array(PRIMES_B, dtype=np.int64)
pidx = {p: i for i, p in enumerate(PRIMES_B)}
QR = [p for p in PRIMES_ALL if 1000 < p][:NQR]
SPF = spf_table(BMAX)
def sqfree_list(B):                # all squarefree b <= B together with their prime factors
    out = []
    for b in range(1, B + 1):
        n, ps, ok = b, [], True
        while n > 1:
            p = SPF[n]; n //= p
            if n % p == 0: ok = False; break
            ps.append(p)
        if ok: out.append((b, ps))
    return out
SQF = sqfree_list(BMAX)   # `SQF` = list of (b, prime factors of b)

# `alphas_for_family` = for every prime p <= BMAX the rank of apparition alpha(p) (first k with p | T_k, via the recurrence
# t_{k+1} = 2*T1*t_k - t_{k-1} mod p; 0 if not found)
def alphas_for_family(T1, U1, m):
    T1m = np.array([T1 % p for p in PRIMES_B], dtype=np.int64)
    c = (2 * T1m) % Pn
    a = np.zeros(len(PRIMES_B), dtype=np.int64)
    t_prev = np.ones(len(PRIMES_B), dtype=np.int64)
    t_cur = T1m.copy()
    kmax_iter = int(Pn.max()) + 2
    for k in range(1, kmax_iter + 1):
        z = (t_cur == 0) & (a == 0)
        if z.any(): a[z] = k
        t_prev, t_cur = t_cur, np.mod(c * t_cur - t_prev, Pn)
    return a

def legendre(a, p):                # Euler criterion value a^((p-1)/2) mod p (1, p-1, or 0)
    return pow(a % p, (p - 1) // 2, p)

# `form_test`: can T_k(m) be b^3 * z^2? Returns 'b3_fail' (b^3 does not divide T_k), 'nonsquare' (T_k/b^3 is a quadratic
# non-residue modulo some q in QR), 'nonsquare_exact', 'form_ok_exact' (exact check passed), or 'form_ok_modular_open'
# (T_k too large for the exact check)
def form_test(T1, U1, m, k, b):
    b3 = b ** 3
    if eps_pow_mod(T1, U1, m, k, b3) != 0: return 'b3_fail'
    for q in QR:
        if b % q == 0 or m % q == 0: continue
        t = eps_pow_mod(T1, U1, m, k, q)
        if t == 0: continue
        r = (t * pow(b3 % q, q - 2, q)) % q
        if legendre(r, q) == q - 1: return 'nonsquare'
    if k * math.log10(T1 + 1) < 200000:
        T = T_exact(T1, U1, m, k)
        return 'form_ok_exact' if (T % b3 == 0 and is_square(T // b3)) else 'nonsquare_exact'
    return 'form_ok_modular_open'

# ---------- positive controls (abort on failure) ----------
T1, U1 = fund(7); assert (T1, U1) == (8, 3)
assert T_exact(8, 3, 7, 7) == 2**3 * 29 * 197 * 2857
a7 = alphas_for_family(8, 3, 7)
assert a7[pidx[2]] == 1 and a7[pidx[11]] == 3 and a7[pidx[23]] == 3
assert a7[pidx[29]] == 7 and a7[pidx[197]] == 7 and a7[pidx[2857]] == 7
assert a7[pidx[7]] == 0
assert eps_pow_mod(8, 3, 7, 9, 11) == 0 and eps_pow_mod(8, 3, 7, 6, 11) != 0 and eps_pow_mod(8, 3, 7, 15, 11) == 0
assert form_test(8, 3, 7, 1, 2) == 'form_ok_exact'
assert form_test(8, 3, 7, 7, 2) in ('nonsquare', 'nonsquare_exact')
assert form_test(8, 3, 7, 2, 127) == 'b3_fail'
print("PK ok: T_7(7)=2^3*29*197*2857 (M-W 1986 S.111); alpha(11)=alpha(23)=3, alpha(29)=alpha(197)=alpha(2857)=7; B-W Lemma 3.3 an p=11; Form-Test ok", flush=True)

# ---------- state (new or from checkpoint) ----------
# `fam` = families, `fam_par` = families dead by parity (T1 odd), `pairs` = pairs (m, b) examined, `cnt` = number of
# exclusions per reason, `survivors` = (m, b, k, result) not excluded by the filters, `open_pts` = those with result 'form_ok_*'
fam = fam_par = pairs = 0
cnt = dict(no_alpha=0, twoadic=0, k_even=0, lemmaL=0, b3_fail=0, nonsquare=0, nonsquare_exact=0, form_ok_exact=0, form_ok_modular_open=0)
survivors, open_pts = [], []
elapsed_prev = 0.0; start_m = 7

def summary():                     # print the totals
    print(f"b-Box: m <= {MMAX}, b <= {BMAX} (squarefree: {len(SQF)}), ALLE Hoehen.  Familien: {fam} (davon paritaetstot: {fam_par})  Paare (lebende Familien): {pairs}")
    print("Ausschluss-Gruende:", cnt)
    print(f"Lemma-L-Ueberlebende (Kandidat k = alpha_m(b) ungerade, m' | k): {len(survivors)}   OFFEN: {len(open_pts)}")

def state():                       # the state as a dict (written to the checkpoint and to the result JSON)
    return dict(MMAX=MMAX, BMAX=BMAX, n_sqfree_b=len(SQF), fam=fam, fam_parity_dead=fam_par, pairs=pairs, counts=cnt,
                survivors=survivors, open=open_pts, NQR=NQR)

if os.path.exists(CKPT):
    d = json.load(open(CKPT))
    assert d["MMAX"] == MMAX and d["BMAX"] == BMAX, "Checkpoint passt nicht zu MMAX/BMAX"
    fam, fam_par, pairs, cnt = d["fam"], d["fam_parity_dead"], d["pairs"], d["counts"]
    survivors = [tuple(s) for s in d["survivors"]]; open_pts = [tuple(o) for o in d["open"]]
    elapsed_prev = d["elapsed"]; start_m = d["next_m"]
    if d.get("done"):
        print(f"Checkpoint {CKPT}: bereits fertig (gespeichert {d['saved']})."); summary(); sys.exit(0)
    print(f"RESUME ab m = {start_m}  (Checkpoint {d['saved']}, {elapsed_prev:.0f}s, {len(survivors)} Ueberlebende, {len(open_pts)} offen)", flush=True)

t0 = time.time(); last_ckpt = t0

def save_ckpt(next_m, done=False):  # atomic checkpoint: write a temporary file, then rename
    d = state(); d.update(next_m=next_m, elapsed=elapsed_prev + time.time() - t0, done=done,
                          saved=time.strftime("%Y-%m-%d %H:%M:%S"))
    tmp = CKPT + ".tmp"
    with open(tmp, "w") as f: json.dump(d, f)
    os.replace(tmp, CKPT)

# ---------- run ----------
for m in range(start_m, MMAX + 1, 8):
    if time.time() - last_ckpt >= CKPT_SEC:
        save_ckpt(m); last_ckpt = time.time()
        print(f"  [ckpt] next m = {m}  Familien {fam}  Ueberlebende {len(survivors)}  offen {len(open_pts)}  {elapsed_prev + time.time()-t0:.0f}s", flush=True)
    if STOP_AT and m >= STOP_AT:
        save_ckpt(m); print(f"STOP_AT {STOP_AT} erreicht, Checkpoint bei next m = {m}.", flush=True); sys.exit(0)
    f = factor_small(m, PRIMES_ALL)
    if any(e > 1 for e in f.values()): continue
    fam += 1
    T1, U1 = fund(m)
    if T1 % 2 == 1:
        fam_par += 1; continue
    mp = 1                            # `mp` = m'
    for p in f:
        if U1 % p: mp *= p
    a = alphas_for_family(T1, U1, m)
    for b, ps in SQF:
        pairs += 1
        if b == 1:
            k = 1
        else:
            al = [int(a[pidx[p]]) for p in ps]       # `al` = ranks alpha(p) of the primes p | b
            if any(x == 0 for x in al): cnt['no_alpha'] += 1; continue      # some p | b never divides a T_k
            if len({v2(x) for x in al}) > 1: cnt['twoadic'] += 1; continue  # the alpha(p) have different 2-adic valuations
            k = 1
            for x in al: k = k * x // gcd(k, x)
            if k % 2 == 0: cnt['k_even'] += 1; continue                     # k = lcm of the alpha(p) must be odd
        if k % mp != 0: cnt['lemmaL'] += 1; continue                        # requires m' | k
        r = form_test(T1, U1, m, k, b)
        cnt[r] += 1
        survivors.append((m, b, k, r))
        if r.startswith('form_ok'): open_pts.append((m, b, k, r))
    if fam % 1000 == 0:
        print(f"  ... m = {m}  Familien {fam}  Paare {pairs}  Ueberlebende {len(survivors)}  offen {len(open_pts)}  {elapsed_prev + time.time()-t0:.0f}s", flush=True)

save_ckpt(MMAX + 8, done=True)
summary()
print(f"Zeit gesamt: {elapsed_prev + time.time()-t0:.0f}s")
print("Ueberlebende (m, b, k, Ergebnis):")
for s in survivors: print("  S", s)
print("OFFEN:", open_pts)
res = state(); res.update(elapsed=elapsed_prev + time.time() - t0, finished=time.strftime("%Y-%m-%d %H:%M:%S"),
                          script=os.path.basename(__file__))
with open(RESULT, "w") as fjs: json.dump(res, fjs, indent=1)
print("Result-JSON:", RESULT)
