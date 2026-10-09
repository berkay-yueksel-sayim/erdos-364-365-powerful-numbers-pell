# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w9_hoehenbox_v3_ckpt_2026-09-03.py
# Height box: all lattice points (m, k) with squarefree kernel m = 7 (mod 8), m <= MMAX, and n < 10^H, with checkpoint/resume.
#
# For every such kernel the modular route (T_1 mod m, U_1 mod m and log10 T_1 from the continued fraction) gives m' and the
# largest k for which T_k can still have at most H digits; kernels with m' > kmax are skipped. For the remaining candidate
# kernels the fundamental solution is computed exactly (cross-check against the modular route), and for every odd multiple k
# of m' up to kmax the exact T_k is classified: 'parity' (T_k odd), 'witness' (a prime p <= PMAX with v_p(T_k) = 1), or 'open'.
#
# Reads:    the checkpoint file, if it exists.
# Writes:   w9_v3_ckpt_M<MMAX>_H<H>.json in the working directory (saved atomically: temporary file, then replace),
#           and printed output: summary, the points killed by a witness (m, m', k, digits, witness) and the OPEN points.
# Usage:    python w9_hoehenbox_v3_ckpt_2026-09-03.py MMAX H [CKPT_SEC=900] [STOP_AT=0]
#           (defaults MMAX = 1000000, H = 1000). The same call resumes from the checkpoint; CKPT_SEC = seconds between saves;
#           STOP_AT is only for testing the resume: the run stops after a save once m >= STOP_AT.
# Controls: positive controls at the start (abort on failure): fund/fund_fast agree for m = 7, 2, 15, 23, 87, 319, 4999;
#           the test function classifies known cases correctly; T_7(7) = 2^3 * 29 * 197 * 2857 (Mollin-Walsh 1986, p. 111) has
#           witness 29. Every candidate kernel is cross-checked (T_1 mod m, U_1 mod m, log10 T_1) against the exact computation.
# Framework: Mollin-Walsh 1986 (criterion p. 110, lemma p. 111). Own code.
import sys, time, math, json, os
from math import isqrt

MMAX     = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000   # bound for the kernel m
H        = int(sys.argv[2]) if len(sys.argv) > 2 else 1000      # digit bound: points with n < 10^H
CKPT_SEC = int(sys.argv[3]) if len(sys.argv) > 3 else 900       # seconds between checkpoint saves
# stop after a save once m >= STOP_AT (0 = never; resume test only)
STOP_AT  = int(sys.argv[4]) if len(sys.argv) > 4 else 0
# primes up to PMAX are used for factoring m and for the witness search
PMAX = 200000
CKPT = f"w9_v3_ckpt_M{MMAX}_H{H}.json"                          # `CKPT` = checkpoint file name

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
    # `fund_fast` = modular route: runs the continued fraction up to the end of an even period, keeping the convergents only
    # modulo m and as a rescaled float (offset kept in `s`); returns (T_1 mod m, U_1 mod m, log10 T_1, i - 1) where i - 1 is the
    # index of the convergent used
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    hm1, hm0 = 1 % m, a0 % m
    km1, km0 = 0, 1 % m
    hf1, hf0 = 1.0, float(a0)
    s = 0; i = 0
    while True:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        i += 1
        if Q == 1 and i % 2 == 0:
            return hm0, km0, math.log10(hf0) + s, i - 1
        hm1, hm0 = hm0, (a*hm0 + hm1) % m
        km1, km0 = km0, (a*km0 + km1) % m
        hf1, hf0 = hf0, a*hf0 + hf1
        if hf0 > 1e100:
            hf0 /= 1e100; hf1 /= 1e100; s += 100

def T_exact(T1, U1, m, k):   # exact T_k from (T_1 + U_1 sqrt(m))^k by square-and-multiply
    ra, rb, ba, bb = 1, 0, T1, U1
    while k:
        if k & 1: ra, rb = ra*ba + m*rb*bb, ra*bb + rb*ba
        ba, bb = ba*ba + m*bb*bb, 2*ba*bb
        k >>= 1
    return ra

def witness(T):   # smallest prime p <= PMAX that divides T exactly once (v_p(T) = 1), or None
    for p in PRIMES:
        if T % p == 0:
            if (T // p) % p != 0: return p
    return None

def test_point(T):   # classify a value T_k: 'parity' (odd), ('witness', p), or 'open'
    if T % 2 == 1: return 'parity'
    w = witness(T)
    return ('witness', w) if w else 'open'

# ---------- Positive controls (abort on failure) ----------
assert fund(7) == (8, 3)
tm, um, l10, per = fund_fast(7);  assert (tm, um) == (1, 3) and abs(l10 - math.log10(8)) < 1e-9
tm, um, l10, per = fund_fast(2);  assert (tm, um) == (1, 0) and abs(l10 - math.log10(3)) < 1e-9
for mm in (7, 15, 23, 87, 319, 4999):
    T1, U1 = fund(mm); tm, um, l10, per = fund_fast(mm)
    assert (T1 % mm, U1 % mm) == (tm, um) and abs(math.log10(T1) - l10) < 1e-6, mm
assert test_point(2**3 * 3**2 * 5**2) == 'open'
assert test_point(2**3 * 3**2 * 5**2 * 7) == ('witness', 7)
assert test_point(3**2 * 5**3) == 'parity'
T7 = T_exact(8, 3, 7, 7)
assert T7 == 2**3 * 29 * 197 * 2857 and test_point(T7) == ('witness', 29)
print("PK ok: fund/fund_fast stimmen (7,2,15,23,87,319,4999); T_7(7)=2^3*29*197*2857 (M-W 1986 S.111); Testfunktion ok", flush=True)

# ---------- State (new, or restored from a checkpoint) ----------
# `fam` = kernels (squarefree, 7 mod 8) seen; `pts` = lattice points in the box; `dead_par` / `dead_wit` = points excluded by
# parity / by a witness; `cand_fam` = candidate kernels (computed exactly); `max_period` = largest `per` seen;
# `hits` = points with a witness (m, m', k, digits, p); `open_pts` = open points (m, k, digits);
# `elapsed_prev` = seconds spent before the last resume; `start_m` = next kernel to process
fam = pts = dead_par = dead_wit = cand_fam = max_period = 0
hits, open_pts = [], []
elapsed_prev = 0.0; start_m = 7

def summary():
    print(f"Box: m <= {MMAX}, n < 10^{H}   Familien(squarefree, 7 mod 8): {fam}   Kandidaten-Familien (exakt gerechnet): {cand_fam}   max Periode: {max_period}")
    print(f"Gitterpunkte in der Box: {pts}   tot(Paritaet): {dead_par}   tot(Zeuge): {dead_wit}   OFFEN: {len(open_pts)}")

if os.path.exists(CKPT):
    d = json.load(open(CKPT))
    assert d["MMAX"] == MMAX and d["H"] == H, "Checkpoint passt nicht zu MMAX/H"
    fam, pts, dead_par, dead_wit = d["fam"], d["pts"], d["dead_par"], d["dead_wit"]
    cand_fam, max_period = d["cand_fam"], d["max_period"]
    hits = [tuple(h) for h in d["hits"]]; open_pts = [tuple(o) for o in d["open_pts"]]
    elapsed_prev = d["elapsed"]; start_m = d["next_m"]
    if d.get("done"):
        print(f"Checkpoint {CKPT} sagt: bereits fertig (gespeichert {d['saved']}).")
        summary(); print(f"Zeit gesamt: {elapsed_prev:.1f}s"); sys.exit(0)
    print(f"RESUME ab m = {start_m}  (Checkpoint gespeichert {d['saved']}, {elapsed_prev:.0f}s bisher, {pts} Punkte, {len(open_pts)} offen)", flush=True)

t0 = time.time(); last_ckpt = t0

def save_ckpt(next_m, done=False):   # write the state atomically (temporary file, then os.replace)
    d = dict(MMAX=MMAX, H=H, next_m=next_m, fam=fam, pts=pts, dead_par=dead_par, dead_wit=dead_wit,
             cand_fam=cand_fam, max_period=max_period, hits=hits, open_pts=open_pts,
             elapsed=elapsed_prev + time.time() - t0, done=done, saved=time.strftime("%Y-%m-%d %H:%M:%S"))
    tmp = CKPT + ".tmp"
    with open(tmp, "w") as f: json.dump(d, f)
    os.replace(tmp, CKPT)

# ---------- Run ----------
for m in range(start_m, MMAX + 1, 8):
    if time.time() - last_ckpt >= CKPT_SEC:
        save_ckpt(m); last_ckpt = time.time()
        print(f"  [ckpt] next m = {m}  Familien {fam}  Punkte {pts}  offen {len(open_pts)}  {elapsed_prev + time.time()-t0:.0f}s", flush=True)
    if STOP_AT and m >= STOP_AT:
        save_ckpt(m); print(f"STOP_AT {STOP_AT} erreicht, Checkpoint bei next m = {m} gespeichert.", flush=True); sys.exit(0)
    if m % 5000000 == 7 and m > 7:
        print(f"  ... m = {m}  Familien {fam}  Punkte {pts}  offen {len(open_pts)}  {elapsed_prev + time.time()-t0:.0f}s", flush=True)
    f = factor_small(m)
    if any(e > 1 for e in f.values()): continue   # m not squarefree
    fam += 1
    tm, um, l10, per = fund_fast(m)
    if per > max_period: max_period = per
    mp = 1   # `mp` = m' = product of the primes p | m with p not dividing U_1
    for p in f:
        if um % p: mp *= p
    kmax = int((H + 0.31) / l10) + 1   # largest k for which T_k can still have at most H digits (T_k >= T_1^k / 2)
    if mp > kmax: continue
    cand_fam += 1
    T1, U1 = fund(m)
    assert (T1 % m, U1 % m) == (tm, um) and abs(math.log10(T1) - l10) < 1e-6, ("Cross-Check", m)
    k = mp
    while k <= kmax:
        if k % 2 == 1:
            T = T_exact(T1, U1, m, k)
            if len(str(T)) <= H:
                pts += 1
                r = test_point(T)
                if r == 'parity': dead_par += 1
                elif r == 'open': open_pts.append((m, k, len(str(T))))
                else:
                    dead_wit += 1; hits.append((m, mp, k, len(str(T)), r[1]))
        k += mp

save_ckpt(MMAX + 8, done=True)
summary()
print(f"Zeit gesamt: {elapsed_prev + time.time()-t0:.1f}s   (Checkpoint final: {CKPT})")
print("Familien mit Punkten in der Box (m, m', k, Stellen, Zeuge):")
for h in hits: print("  ", h)
print("OFFEN:", open_pts)
