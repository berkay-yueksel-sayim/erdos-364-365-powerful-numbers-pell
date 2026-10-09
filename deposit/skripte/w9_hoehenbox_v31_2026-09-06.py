# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w9_hoehenbox_v31_2026-09-06.py
# Purpose: height-box search (`hoehenbox`): for every squarefree kernel m = 7 (mod 8), m <= MMAX, list the lattice points (m, k)
#   with k an odd multiple of m' and digits(T_k) <= H (middles n < 10^H), compute T_k(m) exactly and classify each point: `parity`
#   (T_k odd), `witness` (a prime p <= PMAX = 200000 with v_p(T_k) = 1, so T_k is not powerful) or `open`. All lattice points and
#   all candidate families go into the JSON result. Framework: Mollin and Walsh 1986 (criterion p. 110, lemma p. 111).
# Reads: nothing except its own checkpoint (resume). Writes (current directory): w9_v31_ckpt_M<MMAX>_H<H>.json,
#   w9_v31_M<MMAX>_H<H>_result.json.
# Usage: python w9_hoehenbox_v31_2026-09-06.py MMAX H [CKPT_SEC=900] [STOP_AT=0] (defaults 1000000, 1000; CKPT_SEC = seconds
#   between checkpoints; STOP_AT > 0 = save and stop at the first m >= STOP_AT)
# Controls (positive, before the search; abort on failure): `fund` and `fund_fast` agree on small kernels; T_7(7) =
#   2^3 * 29 * 197 * 2857 (Mollin and Walsh 1986, p. 111) with witness 29; `test_point` classifies constructed numbers correctly.
import sys, time, math, json, os, platform
from math import isqrt
sys.set_int_max_str_digits(2000000)

MMAX     = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000   # largest kernel m
H        = int(sys.argv[2]) if len(sys.argv) > 2 else 1000      # height: digits(T_k) <= H
CKPT_SEC = int(sys.argv[3]) if len(sys.argv) > 3 else 900       # seconds between checkpoints
STOP_AT  = int(sys.argv[4]) if len(sys.argv) > 4 else 0         # stop (after a checkpoint) at the first m >= STOP_AT; 0 = never
PMAX = 200000   # primes up to PMAX serve as witness candidates
CKPT   = f"w9_v31_ckpt_M{MMAX}_H{H}.json"
RESULT = f"w9_v31_M{MMAX}_H{H}_result.json"

# primes up to n (sieve of Eratosthenes)
def sieve(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
PRIMES = sieve(PMAX)

# factorization of n by trial division with the primes up to PMAX (a remaining cofactor > 1 is taken as one factor): dict prime ->
#   exponent
def factor_small(n):
    f = {}
    for p in PRIMES:
        if p * p > n: break
        while n % p == 0:
            f[p] = f.get(p, 0) + 1; n //= p
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

# exact fundamental solution (T1, U1) of x^2 - m*y^2 = 1 from the continued fraction of sqrt(m)
def fund(m):
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1
        k1, k0 = k0, a*k0 + k1
    return h0, k0

# fast fundamental solution: returns (T1 mod m, U1 mod m, log10 T1, period length); the convergent is also carried as a float,
#   rescaled by 1e100 to avoid overflow (`s` counts the removed powers of ten)
def fund_fast(m):
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

def digits(n):
    return len(str(n))

# exact T_k(m): the k-th power of T1 + U1*sqrt(m) by binary exponentiation; (ra, rb) = result, (ba, bb) = current square
def T_exact(T1, U1, m, k):
    ra, rb, ba, bb = 1, 0, T1, U1
    while k:
        if k & 1: ra, rb = ra*ba + m*rb*bb, ra*bb + rb*ba
        ba, bb = ba*ba + m*bb*bb, 2*ba*bb
        k >>= 1
    return ra

# smallest prime p <= PMAX that divides T exactly once (a witness that T is not powerful), or None
def witness(T):
    for p in PRIMES:
        if T % p == 0:
            if (T // p) % p != 0: return p
    return None

# classification of a lattice point: 'parity' (T odd), ('witness', p) or 'open'
def test_point(T):
    if T % 2 == 1: return 'parity'
    w = witness(T)
    return ('witness', w) if w else 'open'

# ---------- Positive controls (abort on failure) ----------
# names: `tm` = T1 mod m, `um` = U1 mod m, `l10` = log10 T1, `per` = period length
assert fund(7) == (8, 3)
tm, um, l10, per = fund_fast(7)
assert (tm, um) == (1, 3) and abs(l10 - math.log10(8)) < 1e-9
tm, um, l10, per = fund_fast(2)
assert (tm, um) == (1, 0) and abs(l10 - math.log10(3)) < 1e-9
for mm in (7, 15, 23, 87, 319, 4999):
    T1, U1 = fund(mm); tm, um, l10, per = fund_fast(mm)
    assert (T1 % mm, U1 % mm) == (tm, um) and abs(math.log10(T1) - l10) < 1e-6, mm
assert test_point(2**3 * 3**2 * 5**2) == 'open'
assert test_point(2**3 * 3**2 * 5**2 * 7) == ('witness', 7)
assert test_point(3**2 * 5**3) == 'parity'
T7 = T_exact(8, 3, 7, 7)
assert T7 == 2**3 * 29 * 197 * 2857 and test_point(T7) == ('witness', 29)
print("PK ok: fund/fund_fast (7,2,15,23,87,319,4999); T_7(7)=2^3*29*197*2857 (M-W 1986 S.111); Testfunktion ok", flush=True)

# ---------- State ----------
# counters: `fam` = squarefree kernels visited, `pts` = lattice points in the box, `dead_par` = parity points, `dead_wit` = points
#   with a witness, `cand_fam` = candidate families (m' small enough for at least one point), `max_period` = longest continued
#   fraction period seen; lists: `hits` = witness points (m, m', k, digits, p), `open_pts` = open points (m, k, digits),
#   `points` = all lattice points, `cands` = all candidate families
fam = pts = dead_par = dead_wit = cand_fam = max_period = 0
hits, open_pts, points, cands = [], [], [], []
elapsed_prev = 0.0; start_m = 7   # run time already spent (from a checkpoint); first kernel of the class

def summary():
    print(f"Box: m <= {MMAX}, n < 10^{H}   Familien(squarefree, 7 mod 8): {fam}   Kandidaten-Familien (exakt gerechnet): {cand_fam}   max Periode: {max_period}")
    print(f"Gitterpunkte in der Box: {pts}   tot(Paritaet): {dead_par}   tot(Zeuge): {dead_wit}   OFFEN: {len(open_pts)}")

def state():
    return dict(MMAX=MMAX, H=H, fam=fam, pts=pts, dead_par=dead_par, dead_wit=dead_wit, cand_fam=cand_fam,
                max_period=max_period, hits=hits, open_pts=open_pts, points=points, cands=cands)

# resume from the checkpoint if one exists for the same MMAX and H
if os.path.exists(CKPT):
    d = json.load(open(CKPT))
    assert d["MMAX"] == MMAX and d["H"] == H, "Checkpoint passt nicht zu MMAX/H"
    fam, pts, dead_par, dead_wit = d["fam"], d["pts"], d["dead_par"], d["dead_wit"]
    cand_fam, max_period = d["cand_fam"], d["max_period"]
    hits = [tuple(h) for h in d["hits"]]; open_pts = [tuple(o) for o in d["open_pts"]]
    points = [tuple(p) for p in d["points"]]; cands = [tuple(c) for c in d["cands"]]
    elapsed_prev = d["elapsed"]; start_m = d["next_m"]
    if d.get("done"):
        print(f"Checkpoint {CKPT}: bereits fertig (gespeichert {d['saved']})."); summary(); sys.exit(0)
    print(f"RESUME ab m = {start_m}  (Checkpoint {d['saved']}, {elapsed_prev:.0f}s, {pts} Punkte, {len(open_pts)} offen)", flush=True)

t0 = time.time(); last_ckpt = t0

# write the checkpoint atomically (temporary file, then rename); `next_m` = the next kernel to process
def save_ckpt(next_m, done=False):
    d = state(); d.update(next_m=next_m, elapsed=elapsed_prev + time.time() - t0, done=done,
                          saved=time.strftime("%Y-%m-%d %H:%M:%S"))
    tmp = CKPT + ".tmp"
    with open(tmp, "w") as f: json.dump(d, f)
    os.replace(tmp, CKPT)

# ---------- Run ----------
# per squarefree kernel m: `mp` = m' = product of the primes p | m with p not dividing U1; `kmax` = largest k that can still fit
#   into the box (digits(T_k) <= H); the family is a candidate only if m' <= kmax
for m in range(start_m, MMAX + 1, 8):
    if time.time() - last_ckpt >= CKPT_SEC:
        save_ckpt(m); last_ckpt = time.time()
        print(f"  [ckpt] next m = {m}  Familien {fam}  Punkte {pts}  offen {len(open_pts)}  {elapsed_prev + time.time()-t0:.0f}s", flush=True)
    if STOP_AT and m >= STOP_AT:
        save_ckpt(m); print(f"STOP_AT {STOP_AT} erreicht, Checkpoint bei next m = {m}.", flush=True); sys.exit(0)
    if m % 5000000 == 7 and m > 7:
        print(f"  ... m = {m}  Familien {fam}  Punkte {pts}  offen {len(open_pts)}  {elapsed_prev + time.time()-t0:.0f}s", flush=True)
    f = factor_small(m)
    if any(e > 1 for e in f.values()): continue
    fam += 1
    tm, um, l10, per = fund_fast(m)
    if per > max_period: max_period = per
    mp = 1
    for p in f:
        if um % p: mp *= p
    kmax = int((H + 0.31) / l10) + 1
    if mp > kmax: continue
    cand_fam += 1
    T1, U1 = fund(m)
    assert (T1 % m, U1 % m) == (tm, um) and abs(math.log10(T1) - l10) < 1e-6, ("Cross-Check", m)
    cands.append((m, mp, digits(T1), round(l10, 6), kmax, per, 'T1odd' if T1 % 2 else 'T1even'))
    # lattice points of this family: the odd multiples k of m' up to kmax
    k = mp
    while k <= kmax:
        if k % 2 == 1:
            T = T_exact(T1, U1, m, k)
            dg = digits(T)
            if dg <= H:
                pts += 1
                r = test_point(T)
                if r == 'parity':
                    dead_par += 1; points.append((m, mp, k, dg, 'parity', 0))
                elif r == 'open':
                    open_pts.append((m, k, dg)); points.append((m, mp, k, dg, 'open', 0))
                else:
                    dead_wit += 1
                    hits.append((m, mp, k, dg, r[1])); points.append((m, mp, k, dg, 'witness', r[1]))
        k += mp

save_ckpt(MMAX + 8, done=True)
summary()
print(f"Zeit gesamt: {elapsed_prev + time.time()-t0:.1f}s")
print("Zeugen-Punkte (m, m', k, Stellen, Zeuge):")
for h in hits: print("  ", h)
print("OFFEN:", open_pts)
print("Kandidaten-Familien (m, m', Stellen T1, log10 T1, kmax, Periode, Paritaet T1):")
for c in cands: print("  C", c)
print("Alle Gitterpunkte (m, m', k, Stellen, Status, Zeuge):")
for p in points: print("  P", p)
res = state(); res.update(elapsed=elapsed_prev + time.time() - t0, finished=time.strftime("%Y-%m-%d %H:%M:%S"),
                          python=platform.python_version(), script=os.path.basename(__file__), PMAX=PMAX,
                          schema=dict(cands="(m, m_strich, digits_T1, log10_T1, kmax, cf_period, T1_parity)",
                                      points="(m, m_strich, k, digits_Tk, status, witness_or_0)",
                                      note="status in {parity, witness, open}; box = kern m <= MMAX and digits(T_k) <= H"))
with open(RESULT, "w") as fjs: json.dump(res, fjs, indent=1)
print("Deposit-JSON:", RESULT)
