# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w10_steckbrief_scan_2026-09-06.py
# Purpose: "profile scan" (`steckbrief` = profile): height box (middle n = T_k(m) < 10^H) for a THIN set of large kernels:
#   squarefree m = 7 (mod 8) in the range (MLO, MHI] whose prime factors are ALL <= PSMOOTH and that have >= TMIN prime factors.
#   Motivation: 380 of 445 points of w9 lie in families with all prime factors <= 100; many prime factors => large 2-class group
#   => small unit; and small prime factors divide U_1 with probability ~2/p, which makes m' small.
# Rigorous shortcut: the convergent numerators satisfy A_j >= F_{j+2} (Fibonacci), so T_1 = A_{l-1} >= phi^(l-1); if the period l
#   exceeds the bound L = (H+0.31)/log10(phi) + 2, then T_1 > 10^H and the family has no point in the box (T_k >= T_1 for k >= 1).
#   Such families are skipped without computing the unit (counted as `unit_big`).
# Frame: Mollin-Walsh 1986 (criterion p. 110, lemma p. 111).
# Usage: python w10_steckbrief_scan_2026-09-06.py MLO MHI H PSMOOTH TMIN [CKPT_SEC=900]   (defaults 6, 1e8, 2000, 100, 1, 900)
# Writes (current directory): checkpoint w10_ckpt_<TAG>.json (resumed automatically if present; every CKPT_SEC seconds) and
#   result w10_<TAG>_result.json, with TAG = L<MLO>_U<MHI>_H<H>_P<PSMOOTH>_T<TMIN>.
# Controls: positive controls run first and abort on failure. Control mode: MLO=6, MHI=1e8 must give exactly the 100-smooth
#   subset of ergebnisse/w9_v31_M100000000_H2000_result.json (compared outside this script).

import sys, math, time, json, os
from math import isqrt
sys.set_int_max_str_digits(2000000)

MLO      = int(float(sys.argv[1])) if len(sys.argv) > 1 else 6
MHI      = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10**8
H        = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
PSMOOTH  = int(sys.argv[4]) if len(sys.argv) > 4 else 100
TMIN     = int(sys.argv[5]) if len(sys.argv) > 5 else 1
CKPT_SEC = int(sys.argv[6]) if len(sys.argv) > 6 else 900
PMAX = 200000
TAG = f"L{MLO}_U{MHI}_H{H}_P{PSMOOTH}_T{TMIN}"
CKPT = f"w10_ckpt_{TAG}.json"; RESULT = f"w10_{TAG}_result.json"
LOG10PHI = math.log10((1 + 5 ** 0.5) / 2)
LMAX = int((H + 0.31) / LOG10PHI) + 2  # period bound (rigorous, see above)

# `PRIMES` = primes up to `PMAX` (used for the witness search); `SMALL` = primes <= PSMOOTH
def sieve(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
PRIMES = sieve(PMAX)
SMALL = [p for p in PRIMES if p <= PSMOOTH]

# `fund` = exact fundamental solution (T1, U1) of x^2 - m*y^2 = 1 by continued fraction (big integers)
def fund(m):
    a0 = isqrt(m); P, Q, a = 0, 1, a0
    h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1
        k1, k0 = k0, a*k0 + k1
    return h0, k0

def fund_fast_capped(m, lmax):
    """(T1 mod m, U1 mod m, log10 T1, period) or None if the period exceeds lmax."""
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
        if i > lmax: return None
        hm1, hm0 = hm0, (a*hm0 + hm1) % m
        km1, km0 = km0, (a*km0 + km1) % m
        hf1, hf0 = hf0, a*hf0 + hf1
        if hf0 > 1e100:
            hf0 /= 1e100; hf1 /= 1e100; s += 100

# `T_exact` = exact T_k by binary powering of T1 + U1*sqrt(m)
def T_exact(T1, U1, m, k):
    ra, rb, ba, bb = 1, 0, T1, U1
    while k:
        if k & 1: ra, rb = ra*ba + m*rb*bb, ra*bb + rb*ba
        ba, bb = ba*ba + m*bb*bb, 2*ba*bb
        k >>= 1
    return ra

def witness(T):
    for p in PRIMES:
        if T % p == 0:
            if (T // p) % p != 0: return p
    return None

# `test_point` = status of T: 'parity' if T is odd, ('witness', p) if a prime p <= PMAX divides T exactly once, else 'open'
def test_point(T):
    if T % 2 == 1: return 'parity'
    w = witness(T)
    return ('witness', w) if w else 'open'

def enumerate_kernels():
    """All squarefree m = 7 (mod 8) in (MLO, MHI] with prime factors in SMALL and >= TMIN factors (DFS)."""
    out = []
    def dfs(idx, prod, cnt):
        if prod > MLO and cnt >= TMIN and prod % 8 == 7: out.append(prod)
        for j in range(idx, len(SMALL)):
            p = SMALL[j]
            if prod * p > MHI: break
            dfs(j + 1, prod * p, cnt + 1)
    dfs(0, 1, 0)
    out.sort()
    return out

# ---------- positive controls (abort on failure) ----------
assert fund(7) == (8, 3)
# the capped version gives (T1 mod 7, U1 mod 7) = (1, 3) and log10 T1 = log10(8) at m = 7
r = fund_fast_capped(7, LMAX); assert r is not None and (r[0], r[1]) == (1, 3) and abs(r[2] - math.log10(8)) < 1e-9
assert fund_fast_capped(7, 2) is None, "Cap muss greifen (Periode von sqrt7 ist 4)"
# T_7(7) = 2^3 * 29 * 197 * 2857 (Mollin-Walsh), witness 29
T7 = T_exact(8, 3, 7, 7); assert T7 == 2**3 * 29 * 197 * 2857 and test_point(T7) == ('witness', 29)
# status classification: 2^3 * 3^2 * 5^2 is 'open', 3^2 * 5^3 is 'parity'
assert test_point(2**3 * 3**2 * 5**2) == 'open' and test_point(3**2 * 5**3) == 'parity'
# Fibonacci bound on examples: log10 T1 >= (l-1)*log10(phi)
for mm in (7, 15, 23, 87, 319, 4999, 123879, 16057223):
    T1, U1 = fund(mm); rr = fund_fast_capped(mm, 10**7)
    assert rr is not None and abs(math.log10(T1) - rr[2]) < 1e-6 and (T1 % mm, U1 % mm) == (rr[0], rr[1]), mm
    assert math.log10(T1) >= (rr[3] - 1) * LOG10PHI - 1e-9, ("Fibonacci-Schranke verletzt", mm)
print(f"PK ok: fund/fund_fast_capped stimmen (8 Kerne); Fibonacci-Schranke haelt; Cap greift; T_7(7) = 2^3*29*197*2857 (M-W 1986). LMAX = {LMAX}", flush=True)

kernels = enumerate_kernels()
print(f"Steckbrief-Menge: {len(kernels)} Kerne in ({MLO}, {MHI}], Primfaktoren <= {PSMOOTH}, >= {TMIN} Faktoren", flush=True)

# ---------- state ----------
# `fam` = families processed, `unit_big` = unit too big for the box, `cand` = candidates, `pts` = lattice points,
# `dead_par` / `dead_wit` = dead by parity / by witness; `points`, `cands`, `open_pts` = their lists (`open_pts` = open points)
i0 = 0; fam = unit_big = cand = pts = dead_par = dead_wit = 0
points, cands, open_pts = [], [], []
elapsed_prev = 0.0
# resume from the checkpoint file if it exists
if os.path.exists(CKPT):
    d = json.load(open(CKPT))
    assert d["TAG"] == TAG
    i0, fam, unit_big, cand, pts, dead_par, dead_wit = d["i"], d["fam"], d["unit_big"], d["cand"], d["pts"], d["dead_par"], d["dead_wit"]
    points = [tuple(p) for p in d["points"]]; cands = [tuple(c) for c in d["cands"]]; open_pts = [tuple(o) for o in d["open"]]
    elapsed_prev = d["elapsed"]
    if d.get("done"): print("Checkpoint sagt: bereits fertig."); sys.exit(0)
    print(f"RESUME ab Index {i0} von {len(kernels)} ({elapsed_prev:.0f}s bisher)", flush=True)
t0 = time.time(); last_ckpt = t0

# `state` = snapshot of all counters and lists for the checkpoint and result file; `save` = write the checkpoint
# (temporary file, then rename)
def state(i, done=False):
    return dict(TAG=TAG, MLO=MLO, MHI=MHI, H=H, PSMOOTH=PSMOOTH, TMIN=TMIN, LMAX=LMAX, n_kernels=len(kernels), i=i,
                fam=fam, unit_big=unit_big, cand=cand, pts=pts, dead_par=dead_par, dead_wit=dead_wit,
                points=points, cands=cands, open=open_pts, elapsed=elapsed_prev + time.time() - t0, done=done,
                saved=time.strftime("%Y-%m-%d %H:%M:%S"))
def save(i, done=False):
    tmp = CKPT + ".tmp"
    with open(tmp, "w") as f: json.dump(state(i, done), f)
    os.replace(tmp, CKPT)

# ---------- run ----------
for i in range(i0, len(kernels)):
    if time.time() - last_ckpt >= CKPT_SEC:
        save(i); last_ckpt = time.time()
        print(f"  [ckpt] {i}/{len(kernels)}  Familien {fam}  Einheit>Box {unit_big}  Kandidaten {cand}  Punkte {pts}  offen {len(open_pts)}  {elapsed_prev + time.time()-t0:.0f}s", flush=True)
    m = kernels[i]; fam += 1
    r = fund_fast_capped(m, LMAX)
    if r is None: unit_big += 1; continue
    tm, um, l10, per = r
    # `mp` = m' (product of the primes p | m with p not dividing U1); `kmax` = largest index k that can stay below 10^H
    mp = 1
    for p in SMALL:
        if m % p == 0 and um % p: mp *= p
    kmax = int((H + 0.31) / l10) + 1
    if mp > kmax: continue
    cand += 1
    T1, U1 = fund(m)
    assert (T1 % m, U1 % m) == (tm, um) and abs(math.log10(T1) - l10) < 1e-6, ("Cross-Check", m)
    cands.append((m, mp, len(str(T1)), round(l10, 6), kmax, per, 'T1odd' if T1 % 2 else 'T1even'))
    k = mp
    while k <= kmax:
        if k % 2 == 1:
            T = T_exact(T1, U1, m, k); dg = len(str(T))
            if dg <= H:
                pts += 1; res = test_point(T)
                if res == 'parity': dead_par += 1; points.append((m, mp, k, dg, 'parity', 0))
                elif res == 'open': open_pts.append((m, k, dg)); points.append((m, mp, k, dg, 'open', 0))
                else: dead_wit += 1; points.append((m, mp, k, dg, 'witness', res[1]))
        k += mp
    if fam % 20000 == 0:
        print(f"  ... {i}/{len(kernels)}  Familien {fam}  Einheit>Box {unit_big}  Kandidaten {cand}  Punkte {pts}  offen {len(open_pts)}  {elapsed_prev + time.time()-t0:.0f}s", flush=True)

save(len(kernels), done=True)
print(f"W10-Box: Kerne {len(kernels)} in ({MLO}, {MHI}], Primfaktoren <= {PSMOOTH}, >= {TMIN} Faktoren, n < 10^{H}")
print(f"Familien {fam}   Einheit zu gross (Periode > {LMAX}, rigoros ohne Punkt): {unit_big}   Kandidaten (exakt): {cand}")
print(f"Gitterpunkte: {pts}   tot(Paritaet): {dead_par}   tot(Zeuge): {dead_wit}   OFFEN: {len(open_pts)}   Zeit: {elapsed_prev + time.time()-t0:.0f}s")
print("Kandidaten (m, m', Stellen T1, log10 T1, kmax, Periode, Paritaet T1):")
for c in cands: print("  C", c)
print("Punkte (m, m', k, Stellen, Status, Zeuge):")
for p in points: print("  P", p)
print("OFFEN:", open_pts)
with open(RESULT, "w") as fjs: json.dump(state(len(kernels), True), fjs, indent=1)
print("Result-JSON:", RESULT)
