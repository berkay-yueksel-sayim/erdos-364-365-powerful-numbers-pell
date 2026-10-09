# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w11_a135735_tuerme_2026-09-07.py
# Purpose: towers (`tuerme`) for the four A135735 kernels = 7 mod 8 (Reinhart 2024: condition (C) is open for d = 4099215 and
#   d = 39028039587479; not even known for d = 7). Purely modular: from the continued fraction of sqrt(m) we obtain
#   T_1, U_1 mod p^2 for all primes 3 <= p <= P (vectorized), plus log10 T_1 (scaled float). With these:
#   stage s: k_s = k_{s-1} * prod{ p : v_p(T_{k_{s-1}}) = 1 }   (ladder; M-W 1986 p. 111 / LTE for V-sequences, p odd).
#   Reason: for p | T_{k0} and odd j, v_p(T_{k0 j}) = v_p(T_{k0}) + v_p(j); T_k powerful forces p | j.
#   Lower bound: log10 n = log10 T_k >= k * log10 T_1 - 0.31.
# Only odd witnesses p are used (the 2-adic LTE is different); prime divisors of m never divide T_k.
# Reads: nothing. Writes: nothing (printed output only).
# Usage: python w11_a135735_tuerme_2026-09-07.py [P=3000] [STUFEN=3] [SEL]
#   `STUFEN` = number of stages; SEL = comma-separated list of kernels (default: all four)
# Controls: positive control with m = 7 (Mollin-Walsh tower), see below.

import numpy as np, time, math, sys
from math import isqrt
sys.set_int_max_str_digits(2000000)
P = int(sys.argv[1]) if len(sys.argv) > 1 else 3000; STUFEN = int(sys.argv[2]) if len(sys.argv) > 2 else 3
SEL = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else None
def sieve(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2, isqrt(n)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
# `PR` = odd primes up to P; `P2` = array of their squares
PR = [p for p in sieve(P) if p > 2]
P2 = np.array([p*p for p in PR], dtype=np.int64)

def cf_mod_p2(m, lmax=10**7):
    """T_1, U_1 mod p^2 (all p in PR) and log10 T_1; None if the period exceeds lmax."""
    a0 = isqrt(m); Pq, Q, a = 0, 1, a0
    h1 = np.ones(len(PR), dtype=np.int64); h0 = np.full(len(PR), a0, dtype=np.int64) % P2
    k1 = np.zeros(len(PR), dtype=np.int64); k0 = np.ones(len(PR), dtype=np.int64)
    hf1, hf0, s = 1.0, float(a0), 0; i = 0
    while True:
        Pq = a*Q - Pq; Q = (m - Pq*Pq)//Q; a = (a0 + Pq)//Q; i += 1
        if Q == 1 and i % 2 == 0: return h0, k0, math.log10(hf0) + s, i - 1
        if i > lmax: return None
        h1, h0 = h0, np.mod(a*h0 + h1, P2); k1, k0 = k0, np.mod(a*k0 + k1, P2)
        hf1, hf0 = hf0, a*hf0 + hf1
        if hf0 > 1e100: hf0 /= 1e100; hf1 /= 1e100; s += 100

# `Tk_mod` = T_k mod M by binary powering of T1 + U1*sqrt(m) (T1m, U1m already reduced)
def Tk_mod(T1m, U1m, m, k, M):
    ra, rb = 1 % M, 0; ba, bb = T1m % M, U1m % M
    while k:
        if k & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        k >>= 1
    return ra

def witnesses(h, kk, m, k):
    """Primes p in PR (p not dividing m) with v_p(T_k) = 1, from T_1, U_1 mod p^2."""
    out = []
    for j, p in enumerate(PR):
        if m % p == 0: continue
        t = Tk_mod(int(h[j]), int(kk[j]), m, k, p*p)
        if t % p == 0 and t != 0: out.append(p)
    return out

# `tower` = ladder of stages for the kernel m: `stufen` = list of (k, witnesses, number of witnesses) per stage;
# `lb_digits` = lower bound for the number of digits of the middle n
def tower(m, label):
    r = cf_mod_p2(m)
    if r is None: print(f"{label}: Periode > 1e7"); return
    h, kk, l10, idx = r
    T1_even = None  # parity handled separately (from earlier runs): 209991 odd, the others even
    k = 1; stufen = []
    for s in range(1, STUFEN + 1):
        w = witnesses(h, kk, m, k)
        if not w: stufen.append((k, [])); break
        prod = 1
        for p in w: prod *= p
        stufen.append((k, w[:10] if len(w) > 10 else w, len(w)))
        k *= prod
    lb_digits = k * l10 - 0.31
    print(f"{label}: log10 T_1 = {l10:.3f} ({int(l10)+1} Stellen), Konvergenten-Index {idx}")
    for s, st in enumerate(stufen, 1):
        if len(st) == 3: print(f"  Stufe {s}: k = {st[0]}  -> Zeugen v_p = 1 (p <= {P}): {st[1]}{' ...' if st[2] > 10 else ''} ({st[2]} Stueck)")
        else: print(f"  Stufe {s}: k = {st[0]}  -> kein Zeuge <= {P} (Kette endet hier)")
    print(f"  ==> erster moeglicher Index k >= {k}  ({len(str(k))} Stellen);  Mitte n hat >= {lb_digits:.3g} Stellen  (= 10^{math.log10(max(lb_digits,1)):.2f} Stellen)")
    return k, lb_digits

# PK: m = 7 -> stage 1: k = 1, T_1 = 8 = 2^3 has no odd prime divisor, hence no odd witnesses.
# So the PK is done at k = 7 (M-W): the witnesses of T_7 = 2^3*29*197*2857 must be {29, 197, 2857}.
r7 = cf_mod_p2(7); assert r7 is not None
w7 = witnesses(r7[0], r7[1], 7, 7); assert w7 == [29, 197, 2857], w7
assert witnesses(r7[0], r7[1], 7, 3) == [11, 23], "T_3(7) = 2^3*11*23"
print("PK ok: Zeugen von T_7(7) = {29, 197, 2857}, von T_3(7) = {11, 23} (M-W 1986 S.111)", flush=True)
t0 = time.time()
for m, lab in ((4099215, "m = 4099215 (T_1 gerade)"), (117477414815, "m = 117477414815 (T_1 gerade)"),
               (39028039587479, "m = 39028039587479 (T_1 gerade)"), (209991, "m = 209991 (T_1 UNGERADE: Familie paritaetstot, Turm nur zur Info)")):
    if SEL and m not in SEL: continue
    tower(m, lab); print(f"  [{time.time()-t0:.0f}s]", flush=True)
print("Hinweis: Stufe 1 fuer m = 7 (Vergleich): k >= 7 * 29 * 197 * 2857 = 114254287 (M-W 1986), Satz C' treibt weiter auf > 10^17.6 Stellen.")
