# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w27_gap1_tiefpass_2026-09-10.py
# Deeper witness pass on the OPEN cases of w26 (gap-1 route).
# Starting point: w26 (D <= 20000, height 10^2000, witnesses < 5*10^5) left 1072 of 6063 cases open; "open" means that no
#   witness was found below the bound, NOT that a triple was found.
# Two improvements over w26:
#   (1) SEARCH ORDER. Of the 4991 witnesses found, 76 % satisfy p = +-1 (mod 8) (2670 with p = 7, 1120 with p = 1, against 714
#       + 487 for p = 3, 5 mod 8), although this class makes up only half of all primes. Reason: T_k = 3 (mod p) means eps^k = 3
#       +- 2*sqrt(2) (mod p), and 3 + 2*sqrt(2) lives in F_p exactly if 2 is a quadratic residue, i.e. p = +-1 (mod 8). Otherwise
#       it is possible only over F_p^2 and is rarer. => test this class FIRST. (No restriction: the other primes are searched
#       completely afterwards; 24 % of the witnesses lie there.)
#   (2) BOUND. P is raised from 5*10^5 to the given value (default 5*10^6).
# The script reads the open cases from the w26 result and processes ONLY those.
# Controls (both abort on failure): PK: the case (D=2, k=30) must find its known witness 424577. PK-NEG: the case (D=3, k=24)
#   has a-1 PRIME, so it must find NO witness below P (the script must not invent anything).
# Usage: python w27_gap1_tiefpass_2026-09-10.py [P] [W26-JSON]   (defaults: P = 5*10^6, w26_gap1_D20000_H2000_result.json).
# Reads: the w26 result file (`QUELLE`). Writes: w27_tiefpass_P<P>_result.json in the current working directory. Own code.
import sys, math, json, time
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
from math import isqrt

# `P2` = prime bound P; `QUELLE` = source file (result of w26).
P2  = int(float(sys.argv[1])) if len(sys.argv) > 1 else 5_000_000
QUELLE = sys.argv[2] if len(sys.argv) > 2 else 'w26_gap1_D20000_H2000_result.json'

# `sieb` = sieve: the odd primes up to n.
def sieb(n):
    s = bytearray(b'\x01') * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(3, n + 1, 2) if s[i]]

t0 = time.time()
alle = sieb(P2)
# Improvement (1): p = +-1 (mod 8) first, then the rest -- NO restriction.
# `alle` = all odd primes up to P; `bevorzugt` = preferred class (p = 1, 7 mod 8); `rest` = p = 3, 5 mod 8; `PRIMES` = search
#   order.
bevorzugt = [p for p in alle if p % 8 in (1, 7)]
rest      = [p for p in alle if p % 8 in (3, 5)]
PRIMES = bevorzugt + rest
print(f'W27 -- Tiefpass.  Primzahlen bis {P2:.0e}: {len(alle)} ({len(bevorzugt)} bevorzugt, {len(rest)} danach), Sieb {time.time()-t0:.0f}s')

# `fund_exakt` = exact fundamental solution (T1, U1) of x^2 - D y^2 = 1 by continued fractions;
# `pot_mod(T1, U1, D, k, M)` = T_k mod M (binary exponentiation of T1 + U1*sqrt(D)).
def fund_exakt(D):
    a0 = isqrt(D); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - D*k0*k0 != 1:
        P = a*Q - P; Q = (D - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0

def pot_mod(T1, U1, D, k, M):
    ra, rb = 1 % M, 0; ba, bb = T1 % M, U1 % M
    while k:
        if k & 1: ra, rb = (ra*ba + D*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + D*bb*bb) % M, (2*ba*bb) % M
        k >>= 1
    return ra

# `zeuge` = witness: the first prime p (not dividing D) in `primes` with p | T_k - ziel but p^2 not dividing it (`ziel` =
#   target
# value +-3), i.e. v_p(T_k - ziel) = 1; None if there is none.
def zeuge(T1, U1, D, k, ziel, primes):
    for p in primes:
        if D % p == 0: continue
        pp = p * p
        t = pot_mod(T1, U1, D, k, pp)
        if (t - ziel) % p == 0 and (t - ziel) % pp != 0:
            return p
    return None

# ---------- PK ----------
T1, U1 = fund_exakt(2)
z = zeuge(T1, U1, 2, 30, 3, PRIMES)
assert z == 424577, f'PK: (D=2, k=30) sollte den Zeugen 424577 liefern, lieferte {z}'
T1b, U1b = fund_exakt(3)
z2 = zeuge(T1b, U1b, 3, 24, 3, PRIMES)
assert z2 is None, f'PK-NEG: (D=3, k=24) hat a-1 PRIM (13325427460799), darf keinen Zeugen < {P2} finden, fand aber {z2}'
print(f'  PK ok: (D=2,k=30) -> Zeuge 424577 gefunden.   PK-NEG ok: (D=3,k=24) -> a-1 ist prim, kein Zeuge erfunden.')

# ---------- Deeper pass ----------
d = json.load(open(QUELLE))
offen = d['offen']
print(f'  offene Faelle aus {QUELLE}: {len(offen)}\n')
fund_cache = {}
neu = []; bleibt = []
# Each open case is (D, k, digits, side, ., .): `seite` = side ('a-1' -> target 3, otherwise -3); `neu` = newly settled,
# `bleibt` = still open.
for i, (D, k, stellen, seite, _, tot_gratis) in enumerate(offen):
    if D not in fund_cache: fund_cache[D] = fund_exakt(D)
    T1, U1 = fund_cache[D]
    ziel = 3 if seite == 'a-1' else -3
    z = zeuge(T1, U1, D, k, ziel, PRIMES)
    (neu if z is not None else bleibt).append((D, k, stellen, seite, z))
    if (i + 1) % 100 == 0:
        print(f'    … {i+1}/{len(offen)} bearbeitet, {len(neu)} neu erledigt ({time.time()-t0:.0f}s)', flush=True)

print(f'\n=== Ergebnis ===')
print(f'  neu erledigt: {len(neu)} von {len(offen)}  ({100*len(neu)/len(offen):.1f} %)')
print(f'  weiterhin offen: {len(bleibt)}')
if neu:
    zs = sorted(r[4] for r in neu)
    print(f'  neue Zeugen: kleinster {zs[0]}, Median {zs[len(zs)//2]}, groesster {zs[-1]}')
    from collections import Counter
    print(f'  davon p = +-1 (mod 8): {sum(1 for p in zs if p % 8 in (1,7))} von {len(zs)}')
print(f'  Beispiele weiterhin offen: {[(r[0], r[1], r[2]) for r in bleibt[:8]]}')

res = {'P': P2, 'quelle': QUELLE, 'neu_erledigt': neu, 'bleibt_offen': bleibt, 'sekunden': round(time.time()-t0, 1)}
fn = f'w27_tiefpass_P{P2}_result.json'
json.dump(res, open(fn, 'w'), indent=1)
print(f'\nErgebnis in {fn}  ({time.time()-t0:.0f}s)')
