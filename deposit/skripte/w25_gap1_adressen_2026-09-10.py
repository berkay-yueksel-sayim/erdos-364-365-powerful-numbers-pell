# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w25_gap1_adressen_2026-09-10.py
# Gap-1 route: the address system for pairs (a, a+1) of consecutive powerful numbers, cross-checked against an INDEPENDENT
# direct search.
# Structure: (a, a+1) both powerful  =>  T := 2a+1 is ODD, D := squarefree kernel of T^2-1, T^2 - D U^2 = 1, and D | U. This
#   is the same Pell condition as in our distance-2 map (Lemma 1.3'), but on the ODD branch, which the scans discard as
#   `paritaetstot` (dead by parity). a = (T-1)/2, a+1 = (T+1)/2.
#   CAUTION: D | U is NECESSARY but NOT SUFFICIENT: T = 3 gives (1, 2), and 2 is not powerful (an additional 2-adic
#   condition). The address system delivers CANDIDATES; the powerful test decides.
# What this run does:
#   Route A (independent): generate all powerful numbers <= NMAX directly (n = s^2 * b^3, b squarefree) and pick out the
#   neighboring pairs. Uses no Pell equation, no kernels, no towers.
#   Route B (address system): for squarefree D <= DMAX form the fundamental solution of T^2 - D U^2 = 1, run through the tower
#   T_k, keep T_k odd with D | U_k, and test (T_k -+ 1)/2 for powerful.
# Controls (both abort on failure):
#   PK      Route B must not produce a pair that the direct search misses. A pair found only by route A is allowed only if its
#           kernel D exceeds DMAX (a limit of reach, not an error); otherwise the run aborts. If the sets agree, this is
#           reported. This is the cross-check that our distance-2 map does not yet have for itself.
#   PK-NEG  25/27 (our distance-2 showcase pair) must not appear in any route (neither a = 25 nor a = 26).
# Usage: python w25_gap1_adressen_2026-09-10.py [NMAX] [DMAX]     (defaults 10**12 and 20000). Reads no file.
# Writes: w25_gap1_N<NMAX>_D<DMAX>_result.json in the current working directory. Own code.
import sys, math, json, time
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
from math import isqrt

NMAX = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**12
DMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 20000

# `quadratfrei_liste` = squarefree flags up to n; `ist_powerful` = powerful test; `kern` = squarefree kernel of n;
# `fund` = fundamental solution of x^2 - D y^2 = 1.
def quadratfrei_liste(n):
    """squarefree flags up to n"""
    qf = bytearray(b'\x01') * (n + 1)
    i = 2
    while i * i <= n:
        q = i * i
        for j in range(q, n + 1, q): qf[j] = 0
        i += 1
    return qf

def ist_powerful(n):
    """every prime factor occurs at least squared"""
    if n < 1: return False
    if n == 1: return True
    d = 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0: n //= d; e += 1
            if e == 1: return False
        d += 1 if d == 2 else 2
    return n == 1          # a remaining prime would have exponent 1

def kern(n):
    m = 1; d = 2
    while d * d <= n:
        e = 0
        while n % d == 0: n //= d; e += 1
        if e % 2: m *= d
        d += 1 if d == 2 else 2
    return m * (n if n > 1 else 1)

def fund(D):
    """Fundamental solution of x^2 - D y^2 = 1; None if D is a square."""
    a0 = isqrt(D)
    if a0 * a0 == D: return None
    P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0 * h0 - D * k0 * k0 != 1:
        P = a * Q - P; Q = (D - P * P) // Q; a = (a0 + P) // Q
        h1, h0 = h0, a * h0 + h1; k1, k0 = k0, a * k0 + k1
    return h0, k0

t0 = time.time()
print(f'W25 -- Abstand-1-Adressen.  NMAX = {NMAX:.0e}, DMAX = {DMAX}')

# ---------- Route A: independent direct search ----------
# `qf` = squarefree flags for b up to NMAX^(1/3); `pw` = set of all powerful numbers <= NMAX; `paare_A` = the a with a and a+1
#   powerful.
qf = quadratfrei_liste(int(round(NMAX ** (1/3))) + 2)
pw = set()
b = 1
while b * b * b <= NMAX:
    if qf[b]:
        b3 = b * b * b; s = 1
        while s * s * b3 <= NMAX:
            pw.add(s * s * b3); s += 1
    b += 1
paare_A = sorted(n for n in pw if (n + 1) in pw)
print(f'  Route A (direkt): {len(pw)} powerful numbers <= {NMAX:.0e}, davon {len(paare_A)} mit powerful Nachfolger ({time.time()-t0:.0f}s)')
print(f'    die ersten zehn a: {paare_A[:10]}')

# ---------- Route B: address system ----------
# `kandidaten` = number of candidates (T odd, D | U); `tuerme` = number of towers that contribute a pair; `genutzt` = this
#   tower
# has contributed.
paare_B = set(); kandidaten = 0; tuerme = 0
for D in range(2, DMAX + 1):
    if kern(D) != D: continue                       # only squarefree D
    f = fund(D)
    if f is None: continue
    T1, U1 = f
    T, U = T1, U1
    genutzt = False
    while (T - 1) // 2 <= NMAX:
        if T % 2 == 1 and U % D == 0:               # necessary condition
            kandidaten += 1
            a = (T - 1) // 2
            if a >= 1 and ist_powerful(a) and ist_powerful(a + 1):
                paare_B.add(a); genutzt = True
        T, U = T1 * T + D * U1 * U, T1 * U + U1 * T  # next rung of the tower
    tuerme += genutzt
paare_B = sorted(paare_B)
print(f'  Route B (Adressen): {kandidaten} Kandidaten aus D <= {DMAX}, {len(paare_B)} echte Paare, {tuerme} beitragende Tuerme ({time.time()-t0:.0f}s)')
print(f'    die ersten zehn a: {paare_B[:10]}')

# ---------- PK ----------
# `nurA` / `nurB` = pairs found only by route A / only by route B.
nurA = sorted(set(paare_A) - set(paare_B))
nurB = sorted(set(paare_B) - set(paare_A))
print(f'\nPK (Gegenrechnung): nur in A: {len(nurA)}   nur in B: {len(nurB)}')
if nurB:
    sys.exit(f'ABBRUCH PK: Route B liefert Paare, die die Direktsuche nicht findet: {nurB[:5]} -- Adress-System falsch.')
if nurA:
    print(f'  ⚠️ nur in A (Adressen mit D > {DMAX}, also KEIN Fehler, sondern Reichweite): {nurA[:8]}')
    fehlend_D = sorted({kern((2*a+1)**2 - 1) for a in nurA})
    print(f'     ihre Kerne D: {fehlend_D[:12]}  (alle > DMAX? {all(d > DMAX for d in fehlend_D)})')
    if not all(d > DMAX for d in fehlend_D):
        sys.exit('ABBRUCH PK: ein fehlendes Paar hat D <= DMAX -- Route B uebersieht etwas im abgedeckten Bereich.')
else:
    print('  ✅ Beide Routen liefern EXAKT dieselbe Menge.')

# ---------- PK-NEG ----------
if 25 in paare_A or 25 in paare_B or 26 in paare_A or 26 in paare_B:
    sys.exit('ABBRUCH PK-NEG: das Abstand-2-Paar 25/27 taucht in einer Abstand-1-Liste auf.')
print('  ✅ PK-NEG ok: 25/27 taucht in keiner Route auf (ist Abstand 2, nicht Abstand 1).')

# ---------- Addresses of the pairs found ----------
print(f'\n=== Adressen der {len(paare_A)} Paare (a, D, U, Sprosse) ===')
# Per pair: T = 2a+1, D = kernel of T^2 - 1, U = sqrt((T^2 - 1)/D); `stufe` = rung k at which T occurs in the tower of D.
for a in paare_A[:14]:
    T = 2*a + 1; D = kern(T*T - 1); U2 = (T*T - 1)//D; U = isqrt(U2)
    f = fund(D); stufe = ''
    if f:
        T1, U1 = f; tt, uu = T1, U1; k = 1
        while tt < T: tt, uu = T1*tt + D*U1*uu, T1*uu + U1*tt; k += 1
        stufe = f'k={k}' if tt == T else 'k=?'
    print(f'  a = {a:>12}  D = {D:>8}  U = {U:>14}  {stufe}')

res = {'NMAX': NMAX, 'DMAX': DMAX, 'paare_direkt': paare_A, 'paare_adressen': paare_B,
       'nur_direkt': nurA, 'kandidaten': kandidaten, 'sekunden': round(time.time()-t0, 1)}
fn = f'w25_gap1_N{NMAX}_D{DMAX}_result.json'
json.dump(res, open(fn, 'w'), indent=1)
print(f'\nErgebnis in {fn}  ({time.time()-t0:.0f}s)')
