# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w155_algebraische_teile_offener_bloecke_2026-10-01.py
# Purpose: the open building blocks (`bloecke`: primitive parts M_d of T_d) on square kernels split algebraically: M_d = A * B
#   with A | s - 1 and B | s + 1.
# Reads: ergebnisse/w138_karte_961_daten.json, ergebnisse/w146_tab_offen_result.json.
# Writes: ergebnisse/w155_algebraische_teile_dezimal.txt, ..._result.json and ..._output.txt (same prefix). No ECM is run here.
# Usage: python w155_algebraische_teile_offener_bloecke_2026-10-01.py (no arguments; the results folder is found relative to the
#   script).
# Controls: positive: for m = 7, d = 5: T_5 = 514088, M = 514088/8 = 64261 = 179 * 359, and A, B are 179 and 359; for all 26
#   square kernels and d = 5, 7, 9: A * B = M_d with both parts > 1. Negative: a kernel with T_1 + 1 not a square (m = 15) is
#   recognized as such and not split.

# Background: for kernels with T_1 + 1 = x^2 one has T_d + 1 = (x*w_d)^2 for odd d (see w153), hence T_d = (s - 1)(s + 1)
#   with s = x*w_d. Two of the ten open blocks lie on such kernels: (m, d) = (7, 1449) and (23, 161). Since gcd(s - 1, s + 1) | 2
#   and M_d is odd, M_d = gcd(M_d, s - 1) * gcd(M_d, s + 1) = A * B with coprime A, B. ECM on A and B separately works with
#   smaller numbers per curve.
# What: (1) for the open blocks listed in w146 (fields m, d): test whether T_1 + 1 is a square; if so, form M_d exactly (T_d with
#   all prime factors of the product of T_e over the divisors e < d of d removed), form s = sqrt(T_d + 1) exactly, form A and B,
#   check A * B = M_d and gcd(A, B) = 1, report digit counts and write A and B in decimal; (2) rough cost estimate: the cost per
#   ECM curve grows roughly with the square of the digit count.
# Expectations (stated before the run): exactly the two blocks (7, 1449) and (23, 161) split; A * B = M_d; both parts have more
#   than one digit; the parts are smaller than M_d (sum of the squares of the digit counts below the square of that of M_d).
import sys, math, json, time, pathlib
sys.set_int_max_str_digits(0)
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent / 'ergebnisse'   # script folder; results folder
aus = []   # `aus` = output lines, written to the _output.txt file
def sag(s=''):   # `sag` = say: print a line and record it in `aus`
    print(s, flush=True); aus.append(s)

# returns T_1 of the fundamental solution (T_1, U_1) of x^2 - m*y^2 = 1 (continued fraction of sqrt(m); an odd period length gives
#   a solution of -1, which is squared)
def fund(m):
    a0 = math.isqrt(m); P, Q, a = 0, 1, a0; p1, p = 1, a0; q1, q = 0, 1; l = 0
    while True:
        P = a * Q - P; Q = (m - P * P) // Q; a = (a0 + P) // Q; l += 1
        p1, p = p, a * p + p1; q1, q = q, a * q + q1
        if Q == 1: break
    T, U = p1, q1
    if l % 2 == 1: T, U = T * T + m * U * U, 2 * T * U
    assert T * T - m * U * U == 1
    return T

# `tfolge` = T sequence: list [T_0, ..., T_dmax] with T_0 = 1 and T_j = 2 T1 T_{j-1} - T_{j-2}
def tfolge(T1, dmax):
    T = [1, T1]
    for j in range(2, dmax + 1): T.append(2 * T1 * T[-1] - T[-2])
    return T

# `baustein` = building block M_d: T_d with all common factors with `Lw` (the product of T_e over the divisors e < d of d) divided
#   out
def baustein(d, Tl):
    Lw = 1
    for e in range(1, d):
        if d % e == 0: Lw *= Tl[e]
    M = Tl[d]
    while True:
        g = math.gcd(M, Lw)
        if g == 1: break
        M //= g
    return M

# `zerlege` = decompose: M_d = A * B with A = gcd(M_d, s - 1), B = gcd(M_d, s + 1), s = sqrt(T_d + 1)
def zerlege(m, d, Tl=None):
    """None if T_1 + 1 is not a square; otherwise (M, A, B, s)."""
    T1 = fund(m)
    if math.isqrt(T1 + 1) ** 2 != T1 + 1: return None
    Tl = Tl or tfolge(T1, d)
    M = baustein(d, Tl); s = math.isqrt(Tl[d] + 1)
    assert s * s == Tl[d] + 1, 'T_d + 1 kein Quadrat'
    A = math.gcd(M, s - 1); B = math.gcd(M, s + 1)
    assert A * B == M and math.gcd(A, B) == 1, f'M_d ≠ A·B bei ({m}, {d})'
    return M, A, B, s

sag('=' * 100); sag(f'w155 — algebraische Teile der offenen Bloecke · {time.strftime("%Y-%m-%d %H:%M")}'); sag('=' * 100)
# positive controls
M, A, B, s = zerlege(7, 5); assert M == 64261 and sorted((A, B)) == [179, 359], (M, A, B)
sag('PK+ ✅ (7, 5): M = 64261 = 179 · 359, A und B = 179 und 359')
w138 = json.load(open(ERG / 'w138_karte_961_daten.json', encoding='utf-8'))
kerne = sorted({r['m'] for r in w138['raenge']})   # `kerne` = the kernels occurring among the 961 blocks
quad = [m for m in kerne if math.isqrt(fund(m) + 1) ** 2 == fund(m) + 1]   # `quad` = square kernels (T_1 + 1 a square)
assert len(quad) == 26
for m in quad:
    T1 = fund(m); Tl = tfolge(T1, 9)
    for d in (5, 7, 9):
        r = zerlege(m, d, Tl); assert r is not None and r[0] > 1 and r[1] > 1 and r[2] > 1, (m, d)
sag(f'PK+ ✅ alle {len(quad)} Quadrat-Kerne, d = 5, 7, 9: A · B = M_d, beide Teile > 1')
assert zerlege(15, 7) is None
sag('NK ✅ m = 15 (T_1 + 1 kein Quadrat) wird nicht zerlegt')
# the open blocks
w146 = json.load(open(ERG / 'w146_tab_offen_result.json', encoding='utf-8'))
offen = [(z['m'], z['d'], z['stellen']) for z in w146['zeilen']]   # `offen` = open blocks as (m, d, number of digits `stellen`)
sag(f'\nOffene Bloecke (w146): {len(offen)}: {[(m, d) for m, d, _ in offen]}')
zeilen = []; datei = []   # `zeilen` = result rows (JSON); `datei` = lines of the decimal-expansion file
for m, d, stellen in offen:
    r = zerlege(m, d)
    if r is None: sag(f'  ({m}, {d}) {stellen} Stellen: T_1 + 1 kein Quadrat — keine algebraische Zerlegung'); continue
    M, A, B, s = r
    dM, dA, dB = len(str(M)), len(str(A)), len(str(B))
    assert dM == stellen or abs(dM - stellen) <= 1, (dM, stellen)
    rel = (dA ** 2 + dB ** 2) / dM ** 2   # `rel` = relative cost per ECM curve of the split parts versus the undivided number
    sag(f'  ({m}, {d}): M_d hat {dM} Stellen = A ({dA} Stellen) · B ({dB} Stellen); Kosten je Kurve etwa {rel:.2f}× gegenueber der ungeteilten Zahl')
    # row: `stellen_M`, `stellen_A`, `stellen_B` = digit counts of M_d, A, B; `kosten_rel` = relative cost per ECM curve
    zeilen.append(dict(m=m, d=d, stellen_M=dM, stellen_A=dA, stellen_B=dB, kosten_rel=round(rel, 3)))
    datei += [f'{m}_{d}_A {A}', f'{m}_{d}_B {B}']
assert sorted((z['m'], z['d']) for z in zeilen) == [(7, 1449), (23, 161)], [(z['m'], z['d']) for z in zeilen]
sag('\nGESAMT ✅ genau die zwei erwarteten Bloecke zerfallen; ECM darauf ist nicht gestartet.')
(ERG / 'w155_algebraische_teile_dezimal.txt').write_text('\n'.join(datei) + '\n', encoding='utf-8')
# result keys: `skript` = script name, `datum` = date, `quad_kerne` = number of square kernels, `bloecke` = the rows `zeilen`
(ERG / 'w155_algebraische_teile_result.json').write_text(json.dumps(dict(
    skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), quad_kerne=len(quad), bloecke=zeilen), indent=1, ensure_ascii=False), encoding='utf-8')
(ERG / 'w155_algebraische_teile_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
