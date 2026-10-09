# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w152_paritaet_t1_reinhart_kerne_2026-09-30.py
# Purpose: the parity of T_1 for the kernels of the Reinhart list (OEIS A135735, classes 3 and 7), COMPUTED rather than assumed.
# Reads: ergebnisse/w16_a135735_klasse3_kern_2026-09-08_output.txt (second line: the terms by residue class),
#   ergebnisse/w11_a135735_tuerme_2026-09-07_output.txt (parity labels). Requires SymPy.
# Writes: ergebnisse/w152_paritaet_t1_reinhart_kerne_result.json and ergebnisse/w152_paritaet_t1_reinhart_kerne_output.txt.
# Usage: python w152_paritaet_t1_reinhart_kerne_2026-09-30.py (no arguments; the results folder is found relative to the script).
# Controls: PK+: for all non-squares 2 <= m < 400, T mod 2^16 and the parity from route A agree with SymPy (library function
#   `diop_DN`, used as a black box). PK (error detector): the comparison with the w11 labels must report EXACTLY one deviation,
#   m = 117 477 414 815. NK: m = 7 (T_1 = 8) and m = 209 991 (T_1 odd) give their known values.

# Background: run w11 does not compute the parity. In `tower()` it sets `T1_even = None` with a note that the parity was taken
#   from earlier runs (209991 odd, the others even), and its output labels m = 117 477 414 815 as having T_1 even
#   (`(T_1 gerade)`). That label is wrong: T_1(117 477 414 815) is ODD (T_1 = 9 mod 16). Run w16 computes the parity with
#   (h mod 2) and is not affected.
# What: for each kernel m of classes 3 and 7 of the list (taken from the output of w16, not typed by hand), the fundamental
#   solution T^2 - m U^2 = 1 from the continued fraction of sqrt(m), and from it T_1 mod 2^16, the parity, "m | U_1" and the
#   period. Two routes: ROUTE A (all kernels): convergents modulo M = 2^16 * m; check at the end of the period: p^2 - m q^2 =
#   (-1)^l (mod M), which a wrong recurrence would fail. ROUTE B (period <= 250 000): exact integers; check T^2 - m U^2 = 1
#   exactly; T mod 2^16 and U mod m must agree with route A. For m = 39 028 039 587 479 (period 3.65e6, T_1 with 1.9e6 digits)
#   route B is not feasible; there route A alone carries the result (with the norm check), marked "route A only" in the output.
# Expectations (stated before the run): 209 991 odd (T_1 = 7 mod 16) · 4 099 215 even (= 4) · 117 477 414 815 ODD (= 9) ·
#   39 028 039 587 479 even (= 8) · class 3: 1 752 299 odd, 20 256 129 307 923 even (w16) · m | U_1 for all of them.
import sys, re, ast, json, math, time, pathlib
sys.set_int_max_str_digits(0)
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent / 'ergebnisse'   # script folder; results folder
W16 = ERG / 'w16_a135735_klasse3_kern_2026-09-08_output.txt'   # output of w16 (provides the kernel list)
W11 = ERG / 'w11_a135735_tuerme_2026-09-07_output.txt'          # output of w11 (provides the parity labels)
EXAKT_BIS = 250_000   # `EXAKT_BIS` = largest period length for which the exact route B is run
aus = []   # `aus` = output lines, written to the _output.txt file
def sag(s=''):   # `sag` = say: print a line and record it in `aus`
    print(s, flush=True); aus.append(s)

def cf_mod(m, M, max_schritte=10 ** 8):
    """Continued fraction of sqrt(m); returns (period l, p_{l-1} mod M, q_{l-1} mod M). P, Q, a exact, convergents modulo M (`max_schritte` = step limit)."""
    a0 = math.isqrt(m); P, Q, a = 0, 1, a0
    p1, p = 1, a0 % M; q1, q = 0, 1; l = 0
    while True:
        P = a * Q - P; Q = (m - P * P) // Q; a = (a0 + P) // Q; l += 1
        p1, p = p, (a * p + p1) % M
        q1, q = q, (a * q + q1) % M
        if Q == 1: return l, p1, q1
        if l > max_schritte: raise RuntimeError('Periode zu lang')

# route A: returns `periode` (period), `T_mod` (T_1 mod 2^16) and `U_mod_m` (U_1 mod m), computed from the convergents modulo
#   M = 2^16 * m (an odd period l gives a solution of -1, which is squared)
def route_a(m):
    M = (1 << 16) * m
    l, p, q = cf_mod(m, M)
    norm = (p * p - m * q * q) % M
    assert norm == (1 if l % 2 == 0 else M - 1), f'Norm-Kontrolle faellt durch: m = {m}'
    if l % 2 == 0: T, U = p, q
    else: T, U = (p * p + m * q * q) % M, (2 * p * q) % M
    return dict(periode=l, T_mod=T % 65536, U_mod_m=U % m)

# route B: exact integers (same recurrence without reduction); `T1_stellen` = number of digits of T_1 (only if the period is
#   <= 30000)
def route_b(m):
    a0 = math.isqrt(m); P, Q, a = 0, 1, a0
    p1, p = 1, a0; q1, q = 0, 1; l = 0
    while True:
        P = a * Q - P; Q = (m - P * P) // Q; a = (a0 + P) // Q; l += 1
        p1, p = p, a * p + p1; q1, q = q, a * q + q1
        if Q == 1: break
    T, U = p1, q1
    if l % 2 == 1: T, U = T * T + m * U * U, 2 * T * U
    assert T * T - m * U * U == 1, f'exakte Norm-Kontrolle faellt durch: m = {m}'
    stellen = len(str(T)) if l <= 30000 else None
    return dict(periode=l, T_mod=T % 65536, U_mod_m=U % m, T1_stellen=stellen)

sag('=' * 100); sag(f'w152 — Paritaet von T_1 fuer die Kerne der Reinhart-Liste (A135735), berechnet · {time.strftime("%Y-%m-%d %H:%M")}'); sag('=' * 100)
# the second line of the w16 output is a dict: residue class mod 8 -> odd terms of A135735 in that class
zeile2 = W16.read_text(encoding='utf-8').splitlines()[1]
klassen = ast.literal_eval(re.search(r'\{.*\}', zeile2).group(0))   # `klassen` = classes
# `kerne` = kernels as (m, class)
kerne = sorted([(m, 3) for m in klassen[3]] + [(m, 7) for m in klassen[7]], key=lambda x: x[0])
sag(f'Liste aus w16: Klasse 3 {klassen[3]} · Klasse 7 {klassen[7]}')
assert len(klassen[7]) == 4 and len(klassen[3]) == 2

# ---- PK+: against SymPy
from sympy.solvers.diophantine.diophantine import diop_DN
pk_n = 0; pk_ungerade = 0   # number of non-squares compared; how many of them have T_1 odd
for m in range(2, 400):
    if math.isqrt(m) ** 2 == m: continue
    T_sympy = diop_DN(m, 1)[0][0]
    a = route_a(m)
    assert a['T_mod'] == T_sympy % 65536, f'PK+: m = {m}: Route A {a["T_mod"]} ≠ SymPy {T_sympy % 65536}'
    pk_n += 1; pk_ungerade += (T_sympy % 2 == 1)
sag(f'PK+ ✅ {pk_n} Nichtquadrate 2 ≤ m < 400: T_1 mod 2^16 aus Route A = SymPy (davon {pk_ungerade} mit T_1 ungerade — die Kontrolle kann beide Faelle unterscheiden)')
r7 = route_a(7); assert r7['T_mod'] == 8 and r7['U_mod_m'] == 3 % 7
sag('NK ✅ m = 7: T_1 = 8')

# ---- the kernels
# `zeilen` = result rows; keys: `klasse` = class, `periode` = period, `T1_mod16` = T_1 mod 16, `T1_gerade` = T_1 is even,
#   `m_teilt_U1` = m divides U_1, `route_b` = status of route B, `T1_stellen` = digits of T_1
zeilen = []
for m, kl in kerne:
    t0 = time.time(); a = route_a(m); b = None
    if a['periode'] <= EXAKT_BIS:
        b = route_b(m)
        assert (b['periode'], b['T_mod'], b['U_mod_m']) == (a['periode'], a['T_mod'], a['U_mod_m']), f'Route A ≠ Route B bei m = {m}'
    z = dict(m=m, klasse=kl, periode=a['periode'], T1_mod16=a['T_mod'] % 16, T1_gerade=(a['T_mod'] % 2 == 0), m_teilt_U1=(a['U_mod_m'] == 0),
             route_b=('exakt, stimmt mit A' if b else 'nicht moeglich (Periode > 250000)'), T1_stellen=(b['T1_stellen'] if b else None))
    zeilen.append(z)
    sag(f'  m = {m} (Klasse {kl}): Periode {z["periode"]}, T_1 mod 16 = {z["T1_mod16"]}, T_1 {"GERADE" if z["T1_gerade"] else "UNGERADE"}, m | U_1: {z["m_teilt_U1"]}; '
        f'Route B: {z["route_b"]}' + (f', {z["T1_stellen"]} Stellen' if z['T1_stellen'] else '') + f'  ({time.time() - t0:.0f} s)')
assert all(z['m_teilt_U1'] for z in zeilen), 'ein Kern der Liste teilt sein U_1 nicht'
paritaet = {z['m']: z['T1_gerade'] for z in zeilen}   # `paritaet` = parity: m -> True if T_1 is even
assert paritaet[209991] is False and paritaet[4099215] is True
sag('NK ✅ 209 991 ungerade, 4 099 215 gerade (wie in w11/w16)')

# ---- PK (error detector): comparison with the labels of w11
# `w11` = labels printed by w11: m -> `gerade` (even) or `UNGERADE` (odd); `abweichung` = deviations: kernels where label and
#   computation differ
w11 = {int(m): p for m, p in re.findall(r'^m = (\d+) \(T_1 (gerade|UNGERADE)', W11.read_text(encoding='utf-8'), flags=re.M)}
abweichung = sorted(m for m, p in w11.items() if (p == 'gerade') != paritaet[m])
sag(f'Etiketten von w11: {w11}')
sag(f'Abweichung zwischen w11-Etikett und Rechnung: {abweichung}')
assert abweichung == [117477414815], 'PK: der Fehler-Detektor muss genau den bekannten Fall melden'
sag('PK ✅ der Vergleich meldet genau m = 117 477 414 815 (w11: „gerade", gerechnet: ungerade)')

# kernels per class with T_1 even (`gerade`) or odd (`ungerade`)
gerade7 = [z['m'] for z in zeilen if z['klasse'] == 7 and z['T1_gerade']]
ungerade7 = [z['m'] for z in zeilen if z['klasse'] == 7 and not z['T1_gerade']]
gerade3 = [z['m'] for z in zeilen if z['klasse'] == 3 and z['T1_gerade']]
ungerade3 = [z['m'] for z in zeilen if z['klasse'] == 3 and not z['T1_gerade']]
sag(f'Klasse 7: gerade {gerade7} · ungerade {ungerade7}; Klasse 3: gerade {gerade3} · ungerade {ungerade3}')
sag('GESAMT ✅ (Werte stehen im Ergebnis; kein Text wird hier geaendert)')
# result keys: `skript` = script name, `datum` = date, `kerne` = the rows `zeilen`, `gerade_klasse7`, `ungerade_klasse7`,
#   `gerade_klasse3`, `ungerade_klasse3` = kernels per class with T_1 even / odd, `w11_abweichung` = deviation from the w11
#   labels, `pk_sympy` = number of non-squares compared with SymPy
(ERG / 'w152_paritaet_t1_reinhart_kerne_result.json').write_text(json.dumps(dict(
    skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), kerne=zeilen, gerade_klasse7=gerade7, ungerade_klasse7=ungerade7,
    gerade_klasse3=gerade3, ungerade_klasse3=ungerade3, w11_abweichung=abweichung, pk_sympy=pk_n), indent=1, ensure_ascii=False), encoding='utf-8')
(ERG / 'w152_paritaet_t1_reinhart_kerne_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
