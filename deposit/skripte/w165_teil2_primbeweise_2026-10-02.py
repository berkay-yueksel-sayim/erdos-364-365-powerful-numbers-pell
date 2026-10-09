# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w165_teil2_primbeweise_2026-10-02.py
# Purpose: primality proofs (`primbeweise`) for the 19 witnesses that w144 could only call "probable prime" (Part II, Section 7,
#   distance one, stage 2).
# Reads: ergebnisse/w144_abstand1_stufe2_result.json, skripte/w144_abstand1_stufe2_2026-09-29.py (imported as a module),
#   ergebnisse/w164_zertifikate/cert_pk.txt.
# Writes (ergebnisse/): w165_zertifikate/ (w165_ecpp.gp, w165_auftrag.json, cert_<i>.txt, zeit.txt),
#   w165_teil2_primbeweise_result.json, w165_teil2_primbeweise_output.txt. Needs SymPy and Docker with the image `dkr-pari:1`.
# Usage: python w165_teil2_primbeweise_2026-10-02.py [`--nur-pruefen`] (`--nur-pruefen` = check only: no PARI run, the
#   certificates already present are verified)
# Controls: E2: the w164 certificate `cert_pk.txt` must be accepted and four falsified copies rejected; E1: every case is rebuilt
#   and compared with w144.

# Background: w144 has 19 cases in which the cofactor C of the third number N is only a probable prime (C divides N exactly
#   once, hence C is a witness). Here each of them is proved prime by an ECPP certificate.
# Method: (1) N and C are rebuilt with the functions of w144 (`dritte`, `klein_abteilen`; imported without the main run) and the
#   digit counts are compared with w144; (2) ECPP certificates from PARI/GP 2.17.2 (container image `dkr-pari:1`, as in w164),
#   run synchronously via docker; (3) EVERY certificate is checked by our own checker: the functions below are VERBATIM those of
#   w164 (w164 has no import guard and cannot be imported).
# Expectations (stated before the run):
#   E1 19 cases of class `prim_wahrscheinlich` (probable prime) in w144; per case C has the digit count given there, C | N exactly
#       once, C is a probable prime.
#   E2 PK: the w164 certificate cert_pk.txt (nextprime(10^40)) is accepted; NK: four falsified copies are rejected (as in w164).
#   E3  All 19 certificates are accepted.
import sys, json, time, ast, pathlib, subprocess, importlib.util
from math import gcd, isqrt
import sympy
sys.stdout.reconfigure(encoding='utf-8')
# `Z` = certificate folder
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'; Z = ERG/'w165_zertifikate'; Z.mkdir(exist_ok=True)
aus = []   # `aus` = output lines, written to the _output.txt file
def sag(s=''):   # `sag` = say: print a line and record it in `aus`
    print(s, flush=True); aus.append(s)
_spec = importlib.util.spec_from_file_location('w144', HIER/'w144_abstand1_stufe2_2026-09-29.py')
# loads the functions only; the main run sits behind __main__
w144 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(w144)

# ------------------------------------------------ checker, verbatim from w164 (theorem of Goldwasser-Kilian / Atkin-Morain)
# Certificate = chain of steps (N_i, t, s, a, P): the curve y^2 = x^3 + a*x + b over Z/N_i (b fixed by the point P = (x, y)) has
#   m = N_i + 1 - t points, m = s*q, and q is the next number N_{i+1}; points are in projective coordinates (X, Y, Z).
BASEN = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)   # `BASEN` = the 12 Miller-Rabin bases
# `mr12` = Miller-Rabin test with the 12 bases `BASEN` (deterministic for n < 3.3e24)
def mr12(n):
    if n < 2: return False
    for p in BASEN:
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in BASEN:
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True
O = None   # `O` = point at infinity
def ist_O(P, N): return P is None or P[2] % N == 0   # `ist_O` = is the point at infinity (modulo N)
def addiere(P, Q, a, N):   # `addiere` = add the points P and Q on the curve with coefficient a, modulo N
    if ist_O(P, N): return Q
    if ist_O(Q, N): return P
    X1, Y1, Z1 = P; X2, Y2, Z2 = Q
    u = (Y2 * Z1 - Y1 * Z2) % N; v = (X2 * Z1 - X1 * Z2) % N
    if v == 0:
        return verdopple(P, a, N) if u == 0 else O
    w = Z1 * Z2 % N; v2 = v * v % N; v3 = v2 * v % N
    A = (u * u * w - v3 - 2 * v2 * X1 * Z2) % N
    return (v * A % N, (u * (v2 * X1 * Z2 - A) - v3 * Y1 * Z2) % N, v3 * w % N)
def verdopple(P, a, N):   # `verdopple` = double the point P
    if ist_O(P, N): return O
    X, Y, Zk = P
    if Y % N == 0: return O
    w = (a * Zk * Zk + 3 * X * X) % N; s = Y * Zk % N; B = X * Y * s % N; h = (w * w - 8 * B) % N
    return (2 * h * s % N, (w * (4 * B - h) - 8 * Y * Y * s * s) % N, 8 * s * s * s % N)
def mal(k, P, a, N):   # `mal` = times: the scalar multiple k*P by double-and-add
    R = O
    for bit in bin(k)[2:]:
        R = verdopple(R, a, N)
        if bit == '1': R = addiere(R, P, a, N)
    return R
# `pruefe_schritt` = check one step (Ni, t, s, a, P): gcd(Ni, 6) = 1; t^2 < 4 Ni; s | m with m = Ni + 1 - t; q = m/s exceeds
#   (r + 1)^2 with r = floor(Ni^(1/4)) + 1; the curve is nonsingular modulo Ni; Q = s*P is not O modulo any divisor of Ni;
#   (q - 1)*Q = -Q. Returns q on success and an error text otherwise.
def pruefe_schritt(Ni, t, s, a, P):
    x, y = P
    if Ni < 2 or gcd(Ni, 6) != 1: return 'ggT(N, 6) ≠ 1'
    if t * t >= 4 * Ni: return 't² ≥ 4N'
    m = Ni + 1 - t
    if s <= 0 or m % s: return 's teilt m nicht'
    q = m // s
    r = isqrt(isqrt(Ni)) + 1
    if q <= (r + 1) ** 2: return 'q zu klein'
    b = (y * y - x * x * x - a * x) % Ni
    if gcd((4 * a ** 3 + 27 * b * b) % Ni, Ni) != 1: return 'Kurve singulaer modulo einem Teiler'
    Pp = (x % Ni, y % Ni, 1)
    Q = mal(s, Pp, a, Ni)
    if Q is None or gcd(Q[2], Ni) != 1: return 's·P ist modulo einem Teiler O'
    R = mal(q - 1, Q, a, Ni)
    if R is None or gcd(R[2], Ni) != 1: return '(q−1)·Q entartet'
    X1, Y1, Z1 = R; X2, Y2, Z2 = Q
    if (X1 * Z2 - X2 * Z1) % Ni or (Y1 * Z2 + Y2 * Z1) % Ni: return '(q−1)·Q ≠ −Q'
    return q
# `pruefe_zertifikat` = check a whole certificate: every step valid and each N_{i+1} equal to the previous q_i; the last q must be
#   at most 2^64 and pass `mr12`. Returns True or an error text.
def pruefe_zertifikat(C):
    erwartet = None   # `erwartet` = expected value of the next N (the previous q)
    for Ni, t, s, a, P in C:
        if erwartet is not None and Ni != erwartet: return 'Kette bricht (N_{i+1} ≠ q_i)'
        q = pruefe_schritt(Ni, t, s, a, P)
        if isinstance(q, str): return q
        erwartet = q
    if erwartet is None or erwartet > 2 ** 64 or not mr12(erwartet): return 'letztes q nicht klein und prim'
    return True
# `lies` = read a certificate file
def lies(p): return [tuple(z) for z in ast.literal_eval(p.read_text(encoding='utf-8').strip())]
# ------------------------------------------------ end of the verbatim part

# E2: positive and negative control on the certificate from w164
C0 = lies(ERG/'w164_zertifikate'/'cert_pk.txt')
e2a = pruefe_zertifikat(C0) is True
def kopie(): return [list(z) for z in C0]   # `kopie` = copy of the control certificate
# four falsifications: n1 changes N, n2 changes t, n3 changes the point, n4 changes s
n1 = kopie(); n1[0][0] += 2
n2 = kopie(); n2[0][1] += 1
n3 = kopie(); n3[0][4] = [n3[0][4][0], n3[0][4][1] + 1]
n4 = kopie(); n4[0][2] *= 2
e2b = all(pruefe_zertifikat([tuple(z) for z in n]) is not True for n in (n1, n2, n3, n4))
sag(f'PK/NK E2 {"✅" if e2a and e2b else "❌"}  PK cert_pk.txt angenommen ({e2a}), vier Verfaelschungen abgelehnt ({e2b})'); assert e2a and e2b

# E1: rebuild the cases
a144 = json.load(open(ERG/'w144_abstand1_stufe2_result.json', encoding='utf-8'))
faelle = [e for e in a144['ergebnisse'] if e['klasse'] == 'prim_wahrscheinlich']   # `faelle` = cases
# `zahlen` = the numbers to prove: per case D, k, `seite` (side), digit counts of N (`stellen_N`) and C (`stellen`), and the
#   cofactor C as a string (`zahl`)
zahlen = []
for e in faelle:
    D, k, seite = int(e['D']), int(e['k']), e['seite']
    N = w144.dritte(D, k, seite); C, exp1 = w144.klein_abteilen(N)
    assert len(str(N)) == int(e['stellen']) and len(str(C)) == int(e['kofaktor_stellen']), ('Stellen ≠ w144', D, k, seite)
    assert N % C == 0 and (N // C) % C != 0 and mr12(C) and sympy.isprime(C), ('C kein einfacher wahrscheinlich primer Teiler', D, k, seite)
    zahlen.append(dict(D=D, k=k, seite=seite, stellen_N=len(str(N)), stellen=len(str(C)), zahl=str(C)))
e1 = len(zahlen) == 19
sag(f'PK E1 {"✅" if e1 else "❌"}  {len(zahlen)} Faelle neu gebaut, Stellen wie w144, C teilt N genau einmal; Kofaktoren {sorted(z["stellen"] for z in zahlen)} Stellen'); assert e1

# certificates with PARI (synchronous)
K = '\\' + '\\'   # `K` = the two-character GP comment marker (two backslashes)
# `gp` = lines of the GP script: `primecert` for each cofactor; writes cert_<i>.txt and, per number, a line with index, digits and
#   time to zeit.txt
gp = [K + ' w165 (02.10.2026): ECPP-Zertifikate fuer die 19 Kofaktoren aus w144', 'default(parisizemax, 800000000);',
      'L = [' + ', '.join(z['zahl'] for z in zahlen) + '];',
      'for(i = 1, #L, my(t0 = getwalltime()); my(c = primecert(L[i])); write("/out/cert_" i ".txt", c); write("/out/zeit.txt", i, " ", #Str(L[i]), " ", getwalltime() - t0));',
      'quit']
# `--nur-pruefen` (check only): a run can be cut off by a time limit while the PARI container keeps running and writing its files;
#   then only the finished certificates are checked. The job must be the same (numbers compared with w165_auftrag.json).
if '--nur-pruefen' in sys.argv:
    alt = json.load(open(Z/'w165_auftrag.json', encoding='utf-8'))
    assert [a['zahl'] for a in alt] == [z['zahl'] for z in zahlen], 'Zahlenliste weicht ab — neu rechnen statt nur pruefen'
    fehlend = [i for i in range(1, len(zahlen) + 1) if not (Z/f'cert_{i}.txt').exists()]   # `fehlend` = missing certificates
    # missing certificates are reported as pending (`ausstehend`), never silently skipped
    sag(f'       nur Pruefung: vorhandene Zertifikate, Zahlenliste unveraendert; ausstehend: {fehlend or "keine"}')
else:
    for alt in Z.glob('cert_*.txt'): alt.unlink()
    if (Z/'zeit.txt').exists(): (Z/'zeit.txt').unlink()
    (Z/'w165_ecpp.gp').write_text('\n'.join(gp) + '\n', encoding='utf-8')
    (Z/'w165_auftrag.json').write_text(json.dumps([dict(i=i + 1, **z) for i, z in enumerate(zahlen)], indent=1), encoding='utf-8')
    t0 = time.time()
    r = subprocess.run(['docker', 'run', '--rm', '-v', f'{Z}:/out', '-v', f'{Z}:/in', 'dkr-pari:1', '/in/w165_ecpp.gp'], capture_output=True, text=True,
                       errors='replace', timeout=1800)
    sag(f'       PARI fertig in {time.time() - t0:.0f} s (Exit {r.returncode})')
    assert r.returncode == 0, r.stderr[-500:]

# E3: check every certificate with our own checker
# `ergebnisse` = results per number: `bewiesen` = proved, `schritte` = number of steps, `pruefung` = check status
ergebnisse = []
for i, z in enumerate(zahlen, 1):
    if not (Z/f'cert_{i}.txt').exists():
        ergebnisse.append(dict(z, bewiesen=False, schritte=0, pruefung='ausstehend'))
        sag(f'  (D, k) = ({z["D"]}, {z["k"]}) {z["seite"]}: Kofaktor {z["stellen"]:>3} St. — ⏳ Zertifikat ausstehend'); continue
    C = lies(Z/f'cert_{i}.txt')
    ok_zahl = C[0][0] == int(z['zahl'])   # the certificate must belong to this number
    res = pruefe_zertifikat(C) if ok_zahl else 'Zertifikat gehoert zu einer anderen Zahl'
    ergebnisse.append(dict(z, bewiesen=res is True, schritte=len(C), pruefung='ok' if res is True else res))
    sag(f'  (D, k) = ({z["D"]}, {z["k"]}) {z["seite"]}: Kofaktor {z["stellen"]:>3} St. — ' + (f'✅ bewiesen (ECPP, {len(C)} Schritte), eigener Pruefer' if res is True else f'🔴 {res}'))
n_ok = sum(1 for e in ergebnisse if e['bewiesen'])
sag(f'E3 {"✅" if n_ok == len(ergebnisse) else "❌"}  {n_ok} von {len(ergebnisse)} Zertifikaten angenommen')
# `out` = result dict: `skript` = script name, `datum` = date, `werkzeug` = tool, `n` = number of cases, `bewiesen` = number
#   proved, `max_stellen` = largest digit count, `ergebnisse` = results
out = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), werkzeug='PARI/GP 2.17.2 (dkr-pari:1), Pruefer wortgleich w164',
           n=len(ergebnisse), bewiesen=n_ok, max_stellen=max(z['stellen'] for z in zahlen), ergebnisse=ergebnisse)
(ERG/'w165_teil2_primbeweise_result.json').write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w165_teil2_primbeweise_output.txt').write_text(f'w165 · {time.strftime("%Y-%m-%d %H:%M")}\n' + '\n'.join(aus) + '\n', encoding='utf-8')
sag('Ergebnis: w165_teil2_primbeweise_result.json')
