# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w134_stufe2_rangsieb_79kerne_2026-09-29.py
# Purpose: stage 2a, rank sieve (`rangsieb`) for the 47 new pairs (m, d) of the 79 kernels without a witness <= 3e6 (list from
#   w133): for each pair, the first witness p = +-1 (mod 4d) with p | T_d, p of rank d and v_p(T_d) = 1 is searched in (3e6, 1e8],
#   then in (1e8, 1e9]. No ECM here (w135).
# Reads: ergebnisse/w133_probelauf_79_kerne_result.json (key `neu_ohne`); resumes from the checkpoint file if present. Writes
#   (ergebnisse/): w134_stufe2_rangsieb_checkpoint.json (after every pair), w134_stufe2_rangsieb_output.txt and
#   w134_stufe2_rangsieb_result.json.
# Usage: python w134_stufe2_rangsieb_79kerne_2026-09-29.py [upper limit, default 1e9] [1 = run only the positive controls]
# Controls: positive (before the run): (7, 47) must give 23 971 691 as the first witness in (3e6, 1e8] (as in w71), and (7, 51)
#   must give 582 497 521 in (1e8, 1e9] (as in w73); on failure the script stops without running the pairs.

# Background: the pairs come from extending Prop. 6.10 to all 79 kernels. Same routes as for the earlier 811 pairs, without t40:
#   sieve up to 1e9 over the primes p = +-1 (mod 4d); ECM up to t35 follows in a separate run (w135). This script is the sieve
#   part (the core of w71 and w73 in one script).
# Method: for each pair (m, d), the candidates p = +-1 (mod 4d) in (3e6, 1e8], then (1e8, 1e9]. p | T_d is tested by
#   exponentiating eps mod p; then a primality test (Miller-Rabin with 12 fixed bases, deterministic below 3.3e24), the order of
#   eps must be exactly 4d (so p has rank d), and v_p(T_d) is read off eps^d mod p^3. A pair is SETTLED at the first witness with
#   v_p = 1. Primes of rank d with v_p >= 2 (Wieferich events) are recorded and the search continues.
import json, sys, time, pathlib
from math import isqrt
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'   # script folder; results folder
OBER = int(float(sys.argv[1])) if len(sys.argv) > 1 else 1_000_000_000   # `OBER` = upper limit of the sieve
PK_NUR = len(sys.argv) > 2 and sys.argv[2] == '1'                        # `PK_NUR` = run only the positive controls
# `FENSTER` = search windows (lower, upper]
FENSTER = [(3_000_000, 100_000_000), (100_000_000, OBER)] if OBER > 100_000_000 else [(3_000_000, OBER)]
CKPT = ERG/'w134_stufe2_rangsieb_checkpoint.json'   # `CKPT` = checkpoint file

# fundamental solution (T1, U1) of x^2 - m*y^2 = 1 from the continued fraction of sqrt(m)
def fund(m):
    a0 = isqrt(m); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (m - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0

# (T1 + U1*sqrt(m))^e modulo M by binary exponentiation; returns (rational part, sqrt(m) part)
def eps_pow(T1, U1, m, e, M):
    ra, rb, ba, bb = 1 % M, 0, T1 % M, U1 % M
    while e:
        if e & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        e >>= 1
    return ra, rb

# `ist_prim` = primality test: trial division by the first 12 primes, then Miller-Rabin with the same 12 bases
def ist_prim(n):
    if n < 2: return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0: return n == p
    d, r = n - 1, 0
    while d % 2 == 0: d //= 2; r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = x*x % n
            if x == n - 1: break
        else: return False
    return True

# `primfaktoren` = set of the prime factors of n by trial division
def primfaktoren(n):
    f, x, p = set(), n, 2
    while p*p <= x:
        while x % p == 0: f.add(p); x //= p
        p += 1 if p == 2 else 2
    if x > 1: f.add(x)
    return f

# v_p(T_d), capped at 3, from eps^d modulo p^3
def v_p(T1, U1, m, d, p):
    Td = eps_pow(T1, U1, m, d, p**3)[0]
    return 0 if Td % p else 1 if Td % (p*p) else 2 if Td % p**3 else 3

# `sieb` = sieve: search the primes p in (`unten`, `oben`] (lower, upper bound) with p = +-1 (mod 4d) for the first witness of (m,
#   d); returns (witness or None, `wief` = list of (p, v_p) with v_p >= 2, `kand` = number of candidates tested)
def sieb(T1, U1, m, d, unten, oben):
    # candidates in ascending order (both classes +1 and -1 merged), so that the "first witness" really is the smallest
    M = 4*d; qs = primfaktoren(M); wief, kand = [], 0   # `qs` = prime factors of 4d
    s1 = ((unten - 1)//M + 1)*M + 1; s2 = ((unten + 1)//M + 1)*M - 1   # first numbers = +1 and = -1 (mod M) above `unten`
    a, b = s1, s2
    while min(a, b) <= oben:
        p = min(a, b)
        if p == a: a += M
        else: b += M
        if p <= unten or p % 2 == 0: continue
        kand += 1
        if eps_pow(T1, U1, m, d, p)[0] % p: continue
        if not ist_prim(p): continue
        if eps_pow(T1, U1, m, M, p) != (1 % p, 0) or any(eps_pow(T1, U1, m, M//q, p) == (1 % p, 0) for q in qs): continue
        v = v_p(T1, U1, m, d, p)
        if v == 1: return p, wief, kand
        wief.append((p, v))
    return None, wief, kand

t0 = time.time(); aus = []   # `aus` = output lines, written to the _output.txt file
def sag(s=''): print(s, flush=True); aus.append(s)   # `sag` = say: print a line and record it in `aus`
sag('=' * 100); sag(f'w134 — STUFE 2a: RANG-SIEB FUER DIE NEUEN PAARE (bis {OBER:.0e})   ' + time.strftime('%Y-%m-%d %H:%M')); sag('=' * 100)
# positive controls: `soll` = expected first witness
pk_ok = True
for (m, d, unten, oben, soll) in ((7, 47, 3_000_000, 100_000_000, 23971691), (7, 51, 100_000_000, 1_000_000_000, 582497521)):
    T1, U1 = fund(m); p, w, k = sieb(T1, U1, m, d, unten, oben)
    ok = p == soll; pk_ok &= ok
    sag(f'PK ({m}, {d}) in ({unten:.0e}, {oben:.0e}]: erster Zeuge {p} · Soll {soll} (w71/w73) · {k} Kandidaten · {"✅" if ok else "❌"}')
if not pk_ok: sag('❌ Positivkontrolle gescheitert — KEIN Lauf.'); sys.exit(1)
if PK_NUR: sys.exit(0)

w133 = json.loads((ERG/'w133_probelauf_79_kerne_result.json').read_text(encoding='utf-8'))
paare = [tuple(x) for x in w133['neu_ohne']]   # `paare` = the new pairs (m, d) without a witness
# `ck` = checkpoint: "m_d" -> entry of a finished pair
ck = json.loads(CKPT.read_text(encoding='utf-8')) if CKPT.exists() else {}
sag(f'Paare aus w133 (neu, ohne Zeugen ≤ 3·10⁶): {len(paare)} · schon im Checkpoint: {len(ck)}')
fund_cache = {}   # kernel m -> (T1, U1)
for m, d in paare:
    schl = f'{m}_{d}'   # `schl` = checkpoint key
    if schl in ck: continue
    if m not in fund_cache: fund_cache[m] = fund(m)
    # `eintrag` = entry: `zeuge` = witness, `route` = window in which it was found (`sieb8` up to 1e8, `sieb9` up to `OBER`),
    #   `wieferich` = Wieferich events, `kandidaten` = candidates tested, `sek` = seconds
    T1, U1 = fund_cache[m]; eintrag = dict(m=m, d=d, zeuge=None, route=None, wieferich=[], kandidaten=0)
    t1 = time.time()
    for (u, o), name in zip(FENSTER, ('sieb8', 'sieb9')):
        p, w, k = sieb(T1, U1, m, d, u, o); eintrag['kandidaten'] += k; eintrag['wieferich'] += w
        if p: eintrag.update(zeuge=p, route=name); break
    eintrag['sek'] = round(time.time() - t1, 1); ck[schl] = eintrag
    CKPT.write_text(json.dumps(ck, ensure_ascii=False, indent=1), encoding='utf-8')
    sag(f'  ({m}, {d}): {"Zeuge " + str(eintrag["zeuge"]) + " (" + eintrag["route"] + ")" if eintrag["zeuge"] else "kein Zeuge bis " + f"{OBER:.0e}"}'
        f'{" · Wieferich " + str(eintrag["wieferich"]) if eintrag["wieferich"] else ""} · {eintrag["sek"]} s')
# `erl` = settled (witness found); `offen` = open
erl = [e for e in ck.values() if e['zeuge']]; offen = [e for e in ck.values() if not e['zeuge']]
sag(f'\nERGEBNIS: {len(ck)} Paare · Zeuge gefunden {len(erl)} (sieb8 {sum(e["route"] == "sieb8" for e in erl)}, sieb9 {sum(e["route"] == "sieb9" for e in erl)}) · '
    f'ohne Zeuge bis {OBER:.0e}: {len(offen)} · Wieferich-Ereignisse {sum(len(e["wieferich"]) for e in ck.values())} · {time.time() - t0:.0f} s')
sag('⛔ Nicht hier: ECM auf den Resten (w135). Die Schwelle fuer neue offene Bloecke gilt NACH ECM, nicht nach dem Sieb.')
# result file: `skript` = script name, `datum` = date, `ober` = upper limit, `pk` = controls passed, `paare` = number of pairs,
#   `erledigt` = settled pairs [m, d, witness, route], `offen` = open pairs, `wieferich` = Wieferich events [m, d, p, v_p]
(ERG/'w134_stufe2_rangsieb_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
(ERG/'w134_stufe2_rangsieb_result.json').write_text(json.dumps(dict(
    skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), ober=OBER, pk=pk_ok, paare=len(ck),
    erledigt=[[e['m'], e['d'], e['zeuge'], e['route']] for e in erl], offen=[[e['m'], e['d']] for e in offen],
    wieferich=[[e['m'], e['d'], p, v] for e in ck.values() for p, v in e['wieferich']]), ensure_ascii=False, indent=1), encoding='utf-8')
