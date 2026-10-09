# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w147_offene_bloecke_dezimal_2026-09-30.py
# -*- coding: utf-8 -*-
# The 10 open blocks as decimal numbers for the deposit.
#
# Reads:    w135_stufe2_ecm_79kerne_2026-09-29.py (imported as a module for `fund` and `block`), and in ../ergebnisse:
#           w92_karte_811_daten.json, w135_stufe2_ecm_result.json, w143_nacht_ecm_checkpoint.json.
# Writes:   ../ergebnisse/w147_offene_bloecke_dezimal.txt (per block: a header line "m = .. d = .. digits = ..", then the number)
#           and ../ergebnisse/w147_offene_bloecke_dezimal_result.json (metadata and SHA-256 per block).
# Usage:    python w147_offene_bloecke_dezimal_2026-09-30.py   (no arguments)
# Controls: PK E1 to E5 (listed below); every failed control aborts the run with an assertion.
#
# Part I, Section 6 states that the 10 open blocks are in the deposit with their decimal expansions, and that a prime
#   factor of exponent one in any one of them settles that block. This script makes the statement true: for each open
#   pair (m, d) the block M_d = the part of T_d(m) supported on the primes of rank exactly d is built exactly and
#   written as a decimal number to a text file.
# CONSTRUCTION: Phi_d = prod_{e | d} T_{d/e}^{mu(e)} (`w135.block`, the same function as Phi in w89). Primes of rank < d
#   divide Phi_d only if they divide d (lifting the exponent), and no prime dividing d has rank d (Cor. I.4.8(i)) => M_d =
#   Phi_d without the prime factors of d.
# EXPECTATIONS (set before the run):
#   E1  PK: (m, d) = (7, 7): Phi_7 = T_7/T_1 = 16 322 041 = 29 · 197 · 2857 (w30b), every factor ≡ ±1 (mod 28); M_7 = Phi_7.
#   E2  Exactly 10 open pairs: 7 from w92 (route `offen`), 3 from w135 (open); no duplicates.
#   E3  Digit count of M_d equals the tabulated value: old pairs = int(lg) + 1 from w92, new pairs = `stellen` (digits) in w135;
#       for the new pairs, M_d must also equal the cofactor `rest` in w135/w143.
#   E4  No open block has a prime divisor below 10^4 (otherwise it would not be open): trial division as a plausibility check.
#   E5  Every line of the output is read back and parsed again; the number must be exactly the same (check of the write path).
import json, sys, time, pathlib, hashlib, importlib.util
from math import gcd
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'   # `HIER` = this directory, `ERG` = results directory
# load w135 as a module without letting it see our command-line arguments
_spec = importlib.util.spec_from_file_location('w135', HIER/'w135_stufe2_ecm_79kerne_2026-09-29.py')
_alt = sys.argv; sys.argv = [sys.argv[0]]; w135 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(w135); sys.argv = _alt

def primteiler(n):   # `primteiler` = list of the distinct prime divisors of n (trial division)
    ps, p = [], 2
    while p * p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0: n //= p
        p += 1
    return ps + ([n] if n > 1 else [])

def M(m, d):
    # returns (Phi_d, M_d): the primitive part `phi` and `x` = phi with all prime factors of d divided out
    T1, U1 = w135.fund(m); phi = w135.block(T1, U1, m, d); x = phi
    for q in primteiler(d):
        while x % q == 0: x //= q
    return phi, x

def main():
    phi, x = M(7, 7)
    e1 = phi == 16322041 == 29 * 197 * 2857 and x == phi and all(p % 28 in (1, 27) for p in (29, 197, 2857))
    print(f"PK E1 {'✅' if e1 else '❌'}  Phi_7(7) = {phi}, M_7 = {x}"); assert e1
    w92 = json.loads((ERG/'w92_karte_811_daten.json').read_text(encoding='utf-8'))
    w135r = json.loads((ERG/'w135_stufe2_ecm_result.json').read_text(encoding='utf-8'))
    ck143 = json.loads((ERG/'w143_nacht_ecm_checkpoint.json').read_text(encoding='utf-8'))
    # `paare` = open pairs as (m, d, origin, expected digit count); 'alt' = old (from w92), 'neu' = new (from w135)
    paare = [(r['m'], r['d'], 'alt', int(float(r['lg'])) + 1) for r in w92['raenge'] if r['route'] == 'offen']
    paare += [(o[0], o[1], 'neu', o[2]) for o in w135r['offen']]
    e2 = len(paare) == 10 == len({(m, d) for m, d, *_ in paare})
    print(f"PK E2 {'✅' if e2 else '❌'}  {len(paare)} offene Paare"); assert e2
    # `zeilen` = text entries of the output file; `meta` = metadata per block for the result JSON
    zeilen, meta = [], []
    # `herkunft` = origin of the pair, `soll` = expected digit count, `klein` = prime divisors below 10^4
    for m, d, herkunft, soll in sorted(paare, key=lambda z: z[3]):
        phi, x = M(m, d); s = str(x)
        e3 = len(s) == soll and (herkunft == 'alt' or int(ck143[f'{m}_{d}']['rest']) == x)
        klein = [p for p in range(2, 10_000) if x % p == 0]
        e4 = not klein
        print(f"  ({m}, {d}) {herkunft}: M_d mit {len(s)} Stellen (Tabelle {soll}) E3 {'✅' if e3 else '❌'} · kleiner Teiler {klein[:3] or '—'} E4 {'✅' if e4 else '❌'}"
              f" · Phi_d/M_d = {phi // x}")
        assert e3 and e4
        zeilen.append(f'm = {m}  d = {d}  digits = {len(s)}\n{s}\n')
        meta.append(dict(m=m, d=d, herkunft=herkunft, stellen=len(s), phi_durch_M=phi // x, sha256=hashlib.sha256(s.encode()).hexdigest()))
    # `kopf` = header text of the output file
    kopf = ('# The open blocks of Part I, Proposition I.6.8 and Table of open pairs\n'
            '# For each pair (m, d): the block M_d, the part of T_d(m) supported on the primes of rank exactly d, in decimal.\n'
            '# Here eps_m = T_1 + U_1 sqrt(m) is the fundamental solution of x^2 - m y^2 = 1 and eps_m^d = T_d + U_d sqrt(m).\n'
            '# A prime that divides M_d exactly once shows that M_d is not powerful; a complete factorization settles the block either way.\n'
            f'# Generated by w147_offene_bloecke_dezimal_2026-09-30.py on {time.strftime("%Y-%m-%d")}; one line with m, d, digits, then the number.\n\n')
    ziel = ERG/'w147_offene_bloecke_dezimal.txt'   # `ziel` = target file
    ziel.write_text(kopf + '\n'.join(zeilen), encoding='utf-8')
    # E5: read back
    # `gelesen` = numbers parsed back from the file, keyed by (m, d)
    t = ziel.read_text(encoding='utf-8').splitlines(); gelesen = {}
    for i, z in enumerate(t):
        if z.startswith('m = '):
            f = z.split(); gelesen[(int(f[2]), int(f[5]))] = int(t[i + 1])
    e5 = len(gelesen) == 10 and all(gelesen[(z['m'], z['d'])] == M(z['m'], z['d'])[1] for z in meta)
    print(f"PK E5 {'✅' if e5 else '❌'}  zurueckgelesen: {len(gelesen)} Zahlen, alle identisch"); assert e5
    (ERG/'w147_offene_bloecke_dezimal_result.json').write_text(json.dumps(dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'),
        datei=ziel.name, bloecke=meta), indent=1, ensure_ascii=False), encoding='utf-8')
    print(f'Ergebnis: {ziel.name} ({ziel.stat().st_size} Bytes) + w147_offene_bloecke_dezimal_result.json')

if __name__ == '__main__':
    main()
