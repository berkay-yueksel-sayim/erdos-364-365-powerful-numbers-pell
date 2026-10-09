# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w116_w12_nachbau_2026-09-26.py
# -*- coding: utf-8 -*-
# Purpose: rebuild the three `w12_*` outputs (originally written without a stored script) from our own data.
# Route: everything comes from the stored w9 results, class 7 (w9_v31, kernel <= 10^8, H = 2000) and class 3 (w9_v32_c3,
#   kernel <= 10^6, H = 2000). They list for each pair (m, m', k, digits, witness); this script computes the value itself:
#   middle n = T_k(m) = the x-part of eps_m^k, eps_m = smallest solution of x^2 - m*y^2 = 1. The pair is (n - 1, n + 1).
# OEIS A076445 (13 terms) and Alekseyev's table (33 terms) are NOT read in (OEIS is CC BY-SA; it is only cited). Earlier, w12 had
#   established that our lists contain exactly these values (both directions empty). Here this is checked by SHA-256: our own
#   values up to 29 resp. 67 digits must give the same hash as the lists compared at that time. The hash does not reveal the
#   values. Every value this script prints is therefore computed by us.
# Reads: ergebnisse/w9_v31_M100000000_H2000_result.json and ergebnisse/w9_v32_c3_M1000000_H2000_result.json; the old outputs
#   ergebnisse/w12_kreuz_a076445_klasse3+7_2026-09-07_output.txt, w12_lite_a076445_check_2026-09-06_output.txt and
#   w12_paare_aus_w9_2026-09-06_output.txt for the line-by-line comparison (if a file is missing this is reported).
# Writes: ergebnisse/w116_w12_nachbau_{kreuz,lite,paare}_output.txt and ergebnisse/w116_w12_nachbau_result.json.
# Usage: python w116_w12_nachbau_2026-09-26.py   (the option --hash only prints the two hashes and exits)
# Controls: the two hashes must match (assert at the end); the cross-check text contains a PK (T_3(3) - 1 from class 3,
#   T_7(7) - 1 from class 7); each rebuilt output is compared line by line with the old w12 output (the cross-check output
#   without its dated first line).
# Own code.
import sys, json, time, pathlib, hashlib
from math import isqrt
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'   # `ERG` = `ergebnisse` (results)
# SHA-256 of the 13 values a (smaller member of each pair), ascending, each written as a followed by a newline
SHA_A076445_13 = '5d257dead2bb0e951430a12240ca1353ae2f160d006a0195c1a0b9f3a5494f11'
SHA_ALEKSEYEV_33 = 'a4af793c863e8a4d6114c7bf446cb66ad5fda3eea03be2936cff3426af0e7a9a'          # same for the 33
# `sha(werte)` = SHA-256 of the sorted values (`werte` = values), one per line
def sha(werte): return hashlib.sha256(''.join(f'{a}\n' for a in sorted(werte)).encode()).hexdigest()

# `pell(m)` = fundamental solution (x1, y1) of x^2 - m*y^2 = 1 by continued fraction; `T(m, k)` = T_k(m) exactly
def pell(m):
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
def T(m, k):
    x1, y1 = pell(m); X, Y, bx, by = 1, 0, x1, y1
    while k:
        if k & 1: X, Y = X*bx + m*Y*by, X*by + Y*bx
        bx, by = bx*bx + m*by*by, 2*bx*by; k >>= 1
    return X
# `zeuge` = witness: smallest prime p with v_p(n) = 1 by trial division up to 2*10^5; failing that, the remaining cofactor if
#   it is
# below (2*10^5)^2 (then it is prime), else None
def zeuge(n):
    p = 2
    while p < 200_000 and p*p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0: n //= p; e += 1
            if e == 1: return p
        p += 1 if p == 2 else 2
    return n if n > 1 and n < 200_000 * 200_000 else None
# `ungerade_exp` = odd exponents: primes < 2*10^5 with odd exponent in z, and the remaining cofactor
def ungerade_exp(z):
    d, p = {}, 2
    while p < 200_000:
        if z % p == 0:
            e = 0
            while z % p == 0: z //= p; e += 1
            if e % 2: d[p] = e
        p += 1 if p == 2 else 2
    return d, z

# `w7`, `w3` = stored w9 results of class 7 and class 3; `hits` = their entries (m, m', k, digits, witness)
w7 = json.load(open(ERG/'w9_v31_M100000000_H2000_result.json', encoding='utf-8'))
w3 = json.load(open(ERG/'w9_v32_c3_M1000000_H2000_result.json', encoding='utf-8'))
paare = []                             # `paare` = pairs: (a = n - 1, class, m, m', k, digits, witness)
for cls, d in ((7, w7), (3, w3)):
    for m, mp, k, st, zg in d['hits']:
        if st > 70: continue
        n = T(m, k); assert len(str(n)) == st, (cls, m, k, st)
        paare.append((n - 1, cls, m, mp, k, st, zg))
paare.sort()
# `eigene13`, `eigene33` = our own values a for the 13 resp. 33 terms; `ok13`, `ok33` = hash matches
eigene13 = [p[0] for p in paare if p[5] <= 29]; eigene33 = [p[0] for p in paare if p[5] <= 67]
ok13 = sha(eigene13) == SHA_A076445_13; ok33 = sha(eigene33) == SHA_ALEKSEYEV_33
if '--hash' in sys.argv: print(sha(eigene13), sha(eigene33)); sys.exit()
A13, A33 = set(eigene13), set(eigene33)

# ---------------------------------------------------------------- (1) cross-check (lines `K`, output `kreuz`)
K = [f'Kreuzpruefung A076445 (Alekseyev, 33 Terme) gegen Klasse 3 (W9 v3.2, m<=1e6) + Klasse 7 (W9 v3.1, m<=1e8), H=2000 -- {time.strftime("%Y-%m-%d %H:%M")}',
     f'Alekseyev-Terme: {len(eigene33)} | kleinster {min(eigene33)} | groesster hat {len(str(max(eigene33)))} Stellen',
     f'unsere Paare (n-1) bis zur Alekseyev-Grenze: {len(eigene33)} | Klasse 3: {sum(1 for p in paare if p[5] <= 67 and p[1] == 3)} | Klasse 7: {sum(1 for p in paare if p[5] <= 67 and p[1] == 7)}',
     'Alekseyev-Terme, die wir NICHT haben: ' + ('[]' if ok33 else '⚠️ Hash weicht ab — Vergleich mit der Tabelle selbst noetig'),
     'unsere Paare, die Alekseyev NICHT hat: ' + ('[]' if ok33 else '⚠️ Hash weicht ab'),
     f'PK ok: {T(3, 3) - 1} <- Klasse 3 (m=3,k=3); {T(7, 7) - 1} <- Klasse 7 (m=7,k=7)',
     'Zuordnung (n-1 -> Klasse, m, k):'] + [f'   {a} -> ({c}, {m}, {k})' for a, c, m, mp, k, st, zg in paare if st <= 67]
# ---------------------------------------------------------------- (2) middles of A076445 (lines `L`, output `lite`)
L = ['Paar a, a+2 aus A076445 -> Mitte n, Status n, kleine Primteiler von n^2-1 mit Exponent (Kern-Info)']
for a in eigene13:
    n = a + 1; d, rest = ungerade_exp(n*n - 1)
    L.append(f'  n = {n} ({len(str(n))} Stellen): nicht powerful (Zeuge {zeuge(n)}); ungerade Exponenten in n^2-1 (p<2e5): {d}, Rest {len(str(rest))} Stellen')
L.append(f'Unter 1e12: {[a for a in eigene13 if a < 10**12]} -> {sum(1 for a in eigene13 if a < 10**12)} Paare, wie in der Signaturkarte (03.09.)')
# ---------------------------------------------------------------- (3) class-7 list against both (lines `P`, output `paare`)
k7 = [p for p in paare if p[1] == 7]
b40 = [p for p in k7 if p[0] + 1 <= 10**40]; b34 = [p for p in b40 if p[0] + 1 <= 1.85e34]; b20 = [p for p in b40 if p[0] + 1 <= 3.88e20]
P = [f'Paare im Abstand 2 mit Kern <= 1e8 und n <= 1e40: {len(b40)}   (davon n <= 1.85e34: {len(b34)}, n <= 3.88e20: {len(b20)})',
     "n (Mitte) | Stellen | Kern m | m' | k | Zeuge | in A076445 | in Alekseyev-33"]
P += [f'  {a + 1} | {st} | {m} | {mp} | {k} | {zg} | {"JA" if a in A13 else "nein"} | {"ja" if a in A33 else "nein"}' for a, c, m, mp, k, st, zg in b34]
s40 = {p[0] for p in b40}; s34 = {p[0] for p in b34}
P.append(f'A076445-Terme (13), die NICHT in unserer Liste sind: {[a for a in eigene13 if a not in s40]}')
P.append(f'Alekseyev-Terme ( 33 ), die nicht in unserer Liste bis 1.85e34 sind: {[a for a in eigene33 if a not in s34 and a + 1 <= 1.85e34]}')
P.append(f'Unsere Paare <= 1.85e34, die NICHT bei Alekseyev stehen: {[a for a in s34 if a not in A33]}')
P.append('=== Alekseyev-33: Mitte mod 4 (nur n = 0 mod 4 ist EMW-relevant) ===')
m0 = [a for a in eigene33 if (a + 1) % 4 == 0]
P.append(f'Terme mit Mitte = 0 mod 4: {len(m0)} von 33 | Mitte = 2 mod 4: {sum(1 for a in eigene33 if (a + 1) % 4 == 2)}')
P += [f'   {a + 1} ({len(str(a + 1))} Stellen)' for a in m0]
P.append(f'groesster Alekseyev-Term: {max(eigene33)} ({len(str(max(eigene33)))} Stellen)')
P.append('=== unsere Paare (Kern <= 1e8) bis 60 Stellen, alle mit n = 0 mod 4 ===')
P += [f"  n = {a + 1} ({st} Stellen) Kern {m} m'={mp} k={k} Zeuge {zg} | bei Alekseyev: {'ja' if a in A33 else 'nein'} | n mod 4 = {(a + 1) % 4}" for a, c, m, mp, k, st, zg in k7 if st <= 60]

# ---------------------------------------------------------------- comparison with the old w12 outputs
# `alt`: name -> (old file, rebuilt lines, number of leading lines to skip); `bericht` = report per output;
# loop variables: `datei` = file, `neu` = rebuilt lines, `ohne` = number of lines to skip, `gleich` = number of equal lines
alt = {'kreuz': ('w12_kreuz_a076445_klasse3+7_2026-09-07_output.txt', K, 1), 'lite': ('w12_lite_a076445_check_2026-09-06_output.txt', L, 0),
       'paare': ('w12_paare_aus_w9_2026-09-06_output.txt', P, 0)}
bericht = {}
for name, (datei, neu, ohne) in alt.items():
    f = ERG/datei
    if not f.exists(): bericht[name] = 'alte Datei fehlt'; continue
    a = f.read_text(encoding='utf-8').rstrip('\n').split('\n')
    gleich = sum(1 for x, y in zip(a[ohne:], neu[ohne:]) if x == y); n_ = max(len(a), len(neu)) - ohne
    bericht[name] = f'{gleich}/{n_} Zeilen gleich' + ('' if gleich == n_ and len(a) == len(neu) else ' ⚠️')
    (ERG/f'w116_w12_nachbau_{name}_output.txt').write_text('\n'.join(neu) + '\n', encoding='utf-8')
print(f'w116 — Hash A076445 (13) {"✅" if ok13 else "🔴"} · Hash Alekseyev (33) {"✅" if ok33 else "🔴"} · ' + ' · '.join(f'{k}: {v}' for k, v in bericht.items()))
# result keys: `hash_13_ok` / `hash_33_ok` = hash matches, `vergleich_alt` = comparison with the old outputs,
# `paare_bis_67_stellen` = number of pairs up to 67 digits, `quelle` = source
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), hash_13_ok=ok13, hash_33_ok=ok33, vergleich_alt=bericht,
           paare_bis_67_stellen=len(eigene33), quelle='eigene W9-Listen (w9_v31, w9_v32_c3); OEIS A076445 und Alekseyevs Tabelle nur per Hash verglichen')
(ERG/'w116_w12_nachbau_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
assert ok13 and ok33
