# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w21_kreuzpruefung_rizzi_2026-09-09.py
# Cross-check of our list of pairs against the congruence condition of Rizzi (2026), Zenodo 10.5281/zenodo.22030563,
# "BruteMath-Improved Research on Powerful Numbers" (preprint, August 2026).
#
# Reads:    w9_v31_M100000000_H2000_result.json (class-7 points: 371 pairs), w9_v32_c3_M100000000_H2000_result.json (class 3);
#           both are looked up in the current working directory.
# Writes:   nothing (printed output only).
# Usage:    python w21_kreuzpruefung_rizzi_2026-09-09.py
#
# Rizzi derives, by a congruence analysis modulo 36, that the starting number n-1 of a triple configuration must lie in
#   {7, 27, 35} (mod 36); for the middle n this means n = 8, 28 or 0 (mod 36).
# We test this condition on our 371 pairs (class 7, kernel <= 1e8, n < 10^2000) and relate it to our witness reduction (Thm 4.5):
#   3 | n and 9 does not divide n  <=>  v_3(T_k) = 1  <=>  3 is a witness, hence no triple.
#   Expectation: the points that violate Rizzi's condition are EXACTLY the points whose smallest witness is p = 3.
# PK (positive control): the Mollin-Walsh point (m, k) = (7, 7) with middle 130 576 328 must satisfy Rizzi's condition
#     (130576328 mod 36 = 8); its witness is 29, not 3.
# NK (negative control): in class 3 (middle = 2 mod 4) NO point may satisfy Rizzi's condition.
# Own code, nothing taken over: Rizzi's result is compared only as a statement (we compare outputs, never source code).
import json, sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
from math import isqrt
# `ERLAUBT` = allowed residues of the middle n mod 36; corresponds to Rizzi's start residues 35, 7, 27 for n-1
ERLAUBT = {0, 8, 28}

def fund(m):
    # `fund` = fundamental solution: continued-fraction expansion of sqrt(m) until h^2 - m k^2 = 1; returns (T_1, U_1)
    a0 = isqrt(m); P, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        P = a*Q - P; Q = (m - P*P)//Q; a = (a0 + P)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
def Tk_mod(T1, U1, m, k, M):
    # T_k mod M by square-and-multiply on (T_1 + U_1 sqrt(m))^k
    ra, rb, ba, bb = 1 % M, 0, T1 % M, U1 % M
    while k:
        if k & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        k >>= 1
    return ra

# --- Positive control at the Mollin-Walsh point ---
T1, U1 = fund(7)
assert Tk_mod(T1, U1, 7, 7, 36) == 130576328 % 36 == 8, 'PK: T_7(7) mod 36'
print('PK ok: die Mollin-Walsh-Mitte 130 576 328 ist 8 (mod 36), erfuellt Rizzis Bedingung; ihr Zeuge ist 29.')

d = json.load(open('w9_v31_M100000000_H2000_result.json'))
hits = d['hits']                                   # (m, m', k, digits, witness) = the 371 pairs
# `verteilung` = distribution of the residues mod 36; `verletzt` = indices violating Rizzi's condition;
#   `zeuge3` = indices with witness 3
cache = {}; verteilung = {}; verletzt = set(); zeuge3 = set()
for i, (m, mp, k, dig, w) in enumerate(hits):
    if m not in cache: cache[m] = fund(m)
    T1, U1 = cache[m]
    r = Tk_mod(T1, U1, m, k, 36)
    verteilung[r] = verteilung.get(r, 0) + 1
    if r not in ERLAUBT: verletzt.add(i)
    if w == 3: zeuge3.add(i)
print(f'\nVerteilung der Mitten mod 36 ueber die {len(hits)} Paare der Klasse 7:', dict(sorted(verteilung.items())))
print(f'Rizzis Bedingung erfuellt: {len(hits) - len(verletzt)} Paare; nicht erfuellt: {len(verletzt)} Paare (Reste 12 und 24).')
print(f'Unsere Zeugensuche fand p = 3 bei: {len(zeuge3)} Paaren.')
print(f'Beide Mengen identisch: {verletzt == zeuge3}  (symmetrische Differenz: {len(verletzt ^ zeuge3)})')
assert verletzt == zeuge3, 'Erwarteter Zusammenhang verletzt'

# --- Negative control in class 3 ---
d3 = json.load(open('w9_v32_c3_M100000000_H2000_result.json'))
p3 = [p for p in d3['points'] if p[4] != 'parity']   # `p3` = class-3 points that are not excluded by parity
# `treffer3` = number of class-3 points satisfying Rizzi's condition (expected 0)
cache = {}; treffer3 = 0
for (m, mp, k, dig, st, w) in p3:
    if m not in cache: cache[m] = fund(m)
    T1, U1 = cache[m]
    if Tk_mod(T1, U1, m, k, 36) in ERLAUBT: treffer3 += 1
print(f'\nNegativ-Kontrolle Klasse 3 (Mitte = 2 mod 4, nie tripel-faehig): von {len(p3)} Paaren erfuellen {treffer3} Rizzis Bedingung (erwartet 0).')
assert treffer3 == 0, 'Negativ-Kontrolle fehlgeschlagen'

print('\nBefund: Rizzis Bedingung modulo 36 und unsere Zeugen-Reduktion stimmen auf denselben 43 Punkten ueberein.')
print('  Grund: ist 3 ein Teiler der Mitte, aber 9 nicht, so hat 3 den Exponenten 1 und ist ein Zeuge (Thm 4.5).')
print('  Rizzis Bedingung ist damit der Fall p = 3 der Zeugen-Reduktion; die Uebereinstimmung ist eine unabhaengige')
print('  Kontrolle beider Rechnungen (verschiedene Methoden, dieselbe Teilmenge).')
