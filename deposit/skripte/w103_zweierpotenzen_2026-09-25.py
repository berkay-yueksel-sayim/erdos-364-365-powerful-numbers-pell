# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w103_zweierpotenzen_2026-09-25.py
# -*- coding: utf-8 -*-
# The powers of two themselves, for the binary view (w100): five claims checked, and what lies behind them.
# The claims concern 1, 2, 4, ..., 2^63: differences = the smaller number, factor 2, last digit 2, 4, 8, 6 (cycle of length 4
#   mod 10), 1 + 2 + ... + 128 = 255 = 2^8 - 1 (IPv4 octet), 2^32 - 1 = 4 294 967 295; beyond them: check differences, digit
#   sums, last 2/3 digits and residues mod 3/7/9/11 systematically ("which patterns are new, which are forced").
# Parts:
#   (A) the five claims, each checked for all n = 0 ... 63.
#   (B) cycles that are FORCED (Euler/Fermat): last k digits (k = 1 ... 4), residues mod 3, 7, 9, 11, 13, 17; digit sum mod 9.
#   (C) leading digits of 2^n, n = 1 ... 10 000, against Benford: a THEOREM (log10 2 is irrational), not number mysticism.
#   (D) THE NEIGHBORS 2^n +- 1, n = 1 ... 63: this is where it becomes number theory: prime? (Mersenne/Fermat), powerful?
#       square factors and their ORIGIN: "LTE" (p | n, harmless and forced) or "Wieferich" (p^2 already divides 2^ord - 1).
#       2^n is powerful for n >= 2, so if 2^n - 1 and 2^n + 1 were both powerful, there would be an Erdos-364 triple.
#   (E) THE BASE-2 VERSION OF OUR PROP. 6.11: the first square factor that does NOT come from LTE is a Wieferich event:
#       1093^2 | 2^ord(1093) - 1 and 3511^2 | 2^ord(3511) - 1; primitive parts Phi_n(2) (classical cyclotomic polynomial at 2).
#   (F) The pairs from w100 by origin: which come from x^2 - 8y^2 = 1 (the Pell family), and for each pair the equation
#       b2^3*a2^2 - b1^3*a1^2 = 1 from the representation n = a1^2*b1^3, n + 1 = a2^2*b2^3.
# Reads: ergebnisse/w100_brille_binaer_result.json. Writes: ergebnisse/w103_zweierpotenzen_result.json and
#   ergebnisse/w103_zweierpotenzen_output.txt. Requires SymPy. Usage: python w103_zweierpotenzen_2026-09-25.py (no arguments).
# Expectations (each is an assert; a failed one aborts the run):
#   E1  The five claims hold (all elementary).
#   E2  Cycle lengths: last 1/2/3/4 digits 4 / 20 / 100 / 500 (= 4*5^(k-1), for n >= k); mod 3/7/9/11/13/17: 2/3/6/10/12/8
#       (= order of 2).
#   E3  Benford: leading digit 1 in about 30.1 %.
#   E4  Mersenne primes 2^n - 1 for n <= 63: n = 2, 3, 5, 7, 13, 17, 19, 31, 61; Fermat primes 2^n + 1: n = 1, 2, 4, 8, 16 [known
#       facts, measured here]. Powerful among 2^n +- 1 (n >= 2): only 2^3 + 1 = 9. Triples: none.
#   E5  All square factors of 2^n +- 1 for n <= 63 are of LTE type; Wieferich type first at 1093 (order 364) and 3511.
#   PK  2^6 - 1 = 63 = 3^2*7 (LTE: 3 | 6) · 2^21 - 1 contains 7^2 (LTE: 7 | 21) · 2^8 - 1 = 255 = 3*5*17.
import sys, json, math, pathlib
from sympy import factorint, isprime, n_order, cyclotomic_poly, symbols
sys.stdout.reconfigure(encoding='utf-8')
ERG = pathlib.Path(__file__).resolve().parent.parent/'ergebnisse'   # `ERG` = `ergebnisse` (results)
aus = []   # `aus` = output lines
# `sag` = say: print a line and record it
def sag(s=''):
    print(s, flush=True); aus.append(s)
res = dict(skript=pathlib.Path(__file__).name, datum='2026-09-25')
sag('='*100); sag('w103 — DIE ZWEIERPOTENZEN 2⁰ … 2⁶³: Aussagen, Zyklen, Benford, die Nachbarn 2ⁿ ± 1, Wieferich zur Basis 2'); sag('='*100)
Z = [2**n for n in range(64)]   # `Z` = the powers 2^0 ... 2^63

# (A)
a1 = all(Z[n+1] - Z[n] == Z[n] for n in range(63))
a2 = all(Z[n+1] == 2*Z[n] for n in range(63))
a3 = [Z[n] % 10 for n in range(1, 13)] == [2, 4, 8, 6] * 3
a4 = sum(Z[:8]) == 255 == 2**8 - 1
a5 = 2**32 - 1 == 4294967295
assert a1 and a2 and a3 and a4 and a5, 'E1 VERLETZT'
# keys of `res['A']`: `abstand` = difference, `faktor` = factor, `endziffer_zyklus` = last-digit cycle, `ipv4_255` = 255 = 2^8
# - 1, `u32_max` = 2^32 - 1, `u64_max` = 2^64 - 1
res['A'] = dict(abstand=a1, faktor=a2, endziffer_zyklus=a3, ipv4_255=a4, u32_max=a5, u64_max=str(2**64 - 1))
sag('(A) geprueft ✅: Abstand = kleinere Zahl · Faktor 2 · Endziffer 2,4,8,6 · 1+…+128 = 255 = 2⁸ − 1 · 2³² − 1 = 4 294 967 295')
sag(f'    dazu: allgemein 1 + 2 + … + 2^(n−1) = 2ⁿ − 1 (binaer n Einsen) · 2⁶⁴ − 1 = {2**64 - 1:,}')

# (B) cycles
# `periode` = period: smallest p such that the sequence `folge_fn`(n), n = start ... start + `lang` - 1, repeats with period p
def periode(folge_fn, start, lang=3000):
    xs = [folge_fn(n) for n in range(start, start + lang)]
    for p in range(1, lang // 2):
        if all(xs[i] == xs[i + p] for i in range(lang - p)): return p
    return None
zyk = {}   # `zyk` = cycle lengths found
for k in (1, 2, 3, 4):
    zyk[f'letzte_{k}'] = periode(lambda n, k=k: pow(2, n, 10**k), k)
for m in (3, 7, 9, 11, 13, 17):
    zyk[f'mod_{m}'] = periode(lambda n, m=m: pow(2, n, m), 0)
    assert zyk[f'mod_{m}'] == n_order(2, m)
zyk['quersumme_mod9'] = periode(lambda n: sum(map(int, str(2**n))) % 9, 0, 300)
assert [zyk[f'letzte_{k}'] for k in (1, 2, 3, 4)] == [4, 20, 100, 500], ('E2 VERLETZT', zyk)
assert [zyk[f'mod_{m}'] for m in (3, 7, 9, 11, 13, 17)] == [2, 3, 6, 10, 12, 8], ('E2 VERLETZT', zyk)
# keys of `zyk`: `letzte_k` = period of the last k digits, `mod_m` = period mod m, `quersumme_mod9` = period of the digit sum
# mod 9; `B_folgen` = the sequences themselves (`quersumme` = digit sum)
res['B'] = zyk
res['B_folgen'] = {'mod_3': [pow(2, n, 3) for n in range(12)], 'mod_7': [pow(2, n, 7) for n in range(12)],
                   'mod_9': [pow(2, n, 9) for n in range(12)], 'mod_11': [pow(2, n, 11) for n in range(12)],
                   'quersumme': [sum(map(int, str(2**n))) for n in range(64)], 'letzte_2': [pow(2, n, 100) for n in range(2, 24)]}
sag(f'(B) Zyklen ✅ (zwangslaeufig): letzte 1/2/3/4 Ziffern {zyk["letzte_1"]}/{zyk["letzte_2"]}/{zyk["letzte_3"]}/{zyk["letzte_4"]} · '
    f'mod 3/7/9/11/13/17: ' + '/'.join(str(zyk[f"mod_{m}"]) for m in (3, 7, 9, 11, 13, 17)) + f' · Quersumme mod 9: {zyk["quersumme_mod9"]}')

# (C) Benford (`lead` = leading digits, `benf` = per digit d: `gemessen` = measured frequency, `benford` = Benford frequency)
lead = [int(str(2**n)[0]) for n in range(1, 10001)]
benf = [dict(d=d, gemessen=round(lead.count(d) / len(lead), 4), benford=round(math.log10(1 + 1/d), 4)) for d in range(1, 10)]
assert abs(benf[0]['gemessen'] - 0.30103) < 0.005, 'E3 VERLETZT'
res['C'] = benf
sag('(C) Benford ✅: Leitziffern von 2¹ … 2¹⁰⁰⁰⁰: ' + ' · '.join(f"{b['d']}: {100*b['gemessen']:.2f} % ({100*b['benford']:.2f})" for b in benf))

# (D) neighbors
def herkunft(p, e, n, vorz):   # `herkunft` = origin
    # classifies a square factor p^e (e >= 2) of 2^n + `vorz` (`vorz` = -1 or +1): 'Wieferich' if p^2 divides 2^o - 1 (o =
    # order of 2 mod p), otherwise LTE type (p divides n); returns (type, o)
    o = n_order(2, p)
    basis = o if vorz == -1 else o // 2 if o % 2 == 0 else None   # 2^n = -1 (mod p) <=> o even and n = o/2 (mod o)
    wief = pow(2, o, p * p) == 1
    return ('Wieferich' if wief else 'LTE (p teilt n)'), o
D = []        # `D` = one entry per n: value, primality, factorization, square factors of 2^n - 1 and 2^n + 1
triple = []   # `triple` = triples 2^n - 1, 2^n, 2^n + 1 (n >= 2) with both neighbors powerful
# row `zeile` per n; keys `minus` / `plus` = 2^n - 1 / 2^n + 1, each with `wert` = value, `prim` = prime?, `powerful`,
# `zerlegung` = factorization, `quadrate` = square factors (`typ` = type)
for n in range(1, 64):
    zeile = dict(n=n)
    for name, x, vz in (('minus', 2**n - 1, -1), ('plus', 2**n + 1, +1)):
        if x == 1:
            zeile[name] = dict(wert='1', prim=False, powerful=True, zerlegung={}, quadrate=[]); continue
        f = factorint(x)
        q = [dict(p=p, e=e, typ=herkunft(p, e, n, vz)[0]) for p, e in f.items() if e >= 2]
        zeile[name] = dict(wert=str(x), prim=isprime(x), powerful=all(e >= 2 for e in f.values()),
                           zerlegung={str(p): e for p, e in f.items()}, quadrate=q)
    if n >= 2 and zeile['minus']['powerful'] and zeile['plus']['powerful']: triple.append(n)
    D.append(zeile)
res['D'] = D
# `mers` = n with 2^n - 1 prime, `ferm` = n with 2^n + 1 prime, `pow_nb` = powerful neighbors, `quadr` = square factors found
mers = [z['n'] for z in D if z['minus']['prim']]; ferm = [z['n'] for z in D if z['plus']['prim']]
pow_nb = [(z['n'], s) for z in D for s in ('minus', 'plus') if z['n'] >= 2 and z[s]['powerful']]
quadr = [(z['n'], s, q['p'], q['e'], q['typ']) for z in D for s in ('minus', 'plus') for q in z[s]['quadrate']]
assert mers == [2, 3, 5, 7, 13, 17, 19, 31, 61] and ferm == [1, 2, 4, 8, 16], ('E4 VERLETZT', mers, ferm)
assert pow_nb == [(3, 'plus')] and not triple, ('E4 VERLETZT', pow_nb, triple)
assert all(t[4].startswith('LTE') for t in quadr), ('E5 VERLETZT', [t for t in quadr if not t[4].startswith('LTE')])
assert D[5]['minus']['zerlegung'] == {'3': 2, '7': 1} and D[20]['minus']['zerlegung'].get('7') == 2 and D[7]['minus']['zerlegung'] == {'3': 1, '5': 1, '17': 1}
sag(f'(D) Nachbarn 2ⁿ ± 1, n = 1 … 63: Mersenne-prim bei n = {mers} · Fermat-prim bei n = {ferm}')
sag(f'    powerful (n ≥ 2): {pow_nb} ⇒ nur 2³ + 1 = 9 (das Paar 8 | 9) · Tripel 2ⁿ − 1, 2ⁿ, 2ⁿ + 1: {"keins" if not triple else triple}')
sag(f'    Quadratfaktoren: {len(quadr)}, alle vom Typ LTE (p teilt n) — z. B. ' + ', '.join(f'{p}^{e} in 2^{n}{"−" if s == "minus" else "+"}1' for n, s, p, e, _ in quadr[:8]))
res['D_zusammen'] = dict(mersenne=mers, fermat=ferm, powerful_nachbarn=pow_nb, tripel=triple, quadratfaktoren=len(quadr))

# (E) Wieferich primes to base 2 -- the base-2 version of our Prop. 6.11 (`E` = one entry per p)
# keys: `ordnung` = order of 2 mod p, `p_quadrat_teilt_2_hoch_ord_minus_1` = p^2 divides 2^ord - 1, `v_p_in_Phi_ord_2` =
# v_p(Phi_ord(2)), `stellen_Phi` = digits of Phi_ord(2)
X = symbols('x'); E = []
for p in (1093, 3511):
    o = n_order(2, p)
    phi_o = int(cyclotomic_poly(o, X).subs(X, 2))
    v = 0; t = phi_o
    while t % p == 0: t //= p; v += 1
    E.append(dict(p=p, ordnung=o, p_quadrat_teilt_2_hoch_ord_minus_1=pow(2, o, p*p) == 1, v_p_in_Phi_ord_2=v, stellen_Phi=len(str(phi_o))))
    sag(f'(E) p = {p}: Ordnung von 2 mod p = {o} · p² | 2^{o} − 1: {pow(2, o, p*p) == 1} · im primitiven Teil Φ_{o}(2) ({len(str(phi_o))} St.) steckt p genau {v}-mal')
assert all(e['p_quadrat_teilt_2_hoch_ord_minus_1'] and e['v_p_in_Phi_ord_2'] == 2 for e in E), 'E5 (Wieferich) VERLETZT'
res['E'] = E

# (F) the pairs from w100 by origin
# `a2b3(n)` = (a, b) with n = a^2*b^3 for a powerful n (b = product of the primes with odd exponent);
# `fam8` = the values 8*y^2 for the solutions (x, y) of x^2 - 8y^2 = 1;
# `pell8` (below) = the pair comes from this family
w100 = json.load(open(ERG/'w100_brille_binaer_result.json', encoding='utf-8'))
def a2b3(n):
    f = factorint(n); b = 1; a = 1
    for p, e in f.items():
        if e % 2: b *= p; a *= p**((e - 3)//2)
        else: a *= p**(e//2)
    assert a*a*b**3 == n
    return a, b
fam8 = set(); x, y = 3, 1
while 8*y*y < 2**47:
    fam8.add(8*y*y); x, y = 3*x + 8*y, x + 3*y
F = []
for s in w100['paare']:
    n = int(s); a1, b1 = a2b3(n); a2, b2 = a2b3(n + 1)
    assert b2**3 * a2**2 - b1**3 * a1**2 == 1
    F.append(dict(n=s, a1=a1, b1=b1, a2=a2, b2=b2, pell8=n in fam8))
res['F'] = F
sag(f'(F) Paare aus w100: {sum(f["pell8"] for f in F)} von {len(F)} aus x² − 8y² = 1 (' + ', '.join(f['n'] for f in F if f['pell8']) +
    '); jedes Paar loest b₂³·a₂² − b₁³·a₁² = 1.')
sag('E1 ✅ E2 ✅ E3 ✅ E4 ✅ E5 ✅')
# SymPy returns mpz values, hence default=int
(ERG/'w103_zweierpotenzen_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False, default=int), encoding='utf-8')
(ERG/'w103_zweierpotenzen_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
