# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w84_ist_es_dasselbe_objekt_2026-09-22.py
# Purpose: is our "Wieferich prime for eps_m" the same object as a "Lucas-Wieferich prime associated to the pair (P, Q)"?
# Background: the established name (OEIS A001220, comment by Felix Froehlich, 2016, citing McIntosh–Roettger 2007, p. 2088):
#   for a Lucas sequence U_n(P, Q) of the first kind, D = P² − 4Q, e the Legendre symbol (D/p) and a prime p not dividing 2QD,
#   the prime p is called a Lucas-Wieferich prime associated to (P, Q) if U_(p−e) ≡ 0 (mod p²).
# Why this is not obviously the same, and why it is measured rather than asserted: their condition lives on the sequence of
#   the FIRST kind U_n. Ours lives on T_k, which is the COMPANION sequence: with eps = T1 + U1*sqrt(m), of norm 1,
#       T_n = (eps^n + eps^-n)/2 = V_n/2     and     U^Lucas_n = (eps^n - eps^-n)/(2*U1*sqrt(m)).
#   p | U_n and p | V_n hit DIFFERENT index sets, so "rank" does not mean the same thing here. Same name, possibly a different
#   object: exactly the situation in which one measures flawlessly on the wrong subject.
# What is measured: for each of the 16 Wieferich events of w77 (there: p² | T_d, our condition, v_p = 2) the script checks
#   whether the Lucas condition also holds: p² | U^Lucas_(p−e), with (P, Q) = (2*T1, 1), D = P² − 4Q = 4*m*U1², and
#   e = (D/p) = (m/p) for p not dividing 2*U1.
# Reads:  ergebnisse/w77_wieferich_ereignisse_result.json, ergebnisse/w71_rangsieb_P100000000_N176_result.json and
#         ergebnisse/w20_primteil_P3000000_N25_NP5_result.json (loaded, but not used further).
# Writes: ergebnisse/w84_ist_es_dasselbe_objekt_result.json and ergebnisse/w84_ist_es_dasselbe_objekt_output.txt.
# Usage:  python w84_ist_es_dasselbe_objekt_2026-09-22.py   (no arguments)
# Controls and expectations, stated before the run:
#   E1  positive control of the translation: for small n, U^Lucas_n and V_n/2 agree with the values computed directly from eps
#       (exact integer arithmetic, no modulus).
#   E2  positive control of the Lucas condition on a KNOWN case outside our family: for the Pell sequence (P, Q) = (2, -1) the
#       primes p = 13, 31 and 1546463 must satisfy the condition (OEIS A238736, "Balancing Wieferich primes" / Pell quotient;
#       Williams 1982, p. 86). If this control fails, the script does not measure the Lucas condition and every verdict below
#       would be worthless.
#   E3  outcome OPEN. Three results are possible and all are recorded:
#         16/16 -> the conditions coincide on our family; the established name applies to us.
#         0/16  -> different objects; the name does NOT belong to our object.
#         in between -> related, not identical; the most interesting case, which needs a statement of its own in the paper.
#   E4  control in the other direction: for a sample of primes that do NOT satisfy our condition (ordinary rank primes), the
#       Lucas condition must not hold either. Otherwise it would be trivially true and the comparison meaningless.

import json, sys, time, pathlib
from math import isqrt
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent
ERG  = HIER.parent/'ergebnisse'
def J(n): return json.load(open(ERG/n, encoding='utf-8'))   # `J` = load the result file `n` from the results folder

def fund(m):   # fundamental solution (h0, k0) = (T1, U1) of x² − m·y² = 1, by the continued fraction of √m
    a0 = isqrt(m); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (m - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0

def eps_pow_exakt(T1, U1, m, n):   # `eps_pow_exakt` = eps to the power n, exact (no modulus)
    """eps^n = A + B*sqrt(m), exactly. Returns (A, B)."""
    ra, rb, ba, bb = 1, 0, T1, U1
    while n:
        if n & 1: ra, rb = ra*ba + m*rb*bb, ra*bb + rb*ba
        ba, bb = ba*ba + m*bb*bb, 2*ba*bb
        n >>= 1
    return ra, rb

def lucas_UV(P, Q, n, M):
    """U_n, V_n of the Lucas sequence for (P, Q), modulo M. Binary method, without division."""
    U, V, Qk = 0, 2, 1
    for bit in bin(n)[2:]:
        U, V, Qk = (U*V) % M, (V*V - 2*Qk) % M, (Qk*Qk) % M
        if bit == '1':
            U, V = ((P*U + V)*pow(2, -1, M)) % M if M % 2 else None, None
            raise SystemExit('unerwartet: M gerade')
    return U, V

def lucas_U(P, Q, n, M):
    """U_n mod M via the recurrence U_{k+1} = P*U_k - Q*U_{k-1}. Linear, hence slow, but it needs no
       inverse and therefore no assumption on M. Used for the control on small n."""
    a, b = 0, 1
    for _ in range(n): a, b = b, (P*b - Q*a) % M
    return a

def lucas_U_schnell(P, Q, n, M):
    """U_n mod M, binary method (doubling formulas). U_{2k} = U_k*V_k, V_{2k} = V_k^2 - 2Q^k,
       U_{2k+1} = (P*U_{2k} + V_{2k})/2, V_{2k+1} = (D*U_{2k} + P*V_{2k})/2. The division by 2
       is avoided by computing modulo 2*M; the result is reduced modulo M at the end."""
    D = P*P - 4*Q
    M2 = 2*M
    U, V, Qk = 0 % M2, 2 % M2, 1 % M2
    for bit in bin(n)[2:]:
        U, V, Qk = (U*V) % M2, (V*V - 2*Qk) % M2, (Qk*Qk) % M2
        if bit == '1':
            U, V = ((P*U + V) // 2 if (P*U + V) % 2 == 0 else (P*U + V + M2) // 2) % M2, \
                   ((D*U + P*V) // 2 if (D*U + P*V) % 2 == 0 else (D*U + P*V + M2) // 2) % M2
            Qk = (Qk*Q) % M2
    return U % M

def legendre(a, p):   # Legendre symbol (a/p) for an odd prime p
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1

t0 = time.time(); aus = []   # `aus` = output lines, written to the output file at the end
def sag(s=''):   # `sag` = "say": print a line and keep it for the output file
    print(s, flush=True); aus.append(s)

sag('='*100)
sag('w84 — IST ES DASSELBE OBJEKT? Unsere T-Bedingung gegen die Lucas-U-Bedingung   ' + time.strftime('%Y-%m-%d %H:%M'))
sag('='*100)
sag('Etablierter Begriff (OEIS A001220, Froehlich 2016, nach McIntosh–Roettger 2007 S. 2088):')
sag('  p ist "Lucas-Wieferich prime associated to the pair (P, Q)"  :⟺  U_(p−e) ≡ 0 (mod p²),  e = (D/p).')
sag('Unsere Fassung (Prop. 6.11): p² | T_{α(p)}, α = Rang in der BEGLEITfolge T = V/2.')
sag()

# ---------------------------------------------------------------- E1: check the translation
sag('--- PK E1: stimmt die Uebersetzung eps <-> Lucas(P, Q) = (2T₁, 1)?')
T1, U1 = fund(7); P, Q = 2*T1, 1
M = 10**40
ok = True
for n in (1, 2, 3, 5, 8, 13):
    A, B = eps_pow_exakt(T1, U1, 7, n)          # eps^n = A + B*sqrt(7); T_n = A, U_n = B
    u_l = lucas_U_schnell(P, Q, n, M)
    u_soll = (B // U1) % M
    assert B % U1 == 0, ('E1: U_n nicht durch U1 teilbar', n)
    if u_l != u_soll: ok = False; sag(f'    🔴 n={n}: Lucas-U {u_l} != U_n/U₁ {u_soll}')
assert ok, 'E1 VERLETZT — die Uebersetzung stimmt nicht, alles weitere waere wertlos'
sag(f'    ✅ m=7: U^Lucas_n = U_n/U₁ fuer n = 1, 2, 3, 5, 8, 13 (exakt, T₁={T1}, U₁={U1}, P={P}).')

# ---------------------------------------------------------------- E2: positive control on a KNOWN case
sag()
sag('--- PK E2: bekannte Lucas-Wieferich-Primzahlen der PELL-Folge (P, Q) = (2, −1)')
sag('    Quelle: OEIS A238736 "Balancing Wieferich primes ... Pell quotient", Werte 13, 31, 1546463')
sag('    (Kommentar dort: Williams 1982, S. 86 — die einzigen unter 10⁸; keine weiteren bis 10¹⁰).')
D_pell = 2*2 - 4*(-1)                              # = 8
treffer_pk = []   # `treffer_pk` = for each control prime: (p, e, whether the Lucas condition holds)
for p in (13, 31, 1546463):
    e = legendre(D_pell, p)
    u = lucas_U_schnell(2, -1, p - e, p*p)
    treffer_pk.append((p, e, u % (p*p) == 0))
    sag(f'    p = {p:>9}: e = (8/p) = {e:>2}, U_(p−e) mod p² = {"0  ✅" if u % (p*p) == 0 else "≠ 0  🔴"}')
assert all(x[2] for x in treffer_pk), ('E2 VERLETZT: bekannte Faelle bestehen den Test nicht — '
                                       'das Skript misst nicht die Lucas-Bedingung')
# Counter-check: a prime that is NOT in A238736 must not satisfy the condition
for p in (7, 11, 17, 41, 1009):
    e = legendre(D_pell, p)
    assert lucas_U_schnell(2, -1, p - e, p*p) % (p*p) != 0, ('E2-Gegenprobe: falsch positiv bei', p)
sag('    ✅ drei bekannte Treffer bestaetigt, fuenf Nicht-Treffer bleiben negativ. Der Test misst, was er soll.')

# ---------------------------------------------------------------- the actual comparison
sag()
sag('--- VERGLEICH auf unseren 16 Wieferich-Ereignissen (w77: p² | T_d)')
# `ev` = the 16 Wieferich events (kernel m, rank d, prime p, valuation v)
ev = J('w77_wieferich_ereignisse_result.json')['ereignisse']
sag(f"{'m':>8} {'d (T-Rang)':>11} {'p':>10} {'v_p(T_d)':>9} {'e=(D/p)':>8} {'Lucas-Bed.':>12}")
beide, nur_unsere = 0, 0   # `beide` = events that satisfy BOTH conditions; `nur_unsere` = events that satisfy only ours
details = []               # `details` = one record per comparable event
for x in ev:
    m, d, p, v = x['m'], x['d'], x['p'], x['v']
    T1, U1 = fund(m); P, Q = 2*T1, 1
    D = P*P - 4*Q                                   # = 4*m*U1^2
    if p % 2 == 0 or D % p == 0 or U1 % p == 0:
        sag(f'{m:>8} {d:>11} {p:>10} {v:>9} {"—":>8} {"ausgeschlossen (p | 2QD)":>12}')
        continue
    e = legendre(D, p)
    lb = lucas_U_schnell(P, Q, p - e, p*p) % (p*p) == 0
    beide += lb; nur_unsere += (not lb)
    details.append(dict(m=m, d=d, p=p, v=v, e=e, lucas=lb))
    sag(f'{m:>8} {d:>11} {p:>10} {v:>9} {e:>8} {("JA" if lb else "nein"):>12}')

# ---------------------------------------------------------------- E4: control in the other direction
sag()
sag('--- PK E4: Kontrolle — erfuellen auch GEWOEHNLICHE Rangprimzahlen die Lucas-Bedingung?')
sag('    (Waere sie trivial wahr, saegte der Vergleich oben nichts.)')
w20 = J('w20_primteil_P3000000_N25_NP5_result.json')
# `probe` = results of the sample; `falsch_positiv` = false positives (Lucas condition holds although v_p = 1)
probe, falsch_positiv = [], 0
for x in w20['P1']['tripel'][:400] if 'tripel' in w20.get('P1', {}) else []:
    pass
# Sample from the sieve hits: there v_p = 1, hence NOT a Wieferich case for us
w71 = J('w71_rangsieb_P100000000_N176_result.json')
for m, d, p, v in w71['treffer'][:25]:
    T1, U1 = fund(m); P, Q = 2*T1, 1; D = P*P - 4*Q
    if p % 2 == 0 or D % p == 0: continue
    e = legendre(D, p)
    lb = lucas_U_schnell(P, Q, p - e, p*p) % (p*p) == 0
    probe.append(lb); falsch_positiv += lb
sag(f'    {len(probe)} Primzahlen mit v_p(T_d) = 1 geprueft: {falsch_positiv} erfuellen die Lucas-Bedingung.')
if falsch_positiv == 0:
    sag('    ✅ E4: die Lucas-Bedingung ist NICHT trivial wahr — der Vergleich oben traegt.')
else:
    sag('    🔴 E4 VERLETZT: die Bedingung feuert auch bei Nicht-Wieferich-Primzahlen. Vergleich wertlos.')

# ---------------------------------------------------------------- verdict
sag()
sag('='*100)
n_ge = len(details)   # `n_ge` = number of comparable events
sag(f'  Von {n_ge} vergleichbaren Ereignissen erfuellen {beide} AUCH die Lucas-Bedingung, {nur_unsere} NICHT.')
if beide == n_ge and n_ge:
    urteil = 'IDENTISCH auf dieser Familie — der etablierte Begriff trifft unser Objekt; McIntosh–Roettger gehoert zitiert.'
elif beide == 0:
    urteil = 'VERSCHIEDENE OBJEKTE — der etablierte Begriff gehoert NICHT auf unser Objekt; die Namensgleichheit taeuscht.'
else:
    urteil = ('VERWANDT, NICHT IDENTISCH — die Begriffe ueberschneiden sich teilweise. Das gehoert als eigene '
              'Aussage in den Satz, nicht als Gleichsetzung.')
sag(f'  ⇒ {urteil}')
sag('='*100)
sag('  ⚠️  Gemessen wurde auf 16 Faellen EINER Familie. Das ist ein Befund ueber diese Faelle,')
sag('      kein Satz ueber die Begriffe. Ein Beweis der Aequivalenz waere eine eigene Aufgabe.')

res = dict(skript=pathlib.Path(__file__).name, datum='2026-09-22',
           quelle_begriff='OEIS A001220, Kommentar Felix Froehlich 27.05.2016, nach McIntosh-Roettger 2007 S. 2088',
           pk_pell=[[p, e, bool(o)] for p, e, o in treffer_pk],
           n_vergleichbar=n_ge, beide=beide, nur_unsere=nur_unsere, details=details,
           e4_geprueft=len(probe), e4_falsch_positiv=falsch_positiv, urteil=urteil,
           laufzeit_s=round(time.time()-t0, 1))
(ERG/'w84_ist_es_dasselbe_objekt_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w84_ist_es_dasselbe_objekt_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
print(f'\nErgebnis: w84_ist_es_dasselbe_objekt_result.json   Zeit {time.time()-t0:.0f}s')
