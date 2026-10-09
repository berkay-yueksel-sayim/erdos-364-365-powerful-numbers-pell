# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w51_ballot_errata_gegenprobe_2026-09-19.py
# Purpose: cross-check of our valuation formula (Prop. 4.3) against the erratum to Ballot's 2019 paper
#   (Fibonacci Quart. 57(3), 265-275; erratum 57(4), 366).
# Background: the erratum corrects exactly the standard formula for REGULAR primes, the one we use in Prop. 4.3:
#     printed (p. 265, eq. 2):     v_p(U_n) = (nu - 1) + v_p(n)     for rho | n
#     corrected (erratum p. 366):  v_p(U_n) = nu + v_p(n/rho)       for rho | n
#   with rho the rank of apparition and nu = v_p(U_rho).
# Relevance: with n = rho*t the corrected version reads v_p(U_{rho t}) = v_p(U_rho) + v_p(t), which is literally our Prop. 4.3.
#   The two versions differ by 1 as soon as p does not divide the rank. The script shows (a) that they differ, (b) that OUR
#   version is the corrected one, (c) that it holds on two independent families of sequences.
# Second observation: Ballot's THEOREM concerns SPECIAL primes, i.e. p | gcd(P, Q). Our sequence T_k has Q = 1 (norm of the
#   fundamental unit), so gcd(P, Q) = 1 and p does not divide Q for EVERY prime. Hence there are no special primes in our setting
#   and Ballot's theorem is not applicable; only the erratum to his introductory formula touches our area, and it confirms ours.
# Reads:  nothing. Writes: ergebnisse/w51_ballot_errata_gegenprobe_output.txt and ..._result.json.
# Usage:  python w51_ballot_errata_gegenprobe_2026-09-19.py   (no arguments)
# Controls: PK+ (Fibonacci, p = 11): rho = 10, U_10 = 55 = 5*11, so nu = 1. The PRINTED version gives 0 there and is therefore
#   demonstrably wrong; the corrected one gives 1. PK-: a deliberately wrong rank must make the check fail (otherwise it checks
#   nothing).

import sys, json, time, pathlib
from math import gcd
sys.stdout.reconfigure(encoding='utf-8')

HIER = pathlib.Path(__file__).resolve().parent
ERG  = HIER.parent / 'ergebnisse'
t0 = time.time(); Z = []   # `Z` = output lines, written to the output file at the end
def out(s=''): Z.append(s); print(s, flush=True)   # `out`: print a line and keep it for the output file

def v_p(n, p):   # p-adic valuation of n (None for n = 0)
    if n == 0: return None
    c = 0
    while n % p == 0: n //= p; c += 1
    return c

def lucas_U(n, P, Q):
    """U_n of the Lucas sequence for (P, Q), exact, iterative."""
    a, b = 0, 1                      # U_0, U_1
    for _ in range(n): a, b = b, P*b - Q*a
    return a
def T_k(k, T1, m):
    """T-coordinate of eps^k with eps = T1 + U1*sqrt(m), N(eps) = 1.  T_{k+1} = 2 T1 T_k - T_{k-1}."""
    a, b = 1, T1                     # T_0 = 1, T_1
    for _ in range(k-1): a, b = b, 2*T1*b - a
    return b if k >= 1 else 1
def rang(p, folge, grenze=4000):   # `rang` = rank of apparition: the smallest n <= `grenze` with p | folge(n), else None
    for n in range(1, grenze+1):
        x = folge(n)
        if x % p == 0: return n
    return None

out('='*100)
out('W51 — UNSERE BEWERTUNGSFORMEL GEGEN BALLOTS ERRATA (Fibonacci Quart. 57(4), 366)')
out()
out('   gedruckt   (Ballot 2019, S. 265, Gl. 2):   v_p(U_n) = (nu - 1) + v_p(n)     fuer rho | n')
out('   korrigiert (Errata 2019, S. 366)      :   v_p(U_n) = nu + v_p(n/rho)       fuer rho | n')
out('   unsere Prop. 4.3                      :   v_p(T_{rt}) = v_p(T_r) + v_p(t)  fuer ungerades t')
out()

# ---------------------------------------------------------------- PK+ : the case where the printed version fails
fib = lambda n: lucas_U(n, 1, -1)   # `fib` = Fibonacci numbers
p = 11; rho = rang(p, fib); nu = v_p(fib(rho), p)
# `gedruckt` = printed formula, `korr` = corrected formula, `echt` = actual valuation
gedruckt  = (nu - 1) + v_p(rho, p)
korr      = nu + v_p(rho // rho, p)
echt      = v_p(fib(rho), p)
out(f'   PK+ Fibonacci, p = {p}: rho = {rho}, U_rho = {fib(rho)} = {fib(rho)//p}*{p}, also nu = {nu}')
out(f'       echt {echt}   |   gedruckte Fassung {gedruckt}   |   korrigierte Fassung {korr}')
assert echt == korr and gedruckt != echt, ('PK+: die gedruckte Fassung muss hier FALSCH sein', echt, gedruckt, korr)
out(f'       ⇒ die gedruckte Fassung ist an dieser Stelle nachweislich falsch, die korrigierte richtig.')
out()

# ---------------------------------------------------------------- (A) our own sequence, Q = 1
out('   (A) UNSERE T-FOLGE (Q = 1, also gibt es KEINE speziellen Primzahlen im Sinne Ballots):')
out(f"      {'m':>5} {'T1':>6} {'p':>6} {'rho':>5} {'nu':>3} | {'t':>4} {'v_p(T_rt) echt':>15} {'unsere Formel':>14} {'gedruckt':>9}")
# `faelle` = all cases tested; `treffer_unser` / `treffer_gedruckt` = cases in which our / the printed formula is correct;
# `gesamt` = number of cases
faelle = []; treffer_unser = 0; treffer_gedruckt = 0; gesamt = 0
KERNE = [(2,3),(3,2),(5,9),(6,5),(7,8),(10,19),(11,10),(13,649)]   # `KERNE` = kernels: (m, T1) from x^2 - m y^2 = 1
for m, T1 in KERNE:
    assert T1*T1 - m*((T1*T1-1)//m if (T1*T1-1)%m==0 else 0) >= 0
    f = lambda n, T1=T1, m=m: T_k(n, T1, m)
    for p in (3,5,7,11,13):
        if T1 % p == 0 and p == 3: pass
        r = rang(p, f, 300)
        if r is None: continue
        nu = v_p(f(r), p)
        if nu is None: continue
        for t in (3,5,7,9):
            n = r*t
            if n > 260: continue
            echt = v_p(f(n), p)
            if echt is None: continue
            unser = nu + v_p(t, p)         # `unser` = our formula (= the corrected one)
            gedr  = (nu - 1) + v_p(n, p)   # `gedr` = the printed formula
            gesamt += 1
            if echt == unser: treffer_unser += 1
            if echt == gedr:  treffer_gedruckt += 1
            if len(faelle) < 8:
                out(f'      {m:>5} {T1:>6} {p:>6} {r:>5} {nu:>3} | {t:>4} {echt:>15} {unser:>14} {gedr:>9}')
            faelle.append((m,T1,p,r,nu,t,echt,unser,gedr))
out(f'      ... insgesamt {gesamt} Faelle')
out(f'      UNSERE Formel korrekt in {treffer_unser} von {gesamt}   ({100*treffer_unser/gesamt:.1f} %)')
out(f'      GEDRUCKTE Fassung korrekt in {treffer_gedruckt} von {gesamt}   ({100*treffer_gedruckt/gesamt:.1f} %)')
assert treffer_unser == gesamt, ('unsere Formel muss ueberall halten', treffer_unser, gesamt)
out()

# ---------------------------------------------------------------- (B) second family: Fibonacci (Q = -1)
out('   (B) ZWEITE, UNABHAENGIGE FAMILIE — Fibonacci (Q = -1, andere Folge, dieselbe Gesetzmaessigkeit):')
g2 = g2g = 0   # `g2` = cases tested; `g2g` = cases in which the corrected formula holds
for p in (3,7,11,13,17,19,23,29,31,37):
    r = rang(p, fib, 200)
    if r is None: continue
    nu = v_p(fib(r), p)
    for t in (2,3,4,5,7):
        n = r*t
        if n > 400: continue
        echt = v_p(fib(n), p)
        if echt is None: continue
        g2 += 1
        if echt == nu + v_p(t, p): g2g += 1
out(f'      korrigierte Fassung korrekt in {g2g} von {g2} Faellen ({100*g2g/g2:.1f} %)')
assert g2g == g2, ('auch bei Fibonacci muss die korrigierte Fassung ueberall halten', g2g, g2)

# ---------------------------------------------------------------- PK-
r = rang(11, fib); nu = v_p(fib(r), 11)
falsch = sum(1 for t in (3,5,7) if v_p(fib((r+1)*t), 11) == nu + v_p(t, 11))
assert falsch < 3, 'PK-: mit falschem Rang darf die Formel NICHT durchgehend stimmen'
out(f'      PK- ok: mit einem absichtlich falschen Rang ({r+1} statt {r}) stimmt sie nur in {falsch} von 3 Faellen.')

out()
out('   ⇒ ERGEBNIS: unsere Prop. 4.3 ist die KORRIGIERTE Fassung. Die Errata bestaetigt uns, sie widerlegt uns nicht.')
out('     Ballots SATZ (spezielle Primzahlen, p | gcd(P,Q)) ist bei uns gegenstandslos, weil Q = 1.')
out(f'\nZeit {time.time()-t0:.0f}s')
(ERG/'w51_ballot_errata_gegenprobe_output.txt').write_text('\n'.join(Z)+'\n', encoding='utf-8')
(ERG/'w51_ballot_errata_gegenprobe_result.json').write_text(json.dumps(dict(
    skript='w51_ballot_errata_gegenprobe_2026-09-19.py',
    unsere_folge=dict(faelle=gesamt, unsere_formel_korrekt=treffer_unser, gedruckte_fassung_korrekt=treffer_gedruckt),
    fibonacci=dict(faelle=g2, korrigierte_fassung_korrekt=g2g),
    laufzeit_s=time.time()-t0, python=sys.version.split()[0]), indent=1), encoding='utf-8')
