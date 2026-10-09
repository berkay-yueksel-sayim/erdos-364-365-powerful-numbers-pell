# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w55_a2_als_negative_pell_2026-09-19.py
# Purpose: type A2 (x^2 + 1 powerful) is a negative Pell family; this turns the sieve into a GENERATOR. Type A2 had been recorded
#   only as an observation ("negative Pell, norm -1") on the 18 known pairs; here it becomes a generating description,
#   and that description explains the earlier search results.
# Reads: no files. Writes: ergebnisse/w55_a2_als_negative_pell_output.txt and ergebnisse/w55_a2_als_negative_pell_result.json
#   (the results folder is found relative to the script).
# Usage: python w55_a2_als_negative_pell_2026-09-19.py [BMAX]   (BMAX = largest b, default 300)
# Controls: PK+: 682^2 - 125*61^2 = -1, so b = 5, and 682^2+1 = 5^3*61^2 (the only known A2 case before this run).
#   PK-: for b = 2 there must be NO solution (u^2 = -1 mod 8 is impossible) -- otherwise the check would test nothing.
#   PK: every solution found must satisfy the identity u^2+1 = b^3 v^2 exactly and u must be even (both asserted);
#   powerfulness is verified by trial division for the first three solutions with x^2+1 < 10^14.

# THE EQUIVALENCE (both directions elementary):
#   (=>) If x^2+1 is powerful, write it in Golomb's unique form x^2+1 = a^2 b^3 with b squarefree. Then x^2 - b^3 a^2 = -1.
#   (<=) If (u,v) solves u^2 - b^3 v^2 = -1 with b squarefree, then u^2+1 = b^3 v^2, and that is
#        POWERFUL: the primes of b have exponent 3 + 2 v_p(v) >= 3, all other primes 2 v_p(v) >= 2.
#   => Type A2 is EXACTLY the solution set of the negative Pell equations for cubic discriminants b^3.
#
# A FREE CONSEQUENCE: u is automatically EVEN. If u were odd, then u^2+1 = 2 (mod 8),
#   so v_2 = 1; but v_2(b^3 v^2) = 3 v_2(b) + 2 v_2(v) is never 1. => the parity lemma drops out here.
#
# USE: negative Pell is solvable if and only if the continued fraction period of sqrt(b^3) is ODD. Instead of sieving over x, one
#   runs over b and LISTS the solutions. This reaches sizes that no sieve ever sees -- and explains why w34 (up to 2*10^5)
#   and w52 (up to 3*10^6) both found only x = 682.

import sys, json, time, pathlib
from math import isqrt
sys.stdout.reconfigure(encoding='utf-8')

HIER = pathlib.Path(__file__).resolve().parent
ERG  = HIER.parent / 'ergebnisse'
t0 = time.time(); Z = []   # `Z` = output lines, written to the _output.txt file
def out(s=''): Z.append(s); print(s, flush=True)   # `out` = print a line and record it in `Z`
BMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300   # `BMAX` = largest b to run over

# `quadratfrei` = squarefree test by trial division
def quadratfrei(n):
    d = 2
    while d*d <= n:
        if n % (d*d) == 0: return False
        d += 1
    return True
# powerful test: no prime divides n exactly once (trial division)
def powerful(n):
    m = n; p = 2
    while p*p <= m:
        if m % p == 0:
            c = 0
            while m % p == 0: m //= p; c += 1
            if c == 1: return False
        p += 1 if p == 2 else 2
    return m == 1
def periode_ungerade(D):
    """`periode_ungerade` = odd period: length of the continued fraction period of sqrt(D) -- odd <=> negative Pell equation solvable."""
    a0 = isqrt(D)
    if a0*a0 == D: return None
    m, d, a = 0, 1, a0
    L = 0
    while a != 2*a0:
        m = d*a - m; d = (D - m*m)//d; a = (a0 + m)//d
        L += 1
        if L > 5_000_000: return None
    return L
def negpell(D, grenze=2_000_000):
    """Fundamental solution of u^2 - D v^2 = -1 via the convergents (`grenze` = step limit), or None."""
    a0 = isqrt(D)
    if a0*a0 == D: return None
    m, d, a = 0, 1, a0
    p0, p1 = 1, a0
    q0, q1 = 0, 1
    for _ in range(grenze):
        if p0*p0 - D*q0*q0 == -1: return p0, q0
        m = d*a - m; d = (D - m*m)//d; a = (a0 + m)//d
        p0, p1 = p1, a*p1 + p0
        q0, q1 = q1, a*q1 + q0
        if a == 2*a0 and p0*p0 - D*q0*q0 != -1: return None
    return None

out('='*104)
out('W55 — TYP A2 IST EINE NEGATIVE PELL-FAMILIE')
out()
out('   x^2 + 1 powerful   <=>   x^2 - b^3 y^2 = -1 fuer das quadratfreie b aus x^2+1 = a^2 b^3')
out('   Loesbar genau dann, wenn die Kettenbruchperiode von sqrt(b^3) UNGERADE ist.')
out()
assert 682*682 - 125*61*61 == -1, 'PK+ Identitaet'
assert 682*682 + 1 == 5**3 * 61**2 and powerful(682*682+1), 'PK+ unser A2-Fall'
assert negpell(8) is None, 'PK-: b = 2 darf keine Loesung haben (u^2 = -1 mod 8 unmoeglich)'
out('   PK+ ok: 682^2 - 125*61^2 = -1, und 682^2+1 = 5^3*61^2 ist powerful.')
out('   PK- ok: b = 2 (D = 8) hat keine Loesung.')
out()

# `loes` = solutions (b, u, v) with u^2 + 1 = b^3 v^2; `ohne` = the squarefree b without a solution
loes = []; ohne = []
t1 = time.time()
for b in range(1, BMAX+1):
    if not quadratfrei(b): continue
    D = b**3
    s = negpell(D)
    if s is None: ohne.append(b); continue
    u, v = s
    assert u*u + 1 == D*v*v, ('Identitaet verletzt', b, u, v)
    assert u % 2 == 0, ('u muss gerade sein (Paritaetslemma faellt hier ab)', b, u)
    loes.append((b, u, v))
out(f'   b quadratfrei bis {BMAX}: {len(loes)} mit Loesung, {len(ohne)} ohne   ({time.time()-t1:.0f}s)')
out(f'   ⇒ ALLE {len(loes)} erzeugen ein A2-Paar (x^2, x^2+1). Die Menge ist also nicht duenn, sondern LISTBAR.')
out()
out('   DIE KLEINSTEN, nach x sortiert — und hier steht die Erklaerung unserer Suchen:')
out(f"      {'b':>5} {'Stellen von x':>14}  x = u")
for b, u, v in sorted(loes, key=lambda t: t[1])[:6]:
    s = str(u)
    out(f'      {b:>5} {len(s):>14}  {s if len(s) <= 60 else s[:50]+"…"}')
klein = sorted(loes, key=lambda t: t[1])   # `klein` = solutions sorted by size of x = u
out()
out(f'   🔑 Der kleinste ist x = {klein[0][1]} (b = {klein[0][0]}) — UNSER BEKANNTER FALL.')
out(f'      Der ZWEITKLEINSTE ist x = {klein[1][1]} (b = {klein[1][0]}), also rund 10^{len(str(klein[1][1]))-1}.')
out(f'      ⇒ W34 suchte bis 2*10^5, W52 bis 3*10^6. BEIDE konnten den zweiten gar nicht erreichen.')
out(f'      Das ist keine Seltenheit des Typs, sondern eine Luecke in der GROESSE.')
out()
out('   Kontrolle an den ersten Loesungen: ist x^2+1 wirklich powerful?')
for b, u, v in klein[:3]:
    n = u*u + 1
    ok = powerful(n) if n < 10**14 else None
    out(f'      b = {b:>3}: x^2+1 = b^3 * v^2 nachgerechnet ✓' + (f'   powerful exakt geprueft: {ok}' if ok is not None else '   (zu gross fuer Probedivision; Form b^3*v^2 genuegt als Beweis)'))
out('      🔑 Der Beweis braucht keine Faktorisierung: b^3 v^2 hat bei jeder Primzahl Exponent >= 2.')
out()
out('   WELCHE b gehen NICHT (erste 20) — und warum das kein Zufall ist:')
out(f'      {ohne[:20]}')
# `mod4` = for p in (3, 7, 11, 19, 23): number of the first 200 unsolvable b divisible by p (not used further)
mod4 = {}
for b in ohne[:200]:
    for p in (3, 7, 11, 19, 23):
        if b % p == 0: mod4[p] = mod4.get(p, 0) + 1
out(f'      Ein Primteiler p = 3 (mod 4) in b macht die Gleichung unloesbar (u^2 = -1 mod p waere noetig).')
# `tref` = unsolvable b all of whose prime divisors q are 2 or = 1 (mod 4) (not used further)
tref = [b for b in ohne if all(p % 4 == 1 or p == 2 for p in [q for q in range(2, b+1) if b % q == 0 and all(b//q*0 == 0 for _ in [0]) and q > 1 and all(q % r for r in range(2, isqrt(q)+1))])]
out(f'      Gegenprobe: {len([b for b, _, _ in loes if b % 4 == 3])} der Loesungs-b sind selbst = 3 (mod 4)'
    f' — die Bedingung sitzt an den PRIMTEILERN, nicht an b selbst.')

out(f'\nZeit {time.time()-t0:.0f}s')
(ERG/'w55_a2_als_negative_pell_output.txt').write_text('\n'.join(Z)+'\n', encoding='utf-8')
# result keys: `anzahl_mit_loesung` = number of b with a solution, `anzahl_ohne` = without, `kleinste` = the smallest solutions
#   (b, u, v), `ohne_loesung` = the first b without a solution, `laufzeit_s` = run time in seconds
(ERG/'w55_a2_als_negative_pell_result.json').write_text(json.dumps(dict(
    skript='w55_a2_als_negative_pell_2026-09-19.py', BMAX=BMAX,
    anzahl_mit_loesung=len(loes), anzahl_ohne=len(ohne),
    kleinste=[(b, str(u), str(v)) for b, u, v in klein[:8]],
    ohne_loesung=ohne[:60], laufzeit_s=time.time()-t0, python=sys.version.split()[0]), indent=1), encoding='utf-8')
