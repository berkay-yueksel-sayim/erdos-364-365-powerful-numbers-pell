# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w58_moduln_und_quadratbaum_2026-09-19.py
# Purpose: two side questions. (1) Can congruence conditions go to higher moduli, and would that help? (2) A square consists of
#   4 triangles; Thales circles on the 4 sides give 4 new squares; taking the areas into account gives formulas for infinite
#   structures. Does that help, and is it related to our tower?
# Reads: nothing. Writes: ergebnisse/w58_moduln_und_quadratbaum_output.txt and ergebnisse/w58_moduln_und_quadratbaum_result.json.
# Usage: python w58_moduln_und_quadratbaum_2026-09-19.py [PMAX]     (default 200000; `PMAX` = prime limit for the density table)
#
# ---------------------------------------------------------------------------------------------------------
# PART A: HOW FAR DO CONGRUENCES CARRY? This is decidable, and the answer is a LIMIT.
#
#   A number n is LOCALLY excluded at p if v_p(n) = 1 (then n is not powerful, Thm 4.5).
#   This can be read off n mod p^2: excluded <=> n = 0 (mod p) and n /= 0 (mod p^2).
#   A triple requires that n-1, n, n+1 are ALL THREE not excluded.
#   => For every modulus M = prod p^2 the set of admissible residues can be computed EXACTLY.
#
#   THE THEORETICAL DENSITY per prime (it is recomputed below):
#     p = 2 : n must be = 0 (mod 4) (if n is odd, n-1 and n+1 are both even and cannot both be divisible by 4,
#             since they differ by 2).                                                  density 1/4
#     p = 3 : exactly one of the three is = 0 (mod 3), and it must be = 0 (mod 9).    density 1/3
#     p >= 5: either none is divisible by p, or exactly one -- and it must be divisible by p^2.
#             density 1 - 3/p + 3/p^2.
#   => total density = 1/4 * 1/3 * prod_{p>=5} (1 - 3/p + 3/p^2).
#   Every factor is STRICTLY POSITIVE, so for EVERY finite modulus the set of residues is NOT EMPTY.
#   At the same time sum 3/p diverges, so the product tends to 0.
#   => ANSWER: higher moduli filter better and better, but they NEVER reach zero.
#      A pure congruence argument can never prove the conjecture.
#
# PK+ : at M = 36 EXACTLY THREE residues must survive -- this is Rizzi's published condition (middle = 0, 8 or 28 mod 36),
#       which our witness reduction contains as the case p = 3 (w21).
# PK- : at M = 4 ONLY the residue 0 may survive (the classical statement "middle = 0 mod 4").
#
# ---------------------------------------------------------------------------------------------------------
# PART B: THE SQUARE TREE, and why it is NOT our tower (but a relative of it).
#
#   Square of side s -> 4 triangles -> Thales circle over each side (diameter s) -> one new square each of side s/sqrt(2),
#   i.e. area s^2/2. Four of them: total area 2 s^2.
#   => Generation n: 4^n squares of area s^2/2^n, total area 2^n s^2. The numbers are 2^n, 4^n, sqrt(2)^n.
#   These are POWERS OF sqrt(2), and sqrt(2) is NOT a unit in Z[sqrt 2] (norm -2).
#   Our tower, by contrast, uses the UNIT 1 + sqrt(2) (norm -1). => different subgroups of the same ring.
#   The geometric figure that really produces our tower is the 1 : sqrt(2) rectangle (the A-series paper format): halving it
#   leaves a similar rectangle. Its continued fraction sqrt(2) = [1;2,2,2,...] gives the convergents p_k/q_k, and those with
#   norm p^2 - 2q^2 = +1 are exactly our tower steps. (The silver rectangle 1 : (1+sqrt 2) = [2;2,2,...] has the same tail but a
#   different head: its denominators match ours (2, 12, 70, 408), its numerators do not.)
#
# PK+ : the convergents of sqrt(2) = [1;2,2,2,...] must contain (3,2), (17,12), (99,70), (577,408) (our tower steps from w56);
#       at least four of them must appear.
# Own code.
import sys, json, time, pathlib
from math import isqrt, gcd
sys.stdout.reconfigure(encoding='utf-8')

HIER = pathlib.Path(__file__).resolve().parent
ERG  = HIER.parent / 'ergebnisse'
t0 = time.time(); Z = []   # `Z` = output lines
# `out` = print a line and record it in `Z`
def out(s=''): Z.append(s); print(s, flush=True)
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200000

# `sieb(n)` = sieve of Eratosthenes: list of the primes up to n
def sieb(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2, isqrt(n)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]

out('='*104)
out('W58 — BERKS ZWEI FRAGEN: hoehere Moduln, und der Quadrat-Baum')
out()
out('   TEIL A — WIE WEIT TRAGEN KONGRUENZEN?')
out('   n ist lokal bei p ausgeschlossen  <=>  v_p(n) = 1  <=>  n = 0 (mod p) und n /= 0 (mod p^2).')
out('   Ein Tripel verlangt: n-1, n, n+1 alle drei NICHT ausgeschlossen. Exakt ausrechenbar je Modul.')
out()

# `zulaessig(primes)` = admissible: returns (modulus M, list of the admissible residues) for the given primes
def zulaessig(primes):
    """exact set of the admissible residues mod prod p^2 for the middle n of a triple."""
    M = 1
    for p in primes: M *= p*p
    ok = []
    for r in range(M):
        gut = True
        for p in primes:
            pp = p*p
            for m in (r-1, r, r+1):
                if m % p == 0 and m % pp != 0: gut = False; break
            if not gut: break
        if gut: ok.append(r)
    return M, ok

out(f"      {'Primzahlen':>16} {'Modul M':>12} {'zulaessige Reste':>17} {'Dichte':>12}  Bemerkung")
# `tabelle` = table rows (primes, modulus M, number of admissible residues, density, first residues); `bem` = remark (printed)
tabelle = []
for ps, bem in (([2], 'klassisch: Mitte = 0 (mod 4)'),
                ([3], ''),
                ([2,3], 'M = 36 -> Rizzis Bedingung'),
                ([2,3,5], ''),
                ([2,3,5,7], '')):
    M, ok = zulaessig(ps)
    d = len(ok)/M
    tabelle.append((ps, M, len(ok), d, ok[:8]))
    out(f'      {str(ps):>16} {M:>12} {len(ok):>17} {d:>12.6f}  {bem}')
M4, ok4 = zulaessig([2])   # PK-: modulus 4
assert ok4 == [0], ('PK-: mod 4 darf nur der Rest 0 ueberleben', ok4)
M36, ok36 = zulaessig([2,3])   # PK+: modulus 36
assert len(ok36) == 3, ('PK+: mod 36 muessen genau 3 Reste ueberleben', ok36)
out()
out(f'      🔑 PK- ok: mod 4 ueberlebt NUR der Rest 0 — die klassische Aussage „Mitte = 0 (mod 4)".')
out(f'      🔑 PK+ ok: mod 36 ueberleben GENAU DREI Reste: {sorted(ok36)}')
out(f'         Das ist **Rizzis publizierte Bedingung** (Mitte = 0, 8 oder 28 mod 36), die wir in W21')
out(f'         schon als Fall p = 3 unserer Zeugenreduktion erkannt haben. Hier faellt sie von selbst heraus.')
out()
out('   DIE DICHTE, wenn man immer mehr Primzahlen mitnimmt:')
out(f"      {'Primzahlen bis':>16} {'Dichte':>14} {'Faktor gegen vorher':>20}")
primes = sieb(PMAX)
d = 0.25 * (1/3)
vor = d
# `marken` = marks: prime thresholds at which the density is printed; `vor` = previous density; `zeiger` = index of the next mark
marken = [5, 10, 100, 1000, 10000, 100000, PMAX]
dichten = []   # `dichten` = densities at the marks
zeiger = 0
for p in primes:
    if p < 5: continue
    d *= (1 - 3/p + 3/(p*p))
    while zeiger < len(marken) and p >= marken[zeiger]:
        out(f'      {marken[zeiger]:>16} {d:>14.3e} {d/vor:>20.4f}')
        dichten.append((marken[zeiger], d)); vor = d; zeiger += 1
out()
out('      🔑 JEDER Faktor (1 - 3/p + 3/p²) ist STRIKT POSITIV. Also ist die Restmenge fuer JEDEN')
out('        endlichen Modul NICHT LEER — es ueberlebt immer mindestens ein Rest.')
out('        Zugleich divergiert Summe 3/p, also geht das Produkt gegen NULL.')
out('      ⇒ ANTWORT AUF FRAGE 1: hoehere Moduln filtern immer besser, aber sie erreichen NIE null.')
out('        **Ein reines Kongruenzargument kann die Vermutung niemals beweisen.**')
out('        Als SUCHFILTER sind sie sehr wohl etwas wert — und wir benutzen schon hohe:')
out('        Cor. 4.8 gibt p = ±1 (mod 4d), und in W30b war d = 43, also **mod 172**; das duennte die')
out('        Probedivision um den Faktor ~86 aus. Das ist die praktische Antwort auf „wuerde das was bringen".')
out()

# ---------------------------------------------------------------------------------------------- PART B
out('   TEIL B — DER QUADRAT-BAUM: was er erzeugt, und was nicht.')
out('      Quadrat Seite s -> Thaleskreis ueber jeder Seite (Durchmesser s) -> neues Quadrat Seite s/√2.')
out(f"      {'Generation':>11} {'Anzahl':>10} {'Seite':>14} {'Flaeche je':>14} {'Gesamtflaeche':>15}")
for n in range(0, 7):
    out(f'      {n:>11} {4**n:>10} {"s/√2^"+str(n):>14} {"s²/2^"+str(n):>14} {"2^"+str(n)+" · s²":>15}')
out('      ⇒ Die Zahlen sind 2^n, 4^n, √2^n. **Reine Potenzen von √2.**')
out('      🔴 Und √2 ist KEINE Einheit in Z[√2] — die Norm ist −2, nicht ±1.')
out('        Unser Turm benutzt die EINHEIT 1 + √2 (Norm −1). ⇒ zwei VERSCHIEDENE Untergruppen desselben Rings.')
out('      ⇒ ANTWORT AUF FRAGE 2, erster Teil: der Quadrat-Baum ist eine unendliche Struktur, ja — aber eine')
out('        GEOMETRISCHE Folge (2^n). Unsere Tuerme sind keine geometrische Folge, sondern eine Einheitengruppe.')
out()
# Why the 1 : sqrt(2) rectangle and not the silver rectangle 1 : (1+sqrt 2): the convergents of the latter, [2;2,2,...], are
#   2/1, 5/2, 12/5, 29/12, ...; their DENOMINATORS agree with ours (2, 12, 70, 408), their NUMERATORS do not. Our tower steps are
#   the convergents of sqrt(2) = [1;2,2,2,...], not of 1+sqrt(2): same tail, different head. The right figure is therefore the
#   1 : sqrt(2) rectangle (A-series format): halving it gives a SIMILAR rectangle, and this self-similarity produces our tower.
out('      ✅ DIE FIGUR, DIE WIRKLICH UNSEREN TURM ERZEUGT: das 1 : √2 - RECHTECK — das DIN-A-Format.')
out('        Halbiert man es, entsteht ein AEHNLICHES Rechteck. Kettenbruch √2 = [1;2,2,2,…].')
out('      🔴 Gegenprobe: das SILBERNE Rechteck 1:(1+√2) = [2;2,2,…] faellt durch die PK:')
out('        dessen Nenner stimmen (2, 12, 70, 408), seine ZAEHLER nicht. Gleicher Schwanz, anderer Kopf.')
out(f"      {'k':>3} {'Naeherungsbruch':>18} {'p_k':>10} {'q_k':>10} {'p²−2q²':>8}  unsere Turmstufe?")
p0, p1 = 1, 1          # √2 = [1; 2,2,2,...]
q0, q1 = 0, 1
# `turm` = known tower steps (x, y) (solutions of x^2 - 2y^2 = 1); `tref` = those found among the convergents (p_k, q_k)
turm = {(3,2), (17,12), (99,70), (577,408), (3363,2378)}
tref = []
for k in range(1, 12):
    n2 = p1*p1 - 2*q1*q1
    ist = (p1, q1) in turm
    if ist: tref.append((p1, q1))
    out(f'      {k:>3} {str(p1)+"/"+str(q1):>18} {p1:>10} {q1:>10} {n2:>8}  {"JA" if ist else ""}')
    p0, p1 = p1, 2*p1 + p0
    q0, q1 = q1, 2*q1 + q0
assert len(tref) >= 4, ('PK+: die Turmstufen (3,2), (17,12), (99,70), (577,408) muessen auftauchen', tref)
out(f'      ✅ PK+ ok: {len(tref)} unserer Turmstufen tauchen als Naeherungsbrueche auf: {tref}')
out('      🔑 Die Normen wechseln −1, +1, −1, +1, … — GENAU die mit Norm +1 sind unsere Turmstufen,')
out('        die mit Norm −1 sind die „halben" Schritte dazwischen (und die negative Pell-Seite aus W55).')
out('      ⇒ ANTWORT AUF FRAGE 2, zweiter Teil: **ja, die Idee traegt — aber nicht mit Quadraten, sondern mit')
out('        dem 1:√2-Rechteck.** Seine unendliche Selbstaehnlichkeit IST unser Turm, Stufe fuer Stufe.')
out('      ⛔ HILFT ES BEIM BEWEIS? Nein. Es ist ein BILD desselben Objekts, kein neues Argument —')
out('        die Kettenbruchentwicklung ist genau das Werkzeug, das wir seit W9 benutzen.')

out(f'\nZeit {time.time()-t0:.0f}s')
(ERG/'w58_moduln_und_quadratbaum_output.txt').write_text('\n'.join(Z)+'\n', encoding='utf-8')
(ERG/'w58_moduln_und_quadratbaum_result.json').write_text(json.dumps(dict(
    skript='w58_moduln_und_quadratbaum_2026-09-19.py', PMAX=PMAX,
    moduln=[(ps, M, n, dd, o) for ps, M, n, dd, o in tabelle],
    dichten=dichten, turmstufen_gefunden=[[int(a), int(b)] for a, b in tref],
    laufzeit_s=time.time()-t0, python=sys.version.split()[0]), indent=1), encoding='utf-8')
