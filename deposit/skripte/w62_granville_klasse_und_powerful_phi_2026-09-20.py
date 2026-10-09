# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w62_granville_klasse_und_powerful_phi_2026-09-20.py
# Purpose: Granville's class boundary (Granville 2012, arXiv:1212.6306, Theorem 3 and Corollary 3) tested against our witness
#   question, and a systematic table of powerful values Φ_d(x).
# Background: Theorem 3: for b, c coprime with c ≡ 2 (mod 4) and Δ = b²+4c > 0, x_n has for n ≠ 1, 2, 3, 6 a primitive prime
#   divisor that divides x_n to EXACTLY an ODD power. Corollary 3 ties this to the primitive number: p divides x_{n_p} and
#   φ_{n_p} to the same power.
# Two design points that explain the code:
#   (a) Bounds are placed on the VALUE, not on x. deg Φ_d = φ(d), and φ(11) = 10: the same x-bound means a value around 10^6
#       at d = 3 and around 10^30 at d = 11. A common x-bound across several degrees is not a bound.
#   (b) "has a prime with odd exponent" is tested through its negation: all exponents even ⟺ the number is a SQUARE, so
#       `isqrt(v)**2 == v` replaces factorization (which timed out). This checks a NECESSARY consequence of Granville's
#       statement, not the statement itself: he asserts a *primitive* prime divisor with odd power, which is stronger. A square
#       refutes it; a non-square does not confirm it.
# No foreign code: Granville was read, the statements were recomputed, nothing was taken over.
# Reads:  nothing. Writes: ergebnisse/w62_granville_klasse_output.txt and ..._result.json.
# Usage:  python w62_granville_klasse_und_powerful_phi_2026-09-20.py   (no arguments; needs numpy)
# Controls: PK+: Granville's exception n = 6 at r = 2, s = 1 (r²−rs+s² = 3 = 3·1²), n = 2 is none (r+s = 3), x₃·x₆ = 441 = 21²;
#   Φ₃(18) = 343 and Φ₅(3) = 121 must appear in the table. PK-: the table must be THIN (< 1 % of the pairs tested), otherwise the
#   powerful test is broken.

import sys, json, time, pathlib
from math import isqrt
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent
ERG = HIER.parent/'ergebnisse'
t0 = time.time(); Z = []   # `Z` = output lines, written to the output file at the end
def out(s=''): Z.append(s); print(s, flush=True)   # `out`: print a line and keep it for the output file

def fac(n):   # factorization of |n| as {prime: exponent}, by trial division
    f = {}; m = abs(n); p = 2
    while p*p <= m:
        while m % p == 0: f[p] = f.get(p,0)+1; m //= p
        p += 1 if p == 2 else 2
    if m > 1: f[m] = f.get(m,0)+1
    return f
def zeig(n):   # `zeig` = show: factorization of n as text, e.g. '2^2·7^3'
    return '·'.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(fac(n).items())) or '1'
def quadrat(n): return isqrt(n)**2 == n   # `quadrat` = n is a perfect square
def squarefree(n): return all(e == 1 for e in fac(n).values())
def powerful_menge(N):   # `powerful_menge` = set of all powerful numbers <= N, enumerated as a²·b³ with b squarefree
    teile = []; b = 1
    while b**3 <= N:
        if squarefree(b):
            amax = isqrt(N//b**3)
            if amax:
                a = np.arange(1, amax+1, dtype=np.int64); teile.append(a*a*(b**3))
        b += 1
    return set(np.unique(np.concatenate(teile)).tolist())
TEILER = {d: [e for e in range(1, d) if d % e == 0] for d in range(1, 16)}   # `TEILER` = proper divisors e < d of d, for d < 16
def phi(d, x):   # cyclotomic polynomial value Φ_d(x): (x^d − 1) divided by Φ_e(x) for all proper divisors e of d
    z = x**d - 1
    for e in TEILER[d]: z //= phi(e, x)
    return z

out('='*100)
out('W62 — GRANVILLES KLASSENGRENZE und powerful Φ_d(x)')
out()

# ---------------------------------------------------------------- 1
out('1. IST DIE KLASSE RICHTIG GELESEN?  Gegenprobe an Granvilles EIGENEM Beispiel r = 2, s = 1')
r, s = 2, 1   # Granville's own example
assert (r*s) % 2 == 0 and (r*s) % 4 != 0, 'PK+ Klassenbedingung 2|rs, 4∤rs'
q2, q6 = r+s, r*r - r*s + s*s
assert not quadrat(q2), 'PK+ n=2 darf KEINE Ausnahme sein (r+s = 3 ist kein Quadrat)'
assert q6 % 3 == 0 and quadrat(q6//3), 'PK+ n=6 MUSS Ausnahme sein (r²−rs+s² = 3 = 3·1²)'
xs = {n: (r**n - s**n)//(r-s) for n in range(1, 8)}
pr = xs[3]*xs[6]
assert pr == 441 and quadrat(pr), ('PK+ x₃·x₆ = 441 = 21²', pr)
out(f'   Klasse erfuellt. n = 2 keine Ausnahme (r+s = {q2}); n = 6 Ausnahme (r²−rs+s² = {q6} = 3·1²).')
out(f'   x₆ = {xs[6]} = {zeig(xs[6])} — beide Primteiler kommen frueher vor ⇒ kein charakteristischer Primteiler.')
out(f'   x₃·x₆ = {xs[3]}·{xs[6]} = {pr} = {isqrt(pr)}²   ✓ genau das Beispiel aus seiner Einleitung.')
out()

# ---------------------------------------------------------------- 2
out('2. DIE KLASSE IN UNSERER SCHREIBWEISE:  x_n = (xⁿ−1)/(x−1), also r = x, s = 1')
out('   "2 | rs und 4 ∤ rs"  ⟺  **x ≡ 2 (mod 4)**')
out('   🔑 Unsere beiden powerful Φ₃-Funde liegen auf VERSCHIEDENEN Seiten dieser Grenze:')
out('      x = 18      (mod 4 = 2)  →  IN der Klasse')
out('      x = 88 916  (mod 4 = 0)  →  AUSSERHALB')
out()

# ---------------------------------------------------------------- 3
out('3. IST DIE KLASSE WESENTLICH ODER DEKORATIV?')
out('   Test: ist Φ_d(x) ein QUADRAT? Dann hat es keine Primzahl mit ungeradem Exponenten,')
out('   und Granvilles Schluss kann dort nicht gelten. (Notwendige Folge, nicht die Aussage selbst.)')
DS = [4, 5, 7, 8, 9, 10, 11, 12]          # `DS` = degrees d tested (Granville's exceptions n = 1, 2, 3, 6 excluded)
VCAP = 10**14                             # `VCAP` = bound on the value Φ_d(x)
innen, aussen = [], []                    # squares inside / outside the class (x ≡ 2 mod 4), as (x, d, value)
for d in DS:
    x = 2
    while True:
        v = phi(d, x)
        if v > VCAP: break
        if v > 1 and quadrat(v):
            (innen if x % 4 == 2 else aussen).append((x, d, v))
        x += 1
out(f'   Wertschranke Φ_d(x) ≤ 1e14, d ∈ {DS}')
out(f'      INNERHALB der Klasse (x ≡ 2 mod 4): {len(innen)} Quadrate')
for x, d, v in innen[:5]: out(f'         x = {x}, d = {d}:  Φ_d = {v} = {zeig(v)}')
out(f'      AUSSERHALB:                         {len(aussen)} Quadrate')
for x, d, v in aussen[:5]: out(f'         x = {x}, d = {d}:  Φ_d = {v} = {zeig(v)}')
assert not innen, ('Innerhalb der Klasse darf keines auftreten', innen[:3])
assert aussen, 'Ausserhalb wird mindestens eines erwartet, sonst waere die Klasse dekorativ'
out('   ⇒ **Die Klassenbedingung ist WESENTLICH, nicht dekorativ.** Innen haelt es, aussen bricht es.')
out(f'   ⇒ Kleinstes Gegenbeispiel: x = {aussen[0][0]}, d = {aussen[0][1]}, Φ_d = {aussen[0][2]} = {zeig(aussen[0][2])}')
out()

# ---------------------------------------------------------------- 4
out('4. TAFEL: powerful Werte von Φ_d(x)')
VMAX = 10**12                  # `VMAX` = bound on the value Φ_d(x) for the table
PW = powerful_menge(VMAX)      # `PW` = all powerful numbers up to VMAX
out(f'   ⚠️ Begrenzt nach dem WERT (Φ_d(x) ≤ 1e12), nicht nach x — deg Φ_d = φ(d), eine gemeinsame')
out(f'      x-Schranke waere ueber die Grade hinweg keine Schranke.')
out(f'   powerful-Menge bis {VMAX:,}: {len(PW):,} Elemente.')
# `tafel` = table of hits (d, x, value); `gepr` = number of pairs (d, x) tested; `spanne` = largest x reached per degree d
tafel = []; gepr = 0; spanne = {}
for d in range(3, 13):
    x = 2
    while True:
        v = phi(d, x)
        if v > VMAX: break
        gepr += 1
        if v > 1 and v in PW: tafel.append((d, x, v))
        x += 1
    spanne[d] = x - 1
out('   x-Reichweite je Grad: ' + ' · '.join(f'd={d}: {spanne[d]:,}' for d in sorted(spanne)))
out(f'   Geprueft: {gepr:,} Paare.   Treffer: {len(tafel)}')
for d, x, v in tafel: out(f'      d = {d:>2},  x = {x:>7}:   Φ_d(x) = {v:>14,} = {zeig(v)}')
assert len(tafel) < gepr/100, 'PK- die Tafel muss duenn sein'
assert (3, 18, 343) in tafel, 'PK+ Φ₃(18) = 343'
assert any(d == 5 and x == 3 for d, x, v in tafel), 'PK+ Φ₅(3) = 121'
gr = sorted({d for d, x, v in tafel})   # `gr` = degrees d that occur in the table
out(f'   ⇒ {len(tafel)} von {gepr:,} ({100*len(tafel)/gepr:.4f} %). **Die Tafel ist duenn.**')
out(f'   🔑 **Betroffene Grade: d ∈ {gr}** — also NICHT nur die quadratischen Faelle d ∈ {{3, 4, 6}}.')
out()
out('5. WAS DAS FUER DEN RECORD HEISST')
out('   `T12` (*„137 von 137 Raengen, 0 powerful"*) ist eine Aussage ueber RAENGE und bleibt richtig.')
out('   Diese Tafel ist eine Aussage ueber ARGUMENTE. ⇒ **Zwei verschiedene Fragen, und nur die erste')
out('   haben wir beantwortet.** Die Tafel steht hier, damit die zweite nicht stillschweigend')
out('   mitbeantwortet scheint.')

out()
out(f'Zeit {time.time()-t0:.1f}s')
(ERG/'w62_granville_klasse_output.txt').write_text('\n'.join(Z)+'\n', encoding='utf-8')
(ERG/'w62_granville_klasse_result.json').write_text(json.dumps(dict(
    skript='w62_granville_klasse_und_powerful_phi_2026-09-20.py',
    granville=dict(quelle='arXiv:1212.6306, Thm 3 + Kor 3', klasse_original='2|rs, 4∤rs',
                   klasse_unsere_schreibweise='x ≡ 2 (mod 4)', ausnahmen_thm3=[1, 2, 3, 6],
                   eigenbeispiel=dict(r=2, s=1, x3_mal_x6=441, wurzel=21)),
    klassengrenze=dict(wert_schranke=VCAP, grade=DS, quadrate_innerhalb=len(innen),
                       quadrate_ausserhalb=len(aussen),
                       kleinstes_ausserhalb=dict(x=aussen[0][0], d=aussen[0][1], wert=aussen[0][2]),
                       test='Quadrat-Test als notwendige Folge, nicht die Aussage selbst'),
    tafel=dict(wert_schranke=VMAX, grade='3..12', x_reichweite=spanne, geprueft=gepr,
               treffer=len(tafel), betroffene_grade=gr,
               eintraege=[dict(d=d, x=x, wert=v, faktorisierung=zeig(v)) for d, x, v in tafel]),
    abgrenzung='Tafel ueber ARGUMENTE; T12 ist eine Aussage ueber RAENGE.',
    laufzeit_s=time.time()-t0, python=sys.version.split()[0]), indent=1), encoding='utf-8')
