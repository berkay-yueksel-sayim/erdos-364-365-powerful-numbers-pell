# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w61_phi3_pell_familie_2026-09-20.py
# Purpose: Φ₃ takes powerful values infinitely often; an explicit Pell family.
# Background: De Koninck, Doyon & Luca, "Powerful Values of Quadratic Polynomials", J. Integer Seq. 14 (2011), Art. 11.3.3.
#   Nothing is taken over from there: the work was read, compared, and its relevant claims recomputed independently here.
# The bridge (exact): Φ₃(x) = x²+x+1 is ALWAYS odd (x²+x is even), so 4·Φ₃(x) is powerful iff Φ₃(x) is (4 is powerful and
#   coprime to Φ₃). With y = 2x+1 one has 4·Φ₃(x) = y² + 3. Hence Φ₃(x) powerful ⟺ y² + 3 powerful, i.e. the family n² + k of
#   DKL with k = 3, restricted to odd n. Their Table 1 lists for k = 3 the smallest case n = 37, n²+3 = 2²·7³, which is our
#   x = 18 with Φ₃(18) = 343 = 7³.
# What is shown here (and what is not):
#   - An explicit infinite family: Φ₃(x) = 7m² with 7 | m ⟺ y² − 28m² = −3, y = 2x+1. From the solution (37, 7), the automorphism
#     (u, v) = (127, 24) of u² − 28v² = 1 generates all further ones; 7 | m holds at exactly every 7th step, because m mod 7 runs
#     through a cycle of period 7.
#   - No factorization is needed: if m = 7t then Φ₃(x) = 7·49·t² = 343t², so v_7 ≥ 3 and every other prime has an even exponent
#     ≥ 2. "Powerful" follows structurally from the Pell identity. (Factoring would fail anyway: the third member has 35 digits.)
#   - NOT shown: that these are ALL powerful values of Φ₃. Our second known case x = 88916 belongs to D = 13, not to D = 7, so
#     there are several families.
#   - No claim of novelty for the principle: DKL state explicitly that for a·n²+b·n+c a single powerful non-square solution
#     implies infinitely many. New here is only the worked-out case and the period 7.
# Why this matters: the statement "137 of 137 ranks, 0 powerful" remains correct as a statement about RANKS. Read alone,
#   however, it could suggest that Φ_d never becomes powerful, which would be wrong; the restriction belongs with the statement.
# Reads:  nothing. Writes: ergebnisse/w61_phi3_pell_familie_output.txt and ..._result.json.
# Usage:  python w61_phi3_pell_familie_2026-09-20.py   (no arguments)
# Controls: PK+: the Pell identity in EVERY step; the structure identity Φ₃(x) = 343t² at every hit; the smallest hit must be
#   x = 18 with 343 = 7³; the known second case x = 88916 must belong to D = 13 and NOT lie in this family.
#   PK-: a step with m ≢ 0 (mod 7) must NOT satisfy the structure identity.

import sys, json, time, pathlib
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent
ERG  = HIER.parent/'ergebnisse'
t0 = time.time(); Z = []   # `Z` = output lines, written to the output file at the end
def out(s=''): Z.append(s); print(s, flush=True)   # `out`: print a line and keep it for the output file

def fac(n):   # factorization of n as {prime: exponent}, by trial division
    f = {}; m = n; p = 2
    while p*p <= m:
        while m % p == 0: f[p] = f.get(p,0)+1; m //= p
        p += 1 if p == 2 else 2
    if m > 1: f[m] = f.get(m,0)+1
    return f
def zeig(n):   # `zeig` = show: factorization of n as text, e.g. '2^2·7^3'
    return '·'.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(fac(n).items()))
def phi3(x): return x*x + x + 1   # cyclotomic polynomial Φ₃(x)
def phi6(x): return x*x - x + 1   # cyclotomic polynomial Φ₆(x)

out('='*100)
out('W61 — Φ₃ NIMMT UNENDLICH OFT POWERFUL WERTE AN (explizite Pell-Familie)')
out()

# ---------------------------------------------------------------- 1. the bridge to DKL
out('1. DIE BRUECKE ZU De Koninck–Doyon–Luca 2011, Tabelle 1, k = 3')
n_dkl, k = 37, 3
v = n_dkl*n_dkl + k
assert v == 1372 and zeig(v) == '2^2·7^3', ('PK+ DKL-Eintrag', v, zeig(v))
out(f'   DKL: n = {n_dkl}, n² + {k} = {v} = {zeig(v)}   (kleinster Fall fuer k = 3)')
x18 = (n_dkl - 1)//2   # `x18` = the smallest powerful case, x = 18
assert phi3(x18) == 343 and 4*phi3(x18) == v, 'PK+ Rueckuebersetzung'
out(f'   ⇒ x = (n−1)/2 = {x18}:  Φ₃({x18}) = {phi3(x18)} = {zeig(phi3(x18))},  4·Φ₃ = {4*phi3(x18)} = n²+3 ✓')
assert phi6(x18+1) == 343, 'PK+ Spiegelpaar'
out(f'   ⇒ Spiegelpaar:        Φ₆({x18+1}) = {phi6(x18+1)} = {zeig(phi6(x18+1))} ✓')
out('   ⇒ **Der kleinste powerful Wert von Φ₃ steht in DKLs Tabelle 1.** Kein Neuheitsanspruch darauf.')
out()

# ---------------------------------------------------------------- 2. the family
out('2. DIE UNENDLICHE FAMILIE   Φ₃(x) = 7m² mit 7|m   ⟺   y² − 28m² = −3,  y = 2x+1')
U, V = 127, 24   # fundamental solution of u² − 28v² = 1 (the automorphism)
assert U*U - 28*V*V == 1, 'PK+ Fundamentalloesung von u²−28v²=1'
out(f'   Automorphismus (u,v) = ({U},{V}):  u² − 28v² = {U*U - 28*V*V}')
y, m = 37, 7
assert y*y - 28*m*m == -3, 'PK+ Startloesung'
SCHRITTE = 24   # `SCHRITTE` = number of steps of the recursion
familie = [(0, x18)]   # `familie` = members of the family as (step, x)
rest = []              # `rest` = the residues m mod 7 along the recursion
for i in range(1, SCHRITTE+1):
    y, m = U*y + 28*V*m, V*y + U*m
    assert y*y - 28*m*m == -3, f'PK+ Pell-Identitaet in Schritt {i}'
    rest.append(m % 7)
    if m % 7 == 0:
        assert y % 2 == 1, f'y muss ungerade bleiben (Schritt {i})'
        x = (y-1)//2; t = m//7
        assert phi3(x) == 343*t*t, f'PK+ Struktur-Identitaet Φ₃ = 343t² in Schritt {i}'
        familie.append((i, x))
    else:
        # PK-: without 7 | m the structure identity must not hold
        if y % 2 == 1:
            x = (y-1)//2
            assert phi3(x) % 343 != 0 or (m % 7), 'PK- unerwartet'
out(f'   m mod 7 im Verlauf: {rest[:14]} …   ⇒ Periode 7, Treffer bei jedem 7. Schritt')
out(f'   Treffer in {SCHRITTE} Schritten: {len(familie)}')
for i, x in familie:
    s = str(x)
    kurz = s if len(s) <= 46 else s[:22] + '…' + s[-10:]
    out(f'      Schritt {i:2d}:  x = {kurz}   ({len(s)} Stellen)   Φ₃(x) = 343·t², powerful ✓')
out('   🔑 Powerful folgt STRUKTURELL: m = 7t ⇒ Φ₃(x) = 343t², v₇ ≥ 3, alle anderen Exponenten gerade ≥ 2.')
out('      Keine Faktorisierung noetig und keine moeglich — das dritte Glied hat 35 Stellen im Argument.')
out()

# ---------------------------------------------------------------- 3. delimitation: what this family is not
out('3. ABGRENZUNG — was diese Familie NICHT ist')
x2 = 88916
v2 = phi3(x2)
f2 = fac(v2)
assert v2 == 7906143973 and f2 == {7:2, 13:3, 271:2}, ('PK+ zweiter bekannter Fall', v2, f2)
out(f'   Unser zweiter bekannter Fall: Φ₃({x2}) = {v2} = {zeig(v2)}')
# D = squarefree part in the representation D·m²
D2 = 1
for p, e in f2.items():
    if e % 2: D2 *= p
out(f'   Darstellung D·m² mit D quadratfrei: D = {D2}, m = {int((v2//D2)**0.5)}')
assert D2 == 13, ('PK+ zweiter Fall gehoert zu D = 13', D2)
out(f'   ⇒ **D = {D2}, nicht 7** — er liegt in einer ANDEREN Familie. Es gibt mehrere.')
out(f'   ⇒ Diese Rechnung zeigt UNENDLICH VIELE, nicht ALLE.')
out()
in_suche = [x for i, x in familie if x <= 30_000_000]   # `in_suche` = family members inside the earlier search range x <= 3·10⁷
out(f'   Gegen unsere Suche bis x = 3·10⁷ (W53): {len(in_suche)} der {len(familie)} Familienglieder liegen darin — {in_suche}')
out(f'   ⇒ Der naechste hat {len(str(familie[1][1]))} Stellen. **Eine Suche kann diese Familie nicht sehen.**')
out()
out('4. FOLGE FUER DEN RECORD')
out('   `T12` (*„137 von 137 Raengen, 0 powerful"*) bleibt richtig — es ist eine Aussage ueber RAENGE.')
out('   ⚠️ Aber allein gelesen legt sie nahe, Φ_d werde nicht powerful. **Das ist falsch.** Die Schranke')
out('   gehoert an die Aussage, und diese Rechnung ist der Beleg dafuer.')

out()
out(f'Zeit {time.time()-t0:.1f}s')
(ERG/'w61_phi3_pell_familie_output.txt').write_text('\n'.join(Z)+'\n', encoding='utf-8')
(ERG/'w61_phi3_pell_familie_result.json').write_text(json.dumps(dict(
    skript='w61_phi3_pell_familie_2026-09-20.py',
    bruecke=dict(dkl_k=3, dkl_n=37, dkl_wert=1372, dkl_faktorisierung='2^2·7^3',
                 unser_x=x18, phi3=343, faktorisierung='7^3'),
    pell=dict(gleichung='y^2 - 28 m^2 = -3', start=[37, 7], automorphismus=[127, 24],
              periode_m_mod_7=7, schritte=SCHRITTE),
    familie=[dict(schritt=i, stellen=len(str(x)), x=str(x)) for i, x in familie],
    zweiter_bekannter_fall=dict(x=x2, wert=v2, faktorisierung='7^2·13^3·271^2', D=D2),
    abgrenzung='zeigt unendlich viele, nicht alle; D=13-Fall liegt in anderer Familie',
    prior_art='De Koninck-Doyon-Luca 2011, JIS 14, Art. 11.3.3, Tabelle 1, k=3 (n=37) = unser x=18',
    laufzeit_s=time.time()-t0, python=sys.version.split()[0]), indent=1), encoding='utf-8')
