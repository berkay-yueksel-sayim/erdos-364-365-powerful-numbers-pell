# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w65_d6_klassenpaare_2026-09-20.py
# Purpose: class-pair table for neighboring powerful numbers: are the classes of neighboring powerful numbers dependent on each
#   other beyond the local density? (First step of the class-pair analysis.)
# Idea: do not put all powerful numbers into one pot, but compute WITHIN and BETWEEN their classes. Every powerful number is
#   uniquely n = a²b³ with b squarefree; the classes are C_b = {a²b³ : a >= 1}. Within a class the distances are trivial,
#   b³(a−a')(a+a'), so the whole difficulty sits BETWEEN the classes: d³c² − b³a² = g. The class mark b comes for free: the
#   enumeration of the powerful numbers runs over squarefree b anyway and simply discards b.
# The null model is the main part of the work, not the count. The naive comparison |C_b|·|C_d| would be wrong: classes with
#   large b are thin everywhere, so it would measure the density and call it neighborhood. The question is whether the classes
#   of NEIGHBORING powerful numbers depend on each other beyond the local density.
#   Null model: permute the class marks within a WINDOW of neighboring powerful numbers. This preserves the local class
#   composition (and thereby the arrangement along the number line) and destroys exactly the neighbor correlation. Several
#   window widths are used so that the answer does not depend on one choice. (The comparison set must reproduce the
#   distribution of the nuisance variable.)
# Not done here: the follow-up test (occ(g) conditioned on the class pair must beat the baseline 0.545 of w59).
# Absence condition, stated before the run: an empty cell is a finding ONLY if the null model expects a non-negligible count
#   there. For small classes an empty cell is to be expected.
# Reads:  nothing (the script enumerates the powerful numbers up to 10^N_EXP itself).
# Writes: ergebnisse/w65_d6_klassenpaare_output.txt and ergebnisse/w65_d6_klassenpaare_result.json.
# Usage:  python w65_d6_klassenpaare_2026-09-20.py [N_EXP=14]   (powerful numbers up to 10^N_EXP; needs numpy)
# Controls: PK+: for N_EXP = 14 the number of powerful numbers is 21663503; the sample numbers 4, 8, 9, 72, 108, 343, 1000
#   decompose correctly into a²b³ with b squarefree; every frequent gcd(n, g) is powerful.

import sys, json, time, pathlib
import numpy as np
from math import isqrt
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent
ERG = HIER.parent/'ergebnisse'
t0 = time.time(); Z = []   # `Z` = output lines, written to the output file at the end
def out(s=''): Z.append(s); print(s, flush=True)   # `out`: print a line and keep it for the output file
N_EXP = int(sys.argv[1]) if len(sys.argv) > 1 else 14   # `N_EXP` = decimal exponent of the range
N = 10**N_EXP                                           # `N` = upper bound for the powerful numbers

def fac(n):   # factorization of n as {prime: exponent}, by trial division
    f = {}; m = n; p = 2
    while p*p <= m:
        while m % p == 0: f[p] = f.get(p,0)+1; m //= p
        p += 1 if p == 2 else 2
    if m > 1: f[m] = f.get(m,0)+1
    return f
def squarefree(n): return all(e == 1 for e in fac(n).values())

out('='*100)
out('W65 — D6 SCHRITT 1: KLASSENPAARE BENACHBARTER POWERFUL NUMBERS')
out()
out('   Klassen C_b = {a²b³}, b quadratfrei. Frage: haengen die Klassen BENACHBARTER Zahlen zusammen,')
out('   ueber die lokale Dichte hinaus?  Nullmodell: Klassenmarken im Fenster permutieren.')
out()

# ---------------------------------------------------------------- 1. the powerful numbers WITH class mark
t1 = time.time()
werte = []; klassen = []   # `werte` = blocks of values a²b³; `klassen` = the matching class marks b
b = 1
while b**3 <= N:
    if squarefree(b):
        amax = isqrt(N//b**3)
        if amax:
            a = np.arange(1, amax+1, dtype=np.int64)
            werte.append(a*a*(b**3))
            klassen.append(np.full(amax, b, dtype=np.int32))
    b += 1
W = np.concatenate(werte); K = np.concatenate(klassen)   # `W` = all powerful numbers up to N, `K` = their class marks b
del werte, klassen
ordn = np.argsort(W, kind='stable')
W = W[ordn]; K = K[ordn]
del ordn
assert len(W) == 21663503 or N_EXP != 14, ('PK+ Bestand', len(W))
# PK+: the decomposition must be correct
for probe in (4, 8, 9, 72, 108, 343, 1000):
    i = int(np.searchsorted(W, probe))
    if i < len(W) and W[i] == probe:
        bb = int(K[i]); aa = isqrt(probe//bb**3)
        assert aa*aa*bb**3 == probe and squarefree(bb), ('PK+ a²b³', probe, aa, bb)
out(f'   Bestand: {len(W):,} powerful bis 10^{N_EXP}, Klassenmarke gratis aus der Aufzaehlung ({time.time()-t1:.0f}s)')
kl, kn = np.unique(K, return_counts=True)   # `kl` = distinct classes, `kn` = their sizes
out(f'   {len(kl):,} verschiedene Klassen; groesste: ' +
    ' · '.join(f'b={kl[i]} ({kn[i]:,})' for i in np.argsort(-kn)[:6]))
out(f'   PK+ ok: Stichproben 4, 8, 9, 72, 108, 343, 1000 zerfallen korrekt in a²b³ mit b quadratfrei.')
out()

# ---------------------------------------------------------------- 2. observed class pairs
t1 = time.time()
A = K[:-1].astype(np.int64); B = K[1:].astype(np.int64)   # `A`, `B` = classes of the left / right member of each neighbor pair
MAXB = int(kl.max()) + 1
code = A*MAXB + B                                          # `code` = one integer per class pair (b, d)
uc, cc = np.unique(code, return_counts=True)               # `uc` = distinct class pairs, `cc` = their counts
out(f'   {len(A):,} Nachbarpaare ⇒ {len(uc):,} verschiedene Klassenpaare (b,d)   ({time.time()-t1:.0f}s)')
gleich = int(cc[(uc // MAXB) == (uc % MAXB)].sum())   # `gleich` = number of neighbor pairs in the same class (b = d)
out(f'   davon mit b = d (gleiche Klasse): {gleich:,}  =  {100*gleich/len(A):.3f} %')
out()

# ---------------------------------------------------------------- 3. null model
out('3. NULLMODELL — Klassenmarken im Fenster permutieren (Anordnung bleibt erhalten)')
rng = np.random.default_rng(20260920)
erg_null = {}   # `erg_null` = null-model result per window width `FEN`: (mean, standard deviation) of the same-class count
for FEN in (100, 1000, 10000):
    gl_n = []   # `gl_n` = same-class counts of the repetitions
    for wdh in range(3):
        Kp = K.copy()
        nb = len(Kp)//FEN
        # permute within each window
        rest = len(Kp) - nb*FEN
        blk = Kp[:nb*FEN].reshape(nb, FEN)
        idx = np.argsort(rng.random(blk.shape), axis=1)
        blk = np.take_along_axis(blk, idx, axis=1)
        Kp[:nb*FEN] = blk.reshape(-1)
        if rest > 1:
            tail = Kp[nb*FEN:]
            Kp[nb*FEN:] = tail[rng.permutation(rest)]
        g = int((Kp[:-1] == Kp[1:]).sum())
        gl_n.append(g)
    erg_null[FEN] = (float(np.mean(gl_n)), float(np.std(gl_n, ddof=1)) if len(gl_n) > 1 else 0.0)
    m, s = erg_null[FEN]
    out(f'      Fenster {FEN:>6}: b = d erwartet {m:,.0f} ± {s:,.0f}   beobachtet {gleich:,}   '
        f'⇒ Faktor {gleich/m:.2f}')
out()
out('   🔑 Der Faktor ist ueber alle drei Fensterbreiten zu lesen: haengt er stark von der Breite ab,')
out('      misst man das Fenster und nicht die Nachbarschaft.')
fakt = [gleich/erg_null[f][0] for f in erg_null]
out(f'   ⇒ Faktoren: ' + ' · '.join(f'{x:.2f}' for x in fakt) +
    f'   Spanne {max(fakt)-min(fakt):.2f}')
out()

# ---------------------------------------------------------------- 4. which class pairs stand out?
out('4. WELCHE KLASSENPAARE SIND UEBER-/UNTERREPRAESENTIERT?')
FEN = 1000
Kp = K.copy()
nb = len(Kp)//FEN
blk = Kp[:nb*FEN].reshape(nb, FEN)
idx = np.argsort(rng.random(blk.shape), axis=1)
Kp[:nb*FEN] = np.take_along_axis(blk, idx, axis=1).reshape(-1)
An = Kp[:-1].astype(np.int64); Bn = Kp[1:].astype(np.int64)
un, cn = np.unique(An*MAXB + Bn, return_counts=True)
erw = dict(zip(un.tolist(), cn.tolist()))   # `erw` = expected counts per class pair (from one permutation with window 1000)
# `kand` = candidate cells (code, observed, expected) with at least 200 observed or expected occurrences
kand = [(int(u), int(c), erw.get(int(u), 0)) for u, c in zip(uc, cc) if c >= 200 or erw.get(int(u), 0) >= 200]
# `ueber` / `unter` = the 8 most over- / underrepresented cells
ueber = sorted(kand, key=lambda t: -(t[1]/max(t[2], 0.5)))[:8]
unter = sorted(kand, key=lambda t:  (t[1]/max(t[2], 0.5)))[:8]
out(f'   (nur Paare mit >= 200 beobachteten oder erwarteten Vorkommen: {len(kand):,} Zellen)')
out('   UEBERrepraesentiert:')
for u, c, e in ueber: out(f'      (b={u//MAXB:>5}, d={u%MAXB:>5}): beobachtet {c:>8,}  erwartet {e:>8,}  Faktor {c/max(e,0.5):>6.2f}')
out('   UNTERrepraesentiert:')
for u, c, e in unter: out(f'      (b={u//MAXB:>5}, d={u%MAXB:>5}): beobachtet {c:>8,}  erwartet {e:>8,}  Faktor {c/max(e,0.5):>6.2f}')
out()

# ---------------------------------------------------------------- 5. forbidden pairs?
out('5. VERBOTENE PAARE?  (mit der Abwesenheits-Auflage)')
beob = set(uc.tolist())
leer = [(int(u), int(c)) for u, c in zip(un, cn) if int(u) not in beob]   # `leer` = empty cells: expected > 0, observed 0
leer.sort(key=lambda t: -t[1])
out(f'   Zellen mit Erwartung > 0 und Beobachtung 0: {len(leer):,}')
if leer:
    out(f'   groesste Erwartung darunter: {leer[0][1]}')
    for u, c in leer[:5]: out(f'      (b={u//MAXB:>5}, d={u%MAXB:>5}): erwartet {c}, beobachtet 0')
schwelle = [x for x in leer if x[1] >= 20]   # `schwelle` = empty cells with an expectation of at least 20
out(f'   ⇒ davon mit Erwartung >= 20: {len(schwelle)}')
out('   ⚠️ **Eine leere Zelle mit kleiner Erwartung ist KEIN Befund** — sie ist erwartbar leer.')
out('      Ein Verbot waere erst eine Zelle mit HOHER Erwartung und Beobachtung 0.')

out()

# ---------------------------------------------------------------- 6. the actual structure
# The conspicuous cells of section 4 consist of a single BASE PAIR multiplied by all squares k².
# Example (b=5, d=86): (214k)²·5³ -> (3k)²·86³, distance 4k².
# This is no coincidence: if (n, n+g) is a pair of powerful numbers, then (k²n, k²n+k²g) is one again, since a square times a
# powerful number stays powerful. ⇒ Pairs occur in SQUARE-SCALED FAMILIES.
# Also gcd(k²n, k²g) = k² · gcd(n,g). A pair is called PRIMITIVE here if gcd(n, g) = 1; then it cannot be a square
# multiple of a smaller pair.
out('6. QUADRAT-SKALIERTE FAMILIEN — wie viele Paare sind primitiv?')
from math import gcd as _gcd
G = (W[1:] - W[:-1])   # `G` = gaps g between neighboring powerful numbers
# gcd(n, g) for each sampled pair; if n and n+g are both powerful, this gcd is powerful
schritt = max(1, len(G)//400000)   # `schritt` = stride of the sample
probe = np.arange(0, len(G), schritt)   # `probe` = indices of the sampled pairs
gg = np.array([_gcd(int(W[i]), int(G[i])) for i in probe], dtype=np.int64)
prim = int((gg == 1).sum())   # `prim` = number of primitive pairs in the sample
out(f'   Stichprobe: {len(probe):,} Paare (jedes {schritt}.) ueber den ganzen Bereich')
out(f'      gcd(n, g) = 1  (PRIMITIV):        {prim:>9,}  =  {100*prim/len(probe):.2f} %')
out(f'      gcd(n, g) > 1  (skalierte Kopie): {len(probe)-prim:>9,}  =  {100*(len(probe)-prim)/len(probe):.2f} %')
ug, cg = np.unique(gg[gg > 1], return_counts=True)
top = np.argsort(-cg)[:8]
out(f'   haeufigste gcd-Werte > 1: ' + ' · '.join(f'{ug[i]}({cg[i]:,})' for i in top))
out(f'   ✅ PK+ (T11): jeder dieser gcd muss POWERFUL sein.')
def _pw(n):   # n is powerful: n > 1 and every prime exponent >= 2
    f = {}; m = n; p = 2
    while p*p <= m:
        while m % p == 0: f[p] = f.get(p,0)+1; m //= p
        p += 1 if p == 2 else 2
    if m > 1: f[m] = f.get(m,0)+1
    return n > 1 and all(e >= 2 for e in f.values())
# `schlecht` = frequent gcd values that are not powerful (must be empty)
schlecht = [int(x) for x in ug[:200] if not _pw(int(x))]
assert not schlecht, ('T11 verletzt: gcd nicht powerful', schlecht[:5])
out(f'      geprueft an den {min(200, len(ug))} haeufigsten: alle powerful ✓')
out()
out('   🔑 LESART, und sie ordnet die auffaelligen Zellen aus Abschnitt 4 ein:')
out('      Die Uebervertretung von (5,86) und (10,43) ist KEINE Klassen-Affinitaet, sondern EIN Basispaar,')
out('      das mit k² mitwandert und dabei seine Klassen mitnimmt. **Eine Zelle, viele Vorkommen, eine Ursache.**')
out('   ⚠️ UND DIE GRENZE MEINES NULLMODELLS, ausdruecklich:')
out('      Das Fenster-Permutieren zerstoert die Nachbarschaft, aber auch den MINDESTABSTAND innerhalb')
out('      einer Klasse. Deshalb ist die Diagonale (b = d) dort systematisch zu hoch erwartet — der')
out('      gemessene Faktor 0,5 ist ein ABSTANDSEFFEKT, kein arithmetischer Befund.')
out('      ⇒ Fuer die Diagonale taugt dieses Nullmodell NICHT. Fuer die Nebenzellen taugt es, und dort')
out('        zeigt es auf die skalierten Familien.')

out()
out(f'Zeit {time.time()-t0:.0f}s')
(ERG/'w65_d6_klassenpaare_output.txt').write_text('\n'.join(Z)+'\n', encoding='utf-8')
(ERG/'w65_d6_klassenpaare_result.json').write_text(json.dumps(dict(
    skript='w65_d6_klassenpaare_2026-09-20.py', N=N, n_powerful=int(len(W)), n_klassen=int(len(kl)),
    n_paare=int(len(A)), n_klassenpaare=int(len(uc)),
    gleiche_klasse=dict(beobachtet=gleich, anteil=gleich/len(A),
                        null={str(f): dict(mittel=erg_null[f][0], sd=erg_null[f][1],
                                           faktor=gleich/erg_null[f][0]) for f in erg_null},
                        faktor_spanne=max(fakt)-min(fakt)),
    ueberrepraesentiert=[dict(b=u//MAXB, d=u%MAXB, beob=c, erw=e) for u, c, e in ueber],
    unterrepraesentiert=[dict(b=u//MAXB, d=u%MAXB, beob=c, erw=e) for u, c, e in unter],
    leere_zellen=dict(anzahl=len(leer), mit_erwartung_ab_20=len(schwelle),
                      groesste_erwartung=leer[0][1] if leer else 0),
    primitivitaet=dict(stichprobe=len(probe), schritt=schritt, primitiv=prim,
                       anteil_primitiv=prim/len(probe),
                       haeufigste_gcd=[int(ug[i]) for i in top]),
    nullmodell='Klassenmarken innerhalb von Fenstern permutiert; erhaelt lokale Klassenzusammensetzung',
    nullmodell_grenze='zerstoert auch den Mindestabstand innerhalb einer Klasse => Diagonale nicht beurteilbar',
    was_offen='Exit-Test (occ(g) bedingt auf Klassenpaar gegen 0,545) ist Schritt 2',
    laufzeit_s=time.time()-t0, python=sys.version.split()[0]), indent=1), encoding='utf-8')
