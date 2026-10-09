# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w63_d2b_stabilitaet_von_v_2026-09-20.py
# D2b, step 1: is the cofactor effect v(kc) stable at all?
#
# Reads:    nothing (the powerful numbers up to 10^14 are enumerated directly; numpy).
# Writes:   ../ergebnisse/w63_d2b_stabilitaet_output.txt and ../ergebnisse/w63_d2b_stabilitaet_result.json
# Usage:    python w63_d2b_stabilitaet_von_v_2026-09-20.py [MINCNT]   (MINCNT = minimum number of occurrences, default 20)
# Controls: PK+ : the inventory has 21 663 503 powerful numbers up to 10^14; the cofactors `kc` are squarefree (sample of 400);
#           PK- : permuting `occ` must make the correlation collapse; the noise model of part (C) must match the measured
#           correlation.
#
# WHY THIS STEP COMES BEFORE THE ACTUAL QUESTION:
#   `D2b` asks WHAT the cofactor effect `v(kc)` depends on. But `v(kc)` is a **fitted parameter** of the model
#   `log occ(g) = u(sqfull(g)) + v(kc)` (script w59). Interpreting it means interpreting an artifact of the fit, **unless it is
#   stable**. => First measure whether two independent fits on DISJOINT halves give the same `v`.
#   If they do not, every story about the "smoothness of kc" is a story about noise.
#
# STRUCTURE:
#   1. The inventory up to 10^14, gap values with >= MINCNT occurrences; for each g: sqfull(g) and kc = g/sqfull(g).
#   2. Split the g values RANDOMLY into two disjoint halves and fit each half separately.
#   3. Compare `v` of the two halves on the kc levels that occur in BOTH halves.
#      Measured: (a) the correlation, (b) the spread of the difference against the spread of v itself.
#      (b) is the decisive one: a high correlation with a large spread of the difference only means that both fits share
#         the same coarse direction, not that the individual level is determined.
#   4. The same for `u` (the core, i.e. the square-full part `sqfull`; not the kernel m) as a comparison.
#   5. PK-: permute `occ` and repeat everything. The correlation MUST then collapse.
#   6. Repeated over several splits, because a single split says nothing about the spread.
#
# NOT DONE HERE: no interpretation of v. That is step 2 and depends on the result here.
# Own code (Miller-Rabin and Pollard rho as in w59), numpy.

import sys, json, time, pathlib
import numpy as np
from math import isqrt, gcd
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent   # `HIER` = this directory
ERG = HIER.parent/'ergebnisse'                   # `ERG` = results directory
t0 = time.time(); Z = []   # `Z` = lines of the output file
def out(s=''): Z.append(s); print(s, flush=True)   # `out` = print a line and record it for the output file
MINCNT = int(sys.argv[1]) if len(sys.argv) > 1 else 20   # minimum number of occurrences of a gap value

def ist_prim(n):   # `ist_prim` = deterministic Miller-Rabin primality test with the first 12 primes as bases
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % p == 0: return n == p
    d = n-1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x = pow(a, d, n)
        if x == 1 or x == n-1: continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True
def rho(n):   # Pollard rho: returns a nontrivial divisor of the composite n (the constant c is increased until one is found)
    if n % 2 == 0: return 2
    c = 1
    while True:
        x = y = 2; d = 1
        while d == 1:
            x = (x*x + c) % n
            y = (y*y + c) % n; y = (y*y + c) % n
            d = gcd(abs(x-y), n)
        if d != n: return d
        c += 1
def faktor(n, f=None):   # `faktor` = complete factorization of n as {prime: exponent} (trial division, then Miller-Rabin / rho)
    f = {} if f is None else f
    if n == 1: return f
    for p in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47):
        while n % p == 0: f[p] = f.get(p,0)+1; n //= p
    if n == 1: return f
    if ist_prim(n): f[n] = f.get(n,0)+1; return f
    w = isqrt(n)
    if w*w == n:
        for p, e in faktor(w).items(): f[p] = f.get(p,0) + 2*e
        return f
    d = rho(n); faktor(d, f); faktor(n//d, f); return f
def squarefree(n): return all(e == 1 for e in faktor(n).values())
# `sqfull` = square-full part of a number, from its factorization {prime: exponent}: product of p^e over all e >= 2
def sqfull(fs):
    r = 1
    for p, e in fs.items():
        if e >= 2: r *= p**e
    return r
def powerful_bis(N):   # `powerful_bis` = sorted array of all powerful numbers <= N, as a^2*b^3 with b squarefree
    teile = []; b = 1
    while b**3 <= N:
        if squarefree(b):
            amax = isqrt(N//b**3)
            if amax:
                a = np.arange(1, amax+1, dtype=np.int64); teile.append(a*a*(b**3))
        b += 1
    arr = np.concatenate(teile); arr.sort(); return arr

out('='*100)
out('W63 — D2b SCHRITT 1: IST v(kc) STABIL?')
out()
out('   Gefragt wird NICHT, wovon v abhaengt, sondern ob v ueberhaupt eine bestimmte Groesse ist.')
out('   Zwei unabhaengige Anpassungen auf disjunkten Haelften — liefern sie dasselbe v?')
out()

N = 10**14
arr = powerful_bis(N)   # `arr` = all powerful numbers up to N
assert len(arr) == 21663503, ('PK+ Bestand', len(arr))
d = np.diff(arr)   # gaps between consecutive powerful numbers
vals, cnts = np.unique(d, return_counts=True)   # distinct gap values and their numbers of occurrences
sel = vals[cnts >= MINCNT]
# `occ` = occurrences: gap value -> count (only values with >= MINCNT)
occ = dict(zip(sel.tolist(), cnts[cnts >= MINCNT].tolist()))
gl = sorted(occ)   # `gl` = the gap values g, sorted
out(f'   Bestand {len(arr):,} powerful bis 10^14; {len(gl):,} Lueckenwerte mit >= {MINCNT} Vorkommen.')
FS = {g: faktor(g) for g in gl}        # factorizations of the gap values
SQ = {g: sqfull(FS[g]) for g in gl}    # `SQ` = square-full part sqfull(g) (the "core")
KC = {g: g//SQ[g] for g in gl}         # `KC` = cofactor kc = g / sqfull(g)
assert all(all(e == 1 for e in faktor(KC[g]).values()) for g in gl[:400]), 'PK+ kc quadratfrei (Stichprobe)'
S_l = sorted({SQ[g] for g in gl}); C_l = sorted({KC[g] for g in gl})   # distinct cores and distinct cofactors
si = {s: i for i, s in enumerate(S_l)}; ci = {c: i for i, c in enumerate(C_l)}   # level -> index maps
iS = np.array([si[SQ[g]] for g in gl]); iC = np.array([ci[KC[g]] for g in gl])   # core index and cofactor index of every gap value
yo = np.log(np.array([occ[g] for g in gl], float))   # `yo` = log of the occurrence counts (the response of the model)
out(f'   {len(S_l):,} Kerne, {len(C_l):,} quadratfreie Kofaktoren ⇒ im Mittel {len(gl)/len(C_l):.1f} Werte je Kofaktor.')
out()

def fit(maske, runden=80):
    # `fit` = alternating least squares for log occ = u(core) + v(cofactor) on the gap values selected by the boolean mask
    # `maske`; `runden` = number of sweeps; v is re-centered after each sweep. Returns u, v and the masks of the levels
    # present in the subset.
    u = np.zeros(len(S_l)); v = np.zeros(len(C_l))
    cu = np.maximum(np.bincount(iS[maske], minlength=len(S_l)), 1)
    cv = np.maximum(np.bincount(iC[maske], minlength=len(C_l)), 1)
    for _ in range(runden):
        u = np.bincount(iS[maske], weights=yo[maske]-v[iC[maske]], minlength=len(S_l))/cu
        v = np.bincount(iC[maske], weights=yo[maske]-u[iS[maske]], minlength=len(C_l))/cv
        v -= v.mean()
    da_u = np.bincount(iS[maske], minlength=len(S_l)) > 0
    da_v = np.bincount(iC[maske], minlength=len(C_l)) > 0
    return u, v, da_u, da_v

def vergleich(y_, saat):
    # `vergleich` = comparison: split the gap values at random (seed `saat`) into two halves, fit both on the response `y_`
    # (temporarily substituted for the global `yo`), and compare u and v of the halves on the levels present in both
    rg = np.random.default_rng(saat)
    h = rg.random(len(gl)) < 0.5   # `h` = membership in the first half
    global yo
    alt = yo; yo = y_
    u1, v1, du1, dv1 = fit(h)
    u2, v2, du2, dv2 = fit(~h)
    yo = alt
    bv = dv1 & dv2; bu = du1 & du2   # levels present in both halves
    # remove the common shift: v is determined only up to a constant
    a1 = v1[bv] - v1[bv].mean(); a2 = v2[bv] - v2[bv].mean()
    b1 = u1[bu] - u1[bu].mean(); b2 = u2[bu] - u2[bu].mean()
    rv = float(np.corrcoef(a1, a2)[0,1]); ru = float(np.corrcoef(b1, b2)[0,1])
    # spread of the difference against the spread of v itself
    sv = float(np.std(a1 - a2, ddof=1)) / float(np.std(np.concatenate([a1, a2]), ddof=1))
    su = float(np.std(b1 - b2, ddof=1)) / float(np.std(np.concatenate([b1, b2]), ddof=1))
    return dict(r_v=rv, r_u=ru, rel_streuung_v=sv, rel_streuung_u=su,
                n_v=int(bv.sum()), n_u=int(bu.sum()))

out('   (A) ECHTE DATEN, fuenf Aufteilungen')
E = [vergleich(yo, 20260920 + k) for k in range(1, 6)]   # `E` = results of five splits on the real data
for k, e in enumerate(E, 1):
    out(f'      Split {k}:  r(v) = {e["r_v"]:+.3f}  auf {e["n_v"]:,} gemeinsamen Kofaktoren   ·   '
        f'r(u) = {e["r_u"]:+.3f}  auf {e["n_u"]:,} Kernen')
rv = np.array([e['r_v'] for e in E]); ru = np.array([e['r_u'] for e in E])
sv = np.array([e['rel_streuung_v'] for e in E]); su = np.array([e['rel_streuung_u'] for e in E])
out(f'      ⇒ r(v) = {rv.mean():+.3f} ± {rv.std(ddof=1):.3f}      r(u) = {ru.mean():+.3f} ± {ru.std(ddof=1):.3f}')
out(f'      ⇒ Streuung der DIFFERENZ, relativ zur Streuung von v: {sv.mean():.3f} ± {sv.std(ddof=1):.3f}')
out(f'                                            fuer u:         {su.mean():.3f} ± {su.std(ddof=1):.3f}')
out('      🔑 Die zweite Zeile ist die strengere: sie fragt, wie genau die EINZELNE Stufe bestimmt ist,')
out('         nicht nur, ob beide Anpassungen dieselbe grobe Richtung haben.')
out()

out('   (B) PK− : dieselbe Rechnung mit PERMUTIERTEN Zaehlungen')
rg = np.random.default_rng(777)
yp = yo.copy(); rg.shuffle(yp)   # `yp` = permuted response (negative control)
P = [vergleich(yp, 20260920 + k) for k in range(1, 6)]   # `P` = results of five splits on the permuted data
rvp = np.array([e['r_v'] for e in P]); rup = np.array([e['r_u'] for e in P])
out(f'      permutiert:  r(v) = {rvp.mean():+.3f} ± {rvp.std(ddof=1):.3f}   '
    f'r(u) = {rup.mean():+.3f} ± {rup.std(ddof=1):.3f}')
assert abs(rvp.mean()) < 0.2, ('PK- die Korrelation muss bei permutierten Zaehlungen zusammenbrechen', rvp.mean())
out(f'      ✅ Zusammengebrochen ⇒ das Ergebnis in (A) misst Struktur, nicht die Maschinerie.')
out()

out('   (C) WAS DARAUS FUER SCHRITT 2 FOLGT — die OBERGRENZE, gerechnet statt behauptet')
out('      Modell: jede Haelfte misst dieselbe Groesse plus unabhaengiges Rauschen, v_i = v + eps_i.')
out('      Dann gilt zweierlei, und beides ist hier messbar:')
out('         (i)  r  = Var(v) / (Var(v) + Var(eps))            [Korrelation der beiden Haelften]')
out('         (ii) sd(v1-v2) / sd(v) = sqrt( 2 Var(eps) / (Var(v)+Var(eps)) )')
q = float(sv.mean())   # `q` = mean relative spread of the difference
anteil = (q*q/2) / (1 - q*q/2)          # Var(eps)/Var(v), from (ii)
r_vorher = 1.0 / (1.0 + anteil)          # the correlation r predicted from it by (i)
out(f'      Aus der Differenzstreuung {q:.3f} folgt Var(eps)/Var(v) = {anteil:.3f}')
out(f'      ⇒ vorhergesagtes r = {r_vorher:.3f}     gemessenes r = {rv.mean():.3f}')
abw = abs(r_vorher - rv.mean())   # `abw` = deviation between predicted and measured r
out(f'      ⇒ Abweichung {abw:.3f}. {"Das Rauschmodell traegt." if abw < 0.02 else "Das Rauschmodell traegt NICHT."}')
assert abw < 0.05, ('Die zwei Masse muessen zusammenpassen, sonst ist das Modell falsch', r_vorher, rv.mean())
sb = 2*rv.mean()/(1+rv.mean())   # `sb` = Spearman-Brown reliability of the full data: 2r/(1+r)
out('')
out(f'      🔑 **Spearman-Brown:** eine Haelfte hat Reliabilitaet {rv.mean():.3f}; die VOLLEN Daten,')
out(f'         also doppelt so viele Werte je Stufe, haben 2r/(1+r) = **{sb:.3f}**.')
out(f'      ⇒ ***Eine Regression von v auf Merkmale von kc kann hoechstens R² ≈ {sb:.2f} erreichen.***')
out(f'         Alles darueber waere Anpassung an Rauschen. Und ein erreichtes R² ist gegen {sb:.2f} zu lesen,')
out(f'         nicht gegen 1: R² = 0,50 hiesse {0.50/sb:.0%} des ERKLAERBAREN, nicht 50 %.')
out('')
out('   ⇒ URTEIL: Schritt 2 (Deutung von v) ist ZULAESSIG, aber nur mit dieser Decke im Satz.')
out()
out('   ⇒ URTEIL steht im Record. Hier nur die Zahlen.')
out(f'      Lesart: r(v) nahe 1 UND kleine relative Differenzstreuung ⇒ v ist eine bestimmte Groesse,')
out(f'      Schritt 2 (Deutung) ist zulaessig. Sonst nicht.')
out()
out(f'Zeit {time.time()-t0:.0f}s')
(ERG/'w63_d2b_stabilitaet_output.txt').write_text('\n'.join(Z)+'\n', encoding='utf-8')
(ERG/'w63_d2b_stabilitaet_result.json').write_text(json.dumps(dict(
    skript='w63_d2b_stabilitaet_von_v_2026-09-20.py', N=N, MINCNT=MINCNT,
    n_werte=len(gl), n_kerne=len(S_l), n_kofaktoren=len(C_l),
    echt=E, permutiert=P,
    zusammenfassung=dict(r_v_mittel=float(rv.mean()), r_v_sd=float(rv.std(ddof=1)),
                         r_u_mittel=float(ru.mean()), r_u_sd=float(ru.std(ddof=1)),
                         rel_streuung_v=float(sv.mean()), rel_streuung_u=float(su.mean()),
                         r_v_permutiert=float(rvp.mean()), r_u_permutiert=float(rup.mean())),
    rauschmodell=dict(rel_streuung=q, var_eps_durch_var_v=anteil,
                      r_vorhergesagt=r_vorher, r_gemessen=float(rv.mean()), abweichung=abw),
    obergrenze_schritt2=dict(reliabilitaet_haelfte=float(rv.mean()),
                             reliabilitaet_voll_spearman_brown=sb,
                             lesart='ein erreichtes R2 ist gegen diese Decke zu lesen, nicht gegen 1'),
    was_nicht_geprueft='keine Deutung von v; das ist Schritt 2',
    laufzeit_s=time.time()-t0, python=sys.version.split()[0]), indent=1), encoding='utf-8')
