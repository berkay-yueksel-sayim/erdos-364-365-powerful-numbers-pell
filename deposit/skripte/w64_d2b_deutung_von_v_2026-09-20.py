# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w64_d2b_deutung_von_v_2026-09-20.py
# D2b, step 2: what does v(kc) depend on?
#
# Reads:    nothing (the powerful numbers up to 10^14 are enumerated directly; numpy).
# Writes:   ../ergebnisse/w64_d2b_deutung_output.txt and ../ergebnisse/w64_d2b_deutung_result.json
# Usage:    python w64_d2b_deutung_von_v_2026-09-20.py [KMIN]   (KMIN = minimum number of gap values per cofactor, default 4)
# Controls: PK+ : the inventory has 21 663 503 powerful numbers up to 10^14;
#           PK- : with `v` permuted over the levels the regression R² must collapse (< 0.05);
#           in addition, a sign-stability test of the coefficients on disjoint halves (section 4b).
#
# PRECONDITION, and the reason this run is admissible at all:
#   Script w63 showed that `v` is a well-defined quantity: two fits on disjoint halves correlate with r = 0.761, PK- collapses
#   to -0.015, and two independently computed measures agree to three decimals. **Without that, any interpretation here would be
#   a story about a fit artifact.**
#
# THE TRAP TO DEFUSE FIRST:
#   The ceiling `0.864` from w63 was computed UNWEIGHTED over the levels that occur in both halves. Anyone who runs the regression
#   on a DIFFERENT set of levels and compares its R² with this ceiling compares two populations again.
#   => **The ceiling and the regression are computed here on the SAME set of levels, unweighted.**
#   (A comparison needs the same population, the same scale and the same estimator.)
#
# STRUCTURE:
#   1. Full fit `log occ(g) = u(sqfull(g)) + v(kc)` over all gap values.
#   2. ANALYSIS SET: cofactors with at least KMIN associated values (below that, v is too noisy to interpret anything).
#      The number of levels and the coverage are reported.
#   3. CEILING determined anew on exactly this set (split-half + Spearman-Brown), not taken over.
#   4. Regression of `v` on the arithmetic of `kc`: log kc, omega(kc), log of the largest prime factor,
#      log of the smallest prime factor. (The indicator kc = 1 is excluded, see section 2 below.)
#   5. Out of sample (20 % of the levels held back) and PK- (v permuted).
#   6. Ablation: which feature carries the fit, and what does omitting it cost.
#
# NOT CLAIMED HERE: no causality, no theorem, no statement about large gap values.
#    A regression over a finite inventory, read against its own ceiling.
# Own code, numpy.

import sys, json, time, pathlib
import numpy as np
from math import isqrt, gcd
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent   # `HIER` = this directory
ERG = HIER.parent/'ergebnisse'                   # `ERG` = results directory
t0 = time.time(); Z = []   # `Z` = lines of the output file
def out(s=''): Z.append(s); print(s, flush=True)   # `out` = print a line and record it for the output file
KMIN = int(sys.argv[1]) if len(sys.argv) > 1 else 4   # minimum number of gap values per cofactor level

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
out('W64 — D2b SCHRITT 2: WOVON HAENGT v(kc) AB?')
out()
out('   Zulaessig nur wegen W63: v ist eine bestimmte Groesse (r = 0,761; PK- -0,015).')
out('   Decke und Regression werden auf DERSELBEN Stufenmenge gerechnet, ungewichtet.')
out()

N = 10**14
arr = powerful_bis(N)   # `arr` = all powerful numbers up to N
assert len(arr) == 21663503, ('PK+ Bestand', len(arr))
d = np.diff(arr)   # gaps between consecutive powerful numbers
vals, cnts = np.unique(d, return_counts=True)   # distinct gap values and their numbers of occurrences
sel = vals[cnts >= 20]
occ = dict(zip(sel.tolist(), cnts[cnts >= 20].tolist()))   # `occ` = occurrences: gap value -> count (only values with >= 20)
gl = sorted(occ)   # `gl` = the gap values g, sorted
FS = {g: faktor(g) for g in gl}        # factorizations of the gap values
SQ = {g: sqfull(FS[g]) for g in gl}    # `SQ` = square-full part sqfull(g) (the "core"; not the kernel m)
KC = {g: g//SQ[g] for g in gl}         # `KC` = cofactor kc = g / sqfull(g)
S_l = sorted({SQ[g] for g in gl}); C_l = sorted({KC[g] for g in gl})   # distinct cores and distinct cofactors
si = {s: i for i, s in enumerate(S_l)}; ci = {c: i for i, c in enumerate(C_l)}   # level -> index maps
iS = np.array([si[SQ[g]] for g in gl]); iC = np.array([ci[KC[g]] for g in gl])   # core index and cofactor index of every gap value
yo = np.log(np.array([occ[g] for g in gl], float))   # `yo` = log of the occurrence counts (the response of the model)
out(f'   {len(gl):,} Lueckenwerte  ·  {len(S_l):,} Kerne  ·  {len(C_l):,} quadratfreie Kofaktoren')

def fit(maske, y_, runden=80):
    # `fit` = alternating least squares for y = u(core) + v(cofactor) on the gap values selected by the boolean mask `maske`;
    # `runden` = number of sweeps; v is re-centered after each sweep. Returns u and v.
    u = np.zeros(len(S_l)); v = np.zeros(len(C_l))
    cu = np.maximum(np.bincount(iS[maske], minlength=len(S_l)), 1)
    cv = np.maximum(np.bincount(iC[maske], minlength=len(C_l)), 1)
    for _ in range(runden):
        u = np.bincount(iS[maske], weights=y_[maske]-v[iC[maske]], minlength=len(S_l))/cu
        v = np.bincount(iC[maske], weights=y_[maske]-u[iS[maske]], minlength=len(C_l))/cv
        v -= v.mean()
    return u, v

alle = np.ones(len(gl), bool)   # mask selecting all gap values
u_all, v_all = fit(alle, yo)    # the full fit
n_je_kc = np.bincount(iC, minlength=len(C_l))   # `n_je_kc` = number of gap values per cofactor level
M = n_je_kc >= KMIN   # `M` = mask of the cofactor levels in the analysis set
out(f'   ANALYSE-MENGE: Kofaktoren mit mindestens {KMIN} Werten  =>  {int(M.sum()):,} von {len(C_l):,} Stufen '
    f'({100*M.sum()/len(C_l):.0f} %); sie decken {int(n_je_kc[M].sum()):,} der {len(gl):,} Lueckenwerte ab.')
out()

# ---------------------------------------------------------------- 1. ceiling on EXACTLY this set
out('1. DIE DECKE, auf der Analyse-Menge NEU bestimmt (nicht aus W63 uebernommen)')
rs = []; n_da = 0   # `rs` = split-half correlations of the five splits; `n_da` = number of levels compared in the last split
for k in range(1, 6):
    rg = np.random.default_rng(20260920 + k)
    h = rg.random(len(gl)) < 0.5   # random split of the gap values into two halves
    _, v1 = fit(h, yo); _, v2 = fit(~h, yo)
    da = (np.bincount(iC[h], minlength=len(C_l)) > 0) & (np.bincount(iC[~h], minlength=len(C_l)) > 0) & M   # levels in both halves and in the set
    a1 = v1[da] - v1[da].mean(); a2 = v2[da] - v2[da].mean()
    rs.append(float(np.corrcoef(a1, a2)[0, 1])); n_da = int(da.sum())
r_h = float(np.mean(rs))   # mean split-half correlation
DECKE = 2*r_h/(1+r_h)      # `DECKE` = ceiling: Spearman-Brown reliability of the full data
out(f'      Split-Half auf der Analyse-Menge: r = {r_h:.3f} +- {float(np.std(rs, ddof=1)):.3f}   ({n_da:,} Stufen)')
out(f'      => DECKE (Spearman-Brown) = {DECKE:.3f}')
out(f'      Zum Vergleich W63 ueber ALLE Stufen: 0,864. '
    f'{"Hoeher" if DECKE > 0.864 else "Niedriger"}, weil Stufen mit wenigen Werten wegfallen —')
out(f'      genau deshalb wird sie hier neu gerechnet statt uebernommen.')
out()

# ---------------------------------------------------------------- 2. features
out('2. DIE REGRESSION')
# kc = 1 is EXCLUDED from the regression. kc = 1 means "g is itself powerful", a STRUCTURAL special case, not a value of a
#    feature. Used as an indicator in the regression it marks EXACTLY ONE level, so its coefficient would not be an estimate
#    but an exact interpolation of that single row. In the stability test it flipped sign immediately (+0.57 +- 0.60) while
#    all other signs held. => It belongs next to the model, not in it. Its v value is reported separately.
idxM = np.nonzero(M)[0]   # `idxM` = indices of the cofactor levels in the analysis set
i_eins = [i for i, ii in enumerate(idxM) if C_l[ii] == 1]   # position of kc = 1 within the analysis set (if present)
v_kc1 = float(v_all[idxM[i_eins[0]]]) if i_eins else float('nan')   # `v_kc1` = v at kc = 1
behalten = [i for i, ii in enumerate(idxM) if C_l[ii] != 1]   # `behalten` = positions kept for the regression
kcs = [C_l[idxM[i]] for i in behalten]   # the kept cofactors
def merkmale(k):   # `merkmale` = features of kc: [1, log kc, omega(kc), log largest prime factor, log smallest prime factor]
    f = sorted(faktor(k))
    return [1.0, np.log(k), float(len(f)), np.log(f[-1]), np.log(f[0])]
X = np.array([merkmale(k) for k in kcs])   # `X` = feature matrix
y = v_all[idxM][behalten]                  # `y` = the fitted v values (regression target)
# `NAMEN` = feature names (P+ / P- = largest / smallest prime factor)
NAMEN = ['1', 'log kc', 'omega(kc)', 'log P+(kc)', 'log P-(kc)']
out(f'      kc = 1 ausgeschlossen (1 Stufe, struktureller Sonderfall). v(1) = {v_kc1:+.3f}')
out(f'      — das ist der Abstand, um den powerful Lueckenwerte ueber dem Rest liegen; vgl. D19.')
def r2(yv, pv): return float(1 - ((yv-pv)**2).sum()/((yv-yv.mean())**2).sum())   # coefficient of determination
def loese(Xa, ya, Xb):   # `loese` = least-squares fit on (Xa, ya); returns the prediction on Xb and the coefficients
    k, *_ = np.linalg.lstsq(Xa, ya, rcond=None)
    return Xb @ k, k
p, koef = loese(X, y, X)   # in-sample prediction and coefficients (`koef`)
R2 = r2(y, p)
out(f'      {len(y):,} Stufen, {X.shape[1]} Merkmale.   R² (in-sample) = {R2:.3f}')
for nm, kk in zip(NAMEN, koef): out(f'         {nm:<12} {kk:+.4f}')
out()

# ---------------------------------------------------------------- 3. out of sample + PK-
out('3. AUSSERHALB DER ANPASSUNG und PK-')
rg = np.random.default_rng(4242)
m = rg.random(len(y)) < 0.8   # `m` = training mask (80 % of the levels); the rest is held back
pt, _ = loese(X[m], y[m], X[~m])
R2_out = r2(y[~m], pt)   # out-of-sample R²
out(f'      ausserhalb (20 % der Stufen zurueckgehalten, {int((~m).sum()):,}): R² = {R2_out:.3f}')
yp = y.copy(); rg.shuffle(yp)   # `yp` = v permuted over the levels (negative control)
pp, _ = loese(X, yp, X)
R2_perm = r2(yp, pp)
out(f'      PK-, v ueber die Stufen permutiert:                       R² = {R2_perm:.3f}')
assert R2_perm < 0.05, ('PK- muss zusammenbrechen', R2_perm)
out()

# ---------------------------------------------------------------- 4. ablation
out('4. ABLATION')
abl = []   # `abl` = (feature, R² with only that feature and the constant)
for j in range(1, len(NAMEN)):
    Xj = X[:, [0, j]]
    pj, _ = loese(Xj, y, Xj)
    abl.append((NAMEN[j], r2(y, pj)))
out('      einzeln:')
for nm, rr in sorted(abl, key=lambda t: -t[1]): out(f'         nur {nm:<12} R² = {rr:.3f}')
ohne = []   # `ohne` = (feature, loss of R² when that feature is left out)
for j in range(1, len(NAMEN)):
    idx = [c for c in range(X.shape[1]) if c != j]
    pj, _ = loese(X[:, idx], y, X[:, idx])
    ohne.append((NAMEN[j], R2 - r2(y, pj)))
out('      Verlust beim Weglassen:')
for nm, dd in sorted(ohne, key=lambda t: -t[1]): out(f'         ohne {nm:<12} -{dd:.3f}')
out()

# ---------------------------------------------------------------- 4b. are the COEFFICIENTS stable?
# The ablation shows: each feature alone explains a lot, and leaving it out costs almost nothing.
#    That is the hallmark of strong COLLINEARITY: the JOINT R² is reliable, the INDIVIDUAL coefficients are not.
#    Anyone who reads a story from their signs ("at fixed size, smoothness helps") must first show that the signs hold.
#    => The same discipline as in w63, one level up: fit separately on two halves.
out('4b. SIND DIE KOEFFIZIENTEN STABIL?  (Kollinearitaets-Probe)')
korr = np.corrcoef(X[:, 1:].T)   # `korr` = correlation matrix of the features (without the constant)
out('      Korrelation der Merkmale untereinander:')
for i in range(korr.shape[0]):
    for j in range(i+1, korr.shape[1]):
        if abs(korr[i, j]) > 0.5:
            out(f'         {NAMEN[i+1]:<12} vs {NAMEN[j+1]:<12} r = {korr[i, j]:+.3f}')
KO = []   # `KO` = coefficient vectors (first half, second half) for each of 7 random splits
for k in range(1, 8):
    rgk = np.random.default_rng(555 + k)
    hk = rgk.random(len(y)) < 0.5
    _, k1 = loese(X[hk], y[hk], X[hk])
    _, k2 = loese(X[~hk], y[~hk], X[~hk])
    KO.append((k1, k2))
out('      Koeffizienten auf je zwei disjunkten Haelften, 7 Aufteilungen:')
stabil = True   # `stabil` = True if every coefficient keeps its sign over all 14 halves
for j in range(1, len(NAMEN)):
    w = np.array([k1[j] for k1, k2 in KO] + [k2[j] for k1, k2 in KO])
    gleich = (w > 0).all() or (w < 0).all()
    stabil &= gleich
    out(f'         {NAMEN[j]:<12} {w.mean():+.4f} +- {w.std(ddof=1):.4f}   '
        f'Vorzeichen {"STABIL" if gleich else "KIPPT"}')
out(f'      => {"Alle Vorzeichen halten ueber alle 14 Haelften." if stabil else "MINDESTENS EIN Vorzeichen kippt."}')
if stabil:
    out('      ⇒ Die RICHTUNGSAUSSAGEN sind lesbar. Die BETRAEGE nicht — dafuer ist die Kollinearitaet zu stark')
    out('        (omega gegen log P-: r = -0,85; log P+ gegen log P-: r = +0,84).')
else:
    out('      ⇒ Keine Richtungsaussage zulaessig, solange ein Vorzeichen kippt.')
out()

# ---------------------------------------------------------------- 5. against the ceiling
out('5. DAS ERGEBNIS, GEGEN SEINE EIGENE DECKE GELESEN')
out(f'      R² (ausserhalb) = {R2_out:.3f}')
out(f'      DECKE           = {DECKE:.3f}')
out(f'      => {100*R2_out/DECKE:.0f} % des ERKLAERBAREN, nicht {100*R2_out:.0f} % des Ganzen.')
out(f'      => erklaerbarer, aber unerklaerter Rest: {DECKE - R2_out:.3f}')
out('      Kein Satz, keine Kausalitaet, endlicher Bestand.')

out()
out(f'Zeit {time.time()-t0:.0f}s')
(ERG/'w64_d2b_deutung_output.txt').write_text('\n'.join(Z)+'\n', encoding='utf-8')
(ERG/'w64_d2b_deutung_result.json').write_text(json.dumps(dict(
    skript='w64_d2b_deutung_von_v_2026-09-20.py', N=N, KMIN=KMIN,
    n_werte=len(gl), n_kofaktoren=len(C_l), n_analyse_stufen=int(M.sum()),
    abdeckung_werte=int(n_je_kc[M].sum()),
    decke=dict(r_split_half=r_h, spearman_brown=DECKE, w63_ueber_alle_stufen=0.864, n_stufen=n_da),
    kc_gleich_1=dict(v=v_kc1, behandlung='ausgeschlossen, struktureller Sonderfall (D19)'),
    regression=dict(merkmale=NAMEN, koeffizienten=[float(k) for k in koef],
                    r2_in_sample=R2, r2_ausserhalb=R2_out, r2_permutiert=R2_perm),
    ablation=dict(einzeln={nm: rr for nm, rr in abl}, verlust={nm: dd for nm, dd in ohne}),
    koeffizienten_stabilitaet=dict(
        vorzeichen_stabil=bool(stabil),
        je_merkmal={NAMEN[j]: dict(
            mittel=float(np.mean([k1[j] for k1, k2 in KO] + [k2[j] for k1, k2 in KO])),
            sd=float(np.std([k1[j] for k1, k2 in KO] + [k2[j] for k1, k2 in KO], ddof=1)))
            for j in range(1, len(NAMEN))},
        hinweis='starke Kollinearitaet: gemeinsames R2 belastbar, Einzelbetraege nicht'),
    urteil=dict(anteil_des_erklaerbaren=R2_out/DECKE, rest=DECKE - R2_out),
    abgrenzung='Regression ueber endlichen Bestand, gegen die eigene Decke gelesen; kein Satz, keine Kausalitaet',
    laufzeit_s=time.time()-t0, python=sys.version.split()[0]), indent=1), encoding='utf-8')
