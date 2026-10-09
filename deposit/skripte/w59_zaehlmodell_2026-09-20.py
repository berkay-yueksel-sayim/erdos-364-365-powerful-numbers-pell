# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w59_zaehlmodell_2026-09-20.py
# Counting model for the gap spectrum: does the gcd theorem of w49 also determine HOW OFTEN a gap g occurs?
# Question: the theorem of w49 says WHICH gcd values occur (the powerful divisors of g). Does it follow how often a gap g occurs?
# Reads: nothing (the powerful numbers up to 10^14 are generated). Writes: ergebnisse/w59_zaehlmodell_output.txt and
#   ergebnisse/w59_zaehlmodell_result.json. Requires numpy.
# Usage: python w59_zaehlmodell_2026-09-20.py [MINCNT]     (default 20; `MINCNT` = minimum number of occurrences of a gap value)
#
# THE MODEL, WRITTEN DOWN BEFORE THE FIT:
#   For a gap value g and a powerful divisor h | g let
#       occ_h(g) = number of occurrences of g with gcd(n, g) = h  (exact),   k = g/h.
#   From n = h*a, g = h*k it follows that n and n+g are both powerful exactly when a and a+k are both
#   "powerful outside the primes of h". The KERNEL h determines which numbers are admitted; the COFACTOR k determines the
#   distance within this set.
#   => HYPOTHESIS (log-separable):   log occ_h(g)  ~  alpha(h) + beta(k)
#   => and from it, A PRIORI from g alone:
#         occ(g) = sum_{h powerful, h | g} exp(alpha(h) + beta(g/h))
#      Warning: the sum runs over ALL powerful divisors of g, not over the OBSERVED ones. Which cell is non-empty is itself a
#      result; summing over those would be selection on the target quantity.
#
# WHAT REFUTES THE MODEL: a systematic INTERACTION, i.e. residuals that depend jointly on h and k. Then the separation into
#   kernel and cofactor is wrong.
#
# TWO SEPARATE MEASUREMENTS, because they are two different questions:
#   (A) split by CELLS -- does the separable structure hold? (prediction of a cell for a known g)
#   (B) split by g -- is occ(g) predictable for a NEW, never seen gap value?
# BASELINE: an earlier exploratory run (raw predictors, R^2 = 0.333) used 4000 values g <= 10^6 and is NOT comparable with (B),
#   which uses a different set. Here it is recomputed on the SAME set and the SAME split, in the better of its two versions
#   (count resp. log count).
# Further parts: (C) negative control with permuted counts, (C2) why the fallback version fails, (D) the strongest version
#   (alpha and beta as smooth functions of the arithmetic of h and k), (E) spread over five splits, (F) the pre-registered
#   predictions P1/P2, (G) decision of (F) against raw features, (H) where the predictive power sits.
#
# PK+ : the number of powerful numbers up to 10^14 must be 21 663 503 (assert).
# PK+ : the decomposition of g = 900 must give 900:57, 225:14, 100:19, 36:26, 25:4, 9:3, 4:4 (w49, part D).
# PK+ : factorization against 44100 = 2^2*3^2*5^2*7^2 and against 63900 (sqfull = 900).
# PK- : with PERMUTED counts, (B) must collapse.
# Own code (Miller-Rabin and Brent rho written ourselves), numpy.
import sys, json, time, pathlib
import numpy as np
from math import isqrt, gcd
sys.stdout.reconfigure(encoding='utf-8')
sys.setrecursionlimit(10000)

HIER = pathlib.Path(__file__).resolve().parent
ERG  = HIER.parent / 'ergebnisse'
t0 = time.time(); Z = []   # `Z` = output lines
# `out` = print a line and record it in `Z`
def out(s=''): Z.append(s); print(s, flush=True)
MINCNT = int(sys.argv[1]) if len(sys.argv) > 1 else 20

# ------------------------------------------------------------------ number-theory tools
# `ist_prim` = primality test (Miller-Rabin)
def ist_prim(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % p == 0: return n == p
    d = n-1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):   # deterministic for n < 3.3e24
        x = pow(a, d, n)
        if x == 1 or x == n-1: continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def rho(n):
    """Brent. Returns a proper divisor of a composite, odd n > 1."""
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

# `faktor(n)` = factorization of n as {prime: exponent}: small primes, Miller-Rabin, perfect squares, Brent rho; `f` = accumulator
def faktor(n, f=None):
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

def squarefree(n): return all(e==1 for e in faktor(n).values())
def sqfull(fs):    # `sqfull` = squarefull part: product of the p^e with e >= 2 (`fs` = factorization)
    r = 1
    for p, e in fs.items():
        if e >= 2: r *= p**e
    return r
# `powerful_teiler(fs)` = all powerful divisors, computed from the factorization `fs`
def powerful_teiler(fs):
    """all powerful divisors: v_p(h) in {0} u {2..v_p(g)}."""
    tl = [1]
    for p, e in fs.items():
        if e < 2: continue
        neu = []
        for t in tl:
            neu.append(t)
            for j in range(2, e+1): neu.append(t * p**j)
        tl = neu
    return tl

# `powerful_bis(N)` = sorted array of all powerful numbers up to N, built as a^2*b^3 with b squarefree
# (`teile` = the arrays per b, `amax` = largest a)
def powerful_bis(N):
    teile=[]; b=1
    while b**3<=N:
        if squarefree(b):
            amax=isqrt(N//b**3)
            if amax:
                a=np.arange(1,amax+1,dtype=np.int64); teile.append(a*a*(b**3))
        b+=1
    arr=np.concatenate(teile); arr.sort(); return arr

assert faktor(44100) == {2:2,3:2,5:2,7:2}, 'PK+ Faktorisierung 44100'
assert sqfull(faktor(63900)) == 900, 'PK+ sqfull(63900) = 900'
assert sorted(powerful_teiler(faktor(900))) == [1,4,9,25,36,100,225,900], 'PK+ powerful Teiler von 900'

out('='*104)
out('W59 — D2: DAS ZAEHLMODELL')
out()
out('   Modell (vor dem Fit festgelegt):   log occ_h(g)  ≈  alpha(h) + beta(k),   k = g/h')
out('   A priori:  occ(g) = Summe ueber ALLE powerful Teiler h von g, nicht ueber die beobachteten.')
out('   (A) Teilung nach Zellen: haelt die Struktur?   (B) Teilung nach g: ist ein NEUES g vorhersagbar?')
out('   PK+ ok: Faktorisierung, sqfull, powerful Teiler gegen 44100 / 63900 / 900 geprueft.')
out()

# `arr` = all powerful numbers up to N, `d` = gaps between consecutive ones, `vals`/`cnts` = gap values and their frequencies,
# `sel` = gap values occurring at least MINCNT times
N = 10**14
t1=time.time(); arr = powerful_bis(N)
assert len(arr) == 21663503, ('PK+ Bestand', len(arr))
d = np.diff(arr)
vals, cnts = np.unique(d, return_counts=True)
sel = vals[cnts >= MINCNT]
out(f'   Bestand: {len(arr):,} powerful bis 10^14, {len(vals):,} Lueckenwerte  ({time.time()-t1:.0f}s)')
out(f'   Untersucht: {len(sel):,} Werte mit mindestens {MINCNT} Vorkommen; groesster {int(sel.max()):,}')

t1=time.time()
# `maske` = mask of the gaps with a selected value; `gs` = their values g, `ns` = their lower ends n
maske = np.isin(d, sel)
gs = d[maske]; ns = arr[:-1][maske]
# `zelle` = cells: (g, h) -> number of occurrences of the gap g whose lower end n has gcd(n, g) = h
zelle = {}
for g, n in zip(gs.tolist(), ns.tolist()):
    h = gcd(n, g)
    zelle[(g, h)] = zelle.get((g, h), 0) + 1
out(f'   {len(gs):,} Vorkommen ⇒ {len(zelle):,} Zellen (g, h)   ({time.time()-t1:.0f}s)')
pk = {h: c for (g, h), c in zelle.items() if g == 900}
assert pk == {900:57, 225:14, 100:19, 36:26, 25:4, 9:3, 4:4}, ('PK+ Zerlegung g = 900', pk)
out(f'   PK+ ok: g = 900 zerfaellt in {dict(sorted(pk.items(), reverse=True))} — wie in W49 Teil D.')
out()

# ------------------------------------------------------------------ indices
# `schl` = cell keys (g, h); `H` = distinct kernels h; `K` = distinct cofactors k = g/h; `hi`/`ki` = index maps;
# `ih`/`ik` = kernel / cofactor index per cell; `cc` = counts; `y` = log counts; `gg` = g per cell;
# `nH`, `nK`, `nZ` = numbers of kernels, cofactors, cells
schl = list(zelle.keys())
H = sorted({h for (g,h) in schl}); K = sorted({g//h for (g,h) in schl})
hi = {h:i for i,h in enumerate(H)}; ki = {k:i for i,k in enumerate(K)}
ih = np.fromiter((hi[h]    for (g,h) in schl), np.int64, len(schl))
ik = np.fromiter((ki[g//h] for (g,h) in schl), np.int64, len(schl))
cc = np.fromiter(zelle.values(), np.float64, len(schl)); y = np.log(cc)
gg = np.fromiter((g for (g,h) in schl), np.int64, len(schl))
nH, nK, nZ = len(H), len(K), len(schl)
out(f'   Parameter: {nH:,} Kerne + {nK:,} Kofaktoren - 1 = {nH+nK-1:,} bei {nZ:,} Zellen'
    f'   ⇒ {(nH+nK-1)/nZ:.2f} je Datenpunkt')
out(f'   (im Mittel {nZ/nK:.1f} Zellen pro Kofaktor, {nZ/nH:.0f} pro Kern — das Modell ist NICHT gesaettigt,')
out(f'    deshalb ist auch das In-sample-R² hier aussagekraeftig. Trotzdem entscheidet (B).)')
out()

# `fit` = backfitting (alternating least squares) for y ~ a[h] + b[k] with b centered; `runden` = rounds;
# `ch`/`ck` = cell counts per stage
def fit(ih_, ik_, y_, runden=80):
    a = np.zeros(nH); b = np.zeros(nK)
    ch = np.maximum(np.bincount(ih_, minlength=nH), 1)
    ck = np.maximum(np.bincount(ik_, minlength=nK), 1)
    for _ in range(runden):
        a = np.bincount(ih_, weights=y_-b[ik_], minlength=nH)/ch
        b = np.bincount(ik_, weights=y_-a[ih_], minlength=nK)/ck
        b -= b.mean()
    return a, b
def r2(yv, pv):    # coefficient of determination of the predictions `pv` for the values `yv`
    yv = np.asarray(yv, float); pv = np.asarray(pv, float)
    return float(1 - ((yv-pv)**2).sum()/((yv-yv.mean())**2).sum())

# ------------------------------------------------------------------ (A) cell split
# `mtr` = training mask over the cells (80 %); `a_all`, `b_all` = fit on all cells; `r2_in` = in-sample R^2;
# `kh`/`kk` = stages seen in training; `tc` = test cells whose stages are known;
# `kor_x` = correlation of the residuals with alpha*beta (interaction)
rng = np.random.default_rng(20260920)
mtr = rng.random(nZ) < 0.80
a_all, b_all = fit(ih, ik, y)
r2_in = r2(y, a_all[ih]+b_all[ik])
a_c, b_c = fit(ih[mtr], ik[mtr], y[mtr])
kh = np.bincount(ih[mtr], minlength=nH) > 0; kk = np.bincount(ik[mtr], minlength=nK) > 0
tc = (~mtr) & kh[ih] & kk[ik]
r2_zell_out = r2(y[tc], a_c[ih[tc]]+b_c[ik[tc]])
res = y - (a_all[ih]+b_all[ik])
kor_x = float(np.corrcoef(a_all[ih]*b_all[ik], res)[0,1])
out('   (A) ZELLEN-SPLIT — haelt log occ_h(g) = alpha(h) + beta(k)?')
out(f'      in-sample R² = {r2_in:.3f};  ausserhalb R² = {r2_zell_out:.3f}  ({int(tc.sum()):,} Testzellen)')
out(f'      Wechselwirkung (Residuen gegen alpha·beta): r = {kor_x:+.3f}')
out()

# ------------------------------------------------------------------ factorization of all g
# `gl` = studied gap values; `FS` = factorizations; `PT` = powerful divisors; `SQ` = squarefull parts; `occ` = occurrences per g
t1=time.time()
gl = sorted({int(v) for v in sel})
FS = {g: faktor(g) for g in gl}
PT = {g: sorted(powerful_teiler(FS[g])) for g in gl}
SQ = {g: sqfull(FS[g]) for g in gl}
occ = {g: 0 for g in gl}
for (g,h), c in zelle.items(): occ[g] += c
out(f'   Faktorisierung von {len(gl):,} Lueckenwerten ({time.time()-t1:.0f}s); '
    f'im Mittel {np.mean([len(PT[g]) for g in gl]):.1f} powerful Teiler je g')
# control: every observed gcd MUST be a powerful divisor (the theorem)
# `schlecht` = bad: observed cells whose gcd is not a powerful divisor
schlecht = [(g,h) for (g,h) in schl if h not in set(PT[g])]
assert not schlecht, ('Satz-Kontrolle: gcd ist kein powerful Teiler', schlecht[:5])
out(f'   ✅ Satz-Kontrolle: alle {nZ:,} beobachteten gcd-Werte sind powerful Teiler ihres g.')
out()

# ------------------------------------------------------------------ (B) g split -- the actual question
# `gi` = index of each g; `gtrain` = training mask over the g values (80 %); `zell_tr` = cells whose g is in the training part;
# `a_g`/`b_g` = fit on the training cells; `bh`/`bk` = kernels / cofactors seen in training;
# `vorhersage` = prediction of occ(g)
out('   (B) g-SPLIT — ist occ(g) fuer einen NEUEN Lueckenwert vorhersagbar?  ⟵ das ist D2')
gi = {g:i for i,g in enumerate(gl)}
gtrain = rng.random(len(gl)) < 0.80
zell_tr = np.fromiter((gtrain[gi[g]] for g in gg.tolist()), bool, nZ)
a_g, b_g = fit(ih[zell_tr], ik[zell_tr], y[zell_tr])
bh = np.bincount(ih[zell_tr], minlength=nH) > 0; bk = np.bincount(ik[zell_tr], minlength=nK) > 0

def vorhersage(g):
    """a priori: sum over ALL powerful divisors. None if a stage is unknown."""
    s = 0.0
    for h in PT[g]:
        k = g//h
        if h not in hi or k not in ki: return None
        iH, iK = hi[h], ki[k]
        if not (bh[iH] and bk[iK]): return None
        s += float(np.exp(a_g[iH] + b_g[iK]))
    return s

# `test_g` = held-out g; `bewert` = the evaluable ones (all stages known in training); `e_sep`/`v_sep` = observed / predicted occ
test_g = [g for g in gl if not gtrain[gi[g]]]
paare = [(g, vorhersage(g)) for g in test_g]
bewert = [(g,v) for g,v in paare if v is not None]
e_sep = np.array([occ[g] for g,v in bewert], float)
v_sep = np.array([v      for g,v in bewert], float)
r2_sep = r2(e_sep, v_sep)
out(f'      {len(bewert):,} von {len(test_g):,} zurueckgehaltenen g bewertbar (alle Stufen im Training bekannt)')
out(f'      SEPARABLES MODELL:   R² = {r2_sep:.3f}   MAE = {float(np.abs(e_sep-v_sep).mean()):.2f}'
    f'   bei Mittel {e_sep.mean():.1f}')

#   The strict version covers only part of the held-out g, namely the EASY ones (frequent cofactors). Hence a second
#   version with FALLBACK: unknown stages receive the mean effect from training. All held-out g are then evaluable, and
#   selection no longer works as an objection.
a_fb = a_g.copy(); a_fb[~bh] = float(a_g[bh].mean()); am = float(a_g[bh].mean())
b_fb = b_g.copy(); b_fb[~bk] = float(b_g[bk].mean()); bm = float(b_g[bk].mean())
# `vorhersage_fb` = prediction with fallback (`fb`): `am`/`bm` = mean effects used for unknown stages
def vorhersage_fb(g, a_=None, b_=None):
    a_ = a_fb if a_ is None else a_; b_ = b_fb if b_ is None else b_
    s = 0.0
    for h in PT[g]:
        k = g//h
        s += float(np.exp((a_[hi[h]] if h in hi else am) + (b_[ki[k]] if k in ki else bm)))
    return s
e_full = np.array([occ[g] for g in test_g], float)
v_full = np.array([vorhersage_fb(g) for g in test_g], float)
r2_sep_full = r2(e_full, v_full)
out(f'      mit RUECKFALL auf den mittleren Effekt, alle {len(test_g):,} zurueckgehaltenen g:')
out(f'         R² = {r2_sep_full:.3f}   MAE = {float(np.abs(e_full-v_full).mean()):.2f}'
    f'   bei Mittel {e_full.mean():.1f}')
out()

# Baseline and floor on the SAME set (all held-out g) and the SAME split
# `merk(g)` = raw features of g (`merkmale`): 1, log g, log(number of powerful divisors), log sqfull, omega, Omega
def merk(g):
    fs = FS[g]
    return [1.0, np.log(g), np.log(len(PT[g])), np.log(SQ[g]) if SQ[g] > 1 else 0.0,
            float(len(fs)), float(sum(fs.values()))]
# `trg` = training g; `Xtr`/`ytr` = training features / targets, `Xte` = test features
trg = [g for g in gl if gtrain[gi[g]]]
Xtr = np.array([merk(g) for g in trg]); ytr = np.array([occ[g] for g in trg], float)
Xte = np.array([merk(g) for g in test_g])
# `lstsq_r2` = least-squares fit on the training data (once on the counts, once on the log counts); returns both test R^2
def lstsq_r2(Xa, ya, Xb, yb):
    kl, *_ = np.linalg.lstsq(Xa, ya, rcond=None)
    kL, *_ = np.linalg.lstsq(Xa, np.log(ya), rcond=None)
    return r2(yb, Xb @ kl), r2(yb, np.exp(Xb @ kL))
r2_roh_lin, r2_roh_log = lstsq_r2(Xtr, ytr, Xte, e_full)
r2_roh = max(r2_roh_lin, r2_roh_log)
out(f'      GRUNDLINIE (rohe Praediktoren, beste von zwei Fassungen), dieselben {len(test_g):,} Werte:')
out(f'         auf Zaehlung R² = {r2_roh_lin:.3f} · auf log-Zaehlung R² = {r2_roh_log:.3f}  ⇒ {r2_roh:.3f}')
out(f'      ⇒ Abstand des separablen Modells: {r2_sep_full - r2_roh:+.3f}')
out(f'      (die 0,333 vom 18.09. gilt fuer g <= 10^6 und ist mit dieser Menge nicht vergleichbar —')
out(f'       deshalb hier neu gerechnet.)')
out()

# THE FLOOR, same machinery: only the number of powerful divisors
Btr = np.array([[1.0, float(len(PT[g]))] for g in trg])
Bte = np.array([[1.0, float(len(PT[g]))] for g in test_g])
r2_bod_lin, r2_bod_log = lstsq_r2(Btr, ytr, Bte, e_full)
r2_bod = max(r2_bod_lin, r2_bod_log)
out(f'      BODEN — nur die Anzahl der powerful Teiler von g, sonst nichts:')
out(f'         auf Zaehlung R² = {r2_bod_lin:.3f} · auf log-Zaehlung R² = {r2_bod_log:.3f}  ⇒ {r2_bod:.3f}')
out(f'      ⇒ Zugewinn des SATZES ueber die bloße Teilerzahl: {r2_sep_full - r2_bod:+.3f}')
out()

# ------------------------------------------------------------------ (C) PK- (negative control)
# `wp` = permuted counts, `yp` = their logs, `occ_p` = permuted occurrences per g,
# `a_p`/`b_p` = fit on the permuted training cells
out('   (C) PK− — permutierte Zaehlungen, derselbe g-Split')
wp = cc.copy(); rng.shuffle(wp); yp = np.log(wp)
occ_p = {g: 0.0 for g in gl}
for i, g in enumerate(gg.tolist()): occ_p[g] += wp[i]
a_p, b_p = fit(ih[zell_tr], ik[zell_tr], yp[zell_tr])
# `vorhersage_p` = strict prediction from the fit on the permuted counts
def vorhersage_p(g):
    s = 0.0
    for h in PT[g]:
        k = g//h
        if h not in hi or k not in ki: return None
        iH, iK = hi[h], ki[k]
        if not (bh[iH] and bk[iK]): return None
        s += float(np.exp(a_p[iH] + b_p[iK]))
    return s
a_pf = a_p.copy(); a_pf[~bh] = float(a_p[bh].mean())
b_pf = b_p.copy(); b_pf[~bk] = float(b_p[bk].mean())
e_fp = np.array([occ_p[g] for g in test_g], float)
v_fp = np.array([vorhersage_fb(g, a_pf, b_pf) for g in test_g], float)
r2_sep_p = r2(e_fp, v_fp)
out(f'      permutiert R² = {r2_sep_p:.3f}   gegen echt {r2_sep_full:.3f}')
out('      Die Permutation mischt die ZAEHLUNGEN, laesst aber die Traegerstruktur intakt —')
out('      sie ist gegen den BODEN zu lesen, nicht gegen 0.')
out()

# --------------------------------------------- WHY the fallback fails: the cofactors are unique
# `kg` = cofactor -> set of the g containing it; `einmal` = number of cofactors that occur in exactly one g
kg = {}
for (g,h) in schl: kg.setdefault(g//h, set()).add(g)
einmal = sum(1 for k, s in kg.items() if len(s) == 1)
out('   (C2) DIE URSACHE — wie einmalig sind die Kofaktoren?')
out(f'      {nK:,} verschiedene Kofaktoren fuer {len(gl):,} Lueckenwerte  ⇒ {nK/len(gl):.2f} je Wert')
out(f'      davon in genau EINEM g: {einmal:,} ({einmal/nK:.0%});  Median der g je Kofaktor: '
    f'{int(np.median([len(s) for s in kg.values()]))}')
out(f'      ⇒ {1-len(bewert)/len(test_g):.0%} der zurueckgehaltenen g bringen mindestens einen nie')
out('        gesehenen Kofaktor mit. Ein freier Parameter je Kofaktor hat dort nichts zu uebertragen.')
out()

# --------------------------------------------- (D) the STRONGEST version: beta as a smooth function of k
out('   (D) DIE STAERKSTE FASSUNG — alpha und beta als glatte Funktionen der Arithmetik')
out('      (nicht ein freier Parameter je Stufe, sondern aus log, omega, Omega, quadratfrei, log sqfull)')
t1=time.time()
# `arith(n)` = arithmetic features of n: 1, log n, omega, Omega, squarefree indicator, log of the squarefull part
def arith(n):
    fs = FS.get(n) or faktor(n)
    sf = sqfull(fs)
    return [1.0, np.log(n) if n > 1 else 0.0, float(len(fs)), float(sum(fs.values())),
            1.0 if all(e == 1 for e in fs.values()) else 0.0, np.log(sf) if sf > 1 else 0.0]
# regress alpha(h) and beta(k) on their arithmetic, using only the stages seen in training
#   stages WEIGHTED by their number of cells (a stage with 80 cells is more reliable than one with a single cell).
# `nzh`/`nzk` = number of training cells per kernel / cofactor
nzh = np.bincount(ih[zell_tr], minlength=nH).astype(float)
nzk = np.bincount(ik[zell_tr], minlength=nK).astype(float)
def wlstsq(X, yv, w):    # weighted least squares
    s = np.sqrt(w)[:, None]
    k, *_ = np.linalg.lstsq(X*s, yv*np.sqrt(w), rcond=None)
    return k
# `Ha`/`Kb` = features of the known kernels / cofactors, `ca`/`cb` = fitted coefficients for alpha / beta
Ha = np.array([arith(H[i]) for i in np.nonzero(bh)[0]]); ya_ = a_g[bh]; wa = nzh[bh]
Kb = np.array([arith(K[i]) for i in np.nonzero(bk)[0]]); yb_ = b_g[bk]; wb = nzk[bk]
ca = wlstsq(Ha, ya_, wa); cb = wlstsq(Kb, yb_, wb)
r2_alpha = r2(ya_, Ha @ ca); r2_beta = r2(yb_, Kb @ cb)
out(f'      alpha(h) durch Arithmetik von h erklaert: R² = {r2_alpha:.3f}  ({int(bh.sum()):,} Kerne)')
out(f'      beta(k)  durch Arithmetik von k erklaert: R² = {r2_beta:.3f}  ({int(bk.sum()):,} Kofaktoren)')
def vorhersage_glatt(g):    # `glatt` = smooth: prediction from the smooth alpha and beta
    s = 0.0
    for h in PT[g]:
        s += float(np.exp(float(np.dot(arith(h), ca)) + float(np.dot(arith(g//h), cb))))
    return s
v_glatt = np.array([vorhersage_glatt(g) for g in test_g], float)
r2_glatt_roh = r2(e_full, v_glatt)
out(f'      ⇒ occ(g) roh, alle {len(test_g):,} zurueckgehaltenen g:  R² = {r2_glatt_roh:.3f}'
    f'   MAE = {float(np.abs(e_full-v_glatt).mean()):.2f}')
#   RECALIBRATION: a sum fitted in the log and then exponentiated is biased as a count value.
#   Affine correction, determined exclusively on the TRAINING part.
vp_tr  = np.array([vorhersage_glatt(g) for g in trg], float)
vfb_tr = np.array([vorhersage_fb(g)    for g in trg], float)
def rekal(vtr, vte):    # `rekal` = recalibration: affine fit of the training counts on the training predictions, applied to `vte`
    A = np.column_stack([np.ones(len(vtr)), vtr])
    k, *_ = np.linalg.lstsq(A, ytr, rcond=None)
    return k[0] + k[1]*vte, k
vg_k, kg_ = rekal(vp_tr, v_glatt)
vf_k, kf_ = rekal(vfb_tr, v_full)
r2_glatt = r2(e_full, vg_k); r2_fb_k = r2(e_full, vf_k)
out(f'      ⇒ nach affiner Rueckkalibrierung (nur Trainingsteil, {kg_[0]:+.2f} {kg_[1]:+.3f}·x):')
out(f'         glatte Stufen  R² = {r2_glatt:.3f}   MAE = {float(np.abs(e_full-vg_k).mean()):.2f}')
out(f'         freie Stufen   R² = {r2_fb_k:.3f}   MAE = {float(np.abs(e_full-vf_k).mean()):.2f}'
    f'   ({time.time()-t1:.0f}s)')
out()

out('   ⇒ DIE BILANZ, an der D2 haengt   (alle {:,} zurueckgehaltenen g, ein Split):'.format(len(test_g)))
out(f'      Boden — nur Anzahl powerful Teiler                     {r2_bod:+.3f}')
out(f'      permutierte Zaehlungen                                 {r2_sep_p:+.3f}')
out(f'      rohe Praediktoren (6 Merkmale von g)                   {r2_roh:+.3f}')
out(f'      separables Modell, freie Stufen, roh                    {r2_sep_full:+.3f}')
out(f'      separables Modell, freie Stufen, rueckkalibriert        {r2_fb_k:+.3f}')
out(f'      separables Modell, glatte Stufen, roh                   {r2_glatt_roh:+.3f}')
out(f'      separables Modell, glatte Stufen, rueckkalibriert  ⟵ staerkste Fassung  {r2_glatt:+.3f}')
out(f'      dasselbe, freie Stufen, strenge Teilmenge ({len(bewert)/len(test_g):.0%})        {r2_sep:+.3f}')
out(f'      ⇒ Zugewinn der staerksten Fassung ueber die rohen Praediktoren: {r2_glatt - r2_roh:+.3f}')

out()

# ------------------------------------------------------------------ (E) SPREAD over several splits
# Two numbers that differ by 0.04 are no difference without a spread. Five splits, same machinery, only the random
# partition changes. The two contenders compared are:
#   raw arithmetic features of g   versus   separable model (free stages, recalibrated).
# `serie` = per split (best R^2 of the raw features, R^2 of the separable model); `saat` = seed;
# `vh` = prediction of the separable model; `MERK` = raw features per g
out('   (E) STREUUNG — fuenf unabhaengige Aufteilungen, dieselbe Maschinerie')
MERK = {g: merk(g) for g in gl}
serie = []
for saat in (1, 2, 3, 4, 5):
    rg = np.random.default_rng(20260920 + saat)
    gtr2 = rg.random(len(gl)) < 0.80
    ztr2 = np.fromiter((gtr2[gi[g]] for g in gg.tolist()), bool, nZ)
    a2, b2 = fit(ih[ztr2], ik[ztr2], y[ztr2])
    s2h = np.bincount(ih[ztr2], minlength=nH) > 0; s2k = np.bincount(ik[ztr2], minlength=nK) > 0
    a2f = a2.copy(); a2f[~s2h] = float(a2[s2h].mean()); am2 = float(a2[s2h].mean())
    b2f = b2.copy(); b2f[~s2k] = float(b2[s2k].mean()); bm2 = float(b2[s2k].mean())
    ea = np.exp(a2f); eb = np.exp(b2f)
    def vh(g):
        s = 0.0
        for h in PT[g]:
            k = g//h
            s += (ea[hi[h]] if h in hi else np.exp(am2)) * (eb[ki[k]] if k in ki else np.exp(bm2))
        return float(s)
    tr2 = [g for g in gl if gtr2[gi[g]]]; te2 = [g for g in gl if not gtr2[gi[g]]]
    yt2 = np.array([occ[g] for g in tr2], float); ye2 = np.array([occ[g] for g in te2], float)
    vt2 = np.array([vh(g) for g in tr2]); ve2 = np.array([vh(g) for g in te2])
    A = np.column_stack([np.ones(len(vt2)), vt2]); kk, *_ = np.linalg.lstsq(A, yt2, rcond=None)
    r_sep = r2(ye2, kk[0] + kk[1]*ve2)
    X2 = np.array([MERK[g] for g in tr2]); Xe2 = np.array([MERK[g] for g in te2])
    rl, rlog = lstsq_r2(X2, yt2, Xe2, ye2)
    serie.append((max(rl, rlog), r_sep))
    out(f'      Split {saat}:  rohe Merkmale {max(rl,rlog):+.3f}   separables Modell {r_sep:+.3f}'
        f'   Abstand {r_sep-max(rl,rlog):+.3f}')
ser = np.array(serie)
out(f'      ⇒ rohe Merkmale     {ser[:,0].mean():+.3f} ± {ser[:,0].std(ddof=1):.3f}')
out(f'      ⇒ separables Modell {ser[:,1].mean():+.3f} ± {ser[:,1].std(ddof=1):.3f}')
ab = ser[:,1]-ser[:,0]
out(f'      ⇒ ABSTAND           {ab.mean():+.3f} ± {ab.std(ddof=1):.3f}'
    f'   (kleinster {ab.min():+.3f}, groesster {ab.max():+.3f})')
out()

# ------------------------------------------------------------------ (F) THE PRE-REGISTERED P1 / P2
# Two predictions were stated in advance, in ONE specific decomposition: g = sqfull(g) * kc with kc squarefree.
# P1: for fixed sqfull, occ depends only on kc.  P2: for fixed kc, occ scales with the kernel.
# Together these are exactly:   log occ(g) = u(sqfull(g)) + v(kc).  This is a DIFFERENT, narrower claim than the cell
# model (A) -- and it is the one that was written down. So it is checked directly, at the level of g, without cells.
# `KC` = kc per g; `S_l`/`C_l` = sorted distinct kernels sqfull / cofactors kc; `si`/`ci` = index maps; `iS`/`iC` = index per g;
# `yo` = log occ per g
out('   (F) DIE VORREGISTRIERTEN VORHERSAGEN P1 und P2  (Roadmap-Wortlaut, Zerlegung g = sqfull(g)·kc)')
KC = {g: g//SQ[g] for g in gl}
# `schlimm` = offending g: those whose cofactor kc is not squarefree
schlimm = [g for g in gl if not all(e == 1 for e in faktor(KC[g]).values())]
assert not schlimm, ('kc muss quadratfrei sein', schlimm[:5])
S_l = sorted({SQ[g] for g in gl}); C_l = sorted({KC[g] for g in gl})
si = {s:i for i,s in enumerate(S_l)}; ci = {c:i for i,c in enumerate(C_l)}
out(f'      {len(S_l):,} verschiedene Kerne sqfull(g), {len(C_l):,} verschiedene quadratfreie Kofaktoren')
out(f'      ✅ alle {len(gl):,} Kofaktoren kc sind quadratfrei — die Zerlegung ist die behauptete.')
iS = np.array([si[SQ[g]] for g in gl]); iC = np.array([ci[KC[g]] for g in gl])
yo = np.log(np.array([occ[g] for g in gl], float))
def fit2(m, runden=80):    # backfitting of log occ(g) = u(sqfull) + v(kc) on the g selected by the mask m
    u = np.zeros(len(S_l)); v = np.zeros(len(C_l))
    cu = np.maximum(np.bincount(iS[m], minlength=len(S_l)), 1)
    cv = np.maximum(np.bincount(iC[m], minlength=len(C_l)), 1)
    for _ in range(runden):
        u = np.bincount(iS[m], weights=yo[m]-v[iC[m]], minlength=len(S_l))/cu
        v = np.bincount(iC[m], weights=yo[m]-u[iS[m]], minlength=len(C_l))/cv
        v -= v.mean()
    return u, v
# `mg` = training mask over the g (same split as in (B)); `mt` (below) = test g with known stages
mg = np.array([gtrain[gi[g]] for g in gl])
u_a, v_a = fit2(np.ones(len(gl), bool))
r2_F_in = r2(yo, u_a[iS]+v_a[iC])
u_t, v_t = fit2(mg)
sS = np.bincount(iS[mg], minlength=len(S_l)) > 0; sC = np.bincount(iC[mg], minlength=len(C_l)) > 0
mt = (~mg) & sS[iS] & sC[iC]
r2_F_out = r2(yo[mt], u_t[iS[mt]]+v_t[iC[mt]]) if mt.sum() > 30 else float('nan')
rF = yo - (u_a[iS]+v_a[iC])
kor_F = float(np.corrcoef(u_a[iS]*v_a[iC], rF)[0,1])
out(f'      log occ(g) = u(sqfull) + v(kc):  in-sample R² = {r2_F_in:.3f}')
out(f'         ausserhalb, {int(mt.sum()):,} von {int((~mg).sum()):,} zurueckgehaltenen g bewertbar '
    f'({mt.sum()/max(1,(~mg).sum()):.0%}):  R² = {r2_F_out:.3f}')
out(f'      Wechselwirkung u·v gegen Residuen: r = {kor_F:+.3f}'
    f'   ⇒ {"P1/P2 vertraeglich" if abs(kor_F) < 0.1 else "P1/P2 verletzt"}')
# P2 on its own: for fixed kc, is the kernel effect the same? Only kc with >= 3 different kernels (`grp`, `mehr`, `sp` = spreads).
grp = {}
for g in gl: grp.setdefault(KC[g], []).append(g)
mehr = [c for c, l in grp.items() if len(l) >= 3]
sp = [float(np.std([np.log(occ[g]) - u_a[si[SQ[g]]] for g in grp[c]], ddof=1)) for c in mehr[:4000]]
out(f'      P2 einzeln: {len(mehr):,} Kofaktoren mit >= 3 Kernen; Streuung von log occ - u(Kern)')
out(f'         Median {float(np.median(sp)):.3f}, gegen Gesamtstreuung von log occ {float(yo.std(ddof=1)):.3f}')
out()

# ------------------------------------------------------------------ (G) THE DECISIVE COMPARISON
# The pre-registered model (F) against the raw features -- on BOTH scales and over five splits.
# R^2 on log occ and R^2 on occ are different quantities; holding them against each other would be the same error as
# comparing different populations (see the note on the baseline above).
# `occ_a` = occ per g; `MK` = raw features per g; `reihen` = per split (model log, raw log, model count, raw count)
out('   (G) ENTSCHEIDUNG — vorregistriertes Modell gegen rohe Merkmale, beide Skalen, fuenf Splits')
occ_a = np.array([occ[g] for g in gl], float)
MK = np.array([merk(g) for g in gl])
reihen = []
for saat in (1,2,3,4,5):
    rg = np.random.default_rng(20260920 + saat)
    m = rg.random(len(gl)) < 0.80
    # pre-registered model
    u = np.zeros(len(S_l)); v = np.zeros(len(C_l))
    cu = np.maximum(np.bincount(iS[m], minlength=len(S_l)),1); cv = np.maximum(np.bincount(iC[m], minlength=len(C_l)),1)
    for _ in range(80):
        u = np.bincount(iS[m], weights=yo[m]-v[iC[m]], minlength=len(S_l))/cu
        v = np.bincount(iC[m], weights=yo[m]-u[iS[m]], minlength=len(C_l))/cv
        v -= v.mean()
    kS = np.bincount(iS[m], minlength=len(S_l))>0; kC = np.bincount(iC[m], minlength=len(C_l))>0
    uf = u.copy(); uf[~kS] = float(u[kS].mean()); vf = v.copy(); vf[~kC] = float(v[kC].mean())
    pl = uf[iS] + vf[iC]                       # prediction on the log scale, all g
    # log scale
    r2F_log = r2(yo[~m], pl[~m])
    kl, *_ = np.linalg.lstsq(np.column_stack([np.ones(int(m.sum())), MK[m][:,1:] @ np.eye(5)]), yo[m], rcond=None) if False else (None,)
    klog, *_ = np.linalg.lstsq(MK[m], yo[m], rcond=None)
    r2R_log = r2(yo[~m], MK[~m] @ klog)
    # count scale, both recalibrated affinely on the training part
    def rk(p):
        A = np.column_stack([np.ones(int(m.sum())), p[m]])
        k, *_ = np.linalg.lstsq(A, occ_a[m], rcond=None)
        return k[0] + k[1]*p[~m]
    r2F_cnt = r2(occ_a[~m], rk(np.exp(pl)))
    kcnt, *_ = np.linalg.lstsq(MK[m], occ_a[m], rcond=None)
    r2R_cnt = r2(occ_a[~m], MK[~m] @ kcnt)
    reihen.append((r2F_log, r2R_log, r2F_cnt, r2R_cnt))
    out(f'      Split {saat}:  log occ — Modell {r2F_log:+.3f} / roh {r2R_log:+.3f}'
        f'      occ — Modell {r2F_cnt:+.3f} / roh {r2R_cnt:+.3f}')
RE = np.array(reihen)
for j, nm in ((0,'log occ, vorregistriertes Modell'), (1,'log occ, rohe Merkmale'),
              (2,'occ,     vorregistriertes Modell'), (3,'occ,     rohe Merkmale')):
    out(f'      ⇒ {nm:34s} {RE[:,j].mean():+.3f} ± {RE[:,j].std(ddof=1):.3f}')
dl = RE[:,0]-RE[:,1]; dc = RE[:,2]-RE[:,3]
out(f'      ⇒ ABSTAND auf log occ: {dl.mean():+.3f} ± {dl.std(ddof=1):.3f}'
    f'   (kleinster {dl.min():+.3f})')
out(f'      ⇒ ABSTAND auf occ:     {dc.mean():+.3f} ± {dc.std(ddof=1):.3f}'
    f'   (kleinster {dc.min():+.3f})')
out()

# ------------------------------------------------------------------ (H) WHERE does the predictive power sit?
# The raw features already know log sqfull(g). If the gain came from the kernel, it would already be contained there.
# So: kernel alone, cofactor alone, both -- on the count scale, five splits. (`abl` = ablation)
out('   (H) WO sitzt die Kraft — Kern u(sqfull) oder Kofaktor v(kc)?')
abl = []
for saat in (1,2,3,4,5):
    rg = np.random.default_rng(20260920 + saat)
    m = rg.random(len(gl)) < 0.80
    cu = np.maximum(np.bincount(iS[m], minlength=len(S_l)),1); cv = np.maximum(np.bincount(iC[m], minlength=len(C_l)),1)
    # both
    u = np.zeros(len(S_l)); v = np.zeros(len(C_l))
    for _ in range(80):
        u = np.bincount(iS[m], weights=yo[m]-v[iC[m]], minlength=len(S_l))/cu
        v = np.bincount(iC[m], weights=yo[m]-u[iS[m]], minlength=len(C_l))/cv
        v -= v.mean()
    kS = np.bincount(iS[m], minlength=len(S_l))>0; kC = np.bincount(iC[m], minlength=len(C_l))>0
    uf = u.copy(); uf[~kS] = float(u[kS].mean()); vf = v.copy(); vf[~kC] = float(v[kC].mean())
    # kernel alone / cofactor alone (each fitted on its own, not cut out of the joint solution)
    u1 = np.bincount(iS[m], weights=yo[m], minlength=len(S_l))/cu; u1[~kS] = float(u1[kS].mean())
    v1 = np.bincount(iC[m], weights=yo[m], minlength=len(C_l))/cv; v1[~kC] = float(v1[kC].mean())
    def rk(p):
        A = np.column_stack([np.ones(int(m.sum())), p[m]])
        k, *_ = np.linalg.lstsq(A, occ_a[m], rcond=None)
        return r2(occ_a[~m], k[0] + k[1]*p[~m])
    abl.append((rk(np.exp(u1[iS])), rk(np.exp(v1[iC])), rk(np.exp(uf[iS]+vf[iC]))))
AB = np.array(abl)
out(f'      nur Kern u(sqfull(g))        {AB[:,0].mean():+.3f} ± {AB[:,0].std(ddof=1):.3f}')
out(f'      nur Kofaktor v(kc)           {AB[:,1].mean():+.3f} ± {AB[:,1].std(ddof=1):.3f}')
out(f'      beide                        {AB[:,2].mean():+.3f} ± {AB[:,2].std(ddof=1):.3f}')
out(f'      rohe Merkmale (aus G)        {RE[:,3].mean():+.3f} ± {RE[:,3].std(ddof=1):.3f}')
out(f'      ⇒ Zugewinn des Kofaktors ueber den Kern allein: {AB[:,2].mean()-AB[:,0].mean():+.3f}')
out(f'      ⇒ Zugewinn des Kerns ueber den Kofaktor allein: {AB[:,2].mean()-AB[:,1].mean():+.3f}')
out()
out(f'Zeit {time.time()-t0:.0f}s   — Urteil in den Record, hier nur Zahlen.')
# Keys of the result JSON: `n_werte` = number of gap values, `groesstes_g` = largest g, `n_zellen` = number of cells,
# `n_kerne` = number of kernels, `n_kofaktoren` = number of cofactors, `par_pro_punkt` = parameters per data point.
# `A_zellsplit` = part (A): `r2_aussen` = R^2 outside the training cells, `n_test` = number of test cells, `wechselwirkung_r`
# = interaction correlation. `B_gsplit` = part (B): `r2_separabel_alle` = R^2 of the separable model on all held-out g (with
# fallback), `mae` = mean absolute error, `mittel` = mean, `r2_separabel_streng` = strict version, `n_bewertbar_streng` /
# `anteil_bewertbar_streng` = number / share of evaluable g, `r2_grundlinie*` = baseline R^2 (`zaehlung` = on counts, `log` =
# on log counts), `r2_boden_teilerzahl` = floor R^2 (number of divisors only), `abstand_zur_grundlinie` = gap to the baseline,
# `zugewinn_ueber_staerksten_vergleich` = gain over the strongest comparison. `C_permutiert` = part (C), permuted counts.
# `H_ablation` = part (H): `nur_kern` = kernel only, `nur_kofaktor` = cofactor only, `beide` = both, `_mittel` = mean, `_sd` =
# standard deviation. `G_entscheidung` = part (G): `log_modell`/`cnt_modell` = pre-registered model on the log / count scale,
# `log_roh`/`cnt_roh` = raw features, `abstand_*` = difference. `F_vorregistriert` = part (F): `n_kerne_sqfull` = number of
# kernels sqfull, `n_kofaktoren_quadratfrei` = number of squarefree cofactors, `n_bewertbar` = evaluable g,
# `n_zurueckgehalten` = held-out g, `P2_streuung_median` = median spread for P2, `P2_n_kofaktoren_mit_3_kernen` = cofactors
# with 3 kernels, `gesamtstreuung_log_occ` = total spread of log occ. `E_streuung` = part (E): `rohe_merkmale` = raw features,
# `separabel` = separable model, `rohe_*` / `sep_*` = mean and sd of each, `abstand_*` = difference. `C2_kofaktoren` = part
# (C2): `je_wert` = per gap value, `nur_ein_g` = cofactors in only one g, `anteil_*` = share. `D_glatt` = part (D):
# `r2_alpha_durch_arithmetik` / `r2_beta_durch_arithmetik` = R^2 of alpha / beta explained by arithmetic, `r2_occ_roh` = raw,
# `r2_occ_rueckkalibriert` = recalibrated, `r2_freie_stufen_rueckkalibriert` = free stages recalibrated, `rekalibrierung` =
# recalibration coefficients, `zugewinn_ueber_rohe_praediktoren` = gain over raw predictors. `hinweis_18_09` = note on the
# earlier baseline (R^2 = 0.333 for g <= 10^6), `laufzeit_s` = run time in seconds.
(ERG/'w59_zaehlmodell_output.txt').write_text('\n'.join(Z)+'\n', encoding='utf-8')
(ERG/'w59_zaehlmodell_result.json').write_text(json.dumps(dict(
    skript='w59_zaehlmodell_2026-09-20.py', N=N, MINCNT=MINCNT,
    n_werte=len(gl), groesstes_g=int(sel.max()), n_zellen=nZ,
    n_kerne=nH, n_kofaktoren=nK, par_pro_punkt=(nH+nK-1)/nZ,
    A_zellsplit=dict(r2_in_sample=r2_in, r2_aussen=r2_zell_out, n_test=int(tc.sum()),
                     wechselwirkung_r=kor_x),
    B_gsplit=dict(n_test=len(test_g),
                  r2_separabel_alle=r2_sep_full, mae=float(np.abs(e_full-v_full).mean()),
                  mittel=float(e_full.mean()),
                  r2_separabel_streng=r2_sep, n_bewertbar_streng=len(bewert),
                  anteil_bewertbar_streng=len(bewert)/len(test_g),
                  r2_grundlinie=r2_roh, r2_grundlinie_zaehlung=r2_roh_lin, r2_grundlinie_log=r2_roh_log,
                  r2_boden_teilerzahl=r2_bod,
                  abstand_zur_grundlinie=r2_sep_full-r2_roh,
                  zugewinn_ueber_staerksten_vergleich=r2_sep_full-max(r2_bod, r2_roh, r2_sep_p)),
    C_permutiert=dict(r2_separabel=r2_sep_p),
    H_ablation=dict(nur_kern=[float(x) for x in AB[:,0]], nur_kofaktor=[float(x) for x in AB[:,1]],
                    beide=[float(x) for x in AB[:,2]],
                    nur_kern_mittel=float(AB[:,0].mean()), nur_kern_sd=float(AB[:,0].std(ddof=1)),
                    nur_kofaktor_mittel=float(AB[:,1].mean()), nur_kofaktor_sd=float(AB[:,1].std(ddof=1)),
                    beide_mittel=float(AB[:,2].mean()), beide_sd=float(AB[:,2].std(ddof=1))),
    G_entscheidung=dict(n_splits=len(reihen),
        log_modell=[float(x) for x in RE[:,0]], log_roh=[float(x) for x in RE[:,1]],
        cnt_modell=[float(x) for x in RE[:,2]], cnt_roh=[float(x) for x in RE[:,3]],
        log_modell_mittel=float(RE[:,0].mean()), log_modell_sd=float(RE[:,0].std(ddof=1)),
        log_roh_mittel=float(RE[:,1].mean()), log_roh_sd=float(RE[:,1].std(ddof=1)),
        cnt_modell_mittel=float(RE[:,2].mean()), cnt_modell_sd=float(RE[:,2].std(ddof=1)),
        cnt_roh_mittel=float(RE[:,3].mean()), cnt_roh_sd=float(RE[:,3].std(ddof=1)),
        abstand_log_mittel=float(dl.mean()), abstand_log_sd=float(dl.std(ddof=1)),
        abstand_cnt_mittel=float(dc.mean()), abstand_cnt_sd=float(dc.std(ddof=1))),
    F_vorregistriert=dict(n_kerne_sqfull=len(S_l), n_kofaktoren_quadratfrei=len(C_l),
                          r2_in_sample=r2_F_in, r2_aussen=r2_F_out,
                          n_bewertbar=int(mt.sum()), n_zurueckgehalten=int((~mg).sum()),
                          wechselwirkung_r=kor_F,
                          P2_streuung_median=float(np.median(sp)),
                          P2_n_kofaktoren_mit_3_kernen=len(mehr),
                          gesamtstreuung_log_occ=float(yo.std(ddof=1))),
    E_streuung=dict(n_splits=len(serie),
                    rohe_merkmale=[float(v) for v in ser[:,0]],
                    separabel=[float(v) for v in ser[:,1]],
                    rohe_mittel=float(ser[:,0].mean()), rohe_sd=float(ser[:,0].std(ddof=1)),
                    sep_mittel=float(ser[:,1].mean()), sep_sd=float(ser[:,1].std(ddof=1)),
                    abstand_mittel=float(ab.mean()), abstand_sd=float(ab.std(ddof=1)),
                    abstand_min=float(ab.min()), abstand_max=float(ab.max())),
    C2_kofaktoren=dict(n_kofaktoren=nK, je_wert=nK/len(gl), nur_ein_g=einmal,
                       anteil_nur_ein_g=einmal/nK,
                       anteil_g_mit_unbekanntem_kofaktor=1-len(bewert)/len(test_g)),
    D_glatt=dict(r2_alpha_durch_arithmetik=r2_alpha, r2_beta_durch_arithmetik=r2_beta,
                 r2_occ_roh=r2_glatt_roh, r2_occ_rueckkalibriert=r2_glatt,
                 r2_freie_stufen_rueckkalibriert=r2_fb_k,
                 mae_rueckkalibriert=float(np.abs(e_full-vg_k).mean()),
                 rekalibrierung=[float(kg_[0]), float(kg_[1])],
                 zugewinn_ueber_rohe_praediktoren=r2_glatt-r2_roh),
    hinweis_18_09='R²=0,333 galt fuer g <= 10^6, nicht vergleichbar; Grundlinie hier neu gerechnet',
    laufzeit_s=time.time()-t0, python=sys.version.split()[0]), indent=1), encoding='utf-8')
