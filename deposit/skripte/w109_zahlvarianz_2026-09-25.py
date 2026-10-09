# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w109_zahlvarianz_2026-09-25.py
# Purpose: number variance Σ²(L) (`zahlvarianz`) of several point sets, as a test for order over long ranges: the variance of the
#   NUMBER of unfolded points in a randomly placed window of length L (in mean spacings).
#   Sets: K- uniform random points, K+ GUE eigenvalues, Z zeta zeros, W Wieferich quotients (kernel 15), P primes, Q powerful
#   numbers in [2^44, 2^46), Q2 the squares among them, Q3 the powerful numbers that are not squares.
# Reads: ergebnisse/w105_zeta_nullstellen.json (key "imag"). Requires numpy and mpmath.
# Writes: ergebnisse/w109_zahlvarianz_result.json and ergebnisse/w109_zahlvarianz_output.txt.
# Usage: python w109_zahlvarianz_2026-09-25.py (no arguments; the results folder is found relative to the script).
# Controls: the predictions are computed and printed BEFORE any measurement. PK1: K- must match Poisson (deviation < 6 %); PK2: K+
#   must match the GUE formula (deviation < 8 % from L = 0.5 on); PK3: Q2 must match {L}(1 - {L}) (|difference| < 0.01).

# Background: w105 measured only the distance to the NEIGHBOR (short range). Σ²(L) asks: if a window of length L (in mean
#   spacings) is placed at random on the unfolded points, how strongly does the NUMBER of points in it fluctuate? Poisson: Σ² = L
#   (the fluctuation grows along). Random matrices (GUE): Σ² grows only logarithmically. Lattice: Σ² <= 1/4, bounded.
# Method: for each set, the unfolded points x (mean spacing 1); window starts t uniformly distributed; count in [t, t + L);
#   Σ² = variance. Error bars: split the set into blocks, Σ² per block, SE = standard deviation / √(blocks) (for K+ every
#   matrix is a block). L = 0.25 ... 5 (step 0.25), then 6, 8, 10, 12, 15, 20, 25, 30, 40, 50.
#
# Expectations (the predictions are computed before the measurement and printed):
#   PK1 K-  uniform: Σ² = L, deviation < 6 %.
#   PK2 K+ GUE (30 matrices 1000 × 1000, middle half, semicircle radius 2√N): Σ² = L − 2∫₀ᴸ (L − r)(sin πr/πr)² dr, deviation < 8
#       % from L = 0.5 on.
#   PK3 Q₂  squares, unfolded with √X (= the integers a themselves): Σ² = {L}(1 − {L}) exactly, |Δ| < 0.01.
#   E1 Z zeta zeros 1 ... 2000: like GUE at small L, flatter at larger L (Berry saturation) [plausible]; only 2000 points, so
#       large error bars.
#   E2  W   Wieferich quotients, kernel 15: Poisson [plausible].
#   E3  P   primes in [10⁹, 10⁹ + 10⁷): BELOW Poisson. Prediction from the Hardy–Littlewood model, computed here:
#           Var = L(1 − 1/ln N) + (2/ln²N)·Σ_(1≤d<h) (h − d)(S(d) − 1),  h = L·ln N,  S(d) = singular series (0 for odd d,
#           otherwise 2·C₂·∏_(p|d, p>2) (p−1)/(p−2); C₂ = ∏_(p>2) (1 − 1/(p−1)²) is computed as well).
#           [Assumption: the HL model holds here.]
#   E4  Q₃, Q: SUPERPOSITION OF INDEPENDENT LATTICES. Every family a²b³ (b squarefree fixed) is a lattice with spacing b^(3/2)
#           in √X; the √b for different squarefree b are linearly independent over ℚ ⇒ the lattice phases are jointly
#           equidistributed ⇒ the counts are uncorrelated in the window average ⇒ Σ²(L) = Σ_b f(L·w_b),
#           f(μ) = {μ}(1 − {μ}), w_b = share of family b in the window (COUNTED). [plausible; the independence of the roots
#           is classical (Besicovitch 1940); the MEASUREMENT decides.]
#   E5  asymptotically (sum → integral over b, density of squarefree numbers 6/π²):  Σ²(L) ≈ K·L^(2/3),
#           K = (4/π²)·c^(−2/3)·I,  I = ∫₀^∞ {μ}(1 − {μ}) μ^(−5/3) dμ,  c = c₁ (Q) resp. c₁ − 1 (Q₃). An exponent BETWEEN
#           Poisson (1) and lattice (0). [Hypothesis; at L <= 50 only a few families are "full", so the discrete sum (E4) is
#           the sharper prediction.] Possibly related result (not verified): for squarefree numbers the variance in short
#           intervals grows like h^(1/2) (R. R. Hall).
# Side measurement: w105 unfolded K+ with radius √(2N); the correct radius is 2√N (E|H_ij|² = 1). Here the actual radius is
#   measured, and K+ with both unfoldings.

import sys, json, math, time, pathlib
from math import isqrt
import numpy as np
import mpmath
sys.stdout.reconfigure(encoding='utf-8')
ERG = pathlib.Path(__file__).resolve().parent.parent/'ergebnisse'   # results folder
aus = []   # `aus` = output lines, written to the _output.txt file
def sag(s=''):   # `sag` = say: print a line and record it in `aus`
    print(s, flush=True); aus.append(s)
t0 = time.time()
sag('='*100); sag('w109 — ZAHLVARIANZ Σ²(L): Ordnung ueber lange Strecken   ' + time.strftime('%Y-%m-%d %H:%M')); sag('='*100)
LS = [round(0.25*i, 2) for i in range(1, 21)] + [6, 8, 10, 12, 15, 20, 25, 30, 40, 50]   # `LS` = window lengths L
rng = np.random.default_rng(20260925)
M_JE_BLOCK = 25000   # `M_JE_BLOCK` = random window starts per block
trap = getattr(np, 'trapezoid', None) or np.trapz   # trapezoid rule (name differs between numpy versions)

# `teile` = parts: split the sorted array x into B consecutive blocks (`bloecke`) of nearly equal length
def teile(x, B):
    return [x[i*len(x)//B:(i+1)*len(x)//B] for i in range(B)]
# `sigma2` = number variance estimate for window length L: per block, M_JE_BLOCK random windows, variance of the counts;
#   returns (mean over blocks, standard error)
def sigma2(bloecke, L):
    v = []
    for x in bloecke:
        lo, hi = float(x[0]), float(x[-1]) - L
        if hi <= lo: continue
        t = rng.uniform(lo, hi, M_JE_BLOCK)
        n = np.searchsorted(x, t + L, 'left') - np.searchsorted(x, t, 'left')
        v.append(float(np.var(n)))
    v = np.array(v); return float(v.mean()), float(v.std(ddof=1)/math.sqrt(v.size))
# `messen` = measure: Σ² and its standard error for all window lengths in `LS`
def messen(name, bloecke):
    r = [sigma2(bloecke, L) for L in LS]
    return dict(name=name, sigma2=[round(a, 5) for a, _ in r], se=[round(b, 5) for _, b in r])
# f(μ) = {μ}(1 − {μ}) with {μ} the fractional part
def frac_f(mu):
    fr = mu - np.floor(mu); return fr*(1 - fr)
# GUE number variance: L − 2∫₀ᴸ (L − r)(sin πr/πr)² dr, by the trapezoid rule
def gue(L):
    r = np.linspace(1e-12, L, 40001); y = (np.sin(np.pi*r)/(np.pi*r))**2
    return float(L - 2*trap((L - r)*y, r))
# `res` = result dict: `skript` = script name, `datum` = date, `L` = window lengths, `sets` = Σ² per set, `vorhersage` =
#   predictions (`gue`; `gitter` = lattice; `hardy_littlewood`; `gitter_summe_Q`, `gitter_summe_Q3` = lattice sums of E4),
#   `k_plus_entfaltung` = unfolding comparison for K+ (`radius_mittel` = mean radius, `klein_richtig`, `klein_w105` =
#   small-spacing share with the correct and the w105 unfolding), `kontrollen` = controls (`alle_ok` = all passed), `urteil` =
#   verdict, `familien` = number of families, `anteil_quadrate` = share of squares, `laufzeit_s` = run time in seconds
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), L=LS, sets={}, vorhersage={})

# ================= PREDICTIONS (before any measurement)
res['vorhersage']['poisson'] = LS
res['vorhersage']['gue'] = [round(gue(L), 5) for L in LS]
res['vorhersage']['gitter'] = [round(float(frac_f(np.array(L))), 5) for L in LS]
# Hardy–Littlewood prediction for primes: `sieb` = sieve up to NP; `C2` = twin prime constant (product over p <= NP);
#   `lnN` = log of the mean location of the window [10^9, 10^9 + 10^7); `S[d]` = singular series for even d <= HMAX
NP = 10**7
sieb = bytearray([1]) * (NP + 1); sieb[0] = sieb[1] = 0
for p in range(2, isqrt(NP) + 1):
    if sieb[p]: sieb[p*p::p] = bytearray(len(range(p*p, NP + 1, p)))
C2 = math.exp(sum(math.log(1 - 1/(p - 1)**2) for p in range(3, NP + 1) if sieb[p]))
lnN = math.log(10**9 + 5*10**6)
HMAX = int(max(LS)*lnN) + 2
S = np.zeros(HMAX + 1)
for d in range(2, HMAX + 1, 2):
    dd, prod = d, 1.0
    while dd % 2 == 0: dd //= 2
    q = 3
    while q*q <= dd:
        if dd % q == 0:
            prod *= (q - 1)/(q - 2)
            while dd % q == 0: dd //= q
        q += 2
    if dd > 1: prod *= (dd - 1)/(dd - 2)
    S[d] = 2*C2*prod
# `hl` = Hardy–Littlewood number variance for window length L (formula in E3 above)
def hl(L):
    h = L*lnN; hi = min(int(math.floor(h)), HMAX)
    d = np.arange(1, hi + 1)
    return L*(1 - 1/lnN) + 2/lnN**2 * float(np.sum((h - d)*(S[1:hi + 1] - 1)))
res['vorhersage']['hardy_littlewood'] = [round(hl(L), 5) for L in LS]
sag(f'Vorhersagen: C₂ = {C2:.7f} (Produkt ueber p ≤ 10⁷) · Mittel von S(d), d ≤ {HMAX}: {S[1:].mean():.4f} (muss → 1) · ln N = {lnN:.3f}')
sag('   L      Poisson   GUE      Gitter   Hardy–Littlewood(P)')
for i, L in enumerate(LS):
    if L in (0.5, 1, 2, 5, 10, 20, 50):
        sag(f'   {L:<6} {L:<8} {res["vorhersage"]["gue"][i]:<8.4f} {res["vorhersage"]["gitter"][i]:<8.4f} {res["vorhersage"]["hardy_littlewood"][i]:.4f}'
            f'  (HL/L = {res["vorhersage"]["hardy_littlewood"][i]/L:.3f})')

# powerful window [2^44, 2^46) as in w105, plus the family shares.
# `c1`, `c2`: the counting function of the powerful numbers is c1·√X + c2·∛X + ... with c1 = ζ(3/2)/ζ(3) and c2 = ζ(2/3)/ζ(2)
c1 = float(mpmath.zeta(1.5)/mpmath.zeta(3)); c2 = float(mpmath.zeta(mpmath.mpf(2)/3)/mpmath.zeta(2))
LO, HI = 2**44, 2**46 - 1
BMAX = int(round(HI ** (1/3))) + 2
sfb = bytearray([1]) * (BMAX + 1); sfb[0] = 0   # `sfb` = squarefree flags up to BMAX
for q in range(2, isqrt(BMAX) + 1):
    sfb[q*q::q*q] = bytearray(len(range(q*q, BMAX + 1, q*q)))
# `qteile` = parts for squares (b = 1), `nteile` = parts for the other families (b > 1); `fam[b]` = number of members of family b
#   in the window
qteile, nteile, fam = [], [], {}
for b in range(1, BMAX + 1):
    if not sfb[b]: continue
    b3 = b**3
    if b3 > HI: break
    a0 = isqrt((LO - 1)//b3) + 1; a1 = isqrt(HI//b3)
    if a1 >= a0:
        a = np.arange(a0, a1 + 1, dtype=np.int64); (qteile if b == 1 else nteile).append(a*a*np.int64(b3)); fam[b] = a1 - a0 + 1
# `qu` = the squares, `nq` = the powerful numbers that are not squares, `alle` = all powerful numbers in the window
qu = np.sort(np.concatenate(qteile)); nq = np.sort(np.concatenate(nteile)); alle = np.sort(np.concatenate([qu, nq]))
nQ, nQ3 = alle.size, nq.size
# weights w_b of the families
wQ = np.array([fam[b]/nQ for b in sorted(fam)]); wQ3 = np.array([fam[b]/nQ3 for b in sorted(fam) if b > 1])
# lattice-sum predictions (E4)
vQ = [round(float(np.sum(frac_f(L*wQ))), 5) for L in LS]; vQ3 = [round(float(np.sum(frac_f(L*wQ3))), 5) for L in LS]
res['vorhersage']['gitter_summe_Q'] = vQ; res['vorhersage']['gitter_summe_Q3'] = vQ3
# I in three pieces: [0, μ₀] EXACT (there f = μ − μ², and the integrand ~ μ^(−2/3) is singular: the trapezoid rule from 10⁻⁹ on
#   would give ~500 instead of ~0.3) · [μ₀, 4000] numerical · the rest beyond 4000 with the mean of f = 1/6.
mu0 = 1e-3; I0 = 3*mu0**(1/3) - 0.75*mu0**(4/3)
mu = np.linspace(mu0, 4000, 4_000_001); I_num = I0 + float(trap(frac_f(mu)*mu**(-5/3), mu)) + 3/2*(1/6)*4000**(-2/3)
# The trapezoid rule overestimates slightly just behind μ₀ (the integrand ~ μ^(−2/3) is convex; error ≈ h²/12·(2/3)·μ₀^(−5/3) ≈
#   0.0056), so `I_num` is only a cross-check (within 0.5 %). The value used is the CLOSED FORM via the Fourier series {μ}(1 −
#   {μ}) = Σ_(k≥1) (1 − cos 2πkμ)/(π²k²): I = (2π)^(2/3)/π² · ζ(4/3) · T, T = ∫₀^∞ (1 − cos t) t^(−5/3) dt = −Γ(−2/3)·cos(π/3).
#   This affects ONLY the asymptote, not the lattice sums.
I = float((2*mpmath.pi)**(mpmath.mpf(2)/3)/mpmath.pi**2 * mpmath.zeta(mpmath.mpf(4)/3) * (-mpmath.gamma(-mpmath.mpf(2)/3)*mpmath.cos(mpmath.pi/3)))
assert abs(I_num/I - 1) < 0.005, (I_num, I)
KQ = 4/math.pi**2 * c1**(-2/3) * I; KQ3 = 4/math.pi**2 * (c1 - 1)**(-2/3) * I   # constants K of the asymptote (E5) for Q and Q₃
res['vorhersage'].update(I=I, K_Q=KQ, K_Q3=KQ3, familien=len(fam), anteil_quadrate=round(fam[1]/nQ, 5))
sag(f'Gitter-Ueberlagerung: {len(fam)} Familien im Fenster, Quadrate-Anteil {fam[1]/nQ:.4f}; I = {I:.4f} ⇒ K(Q) = {KQ:.4f}, K(Q₃) = {KQ3:.4f}')
for i, L in enumerate(LS):
    if L in (0.5, 1, 2, 5, 10, 20, 50):
        sag(f'   L = {L:<5}: Q Summe {vQ[i]:8.4f} · Q Asymptote {KQ*L**(2/3):8.4f}  |  Q₃ Summe {vQ3[i]:8.4f} · Q₃ Asymptote {KQ3*L**(2/3):8.4f}')
sag(f'   (Vorhersagen fertig, {time.time()-t0:.0f}s — ab hier wird gemessen)')
sag()

# ================= MEASUREMENT
# K-: uniform random points
u = np.sort(rng.random(200000)) * 200000
res['sets']['K−'] = messen('K−  gleichverteilt, 200 000 Punkte', teile(u, 8))
# K+: GUE, 30 × 1000, with the correct radius 2√N; additionally the w105 unfolding (√(2N)) for the spacing comparison.
# `bl` = unfolded middle halves (one block per matrix), `radien` = measured spectral radii / √N,
#   `s_richtig`, `s_w105` = spacings under the correct and the w105 unfolding
bl, radien, s_richtig, s_w105 = [], [], [], []
# `halbkreis` = semicircle-law counting function (unfolding map for x in [-1, 1])
def halbkreis(x, N): return N * (0.5 + (x*np.sqrt(np.clip(1 - x*x, 0, None)) + np.arcsin(np.clip(x, -1, 1)))/math.pi)
for _ in range(30):
    N = 1000; A = rng.normal(size=(N, N)) + 1j*rng.normal(size=(N, N)); H = (A + A.conj().T)/2
    ev = np.sort(np.linalg.eigvalsh(H)); radien.append(float(max(-ev[0], ev[-1])/math.sqrt(N)))
    mitte = slice(N//4, 3*N//4)   # `mitte` = the middle half of the spectrum
    xr = halbkreis(ev/(2*math.sqrt(N)), N)[mitte]; bl.append(xr); s_richtig.extend(np.diff(xr).tolist())
    xw = halbkreis(ev/math.sqrt(2*N), N)[mitte]; s_w105.extend(np.diff(xw).tolist())
res['sets']['K+'] = messen('K+  GUE 30 × 1000, mittlere Haelfte', bl)
# `klein_anteil` = small share: fraction of spacings s < 0.25 after normalizing the mean spacing to 1
def klein_anteil(s):
    s = np.asarray(s); s = s/s.mean(); return float(np.mean(s < 0.25))
res['k_plus_entfaltung'] = dict(radius_mittel=round(float(np.mean(radien)), 4), klein_richtig=round(klein_anteil(s_richtig), 4),
                                klein_w105=round(klein_anteil(s_w105), 4), n=len(s_richtig))
sag(f'K+ Radius: max|Eigenwert|/√N = {np.mean(radien):.3f} (Halbkreis-Rand 2 ⇒ Radius 2√N). Anteil s < 0,25: richtig entfaltet '
    f'{klein_anteil(s_richtig):.4f} · mit der w105-Entfaltung {klein_anteil(s_w105):.4f} (GUE-Theorie 0,0161; w105 meldete 0,0154)')
# Z: zeta zeros
Z = json.load(open(ERG/'w105_zeta_nullstellen.json', encoding='utf-8'))['imag']
assert len(Z) == 2000 and abs(Z[0] - 14.134725141734693) < 1e-9
# `N_glatt` = smooth counting function of the zeta zeros (Riemann–von Mangoldt main terms), used to unfold
def N_glatt(t): return t/(2*math.pi)*math.log(t/(2*math.pi)) - t/(2*math.pi) + 7/8
xz = np.array([N_glatt(t) for t in Z])
res['sets']['Z'] = messen('Z   Zeta-Nullstellen 1 … 2000', teile(xz, 4))
# W: Wieferich quotients, kernel 15
# fundamental solution (T1, U1) of x² − m·y² = 1 from the continued fraction of √m
def fund(m):
    a0 = isqrt(m); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - m*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (m - Pp*Pp)//Q; a = (a0 + Pp)//Q
        h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
# (T1 + U1·√m)^e modulo M by binary exponentiation; returns (rational part, √m part)
def eps_pow(T1, U1, m, e, M):
    ra, rb, ba, bb = 1 % M, 0, T1 % M, U1 % M
    while e:
        if e & 1: ra, rb = (ra*ba + m*rb*bb) % M, (ra*bb + rb*ba) % M
        ba, bb = (ba*ba + m*bb*bb) % M, (2*ba*bb) % M
        e >>= 1
    return ra, rb
# `ys` = Wieferich quotients y_p = ((B // p) % p) / p, where ε^(p − e) = A_ + B_·√m (mod p²), e = (m/p)
m = 15; T1, U1 = fund(m); ys = []
for p in range(3, 3_000_001):
    if not sieb[p] or m % p == 0 or U1 % p == 0: continue
    e = 1 if pow(m % p, (p-1)//2, p) == 1 else -1
    A_, B_ = eps_pow(T1, U1, m, p - e, p*p); ys.append(((B_ // p) % p) / p)
ys = np.sort(np.array(ys)) * len(ys)
res['sets']['W'] = messen('W   Wieferich-Quotienten, Kern 15, p ≤ 3·10⁶', teile(ys, 8))
sag(f'   (K−, K+, Z, W gemessen, {time.time()-t0:.0f}s)')
# P: primes in [10^9, 10^9 + 10^7), by a segmented sieve with the primes up to √(10^9 + 10^7) from `sieb`
lo = 10**9; LL = 10**7
seg = bytearray([1]) * LL
for p in range(2, isqrt(lo + LL) + 1):
    if sieb[p]:
        st = max(p*p, ((lo + p - 1)//p)*p); seg[st - lo::p] = bytearray(len(range(st - lo, LL, p)))
pz = np.nonzero(np.frombuffer(bytes(seg), dtype=np.uint8))[0].astype(np.float64) + lo   # `pz` = the primes in the window
xp = (pz - lo)/np.log((pz + lo)/2)   # unfolded with the local density 1/ln
res['sets']['P'] = messen('P   Primzahlen in [10⁹, 10⁹ + 10⁷)', teile(xp, 8))
# Q, Q2, Q3: `F` and `Fn` = smooth counting functions of all powerful numbers and of the non-squares among them (used to unfold)
F = lambda X: c1*np.sqrt(X) + c2*np.cbrt(X); Fn = lambda X: (c1 - 1)*np.sqrt(X) + c2*np.cbrt(X)
xQ = F(alle.astype(np.float64)); xQ2 = np.sqrt(qu.astype(np.float64)); xQ3 = Fn(nq.astype(np.float64))
sag(f'   Entfaltung: mittlerer Abstand Q {np.diff(xQ).mean():.5f} · Q₂ {np.diff(xQ2).mean():.5f} · Q₃ {np.diff(xQ3).mean():.5f} (soll 1)')
res['sets']['Q'] = messen('Q   powerful in [2⁴⁴, 2⁴⁶)', teile(xQ, 8))
res['sets']['Q2'] = messen('Q₂  nur Quadrate', teile(xQ2, 8))
res['sets']['Q3'] = messen('Q₃  powerful ohne Quadrate', teile(xQ3, 8))
sag(f'   (alles gemessen, {time.time()-t0:.0f}s)')
sag()

# ================= CONTROLS
# `pk1`, `pk2`, `pk3` = largest deviations from the predictions of PK1, PK2, PK3
s = res['sets']; V = res['vorhersage']
pk1 = max(abs(a/L - 1) for a, L in zip(s['K−']['sigma2'], LS))
pk2 = max(abs(a/g - 1) for a, g, L in zip(s['K+']['sigma2'], V['gue'], LS) if L >= 0.5)
pk3 = max(abs(a - g) for a, g in zip(s['Q2']['sigma2'], V['gitter']))
sag(f'PK1 K− gegen Poisson: groesste rel. Abweichung {pk1:.3f} (Grenze 0,06) ' + ('✅' if pk1 < 0.06 else '🔴'))
sag(f'PK2 K+ gegen GUE-Formel (L ≥ 0,5): groesste rel. Abweichung {pk2:.3f} (Grenze 0,08) ' + ('✅' if pk2 < 0.08 else '🔴'))
sag(f'PK3 Q₂ gegen {{L}}(1 − {{L}}): groesste Abweichung {pk3:.4f} (Grenze 0,01) ' + ('✅' if pk3 < 0.01 else '🔴'))
res['kontrollen'] = dict(pk1=round(pk1, 4), pk2=round(pk2, 4), pk3=round(pk3, 5), alle_ok=bool(pk1 < 0.06 and pk2 < 0.08 and pk3 < 0.01))
sag()
sag('   L     |  K−      K+      GUE    |  Z       |  W       |  P      HL      |  Q       Q-Summe  |  Q₃      Q₃-Summe |  Q₂')
for i, L in enumerate(LS):
    if L in (0.5, 1, 2, 3, 5, 8, 12, 20, 30, 50):
        g = lambda k: s[k]['sigma2'][i]
        sag(f'   {L:<5} | {g("K−"):7.3f} {g("K+"):7.3f} {V["gue"][i]:6.3f} | {g("Z"):7.3f}  | {g("W"):7.3f}  | {g("P"):7.3f} {V["hardy_littlewood"][i]:7.3f} | '
            f'{g("Q"):7.3f}  {V["gitter_summe_Q"][i]:7.3f}  | {g("Q3"):7.3f}  {V["gitter_summe_Q3"][i]:7.3f}  | {g("Q2"):.3f}')
# largest deviation of set k from the prediction v, in units of the standard error
def z_max(k, v):
    return max(abs(a - b)/max(se, 1e-9) for a, b, se in zip(s[k]['sigma2'], v, s[k]['se']))
# largest relative deviation of set k from the prediction v, for window lengths L >= ab
def rel_max(k, v, ab=0.0):
    return max(abs(a/b - 1) for a, b, L in zip(s[k]['sigma2'], v, LS) if L >= ab and b > 0)
# `urteil` = verdict: deviations of the measured sets from the predictions, and slopes of log Σ² between L = 10 and L = 50
urteil = {}
urteil['W_poisson'] = dict(rel=round(rel_max('W', LS), 4), z=round(z_max('W', LS), 2))
urteil['P_poisson'] = dict(rel=round(rel_max('P', LS), 4), z=round(z_max('P', LS), 2))
urteil['P_hl'] = dict(rel=round(rel_max('P', V['hardy_littlewood'], 0.5), 4), z=round(z_max('P', V['hardy_littlewood']), 2))
urteil['Q_gittersumme'] = dict(rel=round(rel_max('Q', V['gitter_summe_Q'], 0.5), 4), z=round(z_max('Q', V['gitter_summe_Q']), 2))
urteil['Q3_gittersumme'] = dict(rel=round(rel_max('Q3', V['gitter_summe_Q3'], 0.5), 4), z=round(z_max('Q3', V['gitter_summe_Q3']), 2))
i10, i50 = LS.index(10), LS.index(50)
for k in ('Q', 'Q3', 'K−', 'K+', 'P'):
    a, b = s[k]['sigma2'][i10], s[k]['sigma2'][i50]
    urteil[f'steigung_{k}'] = round(math.log(b/a)/math.log(5), 4)
res['urteil'] = urteil
sag()
sag(f'E2  W gegen Poisson: groesste rel. Abweichung {urteil["W_poisson"]["rel"]:.3f}, groesstes |z| {urteil["W_poisson"]["z"]:.1f}')
sag(f'E3  P gegen Poisson: rel. {urteil["P_poisson"]["rel"]:.3f} (|z| {urteil["P_poisson"]["z"]:.0f}) · gegen Hardy–Littlewood: rel. {urteil["P_hl"]["rel"]:.3f} (|z| {urteil["P_hl"]["z"]:.1f})')
sag(f'E4  Q gegen Gitter-Summe: rel. {urteil["Q_gittersumme"]["rel"]:.3f} (|z| {urteil["Q_gittersumme"]["z"]:.1f}) · '
    f'Q₃ gegen Gitter-Summe: rel. {urteil["Q3_gittersumme"]["rel"]:.3f} (|z| {urteil["Q3_gittersumme"]["z"]:.1f})')
sag(f'E5  Steigung log Σ² zwischen L = 10 und 50:  Q {urteil["steigung_Q"]:.3f} · Q₃ {urteil["steigung_Q3"]:.3f} · (Poisson K− {urteil["steigung_K−"]:.3f} · '
    f'GUE K+ {urteil["steigung_K+"]:.3f} · Primzahlen {urteil["steigung_P"]:.3f}) — Hypothese 2/3 = 0,667')
res['laufzeit_s'] = round(time.time() - t0)
(ERG/'w109_zahlvarianz_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w109_zahlvarianz_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
sag(f'Laufzeit {time.time()-t0:.0f}s')
