# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# zahlen_bauen.py
# Builds series_source/gemeinsam/zahlen.tex (`gemeinsam` = shared) from the result files. In the LaTeX text, `\Z{key}` stands
#   for a computed number (for example the key `bewiesen` instead of the digits 804). Each number has exactly ONE entry here:
#   source file(s), computation rule, format. zahlen.tex carries for each number a comment with the source file and its sha256,
#   so every printed number can be traced back to its file.
# Formats: 'ganz' (integer) · 'runden4' (rounded to 4 decimals) · 'abrunden4' (rounded DOWN, for LOWER BOUNDS: a rounded-up
#   lower bound would be a false statement) · 'abrunden_wiss2' (scientific notation with 2 digits, rounded down; e.g. 8.41e7 ->
#   8.4*10^7); the other formats are defined in `fmt` below.
# Two kinds of source (`quellen` = sources): 'file.json' is a JSON file whose parsed object is passed to the computation rule;
#   T('file', pattern) is a TEXT ANCHOR for scripts that write only a text output: the pattern (with ONE group) must match
#   EXACTLY ONCE in the file, otherwise the build aborts, so two similar-looking lines can never silently deliver the wrong one.
#   TA('file', pattern) returns ALL matches (at least one).
# EVERY source is checked against MANIFEST_sha256.txt before it is read (hash AND size); if it is missing there or differs,
#   the build aborts. Thus every printed number hangs on a file that the manifest vouches for.
# Number format in the text outputs: digits, optional thousands separator (a space), decimal point, exponent. A comma is
#   REJECTED (1,234 is ambiguous).
# The German descriptions (`was`) of the entries are copied into the comments of zahlen.tex and are not translated here.
# Reads: the result files in ergebnisse/ named in the entries, and MANIFEST_sha256.txt. Writes:
#   series_source/gemeinsam/zahlen.tex.
# Usage: python zahlen_bauen.py   (no arguments; works from any directory).
# Controls: positive controls at start (formats, number parser, anchors with 1 / 0 / 2 matches, hash check for missing /
#   differing / matching file), consistency checks inside the entry functions, and the theorem checks `PRUEF_SAETZE` (statements
#   in the text without a printed number that still depend on data). Any failure aborts the build.
import hashlib, json, math, pathlib, re, sys, time
from collections import namedtuple
sys.stdout.reconfigure(encoding='utf-8')
# `W` = folder of this tool, `S` = series_source, `DATEN` = root of the data folder, `ERG` = result files, `MANIFEST` = list
#   of hashes, `ZIEL` = output file.
W = pathlib.Path(__file__).resolve().parent; S = W.parent
DATEN = S.parent
ERG = DATEN/'ergebnisse'
MANIFEST = DATEN/'MANIFEST_sha256.txt'
ZIEL = S/'gemeinsam'/'zahlen.tex'
# Text anchors: `datei` = file, `muster` = regular-expression pattern. `T` = exactly one match; `TA` = all matches.
T = namedtuple('T', 'datei muster')
TA = namedtuple('TA', 'datei muster')   # all matches (at least one) - for checks over table rows

# Names of result files: `W<label>` is a file written by the script w<label>; the trailing comments say what the file is used
#   for.
W92, W113, W115 = 'w92_karte_811_daten.json', 'w113_familien_weit_B1000000_result.json', 'w115_walker_familien_D1000000_result.json'
W110, W112, W139 = 'w110_paare_bis_1e21_result.json', 'w112_quadrat_familien_result.json', 'w139_tab1_familien_result.json'
W111, W107, W103 = 'w111_paar_familien_result.json', 'w107_paar_korrektur_result.json', 'w103_zweierpotenzen_result.json'
W107OUT = 'w107_paar_korrektur_output.txt'
TUER, M9K = 'auswertung_tueren_2026-09-07_output.txt', 'm9_korrelation_result.json'
W116, W116K, W116P = 'w116_w12_nachbau_result.json', 'w116_w12_nachbau_kreuz_output.txt', 'w116_w12_nachbau_paare_output.txt'
W159 = 'w159_alekseyev_klasse3_1e8_result.json'   # comparison with Alekseyev's table, class 3 up to 10^8
W140 = 'w140_tab_reinhart_result.json'
W29 = 'w29_marge_result.json'; W29T = 'w29_out.txt'   # extent of the computation without floating point in Thm 3.1
# `W153`: kernels with T_1 + 1 = square; the building blocks there have at least two prime divisors (`rem:zweiprimitive`)
W153 = 'w153_t1plus1_quadrat_bloecke_result.json'
# `W152`: parity of T_1 COMPUTED; the labels `(T_1 gerade)` in the output of w11 are wrong for 117 477 414 815
W152 = 'w152_paritaet_t1_reinhart_kerne_result.json'
W14 = 'w14_zeugen_statistik_2026-09-08_output.txt'
W137, W141 = 'w137_wieferich_79kerne_result.json', 'w141_tab_wieferich_result.json'
W84 = 'w84_ist_es_dasselbe_objekt_result.json'
W138, W133, W134 = 'w138_karte_961_daten.json', 'w133_probelauf_79_kerne_result.json', 'w134_stufe2_rangsieb_result.json'
W20N79, W20N79OUT = 'w20_primteil_P3000000_N79_NP5_result.json', 'w20_primteil_P3000000_N79_NP5_output.txt'
W142 = 'w142_tab_offen_result.json'
W158 = 'w158_teile_wahrscheinlich_prim_result.json'   # two blocks closed via probable-prime algebraic parts
# `W146`: successor of w142; also counts the curves of the overnight run w143 (two routes in the script)
W146 = 'w146_tab_offen_result.json'
W20 = 'w20_primteil_P3000000_N25_NP5_result.json'
W21 = 'w21_kreuzpruefung_rizzi_2026-09-09_output.txt'
W19 = 'w19_teilfaelle_woerterbuch_2026-09-09_output.txt'
W18, W17C, W17D = 'w18_rang_kongruenz_2026-09-09_output.txt', 'w17c_rang_abdeckung_tief_2026-09-08_output.txt', 'w17d_isoliert_NP1_2026-09-09_output.txt'
M9OUT, M9T = 'm9_korrelation_2026-09-11_output.txt', 'm9_trennschaerfe_2026-09-11_output.txt'
W112OUT, W115OUT = 'w112_quadrat_familien_output.txt', 'w115_walker_familien_D1000000_output.txt'
W10 = 'w10_L100000000_U1000000000000000_H2000_P100_T4_result.json'
LADDER = 'ladder_all_families_2026-09-03_output.txt'
L_ALLE = T(LADDER, r'^Familien m<=2000 \(m=7 mod 8, quadratfrei\): (\d+) \|')
L_TOT = T(LADDER, r'\| 2-adisch tot \(T1 ungerade\): (\d+) \|')
L_GEBAUT = T(LADDER, r'\| Leiter gebaut: (\d+) \|')
L_STECKT = T(LADDER, r'\| steckengeblieben: (\d+)\s*$')
L_STUFE2 = T(LADDER, r'Kleinste Stufe-2-Schranke: (\S+) Stellen')
# ---- Part I, Section 5. Theorem 5.1 rests on TWO outputs: the 14 families that got stuck in the first ladder are solved only
#   by the deep search. The bound "more than 8.4*10^7 digits for all 154" holds only with both.
DEEP = 'ladder_stuck_deep_2026-09-03_output.txt'
L_MINM = T(LADDER, r'Kleinste Stufe-2-Schranke: \S+ Stellen.*\(m=(\d+)\)')
L_SIEBEN = T(LADDER, r'^\s+7\s+7\s+\d+\s+\[[^\]]*\]\s+\d+\s+\S+\s+\[[^\]]*\]\s+(\d+\.\d+)\s*$')
L_SCHRANKE = T(LADDER, r'^Familien m<=(\d+) \(')
# `L_ZEILEN` rows: (m, log10 of the digit count after stage 2)
L_ZEILEN = TA(LADDER, r'^\s+(\d+)\s+\d+\s+\d+\s+\[[^\]]*\]\s+\d+\s+\S+\s+\[[^\]]*\]\s+(\d+\.\d+)\s*$')
D_GELOEST = T(DEEP, r'^geloest: (\d+) von \d+')
D_VON = T(DEEP, r'^geloest: \d+ von (\d+)')
D_STELLEN = TA(DEEP, r'T hat ~10\^(\d+\.\d+) Stellen\s+GELOEST')
W11 = 'w11_a135735_tuerme_2026-09-07_output.txt'
# `w11_zeuge(m)` = text anchor for the stage-1 witness of kernel m in the w11 output (`paritaet` = `'gerade'` (even) or
#   `'UNGERADE'` (odd), as printed there).
def w11_zeuge(m, paritaet='gerade'):
    return T(W11, rf'^m = {m} \(T_1 {paritaet}\)[^\n]*\n\s+Stufe 1: k = 1\s+-> Zeugen v_p = 1 \(p <= 3000\): \[(\d+)')
W127 = 'w127_dichte_mitten_result.json'
# ---- Part I, Section 7, Theorem 7.1 (vii) = Prop. 6.3: ladders of the three known A135735 kernels with T1 even. TWO runs: P
#   = 2*10^5 for 4 099 215 and 117 477 414 815, P = 3000 for 39 028 039 587 479. Kernel and bound are read TOGETHER (one match =
#   the m-line up to the next "==>" line), so the order cannot swap anything.
# ---- Part I, Section 8
W58, W22, W122 = 'w58_moduln_und_quadratbaum_result.json', 'w22_yokoi_M20000_result.json', 'w122_design_a_pruefung_result.json'
# `w122_pruef` = consistency check of the w122 result (unit as expected, m' = m, the two families); `w122_fenster` = kernels
#   whose approximate size is inside the window.
def w122_pruef(a):
    k = a['kerne']
    assert all(x['einheit_wie_erwartet'] and x['m_strich'] == x['m'] for x in k), 'w122: Einheit oder m′'
    assert {x['art'] for x in k} == {'t²−2', 't²−1'}, 'w122: Familien'
    return True
def w122_fenster(a):
    return [x for x in a['kerne'] if x['log10_naeherung'] < a['fenster_stellen']]
# `W11P` = w11 output of the P = 2*10^5 run; `TURM_MUSTER` = pattern for one ladder block (kernel, first possible index k,
#   minimum digit count of the middle); `T_TURM_P` / `T_TURM_3` = all such blocks of the two runs.
W11P = 'w11_a135735_tuerme_P200000_2026-09-07_output.txt'
TURM_MUSTER = r'^m = (\d+) \(T_1 gerade\)[^\n]*\n(?:[^\n]*\n)*?\s+==> erster moeglicher Index k >= (\d+)[^\n]*Mitte n hat >= (\S+) Stellen'
T_TURM_P = TA(W11P, TURM_MUSTER)
T_TURM_3 = TA(W11, TURM_MUSTER)
# `turm(p200000, p3000, m, feld)` = value of the ladder block of kernel m, taken from the P = 2*10^5 run for 4 099 215 and 117
#   477 414 815 and from the P = 3000 run otherwise; `feld` 'k' = first possible index, else the digit count.
def turm(p200000, p3000, m, feld):
    d = {x[0]: x for x in p200000}; d3 = {x[0]: x for x in p3000}
    assert set(d) == {4099215, 117477414815}, sorted(d)
    q = d[m] if m in d else d3[m]
    return q[1] if feld == 'k' else q[2]
def turm_min_stellen(p200000, p3000, w152):
    # only kernels with T_1 even can carry a triple (parity from w152); 117 477 414 815 has T_1 odd
    werte = [turm(p200000, p3000, m, 'n') for m in w152['gerade_klasse7']]
    assert sorted(w152['gerade_klasse7']) == [4099215, 39028039587479], w152['gerade_klasse7']
    return min(werte)
# ---- Part I, Section 4
W8 = 'w8_v11_M1000000_B10000_result.json'
W51 = 'w51_ballot_errata_gegenprobe_result.json'
def w8_pruef(a):
    # `w8_pruef` = check of the w8 result. Thm 4.4: the five exclusion reasons add up to all pairs; pairs = (families with T1
    #   even) x squarefree b; the 925 survivors all fail at b^3
    c = a['counts']
    assert sum(c[k] for k in ('no_alpha', 'twoadic', 'k_even', 'lemmaL', 'b3_fail')) == a['pairs'], 'w8: Summe ≠ Paare'
    assert (a['fam'] - a['fam_parity_dead']) * a['n_sqfree_b'] == a['pairs'], 'w8: Paare ≠ Familien × b'
    assert len(a['survivors']) == c['b3_fail'] and all(s[3] == 'b3_fail' for s in a['survivors']) and not a['open'], 'w8: Ueberlebende'
    assert c['nonsquare'] == c['nonsquare_exact'] == c['form_ok_exact'] == c['form_ok_modular_open'] == 0, 'w8: Restfaelle'
    return True
# `w51_pruef` = check of the w51 result: our formula and the corrected Fibonacci version are correct in all cases.
def w51_pruef(a):
    u, f = a['unsere_folge'], a['fibonacci']
    assert u['unsere_formel_korrekt'] == u['faelle'] and f['korrigierte_fassung_korrekt'] == f['faelle'], 'w51'
    return True

def tuerme(alle, tot, gebaut, steckt):
    # `tuerme` (towers) = Thm 5.1: families with T1 even = all minus the 2-adically dead ones; cross-check: = ladders built +
    #   ladders stuck
    assert alle - tot == gebaut + steckt, f'Streifen: {alle} − {tot} ≠ {gebaut} + {steckt}'
    return alle - tot

W9 = 'w9_v31_M100000000_H2000_result.json'
# Literature constant (Reinhart 2024: OEIS A135735 exhaustive up to here) - NOT a data value; in the text it is cited. From it
#   the build computes the method line of Part I, Section 3: an unknown kernel with m | U_1 lies above it, hence n > m^(3/2)
#   (Corollary 2.2).
REINHART_GRENZE = 5.325e13
# Reinhart 2024 claims NO completeness for [1.5*10^12, 5.325*10^13]: the authors say they do not claim it and consider it
#   likely that all were found (there is no independent cross-check). Completeness as a statement (Remark 5.5, proof by computer
#   search, two algorithms up to 10^10, one up to 1.5*10^12) holds only in Reinhart 2023, Acta Arith. 211. Hence TWO bounds:
REINHART_SICHER = 1.5e12

def w9_pruef(a):
    # Part I, Section 3: the counters in w9 must match the points, otherwise the funnel (table of successive exclusions)
    #   prints contradictions
    p = a['points']
    assert a['pts'] == len(p) and not a['open_pts'], 'w9: pts ≠ Punkte oder offene Punkte'
    assert a['dead_par'] == sum(x[4] == 'parity' for x in p) and a['dead_wit'] == sum(x[4] == 'witness' for x in p), 'w9: Paritaet/Zeugen'
    assert a['dead_par'] + a['dead_wit'] == a['pts'] and a['cand_fam'] == len(a['cands']), 'w9: Summen'
    assert len({x[0] for x in p}) <= a['cand_fam'] <= a['fam'], 'w9: Trichter nicht monoton'
    return True

# `pell(m)` = T_1 of the fundamental solution of x^2 - m y^2 = 1 (continued fractions).
def pell(m):
    a0 = math.isqrt(m); assert a0 * a0 != m
    mm, dd, q = 0, 1, a0; h1, h = 1, a0; k1, k = 0, 1
    while h * h - m * k * k != 1:
        mm = dd * q - mm; dd = (m - mm * mm) // dd; q = (a0 + mm) // dd; h1, h = h, q * h + h1; k1, k = k, q * k + k1
    return h

def beispiel(a):
    # `beispiel` (example) = worked example: the first point of w9. T_k exactly from the Pell solution; digit count checked
    #   against w9; the witness divides T_k EXACTLY once.
    m, mp, k, st, art, q = a['points'][0]
    T1 = pell(m); t0, t1 = 1, T1
    for _ in range(k - 1): t0, t1 = t1, 2 * T1 * t1 - t0
    assert len(str(t1)) == st and art == 'witness' and t1 % q == 0 and t1 % (q * q) != 0, 'Beispiel passt nicht zu w9'
    # the text speaks of the point with the smallest middle; this must hold in the data
    assert st < min(x[3] for x in a['points'][1:]), 'Beispiel ist nicht mehr der Punkt mit der kleinsten Mitte'
    return t1

def faktoren(n):
    # `faktoren` = prime factorization as TeX (only for small numbers such as the worked example), with a check by
    #   multiplication
    f, d, r = [], 2, n
    while d * d <= r:
        e = 0
        while r % d == 0: r //= d; e += 1
        if e: f.append((d, e))
        d += 1
    if r > 1: f.append((r, 1))
    assert math.prod(p ** e for p, e in f) == n
    return '\\cdot '.join(f'{p}^{{{e}}}' if e > 1 else str(p) for p, e in f)

def streifen(stufe2, steckt, geloest, von, tief):
    # `streifen` (strip) = Thm 5.1: the bound holds for ALL towers only if every family stuck in the first ladder is solved in
    #   the deep search and has more digits there than the minimum of the first ladder
    assert geloest == von == steckt, f'Tiefensuche: {geloest} von {von} geloest, steckengeblieben waren {steckt}'
    assert min(10 ** x for x in tief) > stufe2, 'eine Familie der Tiefensuche liegt unter der Stufe-2-Schranke'
    return stufe2

def w10_rest(a):
    # `w10_rest` = Thm 5.2: counters agree, no points, no candidates; a period > LMAX forces T1 >= phi^LMAX > 10^H
    assert a['done'] and a['n_kernels'] == a['fam'] == a['i'] and a['pts'] == 0 and a['cand'] == 0 and not a['points'], 'w10'
    assert a['LMAX'] * math.log10((1 + 5 ** 0.5) / 2) > a['H'], 'w10: LMAX zu klein fuer 10^H'
    return a['fam'] - a['unit_big']

def zeta(s, N=2000):
    # Euler-Maclaurin summation (three correction terms), accurate to ~1e-13 for s > 1 - our own second route for the limits
    #   from w127
    return sum(n ** -s for n in range(1, N)) + N ** (1 - s) / (s - 1) + N ** -s / 2 + s * N ** (-s - 1) / 12 - s * (s + 1) * (s + 2) * N ** (-s - 3) / 720

def satz_a(a):
    # `satz_a` = Thm 5.3, scope: share of squares among the powerful n with 4 | n ->
    #   (2*sqrt2+1)*zeta(3)/(2*(sqrt2+1)*zeta(3/2)); double squares / sqrt2.
    # A tolerance of 1e-9 does not work: w127 computes with hand-typed, rounded constants (C = 2.173243 instead of
    #   zeta(3/2)/zeta(3) = 2.1732543...), so it is accurate only to ~1e-6. Hence what is printed is checked: the same four
    #   decimals.
    q = (2 * 2 ** 0.5 + 1) * zeta(3) / (2 * (2 ** 0.5 + 1) * zeta(1.5))
    for exakt, w in ((q, a['grenzwert_quadrat']), (q / 2 ** 0.5, a['grenzwert_doppelquadrat']), (q * (1 + 2 ** -0.5), a['grenzwert_satz_a'])):
        assert abs(exakt - w) < 1e-5 and fmt('runden4', exakt) == fmt('runden4', w), f'w127 {w} ≠ geschlossene Formel {exakt}'
    return True

def css(zeilen):
    # `css` = Thm 5.1, comparison with Chim, Shorey and Sinha (2019), p. 440: t > 10^(3*10^13) for the middle t, m in {15, 23,
    #   31, 39, 47, 55, 87} (literature constant). Our stage-2 digit count 10^x: larger for 15, 31, 87, smaller for 23, 47; 39
    #   and 55 have T1 odd.
    d = {m: x for m, x in zeilen}; grenze = 3e13
    assert all(10 ** d[m] > grenze for m in (15, 31, 87)) and all(10 ** d[m] < grenze for m in (23, 47)), d
    assert 39 not in d and 55 not in d and pell(39) % 2 == 1 and pell(55) % 2 == 1, '39/55'
    return True

def paar_grenze(w115, w111):
    # x0 = n0 + 1, n0 = the 39th term of the list. As the NUMBER we use our own value: w115 generates it (family (19, 35)).
    #   Cross-check: log10 of n0 + 1 = the logarithm of the 39th term recorded in w111 (there without the OEIS value itself),
    #   and it is the 39th.
    n0 = max(int(s) for s in w115['paare_ohne_quadrat_bis_XMAX']); letzt = w111['paare'][-1]
    assert letzt['i'] == 39 and abs(math.log10(n0 + 1) - letzt['lg_n1']) < 1e-5, (n0, letzt)
    return n0 + 1

def mersenne(a):
    # `mersenne` = w103 (part D): 2^n - 1 for n = 1 ... 63 completely factored; powerful only for n = 1 (the number 1). Needed
    #   for the statement that no pair straddles a section boundary.
    D = a['D']; assert [e['n'] for e in D] == list(range(1, len(D) + 1))
    assert all(not e['minus']['powerful'] for e in D if e['n'] >= 2), '2^n − 1 powerful fuer ein n ≥ 2'
    return len(D)

def exponent_mal(M, z, n):
    # `exponent_mal` = M = 10^e exactly; returns e*z/n, which must be an integer (Prop. 6.2: 10^8 => 36 and 24)
    e = len(str(int(M))) - 1; assert int(M) == 10**e and (e * z) % n == 0, (M, z, n); return e * z // n

def reinhart_t1(zeilen):
    # `reinhart_t1`: w11, the four kernels 4099215, 117477414815, 39028039587479, 209991 with the digit count of T1; returns
    #   the smallest digit count minus 1
    ms = sorted(int(m) for m, _ in zeilen)
    assert ms == [209991, 4099215, 117477414815, 39028039587479], ms
    return min(int(s) for _, s in zeilen) - 1

# `WOERTER` (words) = fixed status words and table-row labels for `fmt('wort', ...)`; `WEGE` (routes) = order of the rows of
#   the table; `ZERT_WEGE` = the rows with a certified witness.
WOERTER = {'zeuge': 'witness', 'voll': 'complete factorization', 'offen': 'not settled',
           'frage_ja': 'yes', 'frage_nein': 'no',   # table of Wieferich events, column "Block in Prop. I.6.8?"
           # Prop. 6.10 (I.6.8), table "how the pair was settled": rows by CLASS and route, with no numbers in the wording
           'w_erste': 'certified witness: search by position, small primes', 'w_sieb': 'certified witness: search by position, congruence filter',
           'w_sympy': 'certified witness: elliptic curves and $p - 1$ (SymPy)', 'w_gmpecm': 'certified witness: elliptic curves (GMP-ECM)',
           'w_prim': 'certified witness: $M_d$ itself prime, with a certificate',
           'w_vollbew': 'certified witness: complete factorization, factors proved prime',
           'w_voll': 'complete factorization into probable primes',
           # w158: one part of the algebraic split is a probable prime and divides the block exactly once
           'w_prpzeuge': 'witness that is a probable prime', 'w_offen': 'open',
           # primality proofs w163/w164: the part of the algebraic split proved prime serves as witness, (7, 1449)
           'w_teile': 'certified witness: a factor of the algebraic split, proved prime'}
WEGE = ['w_erste', 'w_sieb', 'w_sympy', 'w_gmpecm', 'w_prim', 'w_vollbew', 'w_teile', 'w_voll', 'w_prpzeuge', 'w_offen']
# `ZERT_WEGE` = all routes whose label starts with "certified witness" (instead of a fixed slice)
ZERT_WEGE = [w for w in WEGE if WOERTER[w].startswith('certified witness')]
def zeile_von(x, w133z):
    # `zeile_von` = (class, route) -> table row. w92 lists (319, 5) under route `erste` although its witness comes from the
    #   complete factorization (no witness <= 3*10^6). The decision is made on the DATA: an old pair of route `erste` without a
    #   witness in the strict first search w133 belongs in the row "complete factorization, factors proved prime" (`w_vollbew`).
    r, kl = x['route'], x['klasse']
    if kl == 'offen': return 'w_offen'
    if kl == 'voll': return 'w_voll'
    if kl == 'prpzeuge': return 'w_prpzeuge'
    if r == 'erste': return 'w_erste' if f"{x['m']}_{x['d']}" in w133z else 'w_vollbew'
    return {'sieb8': 'w_sieb', 'sieb9': 'w_sieb', 'ecm81': 'w_sympy', 'gmpecm': 'w_gmpecm', 'ecm135': 'w_gmpecm', 'prim': 'w_prim', 'voll': 'w_vollbew',
            # routes `teile_*` = routes of the algebraic split: class zert after the primality proofs
            'teile_voll': 'w_vollbew', 'teile_zeuge': 'w_teile'}[r]
def wege(a, w133):
    # `wege` = table rows: for each row the count for the first 25 kernels, the other 54, and all - with a sum check against
    #   the classes
    r = a['raenge']; z = {w: [0, 0] for w in WEGE}
    for x in r: z[zeile_von(x, w133['zeugen'])][0 if x['kern'] == 'alt' else 1] += 1
    # empty rows are dropped (after the primality proofs the rows `w_voll` and `w_prpzeuge` are empty); the row `w_offen`
    #   (open) always stays
    zeilen = [[('wort', w), ('ganz', z[w][0]), ('ganz', z[w][1]), ('ganz', sum(z[w]))] for w in WEGE if sum(z[w]) or w == 'w_offen']
    summe = sum(sum(v) for v in z.values()); zert = sum(sum(z[w]) for w in ZERT_WEGE)
    assert summe == a['n'] and zert == a['zert'] and sum(z['w_voll']) == a['voll'] and sum(z['w_offen']) == a['offen'] \
        and sum(z['w_prpzeuge']) == a.get('prpzeuge', 0), ('Wege-Tabelle', summe, zert)
    # expected count of the row `w_vollbew`: the one pair (319, 5) plus the complete factorizations proved prime by the
    #   primality proofs - the expectation is derived from the data
    erw = 1 + sum(1 for x in r if x['klasse'] == 'zert' and x['route'] in ('voll', 'teile_voll'))
    assert sum(z['w_vollbew']) == erw, ('Zeile Vollzerlegung bewiesen', z['w_vollbew'], erw)
    return zeilen
def bouchard(w107, feld):
    # `bouchard` = Remark `rem:bouchard`: smallest pair n, n+1 of our own list with an odd prime p not dividing n(n+1) and p^2
    #   | n + 2
    def zerl(n):
        f, p = {}, 2
        while p * p <= n:
            while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
            p += 1
        if n > 1: f[n] = f.get(n, 0) + 1
        return f
    for s in w107['paare']:
        n = int(s); f = zerl(n + 2)
        neu = [p for p, e in f.items() if p % 2 and e >= 2 and n % p and (n + 1) % p]
        if neu:
            assert all(e >= 2 for e in zerl(n).values()) and all(e >= 2 for e in zerl(n + 1).values())
            return {'n': n, 'p': neu[0], 'fakt': n + 2}[feld]
    raise AssertionError('kein Beispiel')

def erste_von_zwei(liste):
    # `erste_von_zwei` (first of two) = w20 prints the line `Heuristik ueber DIESE Paare` twice: first for P1 (question
    #   population), then for P2 - exactly two expected
    assert len(liste) == 2, liste; return liste[0]

def gleich(a, b):
    # `gleich` (equal) = two anchors of the same line that must be equal (e.g. "955/955"): otherwise the printed number does
    #   not mean "all"
    assert a == b, f'{a} ≠ {b}'; return a

def einzig(liste):
    # `einzig` (unique) = exactly one entry expected - otherwise "the block" in the text would be ambiguous
    assert len(liste) == 1, f'erwartet genau einen Eintrag, nicht {len(liste)}: {liste}'; return liste[0]

def pruefe_doppelte_schluessel(quelltext):
    # A key that occurs twice in the large dictionary is silently overwritten by the later entry (this once turned `wfkerne`
    #   from 79 into 25 in Part I). This source text is parsed and every dictionary with more than 50 keys is checked for
    #   repeated keys.
    import ast, collections
    fehler = []
    for knoten in ast.walk(ast.parse(quelltext)):
        if isinstance(knoten, ast.Dict) and len(knoten.keys) > 50:
            zaehl = collections.Counter(k.value for k in knoten.keys if isinstance(k, ast.Constant))
            fehler += [f'Schluessel {k!r} steht {n}-mal' for k, n in sorted(zaehl.items()) if n > 1]
    return fehler
_fehler_doppelt = pruefe_doppelte_schluessel(pathlib.Path(__file__).read_text(encoding='utf-8'))
assert not _fehler_doppelt, _fehler_doppelt

def und_liste(v):
    # `und_liste` (and-list) = '399, 663 and 20303' from integers (at least two); numbers with 5 or more digits are grouped in
    #   threes with a thin space, as in 'ganz_gruppiert'
    def g(n):
        t = str(int(n))
        if len(t) < 5: return t
        teile = []
        while t: teile.append(t[-3:]); t = t[:-3]
        return '\\,'.join(reversed(teile))
    s = [g(x) for x in v]; assert len(s) >= 2, s
    return ', '.join(s[:-1]) + ' and ' + s[-1]

# `streifen_sieben` ('strip seven') = log10 of the digit count after stage 2 for m = 7, recomputed from k_1 and the stage-2
#   witnesses; `streifen_eins` ('strip one') = number of families for which only the first step succeeded in the deep search;
#   `zquad_bl` = number of building blocks of the kernels with T_1 + 1 = square (w153 against w138).
def streifen_sieben(z):
    k1, *w2, gedruckt = z[0]
    k2 = int(k1)
    for p in w2: k2 *= int(p)
    wert = math.log10(k2 * math.log10(8 + 3 * math.sqrt(7)) - 0.31)
    assert round(wert, 1) == float(gedruckt), f'gedruckte Stellenzahl {gedruckt} ≠ gerechnet {wert}'
    return wert
def streifen_eins(a, b):
    assert a and all(10 ** float(lg) > float(b) for _, lg in a), 'Tiefensuche: Schranke aus Stufe 1 nicht groesser als die kleinste'
    return len(a)
def zquad_bl(a, b):
    assert a['geprueft'] == {'mind. zwei': a['n_bloecke_quadrat']} and a['uebersprungen'] == 0, 'w153: nicht jeder Baustein geprueft oder Verstoss'
    assert sum(1 for r in b['raenge'] if r['m'] in set(a['kerne_quadrat'])) == a['n_bloecke_quadrat'], 'w153 gegen w138: Anzahl der Bausteine'
    return a['n_bloecke_quadrat']
def alek_klasse3(a):
    # `alek_klasse3` = the hash comparison with Alekseyev's 33 terms holds with class 3 up to 10^6 (PK) AND up to 10^8, with
    #   no additional values
    assert a['hash_wie_alekseyev_33_klasse3_1e6'] and a['hash_wie_alekseyev_33_klasse3_1e8'] and a['zusaetzliche_werte_1e6_1e8'] == 0, \
        'w159: der Abgleich mit Alekseyevs Tabelle haelt mit Klasse 3 bis 10^8 nicht'
    return gleich(a['kern_klasse3_max'], a['kern_klasse7_max'])
def wf_inpop_zu(a, b):
    # `wf_inpop_zu` = the Wieferich events whose block lies in the set of 961 pairs - all decided (status in w141), and their
    #   number agrees with the circles of the building-block map (w138, field `wieferich_p`)
    drin = [z for z in a['zeilen'] if z['frage']]
    assert all(z['status'] != 'offen' for z in drin), ('Wieferich-Ereignis in der 961-Menge nicht entschieden', drin)
    return gleich(len(drin), sum(1 for x in b['raenge'] if 'wieferich_p' in x))
def wf_nullrang(a, b, c, d):
    # `wf_nullrang` = Part III counts y_p = 0 on the 25 kernels (w91/w96), Part I only primes WITH a rank (order of eps mod p
    #   divisible by 4). Two routes: (1) class hits "order = 0 (mod 4)" from w96; (2) the Part I events of the 25 kernels (w137)
    #   as a set (m, p), a subset of the hits of w91. In addition: the same 25 kernels in w91 and w138.
    assert set(c['kerne']) == {r['m'] for r in d['raenge'] if r['kern'] == 'alt'}, 'w91 und w138: nicht dieselben 25 Kerne'
    n4 = sum(k['treffer'] for k in a['b_klassen'] if '0 (mod 4)' in k['klasse'])
    alt = {(x['m'], x['p']) for x in b['ereignisse'] if x['kern'] == 'alt'}
    assert alt <= {(t['m'], t['p']) for t in c['treffer']}, 'Teil-I-Ereignis fehlt unter den Nullstellen von Teil III'
    return gleich(n4, len(alt))
# `zquad_prim` = returns the count `prim_kerne` of w153 after checking that there are no exceptions (prime kernels without T_1
#   + 1 = square).
def zquad_prim(a):
    assert a['prim_kerne_ausnahmen'] == 0, 'w153: prime Kerne ohne T_1 + 1 = Quadrat'
    return a['prim_kerne']

def pruef_811(a):
    # `pruef_811` = sum check: certified + completely factored + open = n - otherwise the split in Prop. 6.10 (a) is wrong
    r = a['routen']; z = r['erste'] + r['sieb8'] + r['sieb9'] + r['ecm81'] + r['gmpecm']
    assert z + r['voll'] + r['prim'] + r['offen'] == a['n'], f"Prop 6.10: {z} + {r['voll']} + {r['prim']} + {r['offen']} ≠ {a['n']}"
    return True

# Part II, Section 6: distance two - class 3 (Prop. `klassedrei`), tower sum (Prop. `turmsumme`), lower bound (Cor.
#   `untergrenze`)
W9C3 = 'w9_v32_c3_M100000000_H2000_result.json'
W13 = 'w13_turmsumme_gesetz_2026-09-07_output.txt'
W15 = 'w15_t1_teil_und_untere_schranke_2026-09-08_output.txt'
W16C3 = 'w16_a135735_klasse3_kern_2026-09-08_output.txt'
# Text outputs that were NOT written in UTF-8 (Windows code page, redirected without PYTHONIOENCODING) - named explicitly
#   (`KODIERUNG` = encoding per file) instead of a silent fallback: every other file must be UTF-8, otherwise the build aborts
#   (w15 contains an em dash as byte 0x97).
KODIERUNG = {W15: 'cp1252'}
def w9c3_pruef(a):
    assert a['CLASS'] == 3 and a['H'] == 2000, 'w9 Klasse 3: falscher Lauf'
    return w9_pruef(a)
def turm_steigung(a):
    # `turm_steigung` = sum of 1/(2 m' log10 eps) over the candidates with T1 even; eps = T1 + sqrt(T1^2 - 1), computed from
    #   log10 T1 (w9 stores 6 digits). Recomputed with the exact eps from the Pell solution for 117 of the 122 kernels: equal to
    #   8 digits.
    s = 0.0
    for m, mstr, _, lg, _, _, par in (c[:7] for c in a['cands']):
        if par != 'T1even': continue
        s += 1 / (2 * mstr * (lg + math.log10(1 + math.sqrt(1 - 10**(-2*lg))) if lg < 150 else lg + math.log10(2)))
    return s
def tuerme_gerade(a): return sum(1 for c in a['cands'] if c[6] == 'T1even')
# `nahe` (close) = returns x if it agrees with y within tol.
def nahe(x, y, tol=1e-4):
    assert abs(x - float(y)) < tol, (x, y); return x
W23C7, W23C3 = 'w23_fastaac_c7_M10000000_result.json', 'w23_fastaac_c3_M10000000_result.json'
W24C7, W24C3 = 'w24_aacanteil_c7_M10000000_result.json', 'w24_aacanteil_c3_M10000000_result.json'
MODELL_P = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
def c_tab(a7, a3):
    # `c_tab` = growth question: c_M per decade M = 10 ... 10^7, both classes (w23; the entry 10^8 only repeats the value at
    #   MMAX = 10^7)
    assert a7['MMAX'] == a3['MMAX'] == 10**7 and a7['klasse'] == 7 and a3['klasse'] == 3
    return [[('ganz', j), ('runden4', a7['c_dekade'][str(10**j)]), ('runden4', a3['c_dekade'][str(10**j)])] for j in range(1, 8)]
def modell_quoten(a):
    # `modell_quoten` (model ratios) = share of the kernels with p | m for which also p | U1 holds, times p (model: 1); w24, p
    #   <= 31
    return [a['p_treffer'][str(p)][1] / a['p_treffer'][str(p)][0] * p for p in MODELL_P]
# `w_von(m)` = (w, m', T1) for kernel m with w = 1/(2 m' log10 eps), computed exactly from the Pell solution.
def w_von(m):
    t1 = pell(m); lg = math.log10(t1) + math.log10(1 + math.sqrt(1 - 1 / t1**2)) if t1 < 10**150 else math.log10(t1) + math.log10(2)
    mstr = math.prod(p for p in range(2, m + 1) if m % p == 0 and all(p % q for q in range(2, math.isqrt(p) + 1)) and (t1*t1 - 1) // m % (p*p))
    return 1 / (2 * mstr * lg), mstr, t1
def einkern(a3):
    # `einkern` (single kernel) = class 3, decade 10^5 -> 10^6: the only kernel with m' = 3 (w23 `Nj`) and its share of the
    #   increase of c_M
    e = a3['Nj']['3']; assert e['anzahl'] == 1; m = e['kleinste'][0]
    w, mstr, t1 = w_von(m); assert mstr == 3 and t1 % 2 == 0 and 10**5 < m <= 10**6, (m, mstr)
    return m, w, w / (a3['c_dekade']['1000000'] - a3['c_dekade']['100000'])
# Part II, Section 5: comparison model (marked heuristic), square pairs against the model, partial sums of 1/ln F
W121 = 'w121_zuwachs_zerlegung_result.json'
W25, W26, W27 = 'w25_gap1_N1000000000000_D100000_result.json', 'w26_gap1_D20000_H2000_result.json', 'w27_tiefpass_P5000000_result.json'
W144 = 'w144_abstand1_stufe2_result.json'
W78 = 'w78_omega_heuristik_result.json'   # `rem:zweiprimitive`: omega(M_d) over the 163 completely factored blocks
# ---- Part III, Section 2: cyclotomic values - w55 (x^2+1 as a negative Pell equation), w61 (Phi_3 family), w62 (table,
#   class)
W55, W61, W62 = 'w55_a2_als_negative_pell_result.json', 'w61_phi3_pell_familie_result.json', 'w62_granville_klasse_result.json'
W129 = 'w129_x4plus1_aac_result.json'   # `rem:phiacht`: own recomputation of the AAC condition for q = 1 (mod 8) up to 10^5
# ---- Part III, Section 1: gap spectrum - w48 (existing data), w125 (head, tie-proof), w49 (ceiling), w59/w63/w64 (counting
#   model), w65 (primitivity)
W48S, W125, W49M = 'w48_spektrumkopf_result.json', 'w125_spektrumkopf_gleichstaende_result.json', 'w49_mechanismen_und_10920_result.json'
# `W35L`: Part III, Section 3 - the printed variant is named explicitly (not the variant P500000_N12)
W35L = 'w35_leiter_P2000000_N25_result.json'
W91L, W96N, W108G = 'w91_wieferich_landschaft_result.json', 'w96_wieferich_nachmessung_result.json', 'w108_gedaechtnis_entlang_p_result.json'
# `schwelle` (threshold) = the entry of `a_schwellen` with threshold t.
def schwelle(a, t):
    return next(s for s in a['a_schwellen'] if abs(s['t'] - t) < 1e-15)
W59, W63, W64, W65 = 'w59_zaehlmodell_result.json', 'w63_d2b_stabilitaet_result.json', 'w64_d2b_deutung_result.json', 'w65_d6_klassenpaare_result.json'
W160 = 'w160_anteil_powerful_haeufige_luecken_result.json'
W161 = 'w161_gleitkomma_genauigkeit_result.json'
W165 = 'w165_teil2_primbeweise_result.json'   # primality proofs of the 19 cofactors of Part II
W162, W163, W164 = ('w162_inventar_wahrscheinliche_primzahlen_result.json', 'w163_pocklington_beweise_result.json',
                    # w162-w164: primality proofs; w161: measured accuracy of the floating-point path in Theorem I.3.1; w160:
                    #   share of powerful numbers among the frequent gap values
                    'w164_ecpp_pruefer_result.json')
def ls_tab(a, b):
    # `ls_tab` = per bound: limit, powerful numbers, gap values, share powerful (%), most frequent value, all powerful up to
    #   rank, first non-powerful
    z = []
    for s in sorted(b['ergebnis'], key=int):
        e, f = b['ergebnis'][s], a['ergebnis'][s]
        gleich(e['n_powerful'], f['n_powerful']); gleich(e['n_werte'], f['n_werte']); gleich(tuple(e['top1']), tuple(f['top1']))
        z.append([('$zehnerpotenz', int(s)), ('ganz_gruppiert', e['n_powerful']), ('ganz_gruppiert', e['n_werte']),
                  ('runden2', 100 * f['anteil_pw_gesamt']), ('ganz_gruppiert', e['top1'][0]), ('ganz', e['alle_powerful_bis_rang']),
                  ('ganz_gruppiert', e['erster_nicht_pw']['wert']),
                  ('bruch', (e['schwellen']['200']['davon_powerful'], e['schwellen']['200']['gruppengroesse']))])
    return z
W109V = 'w109_zahlvarianz_result.json'   # number variance; statement moved to Section 1
W148 = 'w148_decke_gleichstaende_result.json'   # ceiling table that is robust against ties (w49 cut the tie off at 89)
def ls_decke(a):
    # `ls_decke` (gap ceiling) = the 10 non-powerful values of the group around rank 200: g, sqfull(g), rank of sqfull(g) as a
    #   range (ties), occurrences, share with gcd = sqfull
    z = []
    for e in a['g200']['nicht_powerful']:
        assert e['g'] % e['sqfull'] == 0
        r0, r1 = e['rang_sqfull']
        z.append([('ganz_gruppiert', e['g']), ('ganz_gruppiert', e['sqfull']), ('rangspanne', (r0, r1)), ('ganz', e['occ']),
                  ('runden1', 100 * e['gcd_max'] / e['occ'])])
    return z
def ls_anteil(a, pw):
    xs = [100 * e['gcd_max'] / e['occ'] for e in (a['g20']['werte'] if pw else a['g200']['nicht_powerful'])]
    return f'{min(xs):.1f}', f'{max(xs):.1f}'
def faktor_tex(s):
    # `faktor_tex` = '7^2·13^3·271^2' -> '7^{2}\cdot 13^{3}\cdot 271^{2}'; every prime power is recomputed and the product is
    #   checked (second route)
    teile = s.split('·'); aus, prod = [], 1
    for t in teile:
        p, _, e = t.partition('^'); p, e = int(p), int(e or 1); prod *= p ** e
        aus.append(f'{p}^{{{e}}}' if e > 1 else f'{p}')
    return prod, '\\cdot '.join(aus)
# `zyk_tab` (cyclotomic table) = table rows (d, x, value, factorization) with the product checked against the value.
def zyk_tab(a):
    z = []
    for e in a['tafel']['eintraege']:
        prod, tex = faktor_tex(e['faktorisierung']); gleich(prod, e['wert'])
        z.append([('ganz', e['d']), ('ganz_gruppiert', e['x']), ('ganz_gruppiert', e['wert']), ('$tex', tex)])
    gleich(len(z), a['tafel']['treffer']); return z
def quadratfrei_bis(n):
    return sum(1 for b in range(1, n + 1) if all(b % (p * p) for p in range(2, int(b ** 0.5) + 1)))
def s2(a, w27, *klassen):
    # `s2` = w144: sum check (each of the remaining cases of w27 exactly once, no triple candidates, PK/NK passed), then the
    #   number of cases in the given classes
    assert a['pk_ok'] and a['nk_ok'] and a['faelle'] == len(w27['bleibt_offen']) == sum(a['klassen'].values())
    assert 'TRIPEL_KANDIDAT' not in a['klassen'] and 'kleiner_zeuge' not in a['klassen'], a['klassen']
    return sum(a['klassen'].get(k, 0) for k in klassen)
# `pell_tu(D)` = fundamental solution (T1, U1) of x^2 - D y^2 = 1.
def pell_tu(D):
    a0 = math.isqrt(D); Pp, Q, a = 0, 1, a0; h1, h0, k1, k0 = 1, a0, 0, 1
    while h0*h0 - D*k0*k0 != 1:
        Pp = a*Q - Pp; Q = (D - Pp*Pp)//Q; a = (a0 + Pp)//Q; h1, h0 = h0, a*h0 + h1; k1, k0 = k0, a*k0 + k1
    return h0, k0
def a1_ueber(w26):
    # `a1_ueber` = w26 runs per kernel up to k <= H/log10 eps + 1 - a rung can lie just ABOVE 10^H. Count exactly: T_k from
    #   the Pell solution, for all entries with k*log10 eps >= H - 1 (below that T_k is certainly < 10^H).
    n = 0
    for x in w26['ueberlebt'] + w26['offen']:
        if x[2] < w26['H'] - 1: continue
        D, k = x[0], x[1]; t1, u1 = pell_tu(D); ra, rb, ba, bb, e = 1, 0, t1, u1, k
        while e:
            if e & 1: ra, rb = ra*ba + D*rb*bb, ra*bb + rb*ba
            ba, bb = ba*ba + D*bb*bb, 2*ba*bb; e >>= 1
        if len(str(ra)) > int(w26['H']): n += 1
    return n
def s2min(a, w27):
    s2(a, w27); return min(r['stellen'] for r in a['ergebnisse'] if r['klasse'].startswith('offen'))
def modell_tab(a):
    # `modell_tab` = w110: count against the naive and the corrected comparison function at X = 10^8, 10^10, ..., 10^20
    zeilen = {z['X']: z for z in a['vergleich']}
    return [[('ganz', e), ('ganz', zeilen[str(10**e)]['gemessen']), ('runden1', zeilen[str(10**e)]['korrigiert']), ('runden1', zeilen[str(10**e)]['naiv'])]
            for e in range(8, 21, 2)]
def modell_null(a, w111, feld):
    z = a['vergleich'][-1]; assert z['X'].startswith('3.8878e+21') and z['gemessen'] == len(w111['paare']), z
    return z[feld]
def quad_rest(a):
    # `quad_rest` (square remainder) = tail argument: from L = 1000 on, the exact number of square pairs stays above the
    #   comparison function - every family born up to L_0 yields at least dL/ln F - 1 members on a stretch of length dL; the
    #   rate of born families > the rate of the model (whose derivative stays < c_1^2*C/4)
    z = next(r for r in a['untergrenze']['zeilen'] if r['L'] == 1000); geboren = a['geburten_bis_T']['1000']
    rest = z['NQ'] - geboren
    assert rest > z['korr'], (rest, z['korr']); return dict(z, geboren=geboren, rest=rest)
def quad_ab(a):
    # `quad_ab` = smallest L of the table from which ALL following rows have N_Q* > V (the first crossing at 192.9 does not
    #   hold: L = 200 lies below)
    z = a['untergrenze']['zeilen']
    for i in range(len(z)):
        if all(r['NQ'] > r['korr'] for r in z[i:]):
            assert z[i]['L'] == int(z[i]['L']); return int(z[i]['L'])
    raise AssertionError('N_Q* bleibt in der Tabelle nirgends ueber V')
# `zuwachs_klein` = share (%) of the increase in the decade (10^5, 10^6] carried by families with ln eps <= 100.
def zuwachs_klein(art):
    def f(a):
        d = a['zerlegung'][art][-1]; assert d['dekade'] == '(1e5, 1e6]', d['dekade']
        klein = sum(v for k, v in d['ln_eps'].items() if float(k.split('-')[1]) <= 100)
        return 100 * klein / d['summe']
    return f
def turm_tab(zeilen, w9a, w9b):
    # `turm_tab` = w13: six heights per class, class 7 first, then class 3; counted = tower sum in every row; final values =
    #   pairs of the w9 runs
    assert len(zeilen) == 12, len(zeilen)
    z7, z3 = zeilen[:6], zeilen[6:]
    for z in zeilen: assert z[1] == z[2], z
    assert [z[0] for z in z7] == [z[0] for z in z3], 'w13: Hoehen der Klassen verschieden'
    assert int(z7[-1][1]) == w9a['dead_wit'] and int(z3[-1][1]) == w9b['dead_wit'], 'w13-Endwerte ≠ w9'
    return [[('ganz', int(a[0])), ('ganz', int(a[1])), ('ganz', int(b[1])), ('ganz', int(a[3])), ('ganz', int(b[3]))] for a, b in zip(z7, z3)]

ZAHLEN = {   # `ZAHLEN` = key: (sources `quellen`, function `funktion` of the source values, format, description `was`)
    'nenner':        ((W92,), lambda a: a['n'], 'ganz', 'Raenge im Gitter der Prop. 6.10'),
    'offen':         ((W92,), lambda a: a['routen']['offen'], 'ganz', 'offene Raenge der Prop. 6.10'),
    # `erledigt` (settled) = certified + completely factored, with a sum check against n. A key named `bewiesen` (proved)
    #   would be misleading: 38 of the 804 rest on PROBABLE primes.
    'erledigt':      ((W92,), lambda a: a['n'] - a['routen']['offen'], 'ganz', 'erledigte Raenge der Prop. 6.10 (zertifiziert + vollzerlegt)'),
    'zertifiziert':  ((W92,), lambda a: pruef_811(a) and sum(a['routen'][k] for k in ('erste', 'sieb8', 'sieb9', 'ecm81', 'gmpecm')), 'ganz',
                      'Raenge mit zertifiziertem Zeugen (Prop. 6.10 (a))'),
    'vollzerlegt':   ((W92,), lambda a: a['routen']['voll'] + a['routen']['prim'], 'ganz',
                      'Raenge durch Vollzerlegung in wahrscheinliche Primzahlen (Prop. 6.10 (a))'),
    'csechsquadrat': ((W113,), lambda a: a['teilsummen']['1000000'], 'runden4', 'Beitrag der Familien mit Quadrat, b ≤ 10^6 (Teil II, C6)'),
    'csechswalker':  ((W115,), lambda a: a['summe_1_durch_lnF']['1000000'], 'runden4', 'Beitrag der Familien ohne Quadrat, D ≤ 10^6 (Teil II, C6)'),
    'csechs':        ((W113, W115), lambda a, b: a['teilsummen']['1000000'] + b['summe_1_durch_lnF']['1000000'], 'abrunden4',
                      'Untergrenze liminf N(x)/log x (Teil II, C6) — ABGERUNDET'),
    # Part I, Section 5
    'streifenfam':   ((L_ALLE,), lambda a: a, 'ganz', 'Familien m ≤ 2000, m ≡ 7 (mod 8) quadratfrei (Thm 5.1)'),
    'streifentot':   ((L_TOT,), lambda a: a, 'ganz', 'davon 2-adisch tot, T1 ungerade (Thm 5.1)'),
    'streifentuerme': ((L_ALLE, L_TOT, L_GEBAUT, L_STECKT), tuerme, 'ganz', 'Tuerme = Familien mit T1 gerade (Thm 5.1), mit Gegenprobe'),
    'streifenstellen': ((L_STUFE2, L_STECKT, D_GELOEST, D_VON, D_STELLEN), streifen, 'abrunden_wiss2',
                        'kleinste Stufe-2-Schranke in Stellen (Thm 5.1), gilt fuer alle Tuerme erst mit der Tiefensuche — ABGERUNDET'),
    'streifenminm':  ((L_MINM,), lambda a: a, 'ganz', 'Kern, an dem die kleinste Stufe-2-Schranke liegt (Thm 5.1)'),
    # The output prints "17.6" ROUNDED (true value 17.567); rounding that down would give a false lower bound. Hence the value
    #   is recomputed from k_1, the stage-2 witnesses and eps = 8 + 3 sqrt 7; the printed number is the second route (its
    #   rounding to one decimal must agree).
    'streifensieben': ((TA(LADDER, r'^\s+7\s+7\s+\d+\s+\[[^\]]*\]\s+(\d+)\s+\S+\s+\[(\d+), (\d+), (\d+)\]\s+(\d+\.\d+)\s*$'),), streifen_sieben, 'abrunden1', 'm = 7: log10 der Stellenzahl nach Stufe 2 (Thm 5.1), abgerundet, aus k_1 und den Zeugen der Stufe 2 gerechnet'),
    'streifeneins':  ((TA(DEEP, r'^m=\s*(\d+) Stufe 1 k=\d+: Zeugen \[[^\]]*\] -> naechstes k = \d+\s+-> T hat ~10\^([\d.]+) Stellen'), L_STUFE2), streifen_eins, 'ganz', 'Tiefensuche: Familien, bei denen nur der erste Schritt gelang (Bound aus Stufe 1 groesser als die kleinste Schranke)'),
    'margestich':    ((W29,), lambda a: a['stichprobe_verworfen'], 'ganz', 'w29: Stichprobe verworfener Familien (Kerne bis margemax)'),
    'margemax':      ((W29,), lambda a: int(a['MMAX']), 'zehnerpotenz', 'w29: Kerngrenze der Stichprobe'),
    'margerand':     ((T(W29T, r'Kerne <= (\S+) mit Marge'),), lambda a: int(float(a)), 'wiss_exakt', 'w29: alle verworfenen Familien bis zu dieser Kerngrenze'),
    'margemin':      ((W29,), lambda a: min(x[0] for x in a['rand_bis_3e6']), 'abrunden1', 'w29: kleinste Marge (Dezimalordnungen) am Rand der Schwelle, rigoros nach unten geschaetzt'),
    'streifenschranke': ((L_SCHRANKE,), lambda a: a, 'ganz', 'Kernschranke des Streifens (Thm 5.1)'),
    'streifensteckt': ((L_STECKT,), lambda a: a, 'ganz', 'in der ersten Leiter steckengeblieben, in der Tiefensuche geloest (Thm 5.1)'),
    'glattkerne':    ((W10,), lambda a: a['n_kernels'], 'ganz_gruppiert','glatte Kerne in (10^8, 10^15], Primfaktoren ≤ 100 (Thm 5.2)'),
    'glattkand':     ((W10,), lambda a: a['cand'], 'ganz', 'Kandidaten darunter (Thm 5.2)'),
    'glattunten':    ((W10,), lambda a: a['MLO'], 'zehnerpotenz', 'untere Kerngrenze (Thm 5.2)'),
    'glattoben':     ((W10,), lambda a: a['MHI'], 'zehnerpotenz', 'obere Kerngrenze (Thm 5.2)'),
    'glatthoehe':    ((W10,), lambda a: a['H'], 'ganz', 'Hoehe H (Thm 5.2)'),
    'glattprim':     ((W10,), lambda a: a['PSMOOTH'], 'ganz', 'groesster erlaubter Primfaktor (Thm 5.2)'),
    'glattmin':      ((W10,), lambda a: a['TMIN'], 'ganz', 'Mindestzahl der Primfaktoren (Thm 5.2)'),
    'glattperiode':  ((W10,), lambda a: a['LMAX'], 'ganz', 'Periodenschranke (Thm 5.2)'),
    'glattgross':    ((W10,), lambda a: w10_rest(a) and a['unit_big'], 'ganz_gruppiert', 'Kerne mit Periode > LMAX (Thm 5.2)'),
    'glattrest':     ((W10,), w10_rest, 'ganz_gruppiert', 'uebrige Kerne, ohne Index unter 10^H (Thm 5.2)'),
    # ---- Part I, Section 4: Theorem 4.4 (b-box) and Prop. 4.3 (ballot errata)
    'bboxkern':      ((W8,), lambda a: w8_pruef(a) and a['MMAX'], 'zehnerpotenz', 'Kernschranke der b-Box (Thm 4.4)'),
    'bboxb':         ((W8,), lambda a: a['BMAX'], 'zehnerpotenz', 'Schranke fuer den quadratfreien Teil b (Thm 4.4)'),
    'bboxfam':       ((W8,), lambda a: w8_pruef(a) and a['fam'], 'ganz_gruppiert', 'Familien m ≤ MMAX (Thm 4.4)'),
    'bboxparitaet':  ((W8,), lambda a: a['fam_parity_dead'], 'ganz_gruppiert', 'davon T1 ungerade (Thm 4.4)'),
    'bboxpaare':     ((W8,), lambda a: w8_pruef(a) and a['pairs'], 'ganz_gruppiert', 'Paare (m, b) (Thm 4.4)'),
    'bboxohnealpha': ((W8,), lambda a: a['counts']['no_alpha'], 'ganz_gruppiert', 'ein Rang existiert nicht (Thm 4.4)'),
    'bboxzweiadisch': ((W8,), lambda a: a['counts']['twoadic'], 'ganz_gruppiert', 'Raenge mit verschiedener 2-Bewertung (Thm 4.4)'),
    'bboxkgerade':   ((W8,), lambda a: a['counts']['k_even'], 'ganz_gruppiert', 'k = α(b) gerade (Thm 4.4)'),
    'bboxlemmal':    ((W8,), lambda a: a['counts']['lemmaL'], 'ganz_gruppiert', 'm′ teilt α(b) nicht (Thm 4.4)'),
    'bboxrest':      ((W8,), lambda a: w8_pruef(a) and a['counts']['b3_fail'], 'ganz', 'Rest, alle scheitern an b³ (Thm 4.4)'),
    'ballotfaelle':  ((W51,), lambda a: w51_pruef(a) and a['unsere_folge']['faelle'], 'ganz', 'Faelle unserer Folge (Prop 4.3, Ballot-Errata)'),
    'ballotgedruckt': ((W51,), lambda a: a['unsere_folge']['gedruckte_fassung_korrekt'], 'ganz', 'davon stimmt die gedruckte Fassung (Prop 4.3)'),
    'ballotfib':     ((W51,), lambda a: w51_pruef(a) and a['fibonacci']['faelle'], 'ganz', 'Faelle Fibonacci, korrigierte Fassung stimmt ueberall'),
    # ---- Part I, Section 7: Theorem 7.1 (vii) = Prop. 6.3 (ladders) - lower bounds ROUNDED DOWN (rounding 4.59*10^30 to
    #   4.6*10^30 would overstate a lower bound)
    'turmk1':        ((T_TURM_P, T_TURM_3), lambda a, b: turm(a, b, 4099215, 'k'), 'abrunden_wiss2', 'erster moeglicher Index, m = 4 099 215 (P = 2·10^5)'),
    'turmk2':        ((T_TURM_P, T_TURM_3), lambda a, b: turm(a, b, 117477414815, 'k'), 'abrunden_wiss2', 'erster moeglicher Index, m = 117 477 414 815 (P = 2·10^5)'),
    'turmk3':        ((T_TURM_P, T_TURM_3), lambda a, b: turm(a, b, 39028039587479, 'k'), 'abrunden_wiss2', 'erster moeglicher Index, m = 39 028 039 587 479 (P = 3000)'),
    'turmstellen':   ((T_TURM_P, T_TURM_3, W152), turm_min_stellen, 'abrunden_wiss2', 'kleinste Stellenzahl der Mitte ueber die Tuerme der Kerne mit T1 gerade'),
    # ---- Part I, Section 8: congruence density (w58), Yokoi check (w22), Design A (w122)
    'kongruenzbis':  ((W58,), lambda a: a['dichten'][-1][0], 'zehnerpotenz', 'Primzahlschranke der Dichte (N2)'),
    'kongruenzdichte': ((W58,), lambda a: a['dichten'][-1][1], 'runden_wiss3', 'Dichte zulaessiger Mitten-Restklassen bis dort (N2) — GERUNDET, nur hinter „about" drucken'),
    'yokoimax':      ((W22,), lambda a: a['MMAX'], 'ganz', 'Kernschranke der Yokoi-Pruefung (D12)'),
    'yokoipaare':    ((W22,), lambda a: a['fehl'] == 0 and a['treffer'] == a['paare'] and a['paare'], 'ganz', 'Paare (m, p), alle mit p³ | T1 ∓ 1 (D12)'),
    'designkerne':   ((W122,), lambda a: w122_pruef(a) and len(a['kerne']), 'ganz', 'Kerne t²−2 / t²−1 (Design A)'),
    'designtmax':    ((W122,), lambda a: max(x['t'] for x in a['kerne']), 'ganz', 'groesstes t (Design A)'),
    'designfenster': ((W122,), lambda a: w122_pruef(a) and len(w122_fenster(a)), 'ganz', 'davon mit erster Kandidatin unter 10^fenster (Design A)'),
    'designmaxm':    ((W122,), lambda a: max(x['m'] for x in w122_fenster(a)), 'ganz', 'groesster solcher Kern (Design A)'),
    'designstellen': ((W122,), lambda a: a['fenster_stellen'], 'ganz', 'Fenster in Stellen (Design A)'),
    # ---- Part I, Section 5, Theorem 5.3 (square middles): witnesses that T1 of the three known kernels with T1 even is not
    #   powerful (w11)
    'quadratzeuge1': ((w11_zeuge(4099215),), lambda a: a, 'ganz', 'm = 4 099 215: Primzahl, die T1 genau einmal teilt (Thm 5.3)'),
    'quadratzeuge2': ((w11_zeuge(117477414815),), lambda a: a, 'ganz', 'm = 117 477 414 815: Zeuge fuer T1 (Thm 5.3)'),
    'quadratzeuge3': ((w11_zeuge(39028039587479),), lambda a: a, 'ganz', 'm = 39 028 039 587 479: kleinster Zeuge fuer T1 (Thm 5.3)'),
    # ---- scope of Theorem 5.3: limits from w127, checked against the closed formula with our own zeta
    'mittenquadrat': ((W127,), lambda a: satz_a(a) and a['grenzwert_quadrat'], 'runden4', 'Anteil Quadrate unter powerful n, 4 | n (Grenzwert)'),
    'mittendoppel':  ((W127,), lambda a: satz_a(a) and a['grenzwert_doppelquadrat'], 'runden4', 'Anteil doppelte Quadrate (Grenzwert)'),
    'mittensatza':   ((W127,), lambda a: satz_a(a) and a['grenzwert_satz_a'], 'runden4', 'Summe: Geltungsbereich Satz 5.3 (Grenzwert)'),
    'mittenN':       ((W127,), lambda a: max(int(k) for k in a['ergebnis']), 'zehnerpotenz', 'groesste gemessene Schranke (w127)'),
    'mittengemessen': ((W127,), lambda a: a['ergebnis'][str(max(int(k) for k in a['ergebnis']))]['anteil_satz_a'], 'runden4',
                       'gemessener Anteil bis zur groessten Schranke (w127)'),
    # ---- Part I, Section 3, reduction theorem - all from w9, each with the check `w9_pruef`
    'kerne':         ((W9,), lambda a: w9_pruef(a) and a['fam'], 'ganz_gruppiert', 'quadratfreie m ≡ 7 (mod 8), m ≤ MMAX (Thm 3.1)'),
    'kernschranke':  ((W9,), lambda a: a['MMAX'], 'zehnerpotenz', 'Kernschranke MMAX (Thm 3.1)'),
    'hoehe':         ((W9,), lambda a: a['H'], 'ganz', 'Hoehe H: T_k < 10^H (Thm 3.1)'),
    'zeugenschranke': ((W9,), lambda a: a['PMAX'], 'wiss_exakt', 'Suchgrenze fuer Zeugen-Primzahlen (Thm 3.1)'),
    'filterfamilien': ((W9,), lambda a: w9_pruef(a) and a['cand_fam'], 'ganz', 'Familien nach dem Verwerfungskriterium (Thm 3.1)'),
    'punktfamilien': ((W9,), lambda a: w9_pruef(a) and len({x[0] for x in a['points']}), 'ganz', 'Familien mit mindestens einem Punkt (Thm 3.1)'),
    'punkte':        ((W9,), lambda a: w9_pruef(a) and a['pts'], 'ganz', 'Gitterpunkte (m, k) mit T_k < 10^H (Thm 3.1)'),
    'paritaetspunkte': ((W9,), lambda a: w9_pruef(a) and a['dead_par'], 'ganz', 'davon T_k ungerade (Thm 3.1)'),
    'zeugenpunkte':  ((W9,), lambda a: w9_pruef(a) and a['dead_wit'], 'ganz', 'davon mit Zeugen-Primzahl (Thm 3.1)'),
    'restpunkte':    ((W9,), lambda a: w9_pruef(a) and len(a['open_pts']), 'ganz', 'Punkte ohne Ausschluss (Thm 3.1, Tab. 1)'),
    'beispielm':     ((W9,), lambda a: a['points'][0][0], 'ganz', 'Rechenbeispiel: Kern (A39)'),
    'beispielk':     ((W9,), lambda a: a['points'][0][2], 'ganz', 'Rechenbeispiel: Index (A39)'),
    'beispielzeuge': ((W9,), lambda a: a['points'][0][5], 'ganz', 'Rechenbeispiel: Zeuge (A39)'),
    'beispielwert':  ((W9,), beispiel, 'ganz_gruppiert', 'Rechenbeispiel: T_k exakt aus der Pell-Loesung, gegen w9 geprueft (A39)'),
    'beispielfaktoren': ((W9,), lambda a: faktoren(beispiel(a)), 'tex', 'Rechenbeispiel: Zerlegung (A39)'),
    'ohnetripel':    ((), lambda: REINHART_GRENZE ** 1.5, 'abrunden_wiss3', 'Methodenzeile (A6): (Reinhart-Grenze)^(3/2), ABGERUNDET — keine Datei'),
    'ohnetripelsicher': ((), lambda: REINHART_SICHER ** 1.5, 'abrunden_wiss3', '(Reinhart-2023-Grenze 1.5e12)^(3/2), ABGERUNDET — die Schranke OHNE Annahme'),
    # Part II, Section 2: Pell families, Table 1 from w139 (each row there is cross-checked three ways)
    'famtab':        ((W139,), lambda a: [[('paar', (z['b'], z['d'])), ('ganz', z['periode']), ('runden2', z['log10_n1']), ('runden2', z['log10_F']),
                                           ('ganz', z['anzahl'])] for z in a['zeilen']], 'tabzeilen', 'Tab. II.1: Familie, Periode, log10 n1, log10 F, Paare'),
    'famzahl':       ((W139,), lambda a: a['familien'], 'ganz', 'Familien der Paare bis zur Grenze (Tab. II.1)'),
    'famtabzwei':    ((W139,), lambda a: next(z['anzahl'] for z in a['zeilen'] if (z['b'], z['d']) == (2, 1)), 'ganz', 'Paare der Familie (2, 1)'),
    'fampaare':      ((W139,), lambda a: a['paare'], 'ganz', 'Paare bis zur Grenze (Tab. II.1)'),
    'famquadrat':    ((W139,), lambda a: a['mit_quadrat'], 'ganz', 'davon mit einem Quadrat'),
    'famohne':       ((W139,), lambda a: a['ohne_quadrat'], 'ganz', 'davon ohne Quadrat (Walker Typ II)'),
    # Rounding the bound down to "3.88*10^21" would be WRONG as the end of a counting range: the 39th pair has n + 1 =
    #   3.8878*10^21 > 3.88*10^21, and up to 3.88*10^21 there are only 38. Hence x0 = n0 + 1 is used exactly (n0 = our own value
    #   from w115) and printed with "≈".
    'paargrenze':    ((W115, W111), lambda a, b: paar_grenze(a, b), 'runden_wiss3', 'x0 = n0 + 1 (39. Paar), nur mit „≈" drucken'),
    'paarletzt':     ((W115, W111), lambda a, b: paar_grenze(a, b) - 1, 'ganz_gruppiert', 'n0: das 39. Paar, von w115 selbst erzeugt'),
    'periodepk':     ((T(W112OUT, r'^PK1 m_b direkt .*?, (\d+) quadratfreie b ≤ \d+: ✅'),), lambda a: a, 'ganz_gruppiert',
                      'Periodenformel (b,1) gegen direkte Iteration: Zahl der b (w112 PK1)'),
    'periodepkmax':  ((W112,), lambda a: a['bmax'], 'ganz_gruppiert', 'Schranke fuer b in w112 PK1'),
    'walkerpk':      ((T(W115OUT, r'^PK1 Restklassen-Formel gegen Iteration mod b·d \(D ≤ \d+\): (\d+)/\d+ ✅'),
                       T(W115OUT, r'^PK1 Restklassen-Formel gegen Iteration mod b·d \(D ≤ \d+\): \d+/(\d+) ✅')),
                      gleich, 'ganz', 'Typ II: Formel gegen Iteration, alle gleich (w115 PK1)'),
    'walkerpkmax':   ((T(W115OUT, r'^PK1 Restklassen-Formel gegen Iteration mod b·d \(D ≤ (\d+)\)'),), lambda a: a, 'ganz', 'Schranke D = bd in w115 PK1'),
    # Part II, Section 3: enumeration and completeness within the families
    'paare52':       ((W107, T(W107OUT, r'^\(2\) Paare bis 2\^\d+: (\d+)')), lambda a, b: gleich(len(a['paare']), b), 'ganz',
                      'Paare mit n + 1 ≤ 2^52, eigene vollstaendige Aufzaehlung (w107)'),
    'paare52exp':    ((T(W107OUT, r'^\(2\) Paare bis 2\^(\d+): \d+'),), lambda a: a, 'ganz', 'Exponent der Aufzaehlungsgrenze (w107)'),
    'mersennemax':   ((W103,), mersenne, 'ganz', '2^n − 1 nicht powerful fuer 2 ≤ n ≤ dieser Zahl (w103, Vollzerlegung)'),
    'famquadratfam': ((W139,), lambda a: sum(1 for z in a['zeilen'] if z['quadrat']), 'ganz', 'Familien mit Quadrat unter den Paaren bis x0'),
    'quadratbmax':   ((W112,), lambda a: a['bmax'], 'zehnerpotenz', 'Schranke fuer den quadratfreien Teil, Familien mit Quadrat (C5)'),
    'walkerdmax':    ((W115,), lambda a: a['DMAX'], 'zehnerpotenz', 'Schranke fuer D = bd, Familien ohne Quadrat (C5, C6)'),
    # Part I, Section 6, Prop. 6.1: exponents at the rank of apparition; independence tested (m9), discrimination measured
    'expfam':        ((M9K, T(TUER, r'Familien (\d+), \(m,p\)-Paare')), lambda a, b: gleich(a['fams'], b), 'ganz', 'Kerne der Prop. 6.1'),
    'expmmax':       ((T(TUER, r'lebende Familien m <= (\S+) \(m = 7 mod 8\)'),), lambda a: a, 'zehnerpotenz', 'Schranke fuer m (Prop. 6.1)'),
    'exppmax':       ((T(TUER, r'\(m = 7 mod 8\), Primzahlen 3\.\.(\d+)'),), lambda a: a, 'ganz', 'Schranke fuer p (Prop. 6.1)'),
    'exppaare':      ((M9K, T(TUER, r'Apparitionsrang: (\d+)')), lambda a, b: gleich(a['paare'], b), 'ganz_gruppiert', 'Paare (m, p) mit Rang (Prop. 6.1)'),
    'expwzwei':      ((M9K, T(TUER, r'Apparition\): (\d+)   Heuristik')), lambda a, b: gleich(a['w2'], b), 'ganz', 'v_p ≥ 2 beobachtet'),
    'expwdrei':      ((M9K, T(TUER, r'v_p\(T_alpha\) >= 3: (\d+)')), lambda a, b: gleich(a['w3'], b), 'ganz', 'v_p ≥ 3 beobachtet'),
    'exphzwei':      ((T(TUER, r'Heuristik sum 1/p: (\d+)'),), lambda a: a, 'ganz', 'Summe 1/p (Modell fuer v_p ≥ 2)'),
    'exphdrei':      ((T(TUER, r'Heuristik sum 1/p\^2: (\S+)'),), lambda a: a, 'runden1', 'Summe 1/p² (Modell fuer v_p ≥ 3)'),
    'korrs':         ((M9K,), lambda a: a['S_obs'], 'runden1', 'Statistik S beobachtet'),
    'korrmu':        ((M9K,), lambda a: a['S_null_mu'], 'runden1', 'S unter Unabhaengigkeit, Mittel'),
    'korrsd':        ((M9K,), lambda a: a['S_null_sd'], 'runden1', 'S unter Unabhaengigkeit, Standardabweichung'),
    'korrreps':      ((M9K,), lambda a: a['reps'], 'ganz', 'Simulationen'),
    'korrz':         ((M9K,), lambda a: a['z'], 'runden2', 'z-Wert'),
    'korrp':         ((M9K,), lambda a: a['p_emp'], 'runden2', 'empirisches p, einseitig'),
    'korrvert':      ((M9K,), lambda a: a['treffer_je_kern'], 'liste_ganz', 'Kerne mit 0, 1, 2, … Ereignissen'),
    'korrmaxp':      ((T(M9OUT, r'Abweichungen \(p, beobachtet, erwartet, z\):\r?\n\s+p =\s+(\d+)'),), lambda a: a, 'ganz', 'Primzahl mit groesster Abweichung'),
    'korrmaxz':      ((T(M9OUT, r'Abweichungen \(p, beobachtet, erwartet, z\):\r?\n\s+p =\s+\d+\s+beob\s+\d+\s+erw\s+\S+\s+z = \+(\S+)'),), lambda a: a,
                      'runden1', 'ihr z-Wert'),
    'trennvoll':     ((T(M9T, r'^\s+1\.00\s+\+(\S+)'),), lambda a: a, 'ganz_gerundet', 'mittleres z bei voller Kopplung'),
    'trennzwei':     ((T(M9T, r'^\s+0\.20\s+\+\S+\s+(\d+) %'),), lambda a: a, 'ganz', 'Anteil erkannter Laeufe bei rho = 0.2, Prozent'),
    'trennfuenf':    ((T(M9T, r'^\s+0\.05\s+\+\S+\s+(\d+) %'),), lambda a: a, 'ganz', 'Anteil erkannter Laeufe bei rho = 0.05, Prozent'),
    'trennlaeufe':   ((T(M9T, r'\((\d+) Laeufe je rho\)'),), lambda a: a, 'ganz', 'Laeufe je Kopplungsstaerke'),
    # Part I, Section 6, Prop. 6.2: list of the pairs at distance 2 with middle = 0 (mod 4); comparison with Alekseyev;
    #   completeness
    'alekterme':     ((W116, T(W116K, r'^Alekseyev-Terme: (\d+) \|')), lambda a, b: gleich(a['paare_bis_67_stellen'], b), 'ganz',
                      'Terme in Alekseyevs Tabelle = unsere Paare bis zu ihrem groessten Term (w116, Hash)'),
    'alekstellen':   ((T(W116K, r'groesster hat (\d+) Stellen'),), lambda a: a, 'ganz', 'Stellen des groessten Terms'),
    'klassedreimax': ((W159,), alek_klasse3, 'zehnerpotenz', 'Kernschranke im Alekseyev-Vergleich, beide Klassen (w159; bis 01.10. Klasse 3 nur 10^6 aus w116)'),
    'alekneun':      ((T(W116P, r'^Terme mit Mitte = 0 mod 4: (\d+) von'),), lambda a: a, 'ganz', 'davon mit Mitte ≡ 0 (mod 4)'),
    'vollexpo':      ((W9,), lambda a: exponent_mal(a['MMAX'], 9, 2), 'ganz', 'n < 10^e/27 mit e = 9/2 · log10 MMAX (Kor. 2.3 (i))'),
    'vollmexpo':     ((W9,), lambda a: exponent_mal(a['MMAX'], 3, 1), 'ganz', 'm < 10^e/9 mit e = 3 · log10 MMAX (Kor. 2.3 (ii))'),
    'reinhartt1':    ((TA(W11, r'^m = (\d+) \(T_1 [^\n]*\((\d+) Stellen\), Konvergenten'),), reinhart_t1, 'ganz',
                      'die vier bekannten A135735-Kerne ≡ 7 (mod 8): T1 ≥ 10^dieser Zahl'),
    # Part I, Section 6, Prop. 6.3: ladders of the known A135735 kernels = 7 (mod 8), table from w140
    'reinharttab':   ((W140, W152), lambda a, b: [[('ganz_gruppiert', z['m']), ('ganz_gruppiert', z['t1_stellen']), ('$wiss_exakt', z['P']),
                                            ('stufen', z['stufen_zeugen']), ('$abrunden_wiss2', z['k_min']), ('$abrunden_wiss2', z['mitte_stellen'])]
                                           for z in a['zeilen'] if z['m'] in b['gerade_klasse7']], 'tabzeilen', 'Tab. Prop. 6.3: Kern (T1 gerade), Stellen T1, P, Zeugen je Stufe, k ≥, Stellen der Mitte ≥'),
    'reinhartungerade': ((W152,), lambda a: sorted(a['ungerade_klasse7'])[0], 'ganz_gruppiert', 'der kleinere der beiden Kerne ≡ 7 (mod 8) mit T1 ungerade (209 991)'),
    'reinhartungerade2': ((W152,), lambda a: sorted(a['ungerade_klasse7'])[1], 'ganz_gruppiert', 'der groessere der beiden Kerne ≡ 7 (mod 8) mit T1 ungerade (117 477 414 815)'),
    'reinhartboxk':  ((W9,), lambda a: [p[2] for p in sorted(a['points'], key=lambda p: p[2]) if p[0] == 4099215], 'liste_ganz',
                      'Indizes k der Gitterpunkte (4 099 215, k) in der Box von Thm 3.1'),
    'reinhartboxz':  ((W9,), lambda a: [p[5] for p in sorted(a['points'], key=lambda p: p[2]) if p[0] == 4099215], 'liste_ganz',
                      'ihre Zeugen, in derselben Reihenfolge'),
    'reinhartungeradest': ((W152,), lambda a: [z['T1_stellen'] for z in a['kerne'] if z['m'] == sorted(a['ungerade_klasse7'])[0]][0], 'ganz', 'Stellen seines T1'),
    'reinhartungeradest2': ((W152,), lambda a: [z['T1_stellen'] for z in a['kerne'] if z['m'] == sorted(a['ungerade_klasse7'])[1]][0], 'ganz_gruppiert', 'Stellen des T1 des zweiten'),
    # Part I, Section 6, Prop. 6.6: smallest witnesses of the 371 points (w14)
    'zeugenalle':    ((W9, T(W14, r'^W14 Zeugen-Statistik: (\d+) tote Punkte')), lambda a, b: gleich(w9_pruef(a) and a['dead_wit'], b), 'ganz',
                      'Punkte mit Zeuge (w14 = w9)'),
    'zeugenprim':    ((T(W14, r'^Zeuge primitiv \(alpha\(p\) = k\): (\d+) von'),), lambda a: a, 'ganz', 'Punkte mit primitivem kleinsten Zeugen'),
    'zeugenmedian':  ((T(W14, r'Median: (\d+)'),), lambda a: a, 'ganz', 'Median des kleinsten Zeugen'),
    'zeugenmax':     ((T(W14, r'groesster kleinster Zeuge: p = (\d+)'),), lambda a: a, 'ganz_gruppiert', 'groesster kleinster Zeuge'),
    **{f'zeugen{p}': ((T(W14, rf'Zeugen-Primzahl p: \[.*?\({p}, (\d+)\)'),), (lambda a: a), 'ganz', f'Punkte mit kleinstem Zeugen {p}')
       for p in (29, 19, 11, 3, 5, 7)},
    'zeugennkp':     ((T(W14, r'Negativ-Kontrolle v_(\d+)\(T_\d+\(7\)\) = \d+ exakt'),), lambda a: a, 'ganz', 'NK: Primzahl'),
    'zeugennkk':     ((T(W14, r'Negativ-Kontrolle v_\d+\(T_(\d+)\(7\)\) = \d+ exakt'),), lambda a: a, 'ganz', 'NK: Index'),
    'zeugennkv':     ((T(W14, r'Negativ-Kontrolle v_\d+\(T_\d+\(7\)\) = (\d+) exakt'),), lambda a: a, 'ganz', 'NK: Exponent'),
    # Part I, Section 6, Remark 6.6': condition modulo 36 against the witnesses (w21)
    'rizzierfuellt': ((T(W21, r'Rizzis Bedingung erfuellt: (\d+) Paare'),), lambda a: a, 'ganz', 'Mitten, die die Bedingung mod 36 erfuellen'),
    'rizzinicht':    ((T(W21, r'nicht erfuellt: (\d+) Paare'), T(W14, rf'Zeugen-Primzahl p: \[.*?\(3, (\d+)\)')), gleich, 'ganz',
                      'Mitten, die sie nicht erfuellen = Punkte mit kleinstem Zeugen 3'),
    # Part I, Section 6, Prop. 6.8: congruence, number of conditions, occupation of the ranks (w18, w17c, w17d)
    'kongrp':        ((T(W18, r'^W18: P = (\d+) \('),), lambda a: a, 'zehnerpotenz', 'Primzahlschranke (a)'),
    'kongrkerne':    ((T(W18, r'-Paaren ueber (\d+) Kerne'),), lambda a: a, 'ganz', 'Kerne (a)'),
    'kongrpaare':    ((T(W18, r' an (\d+) \(p, d\)-Paaren'),), lambda a: a, 'ganz_gruppiert', 'Paare (p, d) (a)'),
    'bedpunkte':     ((W9, T(W18, r'ueber alle (\d+) Gitterpunkte der Klasse 7')), lambda a, b: gleich(w9_pruef(a) and a['pts'], b), 'ganz', 'Punkte (b)'),
    'bedteiler':     ((T(W18, r'Teiler d >= 5 insgesamt (\d+)'),), lambda a: a, 'ganz', 'Teiler d ≥ 5 insgesamt (b)'),
    'bedausweg':     ((T(W18, r'Fluchtweg p \| k/d unmoeglich\): (\d+) \('),), lambda a: a, 'ganz', 'davon ohne Fluchtweg (b)'),
    'bedanteil':     ((T(W18, r'Fluchtweg p \| k/d unmoeglich\): \d+ \((\S+) %\)'),), lambda a: a, 'runden1', 'Anteil in Prozent (b)'),
    'bedmedian':     ((T(W18, r'^ausweglose Bedingungen je Punkt: min \d+, Median (\d+),'),), lambda a: a, 'ganz', 'ohne Fluchtweg je Punkt: Median'),
    'bedmittel':     ((T(W18, r'^ausweglose Bedingungen je Punkt: .*Mittel (\S+), max'),), lambda a: a, 'runden1', '… Mittel'),
    'bedmax':        ((T(W18, r'^ausweglose Bedingungen je Punkt: .*max (\d+)'),), lambda a: a, 'ganz', '… Maximum'),
    'bedallemedian': ((T(W18, r'^Bedingungen je Punkt insgesamt \(Teiler d >= 5\): min \d+, Median (\d+)'),), lambda a: a, 'ganz', 'alle je Punkt: Median'),
    'bedallemax':    ((T(W18, r'^Bedingungen je Punkt insgesamt \(Teiler d >= 5\): .*max (\d+)'),), lambda a: a, 'ganz', 'alle je Punkt: Maximum'),
    'besp':          ((T(W17C, r'^W17c: P = (\d+) \('),), lambda a: a, 'wiss_exakt', 'Primzahlschranke (c)'),
    'beskerne':      ((T(W17C, r'NKERNE = (\d+), NP'),), lambda a: a, 'ganz', 'Kerne (c)'),
    'besnp':         ((T(W17C, r'NP = (\d+) Primzahlen je Rang'),), lambda a: a, 'ganz', 'Primzahlen je Rang (c)'),
    'bespunkte':     ((T(W17C, r'^GESAMT: (\d+) Punkte'),), lambda a: a, 'ganz', 'Punkte (c)'),
    'besraenge':     ((T(W17C, r'^GESAMT: \d+ Punkte, (\d+) Raenge'), T(W17D, r'^GESAMT: \d+ Punkte, (\d+) Raenge')), gleich, 'ganz', 'Raenge (c)'),
    'besbesetzt':    ((T(W17C, r'besetzt (\d+) \('), T(W17D, r'besetzt (\d+) \(')), gleich, 'ganz', 'besetzte Raenge (c)'),
    'beszeuge':      ((T(W17C, r'mit Zeuge (\d+) \('),), lambda a: a, 'ganz', 'mit Zeuge (c)'),
    'beszeugeant':   ((T(W17C, r'mit Zeuge \d+ \((\S+) % der besetzten'),), lambda a: a, 'runden2', 'Anteil (c)'),
    'beswief':       ((T(W17C, r'nur Wieferich (\d+),'),), lambda a: a, 'ganz', 'nur Wieferich (c)'),
    'besleer':       ((T(W17C, r'unbesetzt \(Primzahl > P\) (\d+)'), T(W17D, r'unbesetzt \(Primzahl > P\) (\d+)')), gleich, 'ganz', 'unbesetzt (c)'),
    'beseinsnp':     ((T(W17D, r'NP = (\d+) Primzahlen je Rang'),), lambda a: a, 'ganz', 'Vergleichslauf: Primzahlen je Rang'),
    'beseinszeuge':  ((T(W17D, r'mit Zeuge (\d+) \('),), lambda a: a, 'ganz', 'Vergleichslauf: mit Zeuge'),
    'beseinsant':    ((T(W17D, r'mit Zeuge \d+ \((\S+) % der besetzten'),), lambda a: a, 'runden2', 'Vergleichslauf: Anteil'),
    'beseinswief':   ((T(W17D, r'nur Wieferich (\d+),'),), lambda a: a, 'ganz', 'Vergleichslauf: nur Wieferich'),
    # Part I, Section 6, Prop. 6.10 = I.6.8: the 961 pairs (w138 after w135/w136), table of routes, measurement (c) from w20
    #   on 79 kernels
    'rgpaare':       ((W138,), lambda a: a['n'], 'ganz', 'Paare (m, d) ueber allen Kernen'),
    'rgkerne':       ((W138,), lambda a: a['kerne'], 'ganz', 'Kerne mit Gitterpunkten'),
    'rgerledigt':    ((W138,), lambda a: gleich(a['erledigt'], a['zert'] + a['voll'] + a.get('prpzeuge', 0)), 'ganz', 'erledigt (zertifiziert + vollzerlegt + wahrscheinlich primer Zeuge, w158)'),
    'rgzert':        ((W138,), lambda a: a['zert'], 'ganz', 'mit zertifiziertem Zeugen'),
    'rgvoll':        ((W138,), lambda a: a['voll'], 'ganz', 'durch Vollzerlegung in wahrscheinliche Primzahlen'),
    # primality proofs of the probable primes (w162 inventory, w163 Pocklington, w164 ECPP with our own checker)
    'pbbausteine':   ((W138, W162), lambda a, b: gleich(a['zert_durch_primbeweis'], b['bausteine']), 'ganz', 'Bausteine, die durch die Primbeweise zert wurden'),
    'pbzahlen':      ((W162, W163, W164), lambda a, b, c: gleich(a['nur_wahrscheinlich'], b['bewiesen'] + c['bewiesen']), 'ganz', 'bewiesene, vorher nur wahrscheinliche Primzahlen'),
    'pbpock':        ((W163,), lambda a: gleich(a['bewiesen'], a['n']), 'ganz', 'davon mit Pocklington'),
    'pbecpp':        ((W164,), lambda a: gleich(a['bewiesen'], a['n']), 'ganz', 'davon mit ECPP'),
    'pbpockmax':     ((W163,), lambda a: max(e['stellen'] for e in a['ergebnisse']), 'ganz', 'groesste mit Pocklington (Stellen)'),
    'pbecppmax':     ((W164,), lambda a: max(e['stellen'] for e in a['ergebnisse']), 'ganz', 'groesste mit ECPP (Stellen)'),
    'rgoffen':       ((W138,), lambda a: gleich(a['offen'], a['n'] - a['erledigt']), 'ganz', 'offen'),
    # the two blocks that w158 closes via the algebraic parts - number, kernel, rank, digits of the parts
    'rgprpzeuge':    ((W138,), lambda a: gleich(a['prpzeuge'], a['erledigt'] - a['zert'] - a['voll']), 'ganz', 'mit wahrscheinlich primem Zeugen (w158)'),
    # settled under the condition that probable primes are prime (complete factorization + probable-prime witness)
    'rgbedingt':     ((W138,), lambda a: gleich(a['voll'] + a['prpzeuge'], a['erledigt'] - a['zert']), 'ganz',
                      'erledigt, falls die wahrscheinlichen Primzahlen prim sind (Vollzerlegung + wahrscheinlich primer Zeuge)'),
    'teilevollm':    ((W158,), lambda a: einzig(a['klassen']['beide_prp'])[0], 'ganz', 'w158: Kern des Blocks mit zwei wahrscheinlich primen Teilen'),
    'teilevolld':    ((W158,), lambda a: einzig(a['klassen']['beide_prp'])[1], 'ganz', 'w158: Rang dazu'),
    'teilevolla':    ((W158,), lambda a: next(z['stellen_A'] for z in a['bloecke'] if z['klasse'] == 'beide_prp'), 'ganz', 'w158: Stellen von A'),
    'teilevollb':    ((W158,), lambda a: next(z['stellen_B'] for z in a['bloecke'] if z['klasse'] == 'beide_prp'), 'ganz', 'w158: Stellen von B'),
    'teilezeugem':   ((W158,), lambda a: einzig(a['klassen']['zeuge_prp'])[0], 'ganz', 'w158: Kern des Blocks mit wahrscheinlich primem Zeugen'),
    'teilezeuged':   ((W158,), lambda a: einzig(a['klassen']['zeuge_prp'])[1], 'ganz', 'w158: Rang dazu'),
    'teilezeugest':  ((W158,), lambda a: next(z['stellen_B'] if z['B_prp'] else z['stellen_A'] for z in a['bloecke'] if z['klasse'] == 'zeuge_prp'),
                      'ganz', 'w158: Stellen des wahrscheinlich primen Teils'),
    'rgaltkerne':    ((W138,), lambda a: len({x['m'] for x in a['raenge'] if x['kern'] == 'alt'}), 'ganz', 'die ersten Kerne (meiste Punkte)'),
    'rgneukerne':    ((W138,), lambda a: len({x['m'] for x in a['raenge'] if x['kern'] == 'neu'}), 'ganz', 'die uebrigen Kerne'),
    'rgungerade':    ((W133, W138), lambda a, b: len(set(a['kerne_T1_ungerade']) & {x['m'] for x in b['raenge']}), 'ganz', 'Kerne mit T1 ungerade'),
    'wegetab':       ((W138, W133), wege, 'tabzeilen', 'Tab. I.6.8: wie jedes Paar erledigt wurde, erste / uebrige / alle Kerne'),
    # Quality check: keys that replace numbers typed by hand in the text
    'bouchardn':     ((W107,), lambda a: bouchard(a, 'n'), 'ganz', 'Remark rem:bouchard: kleinstes Paar mit neuer ungerader Primzahl im Quadrat in n + 2'),
    'bouchardp':     ((W107,), lambda a: bouchard(a, 'p'), 'ganz', 'diese Primzahl'),
    'bouchardfakt':  ((W107,), lambda a: faktoren(bouchard(a, 'fakt')), 'tex', 'Zerlegung von n + 2'),
    'pellwief':      ((W84,), lambda a: ', '.join(str(x[0]) for x in a['pk_pell']), 'wortliste_ganz', 'die drei Wieferich-Primzahlen der Pell-Folge (w84 PK)'),
    'reinhartza':    ((W140,), lambda a: next(z['stufen_zeugen'][0][0] for z in a['zeilen'] if z['m'] == 4099215), 'ganz', 'erster Zeuge m = 4 099 215'),
    'reinhartzc':    ((W140,), lambda a: next(z['stufen_zeugen'][0][0] for z in a['zeilen'] if z['m'] == 39028039587479), 'ganz',
                      'erster Zeuge m = 39 028 039 587 479'),
    'offentab':      ((W146, W138), lambda a, b: gleich(len(a['zeilen']), b['offen']) and
                      [[('ganz_gruppiert', z['m']), ('ganz', z['d']), ('ganz_gruppiert', z['stellen'])]
                       + [('oderstrich:ganz_gruppiert', z['kurven'].get(b1)) for b1 in ('50000', '1000000', '3000000')]
                       + [('runden2', z['t35']), ('runden2', z['t40'])] for z in a['zeilen']], 'tabzeilen',
                      'Tab. I.6.8 offen: Kern, Rang, Stellen, Kurven je B1, t35, t40 (w146: alte aus w92, neue aus w135 + w143 mit dem B2 der Frage)'),
    'wfinpop':       ((W138,), lambda a: sum(1 for x in a['raenge'] if 'wieferich_p' in x), 'ganz', 'Bausteine der 961 mit Wieferich-Ereignis'),
    'offenmin':      ((W146,), lambda a: min(z['stellen'] for z in a['zeilen']), 'ganz', 'kleinste Stellenzahl eines offenen Bausteins'),
    'erstesucheP':   ((W20N79,), lambda a: a['P'], 'wiss_exakt', 'Suchgrenze der ersten Suche'),
    'siebober':      ((W134,), lambda a: a['ober'], 'zehnerpotenz', 'obere Grenze des Rang-Siebs'),
    'prtripel':      ((T(W20N79OUT, r'P1 [^\n]*ueber (\d+) Paare \(m, d, p\)'),), lambda a: a, 'ganz_gruppiert', 'Tripel (m, d, p) der Frage-Population'),
    'prvauf1':       ((T(W20N79OUT, r'P1 [^\n]*\r?\n\s+v = 1:\s+(\d+)'),), lambda a: a, 'ganz_gruppiert', 'davon v = 1'),
    'prvauf2':       ((T(W20N79OUT, r'P1 [^\n]*\r?\n[^\n]*\r?\n\s+v = 2:\s+(\d+)'),), lambda a: a, 'ganz', 'davon v = 2'),
    'prerwartet':    ((TA(W20N79OUT, r'Heuristik ueber DIESE Paare: v >= 2 erwartet (\S+), beobachtet'),), lambda a: erste_von_zwei(a), 'runden2',
                      'Heuristik: erwartete v ≥ 2 in der Frage-Population (erste der zwei Zeilen = P1)'),
    # Part I, Section 6, Prop. 6.9: middles 2^a (w19, part B)
    'zweip':         ((T(W19, r'Mitte = Zweierpotenz\. Primzahlen <= (\d+) \('),), lambda a: a, 'zehnerpotenz', 'Primzahlschranke'),
    'zweiamax':      ((T(W19, r', a <= (\d+)\. Wieferich'),), lambda a: a, 'ganz_gruppiert', 'groesster Exponent a'),
    'zweiexp':       ((T(W19, r', a <= (\d+)\. Wieferich'), T(W19, r'^\s+2\^a - 1 ohne Zeugen <= \d+: \d+ von (\d+) Exponenten')),
                      lambda a, b: gleich(a - 1, b), 'ganz_gruppiert', 'Exponenten 2 ≤ a ≤ amax'),
    'zweiminus':     ((T(W19, r'^\s+2\^a - 1 ohne Zeugen <= \d+: (\d+) von'),), lambda a: a, 'ganz', '2^a − 1 ohne Zeugen'),
    'zweiminusant':  ((T(W19, r'^\s+2\^a - 1 ohne Zeugen <= \d+: \d+ von \d+ Exponenten \((\S+) %\)'),), lambda a: a, 'runden2', 'Anteil'),
    'zweiplus':      ((T(W19, r'^\s+2\^a \+ 1 ohne Zeugen <= \d+: (\d+) von'),), lambda a: a, 'ganz', '2^a + 1 ohne Zeugen'),
    'zweiplusant':   ((T(W19, r'^\s+2\^a \+ 1 ohne Zeugen <= \d+: \d+ von \d+ \((\S+) %\)'),), lambda a: a, 'runden2', 'Anteil'),
    'zweiminusok':   ((T(W19, r'^\s+2\^a - 1 ohne Zeugen <= \d+: \d+ von (\d+) Exponenten'), T(W19, r'^\s+2\^a - 1 ohne Zeugen <= \d+: (\d+) von')),
                      lambda n, x: n - x, 'ganz_gruppiert', '2^a − 1 beweisbar nicht powerful'),
    'zweiplusok':    ((T(W19, r'^\s+2\^a \+ 1 ohne Zeugen <= \d+: \d+ von (\d+) \('), T(W19, r'^\s+2\^a \+ 1 ohne Zeugen <= \d+: (\d+) von')),
                      lambda n, x: n - x, 'ganz_gruppiert', '2^a + 1 beweisbar nicht powerful'),
    'zweistellen':   ((T(W19, r', a <= (\d+)\. Wieferich'),), lambda a: math.floor(a * math.log10(2)) + 1, 'ganz_gruppiert', 'Stellen von 2^amax'),
    # Part I, Section 6, Prop. 6.11: Wieferich events of the 79 units (w137), table of the odd ranks (w141)
    'wfereignisse':  ((W141, W137), lambda a, b: gleich(a['ereignisse'], len(b['ereignisse'])), 'ganz', 'Ereignisse unter 3·10^6, 79 Einheiten'),
    'wfprimzahlen':  ((W141,), lambda a: a['verschiedene_p'], 'ganz', 'verschiedene Primzahlen darunter'),
    'wfkerne':       ((W137,), lambda a: a['kerne'], 'ganz', 'Einheiten (Kerne)'),
    'wfp':           ((W137,), lambda a: a['P'], 'wiss_exakt', 'Suchgrenze'),
    'wfungerade':    ((W141,), lambda a: a['ungerade_rang'], 'ganz', 'Ereignisse ungeraden Rangs ≥ 5'),
    'wfgerade':      ((W141,), lambda a: a['t1_gerade'], 'ganz', 'davon auf Einheiten mit T1 gerade'),
    'wfzeuge':       ((W141,), lambda a: a['zeuge'], 'ganz', 'davon mit Zeuge unter der Grenze'),
    'wfoffen':       ((W141,), lambda a: a['offen'], 'ganz', 'davon hier nicht entschieden'),
    'wfentschieden': ((W141,), lambda a: gleich(a['zeuge'] + a['voll'], a['t1_gerade'] - a['offen']), 'ganz', 'davon nicht powerful (Zeuge oder Vollzerlegung)'),
    'wfmgross':      ((W141,), lambda a: sorted(int(k) for k, v in a['m319_zerlegung'].items() if v == 1 and int(k) != 5)[0], 'ganz_gruppiert',
                      'M_5(319): kleinerer grosser Faktor (w131)'),
    'wfmgroesser':   ((W141,), lambda a: sorted(int(k) for k, v in a['m319_zerlegung'].items() if v == 1 and int(k) != 5)[1], 'ganz_gruppiert',
                      'M_5(319): groesserer grosser Faktor (w131)'),
    'wfmquadrat':    ((W141,), lambda a: next(int(k) for k, v in a['m319_zerlegung'].items() if v == 2), 'ganz', 'M_5(319): der quadrierte Faktor'),
    'erstesuchekerne': ((W20,), lambda a: a['NKERNE'], 'ganz', 'Kerne der ersten Suche (w20/w77/w84)'),
    'lucasgleich':   ((W84,), lambda a: gleich(a['beide'], a['n_vergleichbar']), 'ganz', 'Ereignisse, bei denen beide Bedingungen dieselben Primzahlen waehlen'),
    'lucaspk':       ((W84,), lambda a: len(a['pk_pell']), 'ganz', 'bekannte Primzahlen der Pell-Folge (Positivkontrolle)'),
    'lucasnk':       ((W84,), lambda a: a['e4_geprueft'], 'ganz', 'Primzahlen mit Exponent 1 (Negativkontrolle)'),
    'wftab':         ((W141,), lambda a: [[('ganz_gruppiert', z['m']), ('ganz', z['d']), ('ganz_gruppiert', z['p']), ('wort', z['status']),
                                           ('wort', 'frage_ja' if z['frage'] else 'frage_nein')]
                                          for z in a['zeilen']], 'tabzeilen', 'Tab. Prop. 6.11: Einheit, Rang, Primzahl, Status, Block in Prop. I.6.8'),
    'wfinpopzu':     ((W141, W138), wf_inpop_zu, 'ganz', 'Wieferich-Ereignisse mit Block in der 961-Menge, alle entschieden (w141 gegen w138)'),
    'wfnullrang':    ((W96N, W137, W91L, W138), wf_nullrang, 'ganz',
                      'Teil III: Nullstellen y_p = 0 mit Rang (Ordnung ≡ 0 mod 4) = Teil-I-Ereignisse auf den 25 Kernen'),
    # Part II, Section 4: size of the two sets of families
    'csechsbmax':    ((W113,), lambda a: a['b2'], 'zehnerpotenz', 'Schranke fuer den quadratfreien Teil, Familien mit Quadrat (C6)'),
    'csechsfamq':    ((W113,), lambda a: a['familien'], 'ganz_gruppiert', 'nichtleere Familien mit Quadrat bis zur Schranke (C6)'),
    'csechsfamw':    ((W115,), lambda a: a['familien'], 'ganz_gruppiert', 'nichtleere Familien ohne Quadrat bis zur Schranke (C6)'),
    # Part II, Section 6
    'kdkerne':       ((W9C3,), lambda a: w9c3_pruef(a) and a['fam'], 'ganz_gruppiert', 'quadratfreie m ≡ 3 (mod 8), m ≤ MMAX (Prop. klassedrei)'),
    'kdschranke':    ((W9C3,), lambda a: a['MMAX'], 'zehnerpotenz', 'Kernschranke Klasse 3'),
    'kdkandidaten':  ((W9C3,), lambda a: w9c3_pruef(a) and a['cand_fam'], 'ganz', 'Klasse 3: Familien nach dem Verwerfungskriterium'),
    'kdpunkte':      ((W9C3,), lambda a: w9c3_pruef(a) and a['pts'], 'ganz', 'Klasse 3: Gitterpunkte mit T_k < 10^H'),
    'kdungerade':    ((W9C3,), lambda a: w9c3_pruef(a) and a['dead_par'], 'ganz', 'Klasse 3: davon T_k ungerade (kein Paar)'),
    'kdpaare':       ((W9C3,), lambda a: w9c3_pruef(a) and a['dead_wit'], 'ganz', 'Klasse 3: Paare (Mitte ≡ 2 mod 4)'),
    'kdkernmax':     ((W9C3,), lambda a: max(c[0] for c in a['cands']), 'ganz_gruppiert', 'Klasse 3: groesster Kern einer Kandidaten-Familie'),
    'kdpaarkerne':   ((W9C3,), lambda a: len({x[0] for x in a['points'] if x[4] == 'witness'}), 'ganz', 'Klasse 3: Kerne mit mindestens einem Paar'),
    'kddrei':        ((W9C3,), lambda a: sum(1 for x in a['points'] if x[0] == 3 and x[4] == 'witness'), 'ganz', 'Klasse 3: Paare mit Kern 3'),
    'kdreinhartst':  ((T(W16C3, r'm = 20256129307923: m \| U1 = True, T1 gerade, log10 T1 = [\d.]+ \((\d+) Stellen\)'),), lambda a: a, 'ganz_gruppiert',
                      'Stellen von T1(20256129307923) (Reinhart-Kern der Klasse 3, T1 gerade)'),
    'kdreinhartungst': ((T(W16C3, r'm = 1752299: m \| U1 = True, T1 ungerade, log10 T1 = [\d.]+ \((\d+) Stellen\)'),), lambda a: a, 'ganz',
                      'Stellen von T1(1752299) (T1 ungerade)'),
    'turmtab':       ((TA(W13, r'bis 10\^\s*(\d+): gezaehlt\s+(\d+) \| Turm-Summe\s+(\d+) OK \| naiv \(log10 T1\)\s+(\d+)'), W9, W9C3), turm_tab,
                      'tabzeilen', 'Tab. Prop. turmsumme: Hoehe, gezaehlt Kl. 7, Kl. 3, naive Formel Kl. 7, Kl. 3'),
    'turmc7':        ((W9,), turm_steigung, 'abrunden4', 'Steigung der Tuerme, Klasse 7 (Paare je Dezimalstelle)'),
    'turmc3':        ((W9C3,), turm_steigung, 'abrunden4', 'Steigung der Tuerme, Klasse 3'),
    'turmc':         ((W9, W9C3, T(W15, r'Beide Klassen: \d+ Tuerme, c = (\S+) Paare')), lambda a, b, c: nahe(turm_steigung(a) + turm_steigung(b), c),
                      'abrunden3', 'Kor. untergrenze: Steigung beider Klassen, nach unten (gegen w15)'),
    'turmtuerme':    ((W9, W9C3, T(W15, r'Beide Klassen: (\d+) Tuerme')), lambda a, b, c: gleich(tuerme_gerade(a) + tuerme_gerade(b), c), 'ganz',
                      'Kor. untergrenze: Tuerme mit T1 gerade, beide Klassen'),
    'turmtuerme7':   ((W9,), tuerme_gerade, 'ganz', 'Tuerme Klasse 7'),
    'turmtuerme3':   ((W9C3,), tuerme_gerade, 'ganz', 'Tuerme Klasse 3'),
    'turmueber':     ((T(W15, r'fuer h bis 10\^6: ([\d.]+)'),), lambda a: a, 'abrunden1', 'Kor. untergrenze: kleinster Ueberschuss fuer h ≤ 10^6 (w15)'),
    'wtab':          ((W23C7, W23C3), c_tab, 'tabzeilen', 'Frage wachstum: log10 M, c_M Klasse 7, c_M Klasse 3 (w23)'),
    'modellp':       ((W24C7,), lambda a: all(str(p) in a['p_treffer'] for p in MODELL_P) and MODELL_P[-1], 'ganz', 'Modellpruefung: Primzahlen bis hier'),
    'modell7min':    ((W24C7,), lambda a: min(modell_quoten(a)), 'runden2', 'Modellpruefung Klasse 7: kleinste Quote'),
    'modell7max':    ((W24C7,), lambda a: max(modell_quoten(a)), 'runden2', 'Modellpruefung Klasse 7: groesste Quote'),
    'modell3min':    ((W24C3,), lambda a: min(modell_quoten(a)), 'runden2', 'Modellpruefung Klasse 3: kleinste Quote'),
    'modell3max':    ((W24C3,), lambda a: max(modell_quoten(a)), 'runden2', 'Modellpruefung Klasse 3: groesste Quote'),
    'beins7':        ((W24C7,), lambda a: 100 * a['b_verteilung']['1'] / a['lebend'], 'runden1', 'Anteil (%) der Kerne mit m′ = m, Klasse 7'),
    'beins3':        ((W24C3,), lambda a: 100 * a['b_verteilung']['1'] / a['lebend'], 'runden1', 'Anteil (%) der Kerne mit m′ = m, Klasse 3'),
    'einkernm':      ((W23C3,), lambda a: einkern(a)[0], 'ganz_gruppiert', 'Klasse 3: der Kern mit m′ = 3 in (10^5, 10^6]'),
    'einkernanteil': ((W23C3,), lambda a: 100 * einkern(a)[2], 'ganz_gerundet', 'sein Anteil (%) am Zuwachs von c_M von 10^5 bis 10^6 (gerundet)'),
    # Part II, Section 5
    'modellc':       ((W110, W107), lambda a, b: nahe(a['C'], b['C'], 1e-9), 'runden4', 'Vergleichsmodell: Korrekturfaktor C = ∏ C_p (w107 = w110)'),
    'modellrate':    ((W110,), lambda a: a['steigung_heuristik'], 'runden4', 'Vergleichsmodell: Rate C·c₁²/4 je Einheit ln x'),
    'naivrate':      ((W110,), lambda a: a['steigung_naiv'], 'runden4', 'unabhaengiges Modell: Rate c₁²/4 je Einheit ln x'),
    'modelltab':     ((W110,), modell_tab, 'tabzeilen', 'Tab. II.5: log10 X, gezaehlt, Vergleichsfunktion, unabhaengiges Modell'),
    'modellnull':    ((W110, W111), lambda a, b: modell_null(a, b, 'gemessen'), 'ganz', 'Paare bis n₀ (gezaehlt)'),
    'modellnullkorr': ((W110, W111), lambda a, b: modell_null(a, b, 'korrigiert'), 'runden2', 'Vergleichsfunktion bei n₀'),
    'modellnullnaiv': ((W110, W111), lambda a, b: modell_null(a, b, 'naiv'), 'runden2', 'unabhaengiges Modell bei n₀'),
    'modellneu':     ((W110,), lambda a: a['neuer_teil']['gemessen'], 'ganz', 'Paare in (2^52, n₀] (das Modell wurde bis 2^52 gemessen)'),
    'modellneukorr': ((W110,), lambda a: a['neuer_teil']['korrigiert'], 'runden2', 'Vergleichsfunktion auf (2^52, n₀]'),
    'modellneunaiv': ((W110,), lambda a: a['neuer_teil']['naiv'], 'runden2', 'unabhaengiges Modell auf (2^52, n₀]'),
    'quadbmax':      ((W112,), lambda a: a['bmax'], 'zehnerpotenz', 'Quadrat-Familien mit b, d bis hier (C7)'),
    'quadtab':       ((W112,), lambda a: [[('runden1', z['L']), ('ganz', z['NQ']), ('runden1', z['korr']), ('runden1', z['naiv']), ('runden4', z['steigung'])]
                                          for z in a['untergrenze']['zeilen']], 'tabzeilen', 'Tab. II.6: L = ln x, Quadrat-Paare exakt, Vergleichsfunktion, unabhaengig, Rate der geborenen Familien'),
    'quadkreuz':     ((W112,), lambda a: a['untergrenze']['ueberholt_korrigiert_ab_L'], 'runden1', 'C7: ab diesem L = ln x liegt die exakte Zahl ueber der Vergleichsfunktion'),
    'quadkreuzlog':  ((W112,), lambda a: a['untergrenze']['ueberholt_korrigiert_ab_L'] / math.log(10), 'runden1', 'dasselbe als log10 x'),
    'quadgeboren':   ((W112,), lambda a: quad_rest(a)['geboren'], 'ganz', 'C7-Schwanz: Familien, die bis L = 1000 geboren sind'),
    'quadrate1000':  ((W112,), lambda a: quad_rest(a)['steigung'], 'runden4', 'C7-Schwanz: ihre Rate Σ 1/ln F'),
    'quadnq1000':    ((W112,), lambda a: quad_rest(a)['NQ'], 'ganz', 'C7-Schwanz: Quadrat-Paare bis e^1000'),
    'quadv1000':     ((W112,), lambda a: quad_rest(a)['korr'], 'runden1', 'C7-Schwanz: Vergleichsfunktion bei e^1000'),
    'quadrest1000':  ((W112,), lambda a: quad_rest(a)['rest'], 'ganz', 'C7-Schwanz: Paare minus geborene Familien (> Vergleichsfunktion, geprueft)'),
    'teilsumq1':     ((W113,), lambda a: a['teilsummen']['10'], 'runden4', 'Σ 1/ln F, Familien mit Quadrat, b ≤ 10'),
    'teilsumq2':     ((W113,), lambda a: a['teilsummen']['100'], 'runden4', 'Σ 1/ln F, Familien mit Quadrat, b ≤ 100'),
    'teilsumtab':    ((W113, W115), lambda a, b: [[('ganz', j), ('runden4', a['teilsummen'][str(10**j)]), ('runden4', b['summe_1_durch_lnF'][str(10**j)]),
                                                   ('runden4', a['teilsummen'][str(10**j)] + b['summe_1_durch_lnF'][str(10**j)])] for j in range(3, 7)],
                      'tabzeilen', 'Tab. II.7: log10 B, Σ 1/ln F mit Quadrat (b ≤ B), ohne Quadrat (bd ≤ B), Summe'),
    'zuwq':          ((W121,), zuwachs_klein('Q'), 'ganz_gerundet', 'Zuwachs (10^5, 10^6] aus Familien mit ln ε ≤ 100, mit Quadrat (%)'),
    'zuww':          ((W121,), zuwachs_klein('W'), 'ganz_gerundet', 'dasselbe ohne Quadrat (%)'),
    'zuwtopq':       ((W121,), lambda a: 100 * a['zerlegung']['Q'][-1]['top1pct'], 'ganz_gerundet', 'Anteil des obersten Prozents der Familien am Zuwachs, mit Quadrat (%)'),
    'zuwtopw':       ((W121,), lambda a: 100 * a['zerlegung']['W'][-1]['top1pct'], 'ganz_gerundet', 'dasselbe ohne Quadrat (%)'),
    # Part II, Section 7: gap-1 route (w25 cross-check, w26 pairs + third number, w27 deeper pass)
    'a1nmax':        ((W25,), lambda a: a['NMAX'], 'zehnerpotenz', 'w25: Paare bis hier, zwei Routen'),
    'a1dmax25':      ((W25,), lambda a: a['DMAX'], 'zehnerpotenz', 'w25: D bis hier'),
    'a1paare25':     ((W25,), lambda a: gleich(len(a['paare_direkt']), len(a['paare_adressen'])) and not a['nur_direkt'] and len(a['paare_direkt']),
                      'ganz', 'w25: Paare, Direktsuche = Pell-Route'),
    'a1dmax':        ((W26,), lambda a: a['DMAX'], 'ganz_gruppiert', 'w26: D bis hier'),
    'a1kerne':       ((W26,), lambda a: a['familien'], 'ganz_gruppiert', 'w26: quadratfreie D bis DMAX'),
    'a1paare':       ((W26,), lambda a: gleich(a['kandidaten'], len(a['ueberlebt']) + len(a['offen'])) and a['kandidaten'], 'ganz_gruppiert',
                      'w26: Paare (a, a+1) mit D ≤ DMAX, T < 10^H'),
    'a1hoehe':       ((W26,), lambda a: int(a['H']), 'ganz', 'w26: Hoehe'),
    'a1ueber':       ((W26,), a1_ueber, 'ganz', 'w26: Paare knapp UEBER 10^H (je Kern die naechste Sprosse), exakt gezaehlt'),
    'a1paareh':      ((W26,), lambda a: a['kandidaten'] - a1_ueber(a), 'ganz_gruppiert', 'w26: Paare mit 2a + 1 < 10^H, exakt'),
    'a1p1':          ((W26,), lambda a: a['PMAX'], 'wiss_exakt', 'w26: Zeugensuche unter'),
    'a1zeuge1':      ((W26,), lambda a: len(a['ueberlebt']), 'ganz_gruppiert', 'w26: dritte Zahl mit Zeugen unter PMAX'),
    'a1p2':          ((W27,), lambda a: a['P'], 'wiss_exakt', 'w27: Zeugensuche unter'),
    'a1zeuge2':      ((W27,), lambda a: len(a['neu_erledigt']), 'ganz', 'w27: dazu mit Zeugen unter P'),
    'a1rest2':       ((W27, W26), lambda a, b: gleich(len(a['bleibt_offen']) + len(a['neu_erledigt']), len(b['offen'])) and len(a['bleibt_offen']), 'ganz',
                      'w27: ohne Zeugen unter P'),
    'a1s2prim':      ((W144, W27), lambda a, b: s2(a, b, 'prim_zertifiziert', 'prim_wahrscheinlich'), 'ganz', 'w144: Kofaktor prim ⇒ Zeuge'),
    'a1s2primzert':  ((W144, W27), lambda a, b: s2(a, b, 'prim_zertifiziert'), 'ganz', 'w144: davon zertifiziert (< ψ12)'),
    'a1s2ecm':       ((W144, W27), lambda a, b: s2(a, b, 'ecm_zeuge_zertifiziert', 'ecm_zeuge_wahrscheinlich'), 'ganz', 'w144: Zeuge durch ECM'),
    'a1s2offen':     ((W144, W27), lambda a, b: s2(a, b, 'offen', 'offen_potenz'), 'ganz', 'w144: weiter offen'),
    'a1s2kleinst':   ((W144, W27), lambda a, b: s2min(a, b), 'ganz', 'w144: Stellen der kleinsten weiter offenen dritten Zahl'),
    # Quality check of Part II: hand-typed numbers -> keys
    'a1s2probable':  ((W144, W27), lambda a, b: s2(a, b, 'prim_wahrscheinlich', 'ecm_zeuge_wahrscheinlich'), 'ganz', 'w144: Zeugen, die nur wahrscheinlich prim sind'),
    # all probable-prime cofactors from w144 proved by ECPP (w165, own checker) - otherwise the build aborts
    'a1s2ecpp':      ((W144, W27, W165), lambda a, b, c: gleich(gleich(c['bewiesen'], c['n']), s2(a, b, 'prim_wahrscheinlich', 'ecm_zeuge_wahrscheinlich')),
                      'ganz', 'w165: per ECPP bewiesene Kofaktoren (= alle vorher nur wahrscheinlichen)'),
    'a1s2ecppmax':   ((W165,), lambda a: a['max_stellen'], 'ganz', 'w165: groesster bewiesener Kofaktor (Stellen)'),
    'a1s2nk':        ((W144,), lambda a: a['nk_ok'] and len(a['nk']), 'ganz', 'w144: Paare der Negativkontrolle'),
    'quadzeilen':    ((W112,), lambda a: len(a['untergrenze']['zeilen']), 'ganz', 'Tab. II.4: Zeilen'),
    'quadab':        ((W112,), lambda a: quad_ab(a), 'ganz', 'Tab. II.4: ab diesem L bleibt N_Q* in der Tabelle ueber V'),
    'quadablog':     ((W112,), lambda a: math.floor(quad_ab(a) / math.log(10)), 'ganz', 'dasselbe als log10 x, abgerundet'),
    'quadslack':     ((W112,), lambda a: math.ceil(0.05 * quad_rest(a)['geboren']), 'ganz', 'C7-Schwanz: 0,05 je geborene Familie, aufgerundet'),
    'wtabmax':       ((W23C7, W23C3), lambda a, b: gleich(a['MMAX'], b['MMAX']), 'zehnerpotenz', 'Tab. II.7 und w24: Kerne bis hier'),
    # the parity comes from w152 (computed), not from the labels `(T_1 gerade)` of w11 (wrong for 117 477 414 815)
    'reinhartgerade': ((W152, T(W16C3, r'm = (20256129307923): m \| U1 = True, T1 gerade')), lambda a, b: gleich(len(a['gerade_klasse7']) + len(a['gerade_klasse3']), len(a['gerade_klasse7']) + (1 if a['gerade_klasse3'] == [int(b)] else 0)), 'ganz',
                      'Reinhart-Kerne mit T1 gerade (Klassen 7 und 3), Paritaet aus w152; zweite Route: w16 nennt dieselbe Klasse-3-Zahl'),
    'reinhartgeradesicher': ((W152,),
                      lambda a: sum(m < REINHART_SICHER for m in a['gerade_klasse7'] + a['gerade_klasse3']), 'ganz',
                      'davon unter 1.5e12 (Reinhart 2023, Remark 5.5 — vollstaendig): E_3-Kerne ohne Annahme'),
    'reinhartsieben': ((TA(W11, r'^m = (\d+) \(T_1'),), lambda a: len(a), 'ganz', 'Reinharts Kerne der Klasse 7 (m | U1): alle vier'),
    # two routes per number - the distribution omega and the list of the prime blocks must give the same count
    'znab':          ((W55,), lambda a: a['BMAX'], 'ganz', 'w55: squarefree b bis hier'),
    'znaq':          ((W55,), lambda a: gleich(a['anzahl_mit_loesung'] + a['anzahl_ohne'], quadratfrei_bis(a['BMAX'])), 'ganz',
                      'w55: so viele quadratfreie b ≤ BMAX (zweite Route: direkt gezaehlt)'),
    # the lists in w55 are shortened (smallest: 8, without: 60) - the sum with/without is checked by `znaq`
    'znaqg':         ((W55,), lambda a: (gleich(a['ohne_loesung'][0], 1), a['anzahl_mit_loesung'] + a['anzahl_ohne'] - 1)[1], 'ganz',
                      'w55: quadratfreie b mit 1 < b ≤ BMAX (b = 1 steht in der Liste ohne Loesung und zaehlt nicht)'),
    'znamit':        ((W55,), lambda a: a['anzahl_mit_loesung'], 'ganz', 'w55: b mit Loesung von x^2 - b^3 y^2 = -1'),
    # second branch - Lavi's term and the bound of the direct completeness search from w157
    'phidreilavij':  (('w157_phidrei_zweiter_zweig_result.json',), lambda a: a['lavi_j'][0] if len(a['lavi_j']) == 1 and a['gesamt'] == 0 else None,
                      'ganz', 'w157: Lavis x ist das Glied j des zweiten Zweigs'),
    'phidreibrute':  (('w157_phidrei_zweiter_zweig_result.json',), lambda a: a['mbrute'] if a['direkt_ausserhalb'] == 0 else None,
                      'zehnerpotenz', 'w157: direkte Suche aller Loesungen mit m bis hier, keine ausserhalb der zwei Zweige'),
    'phidreistellen': ((W61,), lambda a: ', '.join(str(f['stellen']) for f in a['familie'][1:-1]) + ' and ' + str(a['familie'][-1]['stellen'])
                      if all(len(f['x']) == f['stellen'] for f in a['familie']) else None, 'wortliste_und', 'w61: Stellen der naechsten Familienglieder'),
    'phidreiperiode': ((W61,), lambda a: gleich(a['pell']['periode_m_mod_7'], a['familie'][1]['schritt'] - a['familie'][0]['schritt']), 'ganz',
                      'w61: jeder so vielte Schritt der Pell-Kette'),
    'zyktab':        ((W62,), zyk_tab, 'tabzeilen', 'Teil III Tab.: powerful Phi_d(x) bis zur Wertschranke, d = 3..12'),
    'zykn':          ((W62,), lambda a: a['tafel']['geprueft'], 'ganz_gruppiert', 'w62: gepruefte Paare (d, x)'),
    'zyktreffer':    ((W62,), lambda a: a['tafel']['treffer'], 'ganz', 'w62: Treffer'),
    'zykschranke':   ((W62,), lambda a: a['tafel']['wert_schranke'], 'zehnerpotenz', 'w62: Wertschranke der Tafel'),
    # additionally checked - outside the class exactly ONE square, Phi_5(3) (the text says it is the only square outside)
    'zykklasse':     ((W62,), lambda a: (gleich(a['klassengrenze']['quadrate_innerhalb'], 0), gleich(a['klassengrenze']['quadrate_ausserhalb'], 1),
                      gleich((a['klassengrenze']['kleinstes_ausserhalb']['x'], a['klassengrenze']['kleinstes_ausserhalb']['d']), (3, 5)),
                      a['klassengrenze']['wert_schranke'])[3], 'zehnerpotenz',
                      'w62: bis hier kein Quadrat Phi_d(x) mit x ≡ 2 (mod 4), d in {4,5,7,...,12}'),
    'lstab':         ((W48S, W125), ls_tab, 'tabzeilen', 'Teil III Tab. Kopf: je Schranke Bestand, Anteil powerful, Kopf tie-fest (w125) gegen w48'),
    'lsnpw':         ((W125,), lambda a: a['ergebnis']['100000000000000']['n_powerful'], 'ganz_gruppiert', 'powerful Zahlen bis 10^14'),
    'lswerte':       ((W125,), lambda a: a['ergebnis']['100000000000000']['n_werte'], 'ganz_gruppiert', 'verschiedene Lueckenwerte bis 10^14'),
    'lsanteil':      ((W48S,), lambda a: 100 * a['ergebnis']['100000000000000']['anteil_pw_gesamt'], 'runden2', 'Anteil powerful unter den Werten (%)'),
    'lsbis':         ((W125,), lambda a: a['ergebnis']['100000000000000']['alle_powerful_bis_rang'], 'ganz', 'die so vielen haeufigsten Werte: alle powerful'),
    'lsocc':         ((W125,), lambda a: gleich(a['ergebnis']['100000000000000']['erster_nicht_pw']['occ'],
                      a['ergebnis']['100000000000000']['schwellen']['50']['schwelle']), 'ganz', 'Vorkommen des ersten nicht-powerful Werts'),
    'lserster':      ((W125,), lambda a: a['ergebnis']['100000000000000']['erster_nicht_pw']['wert'], 'ganz_gruppiert', 'erster nicht-powerful Wert'),
    'lsg100':        ((W125,), lambda a: a['ergebnis']['100000000000000']['schwellen']['100']['gruppengroesse'], 'ganz', 'Gruppe um Rang 100'),
    'lsg100pw':      ((W125,), lambda a: gleich(a['ergebnis']['100000000000000']['schwellen']['100']['davon_powerful'],
                      a['ergebnis']['100000000000000']['schwellen']['100']['gruppengroesse'] - 1), 'ganz', 'davon powerful (alle bis auf einen)'),
    'lsg200':        ((W125,), lambda a: a['ergebnis']['100000000000000']['schwellen']['200']['gruppengroesse'], 'ganz', 'Gruppe um Rang 200'),
    'lsg200pw':      ((W125,), lambda a: a['ergebnis']['100000000000000']['schwellen']['200']['davon_powerful'], 'ganz', 'davon powerful'),
    'lserwartet':    ((W48S,), lambda a: 50 * a['ergebnis']['100000000000000']['anteil_pw_gesamt'], 'runden2', 'Nullmodell: erwartet powerful unter 50'),
    'lsdecke':       ((W148,), ls_decke, 'tabzeilen', 'Teil III Tab. Decke: die 10 nicht-powerful Werte der Gruppe um Rang 200 (w148, tie-fest)'),
    'lsdeckeanz':    ((W148, W125), lambda a, b: gleich(len(a['g200']['nicht_powerful']),
                      b['ergebnis']['100000000000000']['schwellen']['200']['gruppengroesse'] - b['ergebnis']['100000000000000']['schwellen']['200']['davon_powerful']),
                      'ganz', 'nicht-powerful Werte der Gruppe um Rang 200: w148 = w125'),
    'lsg20':         ((W148,), lambda a: a['g20']['groesse'], 'ganz', 'Gruppe um Rang 20 (Gleichstaende eingeschlossen)'),
    'lsanteilpw':    ((W148,), lambda a: '--'.join(ls_anteil(a, True)), 'tex_spanne', 'Anteil gcd = sqfull in der Gruppe um Rang 20 (min--max, %)'),
    'lsanteilnpw':   ((W148,), lambda a: '--'.join(ls_anteil(a, False)), 'tex_spanne', 'dasselbe bei den 10 nicht-powerful Werten'),
    'nvlmax':        ((W109V,), lambda a: max(a['L']), 'ganz', 'w109: groesstes Fenster (in mittleren Abstaenden)'),
    'nvrel':         ((W109V,), lambda a: 100 * a['urteil']['Q_gittersumme']['rel'], 'runden1', 'w109: relative Abweichung von der Gittersumme (%)'),
    'nvsteig':       ((W109V,), lambda a: a['urteil']['steigung_Q'], 'runden2', 'w109: gemessene Steigung (log-log) bei grossen Fenstern'),
    'lstop1':        ((W125,), lambda a: a['ergebnis']['100000000000000']['top1'][0], 'ganz_gruppiert', 'haeufigster Wert bis 10^14'),
    'zmn':           ((W59,), lambda a: a['n_werte'], 'ganz_gruppiert', 'Zaehlmodell: Lueckenwerte mit mind. MINCNT Vorkommen'),
    'zmmin':         ((W59,), lambda a: a['MINCNT'], 'ganz', 'Mindestzahl der Vorkommen'),
    # same basis as the counting model, second route via w160 (own count, compared with w59)
    'lsbasiswerte':  ((W160, W59), lambda a, b: gleich(a['basis_werte'], b['n_werte']), 'ganz_gruppiert', 'w160: Lueckenwerte mit mind. 20 Vorkommen (= w59)'),
    'lsbasismind':   ((W160, W59), lambda a, b: gleich(a['basis_mind'], b['MINCNT']), 'ganz', 'w160: Mindestzahl der Vorkommen (= w59)'),
    'lsbasisanteil': ((W160,), lambda a: a['basis_anteil_prozent'], 'runden2', 'w160: Anteil powerful unter diesen Werten (%)'),
    'lsbasiserw':    ((W160,), lambda a: a['nullmodell_50_basis'], 'runden1', 'w160: Nullmodell, erwartet powerful unter 50 aus dieser Basis'),
    # the accuracy of about 10^{-12} stated in the text is measured in w161 (no earlier measurement existed)
    'fpn':           ((W161,), lambda a: gleich(a['n_kerne'], a['n_bis_2e4'] + a['n_zufall_bis_1e8']), 'ganz_gruppiert', 'w161: Kerne der Stichprobe'),
    'fpmax':         ((W161,), lambda a: a['max_rel_fehler'], 'aufrunden_wiss1', 'w161: groesster relativer Fehler von log10 T_1 (Gleitkomma), als Schranke'),
    'zmr2':          ((W59,), lambda a: a['G_entscheidung']['cnt_modell_mittel'], 'runden3', 'R^2 ausserhalb der Anpassung, Modell (Mittel ueber 5 Splits)'),
    'zmr2sd':        ((W59,), lambda a: a['G_entscheidung']['cnt_modell_sd'], 'runden3', 'dessen Streuung'),
    'zmroh':         ((W59,), lambda a: a['G_entscheidung']['cnt_roh_mittel'], 'runden3', 'R^2 Grundlinie (rohe Merkmale)'),
    'zmrohsd':       ((W59,), lambda a: a['G_entscheidung']['cnt_roh_sd'], 'runden3', 'dessen Streuung'),
    'zmr2log':       ((W59,), lambda a: a['G_entscheidung']['log_modell_mittel'], 'runden3', 'R^2 auf log occ, Modell (G)'),
    'zmr2logsd':     ((W59,), lambda a: a['G_entscheidung']['log_modell_sd'], 'runden3', 'dessen Streuung'),
    'zmallesep':     ((W59,), lambda a: a['E_streuung']['sep_mittel'], 'runden3', '(E) Fassung ueber alle powerful Teiler, neue g'),
    'zmallesepsd':   ((W59,), lambda a: a['E_streuung']['sep_sd'], 'runden3', 'dessen Streuung'),
    'zmalleroh':     ((W59,), lambda a: a['E_streuung']['rohe_mittel'], 'runden3', '(E) rohe Merkmale'),
    'zmallerohsd':   ((W59,), lambda a: a['E_streuung']['rohe_sd'], 'runden3', 'dessen Streuung'),
    'zmdeutstufen':  ((W64,), lambda a: a['n_analyse_stufen'], 'ganz_gruppiert', 'w64: Kofaktoren mit mindestens KMIN Lueckenwerten'),
    'zmdeutkmin':    ((W64,), lambda a: a['KMIN'], 'ganz', 'w64: KMIN'),
    'zmkern':        ((W59,), lambda a: a['H_ablation']['nur_kern_mittel'], 'runden3', 'Ablation: nur Kern'),
    'zmkof':         ((W59,), lambda a: a['H_ablation']['nur_kofaktor_mittel'], 'runden3', 'Ablation: nur Kofaktor'),
    'zmww':          ((W59,), lambda a: a['F_vorregistriert']['wechselwirkung_r'], 'runden3', 'Wechselwirkung u*v gegen Residuen'),
    'zmkerne':       ((W59,), lambda a: a['F_vorregistriert']['n_kerne_sqfull'], 'ganz_gruppiert', 'Kerne sqfull(g)'),
    'zmkof_n':       ((W59,), lambda a: a['F_vorregistriert']['n_kofaktoren_quadratfrei'], 'ganz_gruppiert', 'quadratfreie Kofaktoren'),
    'zmdecke':       ((W63,), lambda a: a['obergrenze_schritt2']['reliabilitaet_voll_spearman_brown'], 'runden3', 'Reliabilitaet von v (Decke)'),
    'zmdeutung':     ((W64,), lambda a: a['regression']['r2_ausserhalb'], 'runden3', 'Deutung von v: R^2 ausserhalb'),
    'zmdeutdecke':   ((W64,), lambda a: a['decke']['spearman_brown'], 'runden3', 'Decke auf derselben Menge'),
    'pmstich':       ((W65,), lambda a: a['primitivitaet']['stichprobe'], 'ganz_gruppiert', 'Stichprobe: jedes k-te Paar'),
    'pmschritt':     ((W65,), lambda a: a['primitivitaet']['schritt'], 'ganz', 'Schrittweite der Stichprobe'),
    'pmanteil':      ((W65,), lambda a: 100 * gleich(a['primitivitaet']['primitiv'], round(a['primitivitaet']['anteil_primitiv'] * a['primitivitaet']['stichprobe']))
                      / a['primitivitaet']['stichprobe'], 'runden2', 'Anteil primitiver Paare in der Stichprobe (%)'),
    'pmgcd':         ((W65,), lambda a: ', '.join(str(x) for x in a['primitivitaet']['haeufigste_gcd']), 'wortliste_ganz', 'haeufigste gcd'),
    # ---- Part III, Section 3: independence along the divisors of the index, w35 (printed variant: P = 2*10^6, 25 kernels)
    'ukerne':        ((W35L,), lambda a: gleich(a['NKERNE'], len(a['kerne'])), 'ganz', 'w35: Kerne'),
    'up':            ((W35L,), lambda a: a['P'], 'wiss_exakt', 'w35: Primzahlschranke'),
    'upaare':        ((W35L,), lambda a: gleich(a['n_paare'], sum(a['tafel'].values())), 'ganz_gruppiert', 'w35: Paare (d, qd), beide Raenge besetzt'),
    'u11':           ((W35L,), lambda a: a['tafel']['n11'], 'ganz', 'w35: beide Wieferich'),
    'u10':           ((W35L,), lambda a: a['tafel']['n10'], 'ganz', 'w35: nur d'),
    'u01':           ((W35L,), lambda a: a['tafel']['n01'], 'ganz', 'w35: nur qd'),
    'uerw':          ((W35L,), lambda a: a['erwartet_n11'], 'runden3', 'w35: erwartete gemeinsame Ereignisse unter Unabhaengigkeit'),
    'usim':          ((W35L,), lambda a: a['SIM'], 'ganz', 'w35: Simulationen je Kopplungsstaerke'),
    # w35 calls sim() a SECOND time for the JSON (98.75 % / 0.0 %), while the output shows the first run (98.0 % / 0.2 %).
    #   Hence what is printed is rounded down and given as a bound that holds for BOTH runs.
    'utrenn':        ((W35L,), lambda a: math.floor(100 * a['trennschaerfe']['0.1']), 'ganz', 'w35: mindestens so viel % der Simulationen mit Befund bei rho = 0.1'),
    'utrennnull':    ((W35L,), lambda a: (1 if 100 * a['trennschaerfe']['0.0'] < 1 else None), 'ganz', 'w35: Fehlalarm bei rho = 0 unter so viel %'),
    'ukopie':        ((W35L,), lambda a: 0.1 * a['tafel']['n10'], 'runden1', 'erwartete Kopien bei rho = 0.1: rho mal die Faelle mit Status(d) = 1'),
    # ---- Part III, Section 4: Wieferich fractions of the Pell units - w91 (kernels, hits), w96 (re-measurement), w108
    #   (memory)
    'wf3kerne':      ((W91L,), lambda a: len(a['kerne']), 'ganz', 'w91: Kerne (dieselben 25 wie w35) — Teil III; hiess vorher wfkerne und ueberschrieb die Teil-I-Zahl (79)'),
    'wf3p':          ((W96N,), lambda a: a['p_max'], 'wiss_exakt', 'w96: Primzahlschranke — Teil III; hiess vorher wfp'),
    'wfn':           ((W96N,), lambda a: a['n_paare'], 'ganz_gruppiert', 'w96: Paare (m, p)'),
    'wftreffer':     ((W96N, W91L), lambda a, b: gleich(schwelle(a, 1e-7)['beobachtet'], b['n_treffer']), 'ganz', 'Wieferich-Ereignisse (y_p = 0): w96 = w91'),
    'wferw':         ((W96N, W91L), lambda a, b: gleich(round(schwelle(a, 1e-7)['erwartet'], 2), round(b['heuristik_summe'], 2)), 'runden2',
                      'Erwartung Σ 1/p: w96 = w91'),
    'wft4b':         ((W96N,), lambda a: schwelle(a, 1e-4)['beobachtet'], 'ganz', 'w96: y_p < 10^-4 beobachtet'),
    'wft4e':         ((W96N,), lambda a: schwelle(a, 1e-4)['erwartet'], 'runden1', 'w96: y_p < 10^-4 erwartet'),
    'wfzmax':        ((W96N,), lambda a: max(abs(s['z']) for s in a['a_schwellen']), 'runden2', 'w96: groesstes |z| ueber die sechs Schwellen'),
    'wfklpmin':      ((W96N,), lambda a: min(k['p'] for k in a['b_klassen']), 'runden3', 'w96: kleinstes p der Klassentests'),
    'wfklbon':       ((W96N,), lambda a: min(1.0, len(a['b_klassen']) * min(k['p'] for k in a['b_klassen'])), 'runden2', 'w96: Bonferroni der Klassentests'),
    'wfklanz':       ((W96N,), lambda a: len(a['b_klassen']), 'ganz', 'w96: Zahl der Klassentests'),
    'wglags':        ((W108G,), lambda a: a['lags'], 'ganz', 'w108: Verschiebungen'),
    'wgtests':       ((W108G,), lambda a: gleich(a['n_tests'], a['lags'] * len(a['kerne'])), 'ganz', 'w108: Tests = Verschiebungen × Kerne'),
    'wgueber':       ((W108G,), lambda a: a['ueber_2_576'], 'ganz', 'w108: |z| > 2.576'),
    'wgerw':         ((W108G,), lambda a: 0.01 * a['n_tests'], 'ganz_gerundet', 'w108: erwartet bei 1 %'),
    'wgzmax':        ((W108G,), lambda a: a['z_max'][0], 'runden2', 'w108: groesstes |z|'),
    'wgpbon':        ((W108G,), lambda a: a['p_bonferroni'], 'runden2', 'w108: Bonferroni-p'),
    'wgpk':          ((W108G,), lambda a: a['pk']['z_eingebaut'], 'runden1', 'w108: z bei eingebautem Gedaechtnis (PK)'),
    'wkpaare':       ((W96N,), lambda a: a['c']['paare'], 'ganz', 'w96: Kernpaare'),
    'wkzmax':        ((W96N,), lambda a: a['c']['z_max'][0], 'runden2', 'w96: groesstes |z| der Kreuzkorrelation'),
    'wkpbon':        ((W96N,), lambda a: a['c']['p_bonferroni'], 'runden2', 'w96: Bonferroni-p der Kreuzkorrelation'),
    'wkpk':          ((W96N,), lambda a: a['c']['pk_gepflanzt_z'], 'runden1', 'w96: z bei gepflanzter Korrelation (PK)'),
    'achtq':         ((W129,), lambda a: (gleich(a['treffer'], []), gleich(a['pk'], True), a['anzahl_q'])[2], 'ganz_gruppiert',
                      'w129: Primzahlen q ≡ 1 (mod 8) bis zur Grenze, bei keiner teilt q den Koeffizienten der Grundloesung (PK gruen)'),
    'achtgrenze':    ((W129,), lambda a: a['grenze'], 'zehnerpotenz', 'w129: Grenze'),
    'zpn':           ((W78,), lambda a: gleich(a['n'], sum(a['omega_verteilung'].values())), 'ganz', 'w78: vollstaendig faktorisierte Bloecke M_d'),
    'zpprim':        ((W78,), lambda a: gleich(len(a['phi_prim']), a['omega_verteilung']['1']), 'ganz', 'w78: davon prim (omega = 1)'),
    'zpdmax':        ((W78,), lambda a: gleich(a['d_max_omega1'], max(d for _, d in a['phi_prim'])), 'ganz', 'w78: groesstes d mit M_d prim'),
    'zpab16':        ((W78,), lambda a: gleich(len(a['pruefstein_verletzt']), sum(1 for _, d in a['phi_prim'] if 2 * d >= 31)), 'ganz',
                      'w78: prime Bloecke mit 2d >= 31 (dort haette die naive Uebertragung von Juricevic Thm 2.1 zwei primitive Teiler verlangt)'),
    'offenschwach':  ((W146,), lambda a: sum(1 for z in a['zeilen'] if z['kern'] == 'alt' and z['t35'] < 1), 'ganz', 'w146: offene Paare auf den ersten 25 Kernen mit t35 < 1 (kleinerer ECM-Aufwand)'),
    'offenneum':     ((W146,), lambda a: und_liste(sorted(z['m'] for z in a['zeilen'] if z['kern'] == 'neu')), 'wortliste_und', 'w146: Kerne der offenen Paare auf den anderen 54 Kernen'),
    'zquadkerne':    ((W153, W138), lambda a, b: (gleich(a['n_kerne'], b['kerne']), len(a['kerne_quadrat']))[1], 'ganz', 'w153: Kerne unter den 79 mit T_1 + 1 = Quadrat'),
    'zquadbl':       ((W153, W138), zquad_bl, 'ganz', 'w153: Bausteine auf diesen Kernen, jeder zusammengesetzt und keine Potenz (mindestens zwei Primteiler)'),
    'zquadprim':     ((W153,), zquad_prim, 'ganz', 'w153: prime Kerne m ≡ 7 (mod 8) unter der Grenze, alle mit T_1 + 1 = Quadrat'),
    'zquadgrenze':   ((W153,), lambda a: a['prim_grenze'], 'ganz', 'w153: Grenze fuer die primen Kerne'),
    'reinhartsiebensicher': ((TA(W11, r'^m = (\d+) \(T_1'),), lambda a: sum(int(m) < REINHART_SICHER for m in a), 'ganz',
                      'davon unter 1.5e12 (Reinhart 2023, Remark 5.5): 209991, 4099215, 117477414815'),
}

# THEOREM CHECKS (`PRUEF_SAETZE`): statements in the text WITHOUT a printed number of their own that still depend on data. If
#   one fails, the build aborts - exactly as for a number. Name = the statement as it appears in the text.
PRUEF_SAETZE = {
    'Prop. II.5 (C7): Rate der bis L = 1000 geborenen Quadrat-Familien > Modellrate, und Paare minus Familien > Vergleichsfunktion':
        ((W112, W110), lambda a, b: quad_rest(a)['steigung'] > b['steigung_heuristik'] and quad_rest(a)['rest'] > quad_rest(a)['korr']),
    'Thm 5.1: Leiter gibt mehr als CSS fuer 15, 31, 87, weniger fuer 23, 47; 39 und 55 haben T1 ungerade': ((L_ZEILEN,), css),
    'Thm 5.3: m = 209 991 hat T1 ungerade (Paritaet)': ((T(W11, r'^m = (209991) \(T_1 UNGERADE'),), lambda a: a == 209991),
    'Thm 5.2: die Periodenschranke erzwingt T1 > 10^H, keine Punkte, keine Kandidaten': ((W10,), lambda a: w10_rest(a) >= 0),
    'Prop. 6.9: kein Exponent, bei dem BEIDE ohne Zeugen sind; PK die zwei Wieferich-Primzahlen zur Basis 2; NK a = 3 (2^3 + 1 = 9) ohne Zeugen':
        ((T(W19, r'BEIDE ohne Zeugen [^\n]*?: (\d+) ->'), T(W19, r'Wieferich-Primzahlen zur Basis 2 gefunden: \[(\d+), \d+\]'),
          T(W19, r'Wieferich-Primzahlen zur Basis 2 gefunden: \[\d+, (\d+)\]'), T(W19, r'^\s+2\^a \+ 1 ohne Zeugen [^\n]*kleinste: \[(\d+),')),
         lambda beide, w1, w2, erst: beide == 0 and (w1, w2) == (1093, 3511) and erst == 3),
    'Prop. 6.11: auf den 16 Ereignissen waehlen beide Bedingungen dieselben Primzahlen; keine nur bei uns; kein Falsch-Positiv in der NK':
        ((W84,), lambda a: a['beide'] == a['n_vergleichbar'] and a['nur_unsere'] == 0 and a['e4_falsch_positiv'] == 0
         and all(x[2] is True for x in a['pk_pell'])),
    'Prop. 6.8 (a): keine Ausnahme von p ≡ ±1 (mod 4d)': ((T(W18, r'erfuellt \d+, Verstoesse (\d+) \[\]'),), lambda a: a == 0),
    'Prop. 6.8 (c): Summenprobe besetzt = mit Zeuge + nur Wieferich, in beiden Laeufen':
        ((T(W17C, r'besetzt (\d+) \('), T(W17C, r'mit Zeuge (\d+) \('), T(W17C, r'nur Wieferich (\d+),'),
          T(W17D, r'mit Zeuge (\d+) \('), T(W17D, r'nur Wieferich (\d+),')), lambda b, z, w, z1, w1: b == z + w == z1 + w1),
    'Remark 6.6′: die Mitten ausserhalb der Bedingung mod 36 und die Punkte mit Zeuge 3 sind DIESELBE Menge (symmetrische Differenz 0)':
        ((T(W21, r'Beide Mengen identisch: True  \(symmetrische Differenz: (\d+)\)'),), lambda a: a == 0),
    'Prop. 6.6: in allen Punkten v_p(T_alpha) = 1, und kein Punkt mit p | k':
        ((T(W14, r'^v_p\(T_alpha\) = (1) in allen Faellen: True'), T(W14, r'Index-Boost moeglich, aber Exponent trotzdem 1\): (\d+)')),
         lambda a, b: a == 1 and b == 0),
    'Prop. 6.3: die Tabelle (w140) und die Schluessel turmk1–3 aus §7 nennen fuer jeden Kern denselben ersten moeglichen Index':
        ((W140, T_TURM_P, T_TURM_3), lambda a, p, s: all(turm(p, s, z['m'], 'k') == z['k_min'] for z in a['zeilen'])),
    'Prop. 6.2: unsere Paare bis zum groessten Term = Alekseyevs Tabelle, per Hash, und die Nachbauten gleichen den alten Ausgaben':
        ((W116,), lambda a: a['hash_33_ok'] is True and a['hash_13_ok'] is True and a['vergleich_alt']['kreuz'] == '39/39 Zeilen gleich'),
    'Prop. II.3.2 (C5): Familien mit Quadrat, b ≤ 10^4 — vorhergesagt = in den Daten, keine Abweichung in beide Richtungen':
        ((W112,), lambda a: a['pk2']['vorhergesagt'] == a['pk2']['daten'] and a['pk2']['abw_hier'] == '{}' and a['pk2']['abw_daten'] == '{}'),
    'Prop. II.3.2: kein selbst erzeugtes Paar ohne Quadrat (D ≤ 10^6) fehlt in der Liste; die Glieder ohne Quadrat ueber 2^52 selbst erzeugt':
        ((W115,), lambda a: a['nicht_in_liste'] == [] and all(g['selbst_erzeugt'] for g in a['vergleich_glieder_27_39'])),
}

# ------------------------------------------------------------------ Formats
def fmt(art, x):
    # `fmt(art, x)` formats the value x in the format `art`. The test "float(x) == int(x)" failed for a 22-digit INTEGER
    #   (float is inexact there): integers are accepted directly, only floating-point values must be integral.
    if art == 'ganz':
        assert not isinstance(x, bool) and (isinstance(x, int) or float(x) == int(x)), x; return str(int(x))
    if art == 'ganz_gruppiert':        # from 5 digits on, groups of three with a thin space; 4 digits stay together
        assert not isinstance(x, bool) and (isinstance(x, int) or float(x) == int(x)), x; s = str(abs(int(x)))
        if len(s) >= 5:
            teile = []
            while s: teile.append(s[-3:]); s = s[:-3]
            s = '\\,'.join(reversed(teile))
        return ('-' if int(x) < 0 else '') + s
    if art == 'runden4': return f'{x:.4f}'
    if art == 'abrunden4': return f'{math.floor(x * 10**4) / 10**4:.4f}'
    if art == 'abrunden3': return f'{math.floor(x * 10**3) / 10**3:.3f}'   # lower bounds with 3 digits
    if art == 'abrunden_wiss2':
        e = math.floor(math.log10(x)); m = math.floor(x / 10**e * 10) / 10
        return f'{m:.1f}\\cdot 10^{{{e}}}'
    if art == 'abrunden1':             # one decimal, rounded down (lower bounds)
        return f'{math.floor(x * 10 + 1e-9) / 10:.1f}'
    if art == 'abrunden_wiss3':        # 3 digits, rounded down (lower bounds)
        e = math.floor(math.log10(x)); m = math.floor(x / 10**e * 100) / 100
        return f'{m:.2f}\\cdot 10^{{{e}}}'
    if art == 'zehnerpotenz':          # exact powers of ten only
        e = round(math.log10(x)); assert int(x) == 10**e, x; return f'10^{{{e}}}'
    if art == 'wiss_exakt':            # a*10^e with integer a in 1..9, exact
        e = len(str(int(x))) - 1; a = int(x) // 10**e; assert a * 10**e == int(x) and float(x) == int(x), x
        return f'{a}\\cdot 10^{{{e}}}' if a > 1 else f'10^{{{e}}}'
    if art == 'tex':                   # ready-made math expression made of digits, ^{...} and \cdot - nothing else
        assert re.fullmatch(r'[0-9]+(\^\{[0-9]+\})?(\\cdot [0-9]+(\^\{[0-9]+\})?)*', x), x; return x
    if art == 'runden2': return f'{x:.2f}'
    if art == 'aufrunden_wiss1':       # 1 digit, rounded UP - can be printed as an upper bound
        e = math.floor(math.log10(x)); m = math.ceil(x / 10**e - 1e-12)
        if m >= 10: m, e = 1, e + 1
        assert m * 10**e >= x, (x, m, e); return f'{m}\\cdot 10^{{{e}}}' if m > 1 else f'10^{{{e}}}'
    if art == 'runden_wiss3':          # 3 digits, ROUNDED - print only after "≈", never as a bound
        e = math.floor(math.log10(x)); m = round(x / 10**e, 2)
        if m >= 10: m, e = m / 10, e + 1
        return f'{m:.2f}\\cdot 10^{{{e}}}'
    if art.startswith('$'):            # same form in math mode, for table cells
        return '$' + fmt(art[1:], x) + '$'
    if art == 'wort':                  # status word from a FIXED list - no free text from data files
        return WOERTER[x]
    if art == 'stufen':                # witnesses per stage: "5, 7, 71, 1657; 19, 29, ..."
        assert x and all(st and all(isinstance(v, int) and not isinstance(v, bool) for v in st) for st in x), x
        return '; '.join(', '.join(str(v) for v in st) for st in x)
    if art == 'wortliste_ganz':        # "13, 31, 1546463" as a finished string made of integers
        assert re.fullmatch(r'\d+(, \d+)*', x), x; return x.replace(', ', ',\\ ')
    if art == 'runden3': return f'{x:.3f}'   # three digits (R^2 values)
    if art == 'rangspanne':            # "36--41" or "9" - rank in case of ties as a range
        a, b = x; assert int(a) == a and int(b) == b and 1 <= a <= b, x; return f'{int(a)}' if a == b else f'{int(a)}--{int(b)}'
    if art == 'bruch':                 # "191 of 201" from two integers, numerator <= denominator
        a, b = x; assert int(a) == a and int(b) == b and 0 <= a <= b, x; return f'{int(a)} of {int(b)}'
    if art == 'tex_spanne':            # "24.2--44.9"
        assert re.fullmatch(r'\d+\.\d--\d+\.\d', x), x; return x
    if art == 'wortliste_und':         # "19, 35 and 52" from integers
        assert x is not None and re.fullmatch(r'\d+(\\,\d{3})*(, \d+(\\,\d{3})*)* and \d+(\\,\d{3})*', x), x; return x
    if art == 'liste_ganz':            # "396, 255, 83, 15, 1" - every number an integer
        assert x and all(not isinstance(v, bool) and (isinstance(v, int) or float(v) == int(v)) for v in x), x
        return ', '.join(str(int(v)) for v in x)
    if art == 'runden1': return f'{x:.1f}'
    if art.startswith('oderstrich:'):  # missing value (e.g. no curves at this B1) -> "--"; otherwise the format after the colon
        return '--' if x is None else fmt(art.split(':', 1)[1], x)
    if art == 'ganz_gerundet': return str(round(x))
    if art == 'paar':                  # family (b, d) from two integers
        b, d = x; assert int(b) == b and int(d) == d, x; return f'({int(b)}, {int(d)})'
    # table body: list of rows, each row a list of (format, value); every cell goes through `fmt`
    if art == 'tabzeilen':
        assert x and all(z and all(isinstance(c, tuple) and len(c) == 2 and c[0] != 'tabzeilen' for c in z) for z in x), x
        assert len({len(z) for z in x}) == 1, 'Tabelle: Zeilen verschieden lang'
        return '\n'.join(' & '.join(fmt(a, w) for a, w in z) + ' \\\\' for z in x)
    raise ValueError(art)

# ------------------------------------------------------------------ Text anchors and manifest (pure functions - hence
#   testable with positive controls)
TRENN = '[    ]'
ZAHLFORM = re.compile(rf'^[+-]?\d{{1,3}}(?:{TRENN}?\d{{3}})*(?:\.\d+)?(?:[eE][+-]?\d+)?$')
# `zahl` = parses a number in the allowed format (a comma is rejected); `anker` = value of a pattern that must match exactly
#   once; `alle` = all matches; `manifest_lesen` = reads MANIFEST_sha256.txt into {path: (hash, size)}; `pruef_hash` = compares
#   hash and size of a file with the manifest and returns the first 12 characters of the hash.
def zahl(s):
    s = s.strip()
    if ',' in s: raise ValueError(f'Komma in „{s}" ist mehrdeutig')
    if not ZAHLFORM.match(s): raise ValueError(f'„{s}" ist keine Zahl im erlaubten Format')
    t = re.sub(TRENN, '', s)
    return int(t) if re.fullmatch(r'[+-]?\d+', t) else float(t)
def anker(text, muster):
    treffer = re.findall(muster, text, flags=re.M)
    if len(treffer) != 1: raise ValueError(f'Anker {muster!r}: {len(treffer)} Treffer, verlangt ist genau einer')
    return zahl(treffer[0] if isinstance(treffer[0], str) else treffer[0][0])
def alle(text, muster):
    # TA source: all matches, at least one; one group -> number, several -> tuple of numbers
    treffer = re.findall(muster, text, flags=re.M)
    if not treffer: raise ValueError(f'Muster {muster!r}: kein Treffer')
    return [tuple(zahl(x) for x in t) if isinstance(t, tuple) else zahl(t) for t in treffer]
def manifest_lesen(pfad):
    m = {}
    for z in pfad.read_text(encoding='utf-8').splitlines():
        teile = z.split(None, 2)
        if len(teile) == 3 and re.fullmatch(r'[0-9a-f]{64}', teile[0]):
            m[teile[2].strip().replace('\\', '/')] = (teile[0], int(teile[1]))
    return m
def pruef_hash(name, roh, man):
    soll = man.get('ergebnisse/' + name)
    if soll is None: raise ValueError(f'{name} steht nicht im Manifest')
    ist = (hashlib.sha256(roh).hexdigest(), len(roh))
    if ist != soll: raise ValueError(f'{name} weicht vom Manifest ab (erst pruefen, dann manifest_bauen.py)')
    return ist[0][:12]

# ------------------------------------------------------------------ Positive controls at start
# `wirft` (throws) = True if the call raises ValueError or AssertionError.
def wirft(f, *a):
    try: f(*a)
    except (ValueError, AssertionError): return True
    return False
assert fmt('runden4', 1.68226) == '1.6823' and fmt('abrunden4', 1.68226) == '1.6822' and fmt('ganz', 804) == '804', 'PK Formate'
assert fmt('abrunden_wiss2', 8.41e7) == '8.4\\cdot 10^{7}' and fmt('abrunden_wiss2', 8.499e7) == '8.4\\cdot 10^{7}', 'PK wiss2 rundet ab'
assert (fmt('ganz_gruppiert', 889437) == '889\\,437' and fmt('ganz_gruppiert', 10132126) == '10\\,132\\,126'
        and fmt('ganz_gruppiert', 2027) == '2027' and fmt('ganz_gruppiert', 10000) == '10\\,000'), 'PK Zifferngruppen'
assert (fmt('abrunden_wiss3', 3.8858e20) == '3.88\\cdot 10^{20}' and fmt('zehnerpotenz', 10**8) == '10^{8}' and wirft(fmt, 'zehnerpotenz', 2 * 10**8)
        and fmt('wiss_exakt', 200000) == '2\\cdot 10^{5}' and wirft(fmt, 'wiss_exakt', 210000)
        and fmt('tex', '2^{3}\\cdot 29') == '2^{3}\\cdot 29' and wirft(fmt, 'tex', '2+3')), 'PK neue Formate (29.09.)'
assert faktoren(130576328) == '2^{3}\\cdot 29\\cdot 197\\cdot 2857' and pell(7) == 8 and pell(39) == 25, 'PK Zerlegung/Pell (29.09.)'
assert (alle('a: 5\na: 6', r'^a: (\d+)') == [5, 6] and alle('x 1 2.5', r'x (\d+) (\S+)') == [(1, 2.5)] and wirft(alle, 'b: 1', r'^a: (\d+)')
        and fmt('abrunden1', 17.6) == '17.6' and fmt('abrunden1', 17.69) == '17.6'), 'PK TA-Quelle / abrunden1 (29.09.)'
assert fmt('wort', 'offen') == 'not settled' and fmt('wort', 'voll') == 'complete factorization', 'PK wort'
try:
    fmt('wort', 'irgendwas'); raise AssertionError('PK wort: unbekanntes Wort wurde angenommen')
except KeyError:
    pass
assert fmt('oderstrich:ganz', None) == '--' and fmt('oderstrich:ganz_gruppiert', 2753) == '2753' and wirft(fmt, 'oderstrich:ganz', 1.5), 'PK oderstrich'
assert (fmt('$wiss_exakt', 200000) == '$2\\cdot 10^{5}$' and fmt('stufen', [[5, 7], [19]]) == '5, 7; 19' and wirft(fmt, 'stufen', [[]])
        and wirft(fmt, 'stufen', [[1.5]])), 'PK Mathe-Zellen und Stufen (29.09., Prop. 6.3)'
assert (fmt('liste_ganz', [396, 255, 83, 15, 1]) == '396, 255, 83, 15, 1' and wirft(fmt, 'liste_ganz', [1, 2.5]) and wirft(fmt, 'liste_ganz', [])
        and fmt('runden1', 411.688) == '411.7' and fmt('ganz_gerundet', 652.27) == '652'), 'PK Listen/Rundung (29.09., Prop. 6.1)'
assert (fmt('ganz_gruppiert', 3887785221910670811499) == '3\\,887\\,785\\,221\\,910\\,670\\,811\\,499' and wirft(fmt, 'ganz', 2.5)
        and wirft(fmt, 'ganz', True)), 'PK ganze Zahlen ueber 2^53 (29.09.)'
assert (fmt('runden_wiss3', 3887785221910670811500) == '3.89\\cdot 10^{21}' and fmt('runden_wiss3', 9.996e7) == '1.00\\cdot 10^{8}'
        and fmt('abrunden_wiss3', 3887785221910670811500) == '3.88\\cdot 10^{21}'), 'PK runden_wiss3 (29.09.): gerundet, nicht abgerundet'
assert (fmt('paar', (23, 2)) == '(23, 2)' and wirft(fmt, 'paar', (2.5, 1)) and fmt('runden2', 1.531) == '1.53'
        and fmt('tabzeilen', [[('paar', (2, 1)), ('ganz', 14)], [('paar', (3, 1)), ('ganz', 6)]]) == '(2, 1) & 14 \\\\\n(3, 1) & 6 \\\\'
        and wirft(fmt, 'tabzeilen', [[('ganz', 1)], [('ganz', 1), ('ganz', 2)]]) and wirft(fmt, 'tabzeilen', [[('ganz', 1.5)]])
        and wirft(fmt, 'tabzeilen', [[('tabzeilen', [])]])), 'PK Tabellen-Formate (29.09., Teil II Tab. 1)'
assert abs(zeta(2) - math.pi ** 2 / 6) < 1e-11 and abs(zeta(4) - math.pi ** 4 / 90) < 1e-11, 'PK zeta (exakte Werte, keine Gedaechtniskonstante)'
assert zahl('10 132 126') == 10132126 and zahl('8.41e+07') == 8.41e7 and zahl('0.595') == 0.595 and wirft(zahl, '1,234'), 'PK Zahl-Parser'
assert anker('a: 5\nb: 6', r'^a: (\d+)') == 5 and wirft(anker, 'a: 5\na: 6', r'^a: (\d+)') and wirft(anker, 'b: 6', r'^a: (\d+)'), 'PK Anker'
_h = hashlib.sha256(b'abc').hexdigest()
assert pruef_hash('x.json', b'abc', {'ergebnisse/x.json': (_h, 3)}) == _h[:12], 'PK Hash stimmt'
assert wirft(pruef_hash, 'x.json', b'abd', {'ergebnisse/x.json': (_h, 3)}) and wirft(pruef_hash, 'x.json', b'abc', {}), 'PK Hash weicht ab / fehlt'

# ------------------------------------------------------------------ Build
# `bauen` (build) = checks every source against the manifest, computes every entry of `ZAHLEN` and every theorem check of
#   `PRUEF_SAETZE`, and writes zahlen.tex (one `\@namedef{Z@key}{value}` per number, preceded by a comment with description,
#   format and source hashes).
def bauen():
    man = manifest_lesen(MANIFEST)
    zeilen = ['% zahlen.tex — ERZEUGT von werkzeug/zahlen_bauen.py am ' + time.strftime('%Y-%m-%d %H:%M') + '. NICHT von Hand aendern.',
              '% Je Zahl: Quelle(n) mit sha256 (erste 12 Zeichen, gegen MANIFEST_sha256.txt geprueft), Format; T[...] = Textanker.', '\\makeatletter']
    roh_cache = {}
    def roh(name):
        if name not in roh_cache:
            b = (ERG/name).read_bytes(); roh_cache[name] = (b, pruef_hash(name, b, man))
        return roh_cache[name]
    def lade(quellen):
        werte, herkunft = [], []
        for q in quellen:
            if isinstance(q, T):
                b, h = roh(q.datei); werte.append(anker(b.decode(KODIERUNG.get(q.datei, 'utf-8')), q.muster)); herkunft.append(f'T[{q.datei} ({h})]')
            elif isinstance(q, TA):
                b, h = roh(q.datei); werte.append(alle(b.decode(KODIERUNG.get(q.datei, 'utf-8')), q.muster)); herkunft.append(f'TA[{q.datei} ({h})]')
            else:
                b, h = roh(q); werte.append(json.loads(b.decode('utf-8'))); herkunft.append(f'{q} ({h})')
        return werte, herkunft
    for schl, (quellen, f, art, was) in ZAHLEN.items():
        werte, herkunft = lade(quellen)
        wert = fmt(art, f(*werte))
        zeilen.append(f'% {schl}: {was} · {art} · {" + ".join(dict.fromkeys(herkunft)) or "keine Datei (Literaturkonstante)"}')
        zeilen.append(f'\\@namedef{{Z@{schl}}}{{{wert}}}')
    for aussage, (quellen, f) in PRUEF_SAETZE.items():
        werte, herkunft = lade(quellen)
        if f(*werte) is not True: raise AssertionError(f'Satz-Pruefung gescheitert: {aussage}')
        zeilen.append(f'% Satz-Pruefung ✅ {aussage} · {" + ".join(dict.fromkeys(herkunft))}')
    zeilen.append('\\makeatother')
    ZIEL.write_text('\n'.join(zeilen) + '\n', encoding='utf-8')
    return {k: z for k, z in zip(ZAHLEN, [l.split('}{', 1)[1][:-1] for l in zeilen if l.startswith('\\@namedef')])}

if __name__ == '__main__':
    werte = bauen()
    print('zahlen.tex geschrieben:', ' · '.join(f'{k} = {v}' for k, v in werte.items()))
