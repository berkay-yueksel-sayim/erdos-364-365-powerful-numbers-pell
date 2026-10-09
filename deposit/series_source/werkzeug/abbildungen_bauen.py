# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# abbildungen_bauen.py
# Purpose: builds the figures of the series from the same data as the text.
# Rules (as in zahlen_bauen): no number by hand. Numbers that also appear in the text come from gemeinsam/zahlen.tex (so EXACTLY
#   the printed form); data series come from result files whose sha256 is checked against MANIFEST_sha256.txt. Vector PDF, width =
#   text width (Letter paper, 1.2 in margins, hence 6.1 in), font as in the text (Computer Modern: cmr10 + mathtext "cm"),
#   Okabe-Ito colors (color-blind safe).
# Figures: fig_trichter (Part I section 3, Theorem I.3.1) · fig_leiter (Part I section 6, Prop. I.6.3) · fig_karte (Part I, hero
#   figure; needs w138) · fig_turmwald and fig_treppe (Part II) · fig_kartemn (Part I opener) · fig_luecken (Part III opener).
# Reads:  gemeinsam/zahlen.tex, MANIFEST_sha256.txt and the result files named in the figure functions (the w9_v31 and
#         w9_v32_c3 lists, w110, w125, w133, w138, w139, w140, w152).
# Writes: abbildungen/fig_*.pdf and abbildungen/fig_*.png.
# Usage:  called from bauen.py (after zahlen_bauen), or directly: python abbildungen_bauen.py
#         (needs matplotlib; fig_luecken also needs sympy)
# Controls: assertions in each figure function tie the plotted data to the printed numbers (counts, bounds, known kernels).
import json, math, re, sys, pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
W = pathlib.Path(__file__).resolve().parent; S = W.parent
sys.path.insert(0, str(W)); import zahlen_bauen as zb
ZIEL = S/'abbildungen'   # `ZIEL` = target folder of the figures
BREITE = 6.1             # `BREITE` = figure width in inches (= text width)
# `OKABE` = Okabe-Ito palette; German color names: `blau` blue, `hellblau` light blue, `gruen` green, `zinnober` vermilion,
# `rosa` reddish purple, `gelb` yellow, `schwarz` black, `grau` gray
OKABE = {'blau': '#0072B2', 'hellblau': '#56B4E9', 'gruen': '#009E73', 'orange': '#E69F00', 'zinnober': '#D55E00', 'rosa': '#CC79A7',
         'gelb': '#F0E442', 'schwarz': '#000000', 'grau': '#999999'}
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm', 'axes.formatter.use_mathtext': True,
                     'axes.unicode_minus': False, 'font.size': 10, 'axes.titlesize': 10, 'axes.labelsize': 10, 'legend.fontsize': 10,
                     # with 9 pt, superscripts came out at 5.6 pt (measured in the PDF)
                     'xtick.labelsize': 10, 'ytick.labelsize': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42})
# No "bbox tight": the labels on the left would make the image wider than 6.1 in, LaTeX would scale it down to text width
#   and the font would shrink to about 6 pt. Instead each figure has fixed margins: image width = text width, the font stays
#   at 8-9 pt.

def zahl_tex(schl):   # `zahl_tex`: look up the key `schl` in zahlen.tex
    # the printed form from zahlen.tex: "10\,132\,126" → (display "10 132 126", value 10132126)
    t = (S/'gemeinsam'/'zahlen.tex').read_text(encoding='utf-8')
    # read up to the end of the line: values like "10^{8}" contain braces themselves
    #   (stopping at the first "}" would truncate them)
    m = re.search(r'^\\@namedef\{Z@' + re.escape(schl) + r'\}\{(.*)\}\s*$', t, flags=re.M)
    assert m, f'{schl} fehlt in zahlen.tex (erst zahlen_bauen)'
    # \; instead of \, (the thin space almost vanished); math form (10^{8}, 10\,132\,126) in math mode
    roh = m.group(1); anzeige = '$' + roh.replace('\\,', '\\;') + '$' if re.search(r'[\\^]', roh) else roh
    wert = int(roh.replace('\\,', '')) if re.fullmatch(r'[0-9\\,]+', roh) else None
    return anzeige, wert

def daten(name):   # `daten` = data: load a result file after checking its hash against the manifest
    man = zb.manifest_lesen(zb.MANIFEST); b = (zb.ERG/name).read_bytes(); zb.pruef_hash(name, b, man)
    return json.loads(b.decode('utf-8'))

def fig_trichter():
    # Two stages, honestly separated: KERNELS/FAMILIES (top) and LATTICE POINTS (bottom); the axis is logarithmic and labeled.
    # Wording as in the earlier table in the text (Part I section 3), so that text and figure use the same terms.
    # `zeilen` = rows (key in zahlen.tex, label, color); keys: `kerne` kernels, `filterfamilien` kernels passing the discard
    # criterion, `punktfamilien` kernels with at least one point, `punkte` points, `paritaetspunkte` points excluded by parity,
    # `zeugenpunkte` points excluded by a witness prime, `restpunkte` points remaining.
    zeilen = [('kerne', 'squarefree $m \\equiv 7\\ (\\mathrm{mod}\\ 8)$ with $m \\leq$ ' + zahl_tex('kernschranke')[0], OKABE['hellblau']),
              # "family" means something else in Part II, hence the neutral wording
              ('filterfamilien', 'kernels that pass the discard criterion', OKABE['hellblau']),
              ('punktfamilien', 'kernels with at least one point', OKABE['hellblau']),
              ('punkte', 'points $(m, k)$ with $T_k < 10^{' + zahl_tex('hoehe')[0] + '}$', OKABE['orange']),
              ('paritaetspunkte', 'excluded because $T_k$ is odd', OKABE['grau']),
              ('zeugenpunkte', 'excluded by a prime $p <$ ' + zahl_tex('zeugenschranke')[0] + ' with $v_p(T_k) = 1$', OKABE['gruen'])]
    werte = [(zahl_tex(k), txt, c) for k, txt, c in zeilen]
    rest = zahl_tex('restpunkte')
    assert rest[1] == 0 and werte[4][0][1] + werte[5][0][1] == werte[3][0][1], 'Trichter: Paritaet + Zeugen ≠ Punkte, oder Rest ≠ 0'
    fig, ax = plt.subplots(figsize=(BREITE, 3.2))
    fig.subplots_adjust(left=0.53, right=0.91, bottom=0.16, top=0.99)
    ys = [6, 5, 4, 2.6, 1.6, 0.6]   # `ys` = vertical positions of the bars
    for y, ((anz, w), txt, c) in zip(ys, werte):
        ax.barh(y, math.log10(w), height=0.62, color=c, left=0)
        ax.text(math.log10(w) + 0.08, y, anz, va='center', ha='left', fontsize=10)
        ax.text(-0.15, y, txt, va='center', ha='right', fontsize=10)
    ax.text(-0.15, -0.4, 'remaining', va='center', ha='right', fontsize=10)
    ax.text(0.08, -0.4, f'{rest[0]}', va='center', ha='left', fontsize=10, weight='bold')
    ax.axhline(3.3, color=OKABE['grau'], lw=0.6, ls=(0, (2, 2)))
    ax.text(7.9, 3.45, 'kernels', ha='right', va='bottom', fontsize=9, color='#555555')
    ax.text(7.9, 3.15, 'points', ha='right', va='top', fontsize=9, color='#555555')
    ax.set_xlim(0, 8); ax.set_ylim(-0.9, 6.5); ax.set_yticks([])
    ax.set_xticks(range(0, 9)); ax.set_xticklabels([f'$10^{{{e}}}$' for e in range(0, 9)])
    ax.set_xlabel('count (logarithmic scale)'); ax.spines['left'].set_visible(False)
    fig.savefig(ZIEL/'fig_trichter.pdf'); fig.savefig(ZIEL/'fig_trichter.png', dpi=150); plt.close(fig)
    return 'fig_trichter'

def fig_leiter():
    w140 = daten('w140_tab_reinhart_result.json'); w152 = daten('w152_paritaet_t1_reinhart_kerne_result.json')
    # only kernels with T_1 even can carry a triple (w152)
    zeilen = [z for z in w140['zeilen'] if z['m'] in w152['gerade_klasse7']]
    assert [z['m'] for z in zeilen] == [4099215, 39028039587479], [z['m'] for z in zeilen]
    farben = [OKABE['blau'], OKABE['gruen']]; marker = ['o', '^']
    fig, ax = plt.subplots(figsize=(BREITE, 3.1))
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.15, top=0.98)
    for z, c, mk in zip(zeilen, farben, marker):
        # `ks` = log10 of the smallest possible index k after each stage; `k` = running product of the witnesses
        ks, k = [0.0], 1
        for st in z['stufen_zeugen']:   # `stufen_zeugen` = witnesses per stage of the ladder
            for p in st: k *= p
            ks.append(math.log10(k))
        assert k == z['k_min'] or z['k_min'] >= k, (z['m'], k, z['k_min'])
        ks[-1] = math.log10(z['k_min'])
        e = math.floor(math.log10(z['mitte_stellen'])); a = math.floor(z['mitte_stellen'] / 10**e * 10) / 10
        name = f"{z['m']:,}".replace(',', '\\;')
        ax.step(range(len(ks)), ks, where='post', color=c, lw=1.6)
        ax.plot(range(len(ks)), ks, mk, color=c, ms=4.5, label=f'$m = {name}$: a middle would have at least ${a}\\cdot 10^{{{e}}}$ digits')
    ax.set_xticks(range(4)); ax.set_xlabel('stage of the ladder')
    ax.set_ylabel('$\\log_{10}$ of the smallest possible index $k$')
    ax.legend(frameon=False, loc='upper left'); ax.grid(True, axis='y', alpha=0.25, lw=0.5)
    fig.savefig(ZIEL/'fig_leiter.pdf'); fig.savefig(ZIEL/'fig_leiter.png', dpi=150); plt.close(fig)
    return 'fig_leiter'

def fig_karte():
    # Hero figure of Part I: the 961 building blocks (w138). x = rank d, y = number of digits of M_d (both logarithmic),
    # color + shape = how the block was settled.
    # The palette was tested with a color checker: all checks passed; contrast warning for three light colors, hence a second
    #   encoding (own marker shape per class) and a legend with counts. Runs only once w138 is listed in the manifest,
    #   otherwise it is skipped.
    man = zb.manifest_lesen(zb.MANIFEST)
    if 'ergebnisse/w138_karte_961_daten.json' not in man: return None
    w138 = daten('w138_karte_961_daten.json')
    w133 = daten('w133_probelauf_79_kerne_result.json')
    # Classes EXACTLY like the rows of Table I.6.8: the same function zb.zeile_von, the same wording (without the common prefix
    #   "certified witness:"; the caption states that the first six classes are certified witnesses).
    # `STIL` = style per class: (color, marker, size)
    STIL = {'w_erste': (OKABE['hellblau'], 'o', 10), 'w_sieb': (OKABE['gruen'], 's', 13), 'w_sympy': (OKABE['orange'], 'D', 13),
            'w_gmpecm': (OKABE['zinnober'], '^', 15), 'w_prim': (OKABE['blau'], 'P', 22), 'w_vollbew': (OKABE['schwarz'], '*', 45),
            'w_voll': (OKABE['rosa'], 'v', 15),
            'w_prpzeuge': (OKABE['rosa'], 'h', 26),   # same color as 'w_voll': both rest on probable primes
            # 'w_teile': e.g. (7, 1449), a part of the factorization serves as witness, proved by ECPP
            'w_teile': (OKABE['blau'], 'h', 26)}
    # `r` = the building blocks, one per rank (key `raenge`); `zeile` = class of each block;
    #   `gesamt` = running total of classified blocks
    r = w138['raenge']; zeile = {id(x): zb.zeile_von(x, w133['zeugen']) for x in r}; gesamt = 0
    fig, ax = plt.subplots(figsize=(BREITE, 4.9))
    # legend BELOW the axis: points lay beneath it at the upper left
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.36, top=0.98)
    for w in zb.WEGE:
        if w == 'w_offen': continue
        c, mk, s = STIL[w]; pts = [x for x in r if zeile[id(x)] == w]; gesamt += len(pts)
        if not pts: continue                    # empty classes are left out of the legend
        name = zb.WOERTER[w].replace('certified witness: ', '')
        # a single point vanishes in the dense band: black edge and drawn on top, so that it can be found
        rand = dict(edgecolors=OKABE['schwarz'], linewidths=0.7, zorder=4) if len(pts) == 1 else dict(linewidths=0)
        ax.scatter([x['d'] for x in pts], [x['lg'] + 1 for x in pts], s=s, c=c, marker=mk, alpha=0.85, label=f'{name} ({len(pts)})', **rand)
    offen = [x for x in r if zeile[id(x)] == 'w_offen']; gesamt += len(offen)   # `offen` = blocks that are still open
    ax.scatter([x['d'] for x in offen], [x['lg'] + 1 for x in offen], s=50, c=OKABE['schwarz'], marker='x', linewidths=1.5,
               label=f'open ({len(offen)})', zorder=5)
    wf = [x for x in r if 'wieferich_p' in x]   # `wf` = blocks with a Wieferich event
    ax.scatter([x['d'] for x in wf], [x['lg'] + 1 for x in wf], s=120, facecolors='none', edgecolors=OKABE['schwarz'], linewidths=0.9,
               label=f'Wieferich event ({len(wf)})', zorder=6)
    assert gesamt == w138['n'], ('Karte: nicht jeder Baustein hat genau eine Klasse', gesamt, w138['n'])
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('rank $d$'); ax.set_ylabel('number of digits of $M_d$')
    ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.45, -0.14), ncol=2, handletextpad=0.3, columnspacing=1.0, fontsize=9)
    ax.grid(True, which='major', alpha=0.2, lw=0.5)
    fig.savefig(ZIEL/'fig_karte.pdf'); fig.savefig(ZIEL/'fig_karte.png', dpi=150); plt.close(fig)
    return 'fig_karte'

def fig_turmwald():
    # Hero figure of Part II section 6: each point is a pair at distance two (middle < 10^H, kernel <= 10^8), x = kernel m (log),
    #   y = number of digits of the middle. The rungs of ONE kernel stand vertically above each other at a fixed distance
    #   (Prop. `turmsumme`), hence "tower forest" (`turmwald`). Two classes: color AND marker shape (color-blind safe); the counts
    #   in the legend come from zahlen.tex and are checked against the data.
    # `teile` = parts: (result file, key in zahlen.tex, legend label, color, marker)
    teile = [('w9_v31_M100000000_H2000_result.json', 'zeugenpunkte', 'middle $\\equiv 0\\ (\\mathrm{mod}\\ 4)$', OKABE['blau'], 'o'),
             ('w9_v32_c3_M100000000_H2000_result.json', 'kdpaare', 'middle $\\equiv 2\\ (\\mathrm{mod}\\ 4)$', OKABE['zinnober'], 's')]
    fig, ax = plt.subplots(figsize=(BREITE, 3.9))
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.27, top=0.98)
    hmax = zahl_tex('hoehe')[1]   # `hmax` = height bound: maximal number of digits of the middle
    for datei, schl, name, c, mk in teile:
        # `paare` = the pairs, i.e. lattice points with a witness
        a = daten(datei); paare = [x for x in a['points'] if x[4] == 'witness']
        anz, wert = zahl_tex(schl); assert wert == len(paare), (schl, wert, len(paare))
        assert all(x[3] <= hmax for x in paare)
        ax.scatter([x[0] for x in paare], [x[3] for x in paare], s=5, c=c, marker=mk, linewidths=0, alpha=0.8, label=f'{name} ({anz} pairs)')
    ax.set_xscale('log'); ax.set_ylim(0, hmax * 1.02)
    ax.set_xlabel('kernel $m$'); ax.set_ylabel('number of digits of the middle')
    ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.45, -0.17), ncol=2, handletextpad=0.3, markerscale=2.5)
    ax.grid(True, which='major', alpha=0.2, lw=0.5)
    fig.savefig(ZIEL/'fig_turmwald.pdf'); fig.savefig(ZIEL/'fig_turmwald.png', dpi=150); plt.close(fig)
    return 'fig_turmwald'

def fig_treppe():
    # Hero figure of Part II: the staircase of the pairs n, n + 1 up to n₀, each step in the color of its family.
    #   Positions come from OUR OWN computation (w139: log10 n₁ and log10 F per family, cross-checked against w112/w115),
    #   not from OEIS values.
    #   Controls: (1) each family has as many steps as w139 counts, together `\Z{fampaare}`; (2) at 10^8, 10^10, …, 10^20
    #   the staircase has exactly the counts that w110 obtains from the data (comparison of outputs); (3) the comparison
    #   function V (a heuristic, section 5) is drawn only at the support points of w110.
    # `w139['zeilen']` = rows, one per family: `b`, `d` = class pair, `anzahl` = number of pairs,
    #   `log10_n1` = log10 of the first pair, `log10_F` = log10 of the step factor
    w139 = daten('w139_tab1_familien_result.json'); w110 = daten('w110_paare_bis_1e21_result.json')
    paare = []   # `paare` = all pairs as (log10 n, family (b, d)), sorted by size
    for z in w139['zeilen']:
        for j in range(z['anzahl']):
            paare.append((z['log10_n1'] + j * z['log10_F'], (z['b'], z['d'])))
    paare.sort()
    assert len(paare) == zahl_tex('fampaare')[1], (len(paare), zahl_tex('fampaare'))
    for v in w110['vergleich']:
        if not v['X'].isdigit(): continue
        e = len(v['X']) - 1; assert v['X'] == '1' + '0' * e
        assert sum(1 for lg, _ in paare if lg + 1e-12 < e) == v['gemessen'], ('Treppe ≠ w110 bei 10^%d' % e, v['gemessen'])
    # Rule without ties: families with at least THREE pairs get a color (ranked by count, then by position)
    rang = sorted([z for z in w139['zeilen'] if z['anzahl'] >= 3], key=lambda z: (-z['anzahl'], z['log10_n1']))
    farbe = {}   # `farbe` = color: family (b, d) -> (color, number of pairs)
    for z, c in zip(rang, [OKABE['blau'], OKABE['zinnober'], OKABE['gruen'], OKABE['orange'], OKABE['rosa']]):
        farbe[(z['b'], z['d'])] = (c, z['anzahl'])
    fig, ax = plt.subplots(figsize=(BREITE, 4.4))
    # 7 legend entries in 2 columns = 4 rows (room is needed for the 4th row)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.34, top=0.98)
    # `xs` = positions (log10 n); `ys` = running count of pairs
    xs = [lg for lg, _ in paare]; ys = list(range(1, len(paare) + 1))
    ax.step([0] + xs + [xs[-1] + 0.4], [0] + ys + [ys[-1]], where='post', color=OKABE['grau'], lw=1.0, zorder=1)
    for fam, (c, anz) in farbe.items():
        pts = [(x, y) for (x, f), y in zip(paare, ys) if f == fam]
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=26, c=c, linewidths=0, zorder=3, label=f'family $({fam[0]}, {fam[1]})$: {anz}')
    rest = [(x, y) for (x, f), y in zip(paare, ys) if f not in farbe]   # `rest` = pairs of the families without their own color
    ax.scatter([p[0] for p in rest], [p[1] for p in rest], s=22, facecolors='white', edgecolors=OKABE['schwarz'], linewidths=0.8, zorder=3,
               label=f'the other {len(w139["zeilen"]) - len(farbe)} families: {len(rest)}')
    # `vx`, `vy` = support points of the comparison function V (exponent of 10, corrected value from w110)
    vx = [len(v['X']) - 1 for v in w110['vergleich'] if v['X'].isdigit()]; vy = [v['korrigiert'] for v in w110['vergleich'] if v['X'].isdigit()]
    ax.plot(vx, vy, ls=(0, (4, 3)), color=OKABE['schwarz'], lw=1.0, zorder=2, label='comparison function $V$ (a heuristic, Section 5)')
    ax.set_xlim(0, 22.5); ax.set_ylim(0, len(paare) + 2)
    ax.set_xticks(range(0, 23, 2)); ax.set_xticklabels([f'$10^{{{e}}}$' for e in range(0, 23, 2)])
    ax.set_xlabel('$x$ (logarithmic scale)'); ax.set_ylabel('pairs $n, n + 1$ with $n + 1 \\leq x$')
    ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.45, -0.15), ncol=2, handletextpad=0.3, columnspacing=1.0, fontsize=9)
    ax.grid(True, which='major', alpha=0.2, lw=0.5)
    fig.savefig(ZIEL/'fig_treppe.pdf'); fig.savefig(ZIEL/'fig_treppe.png', dpi=150); plt.close(fig)
    return 'fig_treppe'

# Constant from the literature: Reinhart's list of the kernels with m | U_1 is complete up to this bound
# [Reinhart 2023, Rem. 5.5]; it is typed in the text as there.
REINHART = 1.5e12

def fig_kartemn():
    # Opener figure of Part I: where a triple could still lie at all. x = log10 m (kernel), y = log10 n (middle) on a
    #   logarithmic axis. ONLY proved or computed bounds; the statement numbers are in the caption:
    #   Cor. I.2.3 (n > m^{3/2}; n > m^{4.5}/27 if m does not divide U_1) · Reinhart's complete list up to 1.5·10^12 with w152
    #   (class 7, T_1 even: only one kernel below 1.5·10^12, and it lies in the box) · Prop. I.6.2(d) (no middle below
    #   `\Z{ohnetripelsicher}`) · Theorem I.3.1 (box).
    #   Controls: (1.5·10^12)^{3/2} gives the printed bound; every known kernel below 1.5·10^12 lies in the box; the known one
    #   above lies in the open zone. The figure does NOT draw the regions of Part I section 5 (classes of kernels); the
    #   caption says so.
    import numpy as np
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    # `kb_anz` = display string of the kernel bound (`kernschranke`), `xb` = its decimal exponent; `hoehe` = height bound;
    # `ot_anz` = display string of the bound below which no middle exists (`ohnetripelsicher`), `y_ot` = its log10
    kb_anz, _ = zahl_tex('kernschranke'); xb = int(re.search(r'10\^\{(\d+)\}', kb_anz).group(1))
    hoehe = zahl_tex('hoehe')[1]
    ot_anz, _ = zahl_tex('ohnetripelsicher')
    mm = re.search(r'([0-9.]+)\\cdot 10\^\{(\d+)\}', ot_anz); y_ot = math.log10(float(mm.group(1))) + int(mm.group(2))
    xr = math.log10(REINHART)   # `xr` = log10 of Reinhart's bound
    assert abs(1.5 * xr - y_ot) < 0.01, ('PK: (1,5e12)^1,5 ≠ gedruckte Schranke', 1.5 * xr, y_ot)
    w140 = daten('w140_tab_reinhart_result.json'); w152 = daten('w152_paritaet_t1_reinhart_kerne_result.json')
    gerade7 = sorted(w152['gerade_klasse7'])   # `gerade7` = class-7 kernels with T_1 even (key `gerade_klasse7`)
    assert all(math.log10(m) <= xb for m in gerade7 if m <= REINHART), ('PK: ein Kern mit m | U_1, T_1 gerade unter 1,5e12 liegt nicht in der Box', gerade7)
    # `ueber` = such kernels above Reinhart's bound
    ueber = [m for m in gerade7 if m > REINHART]; assert ueber, 'PK: kein bekannter Kern ueber 1,5e12 (erwartet: einer)'
    # `stellen` = number of digits of the smallest possible middle
    stellen = {z['m']: z['mitte_stellen'] for z in w140['zeilen']}
    x = np.linspace(0, 20, 1601); y15 = 1.5 * x; y45 = 4.5 * x - math.log10(27)
    # 4.6 in was too tall together with the caption for the upper page area (the figure slid to the next page)
    fig, ax = plt.subplots(figsize=(BREITE, 4.3))
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.43, top=0.98)
    # gray, light gray, green; green is opaque (a semi-transparent green left a dark stripe where it overlapped the gray)
    GRAU, HELL, GRUEN = '#BDBDBD', '#E3E3E3', '#A6DCC8'
    # impossible for every kernel
    ax.fill_between(x, 1, np.maximum(y15, y_ot), color=GRAU, lw=0, zorder=1)
    zw = (x > xb) & (x <= xr)   # `zw` = zone between the box edge and Reinhart's bound
    ax.fill_between(x, np.maximum(y15, y_ot), np.maximum(y45, y_ot), where=zw, facecolor=HELL, edgecolor='#8C8C8C', hatch='////', lw=0, zorder=1)
    of = x > xr   # `of` = open zone beyond Reinhart's bound
    ax.fill_between(x, y15, y45, where=of, color=OKABE['orange'], alpha=0.55, lw=0, zorder=1)
    ax.fill_between([0, xb], 1, hoehe, color=GRUEN, lw=0, zorder=2)
    ax.plot(x, y15, color='#555555', lw=0.8, zorder=3); ax.plot(x, np.maximum(y45, 1), color='#555555', lw=0.8, ls=(0, (3, 2)), zorder=3)
    ax.text(18.6, 1.5 * 18.6 * 0.80, '$n = m^{3/2}$', ha='center', va='top', fontsize=9, color='#333333')
    ax.text(15.2, (4.5 * 15.2 - math.log10(27)) * 1.18, '$n = m^{4.5}/27$', ha='right', va='bottom', fontsize=9, color='#333333')
    for m in ueber:
        xm = math.log10(m)
        assert stellen[m] > 1e4, ('PK: der Kern liegt nicht „far beyond the picture" (obere Bildkante: 10^10000)', m, stellen[m])
        ax.plot([xm, xm], [1.5 * xm, 1e4], color=OKABE['blau'], lw=1.6, zorder=4)
        ax.annotate('', xy=(xm, 1e4), xytext=(xm, 3e3), arrowprops=dict(arrowstyle='-|>', color=OKABE['blau'], lw=1.2), zorder=4)
    ax.text(4.0, 60, 'computed:\nno point left', ha='center', va='center', fontsize=10, color='black', zorder=5)
    ax.text(15.5, 4.0, 'excluded for every kernel', ha='center', va='center', fontsize=10, zorder=5)
    # Orange means "open: possible", so the sentence must include it
    ax.text(4.0, 4400, 'a triple can only lie in the\nwhite or orange region', ha='center', va='center', fontsize=10, style='italic', zorder=5)
    ax.set_xlim(0, 20); ax.set_yscale('log'); ax.set_ylim(1, 1e4)
    ax.set_xticks(range(0, 21, 2)); ax.set_xticklabels([f'$10^{{{e}}}$' for e in range(0, 21, 2)])
    ax.set_yticks([1, 10, 100, 1000, 10000]); ax.set_yticklabels(['$10^{1}$', '$10^{10}$', '$10^{100}$', '$10^{1000}$', '$10^{10000}$'])
    ax.set_xlabel('kernel $m$'); ax.set_ylabel('middle $n$')
    hand = [Patch(facecolor=GRAU, label='excluded for every kernel: $n \\leq m^{3/2}$, or $n < ' + ot_anz.strip('$') + '$'),
            Patch(facecolor=GRUEN, label='$m \\leq ' + kb_anz.strip('$') + f'$, $n < 10^{{{hoehe}}}$: every point excluded by computation'),
            Patch(facecolor=HELL, edgecolor='#8C8C8C', hatch='////', label='needs $m \\mid U_1(m)$ with $T_1$ even: no such kernel in this range'),
            Patch(facecolor=OKABE['orange'], alpha=0.55, label='open: possible only for kernels with $m \\mid U_1(m)$'),
            Line2D([0], [0], color=OKABE['blau'], lw=1.6, label='the known kernel $m \\equiv 7$ (mod 8) of this kind above $1.5\\cdot 10^{12}$: excluded far beyond the picture')]
    ax.legend(handles=hand, frameon=False, loc='upper left', bbox_to_anchor=(-0.09, -0.15), ncol=1, handletextpad=0.5, fontsize=9)
    fig.savefig(ZIEL/'fig_kartemn.pdf'); fig.savefig(ZIEL/'fig_kartemn.png', dpi=150); plt.close(fig)
    return 'fig_kartemn'

# `ist_powerful` = is powerful: trial division, True if every prime exponent of n is at least 2 (n = 1 counts)
def ist_powerful(n):
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0: n //= p; e += 1
            if e == 1: return False
        p += 1
    return n == 1

def fig_luecken():
    # Opener figure of Part III: the head of the gap spectrum up to 10^14 (w125; `kopf60` = the 60 most frequent values with their
    #   occurrence counts). Powerfulness is tested here (trial division) and held against the printed keys: the first `\Z{lsbis}`
    #   are all powerful, the first other one is `\Z{lserster}` with `\Z{lsocc}` occurrences. Color AND pattern distinguish the
    #   classes (color-blind and print safe).
    from matplotlib.patches import Patch
    w125 = daten('w125_spektrumkopf_gleichstaende_result.json')
    # `e14` = result for the range up to 10^14; `kopf` = head of the spectrum as (gap value, occurrences); `pw` = powerful flags
    e14 = w125['ergebnis']['100000000000000']; kopf = e14['kopf60']
    pw = [ist_powerful(v) for v, _ in kopf]
    bis = zahl_tex('lsbis')[1]; erster = zahl_tex('lserster')[1]; occ = zahl_tex('lsocc')[1]
    # Ties (w125): values with the same frequency share ranks; their order within a tie is by value and means nothing.
    #   Taking rank `lsbis` + 1 as "the first non-powerful value" would be wrong (63 900 shares the ranks 49-53), so the
    #   check follows the statement of the text:
    ueber = [p for (v, c), p in zip(kopf, pw) if c > occ]
    assert len(ueber) == bis and all(ueber), ('PK: Werte mit mehr als \\Z{lsocc} Vorkommen ≠ \\Z{lsbis} powerful', len(ueber), ueber)
    gruppe = [(v, p) for (v, c), p in zip(kopf, pw) if c == occ]
    assert [v for v, p in gruppe if not p] == [erster], ('PK: in der Gruppe mit \\Z{lsocc} Vorkommen ist nicht genau \\Z{lserster} nicht-powerful', gruppe)
    spanne = e14['erster_nicht_pw']['rang_spanne']; assert spanne == [bis + 1, bis + len(gruppe)], (spanne, bis, len(gruppe))
    anteil = zahl_tex('lsanteil')[0]; werte = zahl_tex('lswerte')[0]
    fig, ax = plt.subplots(figsize=(BREITE, 3.6))
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.30, top=0.97)
    for i, ((v, c), p) in enumerate(zip(kopf, pw), start=1):
        if p: ax.bar(i, c, width=0.8, color=OKABE['blau'], lw=0)
        else: ax.bar(i, c, width=0.8, facecolor='white', edgecolor=OKABE['zinnober'], hatch='////', lw=0.9)
    v0, c0 = kopf[0]
    from sympy import factorint   # the factorization is computed, not typed
    fak = ' \\cdot '.join(f'{p}^{{{e}}}' if e > 1 else f'{p}' for p, e in sorted(factorint(v0).items()))
    ax.annotate(f'${v0:,}$'.replace(',', '\\;') + f' $= {fak}$', xy=(1, c0), xytext=(6, c0 * 0.98),
                fontsize=9, va='center', arrowprops=dict(arrowstyle='-', color='#555555', lw=0.6))
    i_erster = next(i for i, (v, c) in enumerate(kopf, start=1) if v == erster)
    ax.annotate(f'first value that is not powerful: ${erster:,}$'.replace(',', '\\;') + f' (ranks {spanne[0]} to {spanne[1]}, a tie)',
                xy=(i_erster, occ), xytext=(len(kopf), occ + 56), ha='right', fontsize=9, va='center',
                arrowprops=dict(arrowstyle='-', color='#555555', lw=0.6))
    ax.set_xlim(0, len(kopf) + 1); ax.set_ylim(0, 205); ax.set_xlabel('the 60 most frequent gap values up to $10^{14}$, by rank')
    ax.set_ylabel('occurrences')
    hand = [Patch(facecolor=OKABE['blau'], label='powerful gap value'),
            Patch(facecolor='white', edgecolor=OKABE['zinnober'], hatch='////', label='gap value that is not powerful')]
    ax.legend(handles=hand, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=9)
    ax.text(len(kopf) / 2 + 0.5, 200, f'among all {werte} gap values up to $10^{{14}}$, only {anteil}% are powerful', ha='center', va='top', fontsize=9.5,
            weight='bold')
    ax.grid(True, axis='y', alpha=0.25, lw=0.5)
    fig.savefig(ZIEL/'fig_luecken.pdf'); fig.savefig(ZIEL/'fig_luecken.png', dpi=150); plt.close(fig)
    return 'fig_luecken'

def bauen():   # `bauen` = build: all figures; returns the names of the figures that were built
    ZIEL.mkdir(exist_ok=True)
    fertig = [fig_trichter(), fig_leiter(), fig_karte(), fig_turmwald(), fig_treppe(), fig_kartemn(), fig_luecken()]
    return [f for f in fertig if f]

if __name__ == '__main__':
    print('Abbildungen:', ', '.join(bauen()))
