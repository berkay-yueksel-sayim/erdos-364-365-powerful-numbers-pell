# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# bauen.py
# Purpose: builds the four PDFs (Parts I, II, III and the volume) from ONE source and checks them.
# Usage:
#   python bauen.py      generate the numbers → build the parts (3 passes, so that references across part boundaries resolve)
#                        → build the volume → run the checks
#   python bauen.py pk   the same, then EVERY check is run against a planted error and must report it (positive controls)
# Reads:  the sources under series_source/ (teil1 to teil3/inhalt.tex, band/main_band.tex, gemeinsam/*.tex, literatur.bib), the
#         .aux tables and the PDFs it builds; runs pdflatex and bibtex; calls zahlen_bauen, abbildungen_bauen, belege_bauen,
#         probleme_bauen (and konkordanz_bauen if present).
# Writes: gemeinsam/zahlen.tex, figures, evidence-type tables, the PDFs in ausgabe/ and the build log ausgabe/BUILD_PROTOKOLL.txt.
# Checks (pure functions over source texts, .aux tables and PDF text, so they can be fed planted errors without a rebuild):
#   Z  every \Z{…} in the text is defined in zahlen.tex (LaTeX itself also aborts; the control checks that as well)
#   V  every \tref{x} has a \tlabel{x} in its own part, every \partref{P}{x} a \tlabel{x} in part P; plus the LaTeX log
#   G  no dash character (— or ---) in the English text (comments excluded)
#   P  no e-mail address other than the allowed contact address, in the sources AND in the text of the PDFs
#   L  every label <Part>:x carries the same number in the part build and in the volume (one number in every cover)
#   N  no direct \Phi / Φ in the text: primitive part = \Prim{d}, cyclotomic polynomial = \Kreis{n} (symbol defined in ONE place)
#   B  every entry in literatur.bib names its first-hand source in the comment lines directly above it
#   A  American English: no British form from the patterns below. Only these patterns are checked (-ise/-isation, -our, -yse,
#      -lled, -re, whilst/amongst, programme); a British word outside them passes. Negative control in `pk`: Paper 1 (v1.3,
#      American English) must stay silent, otherwise an exception is too narrow.
#   F  no figure or table after the heading "References" (a too-tall float had silently pushed ALL floats of Part I to the end;
#      LaTeX does not report this as an error; forward references in the text are normal).
import re, sys, shutil, subprocess, time, pathlib, copy
sys.stdout.reconfigure(encoding='utf-8')
# `W` = folder of this tool; `S` = series_source folder; `AUS` = output folder `ausgabe` (the PDFs and the build log)
W = pathlib.Path(__file__).resolve().parent; S = W.parent; AUS = S/'ausgabe'
TEILE = [('I', 'teil1'), ('II', 'teil2'), ('III', 'teil3')]   # `TEILE` = parts: (part label, folder)
ERLAUBT = {'berksa@tutamail.com'}   # `ERLAUBT` = allowed e-mail addresses (check P)
prot = []   # `prot` = protocol: all output lines, written to the build log at the end
def sag(s=''):   # `sag` = "say": print a line and keep it in the protocol
    print(s, flush=True); prot.append(s)

# ------------------------------------------------------------------ build
def pdflatex(ordner, name):   # run pdflatex on `name`.tex in the folder `ordner`; returns (exit code, log text)
    r = subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error', '-file-line-error', name + '.tex'],
                       cwd=S/ordner, capture_output=True, text=True, errors='replace')
    log = (S/ordner/(name + '.log')).read_text(encoding='utf-8', errors='replace') if (S/ordner/(name + '.log')).exists() else r.stdout
    return r.returncode, log
def bibtex(ordner, name):   # run bibtex on `name` in the folder `ordner`; returns (exit code, output)
    r = subprocess.run(['bibtex', name], cwd=S/ordner, capture_output=True, text=True, errors='replace')
    return r.returncode, r.stdout
def bauen_alles():   # `bauen_alles` = build everything; returns (numbers, logs), or (None, logs) if a LaTeX step fails
    sys.path.insert(0, str(W)); import zahlen_bauen
    # `werte` = the computed numbers
    werte = zahlen_bauen.bauen(); sag('Zahlen: ' + ' · '.join(f'{k} = {v}' for k, v in werte.items()))
    import abbildungen_bauen   # figures from the same data, after the numbers (they read zahlen.tex)
    sag('Abbildungen: ' + ', '.join(abbildungen_bauen.bauen()))
    logs = {}   # `logs` = LaTeX log of each folder
    # evidence-type tables: built before pass 1 from the source text (the numbers exist only after the .aux files), afterwards
    #   rebuilt after every pass, comparing the LaTeX route (.bel) with the source-text route
    import belege_bauen
    belege_bauen.bauen(pruefen=False)
    # probleme_def.tex (list of open problems) must exist before pass 1
    import probleme_bauen; probleme_bauen.bauen(pruefen=False)
    for durchgang in (1, 2, 3):   # `durchgang` = pass
        for t, o in TEILE:
            rc, log = pdflatex(o, f'main_{o}'); logs[o] = log
            if rc: sag(f'🔴 pdflatex {o} (Durchgang {durchgang}) Exit {rc}'); return None, logs
        sag('Belegarten: ' + ' · '.join(belege_bauen.bauen()))
        if durchgang == 1:
            for t, o in TEILE:
                rc, out = bibtex(o, f'main_{o}')
                if rc: sag(f'🔴 bibtex {o} Exit {rc}: {out[-300:]}'); return None, logs
    # list of open problems, now with the numbers from the finished .aux files of the parts
    sag('Offene Probleme: ' + probleme_bauen.bauen())
    for schritt in ('latex', 'bibtex', 'latex', 'latex'):   # `schritt` = step
        rc, log = pdflatex('band', 'main_band') if schritt == 'latex' else bibtex('band', 'main_band')
        if schritt == 'latex': logs['band'] = log
        if rc: sag(f'🔴 {schritt} band Exit {rc}'); return None, logs
    AUS.mkdir(exist_ok=True)
    for o in [o for _, o in TEILE] + ['band']:
        shutil.copyfile(S/o/f'main_{o}.pdf', AUS/f'{o}.pdf')
    return werte, logs

# ------------------------------------------------------------------ inputs of the checks
def quellen():   # `quellen` = sources: {file: (text, part label or None)}
    q = {f'{o}/inhalt.tex': ((S/o/'inhalt.tex').read_text(encoding='utf-8'), t) for t, o in TEILE}
    q['band/main_band.tex'] = ((S/'band'/'main_band.tex').read_text(encoding='utf-8'), None)
    q['gemeinsam/serie.tex'] = ((S/'gemeinsam'/'serie.tex').read_text(encoding='utf-8'), None)
    return q
def aux_tabelle(pfad):   # label -> number table from an .aux file (empty if the file does not exist)
    if not pfad.exists(): return {}
    t = pfad.read_text(encoding='utf-8', errors='replace')
    return {m.group(1): m.group(2) for m in re.finditer(r'\\newlabel\{([^}]+)\}\{\{([^}]*)\}', t)}
def pdf_seiten():   # `pdf_seiten` = PDF pages: {file: [text of each page]} for all PDFs in ausgabe/
    import pypdf
    return {f'ausgabe/{p.name}': [(s.extract_text() or '') for s in pypdf.PdfReader(str(p)).pages] for p in sorted(AUS.glob('*.pdf'))}
def pdf_texte(seiten=None):   # `pdf_texte` = PDF texts: the pages of each PDF joined into one text
    return {d: '\n'.join(s) for d, s in (seiten or pdf_seiten()).items()}
def pruef_F(seiten):
    """The signature of a float jam: a caption "Figure k:" / "Table k:" on a page AFTER the heading "References"
    (forward references in the text are normal, a float behind the bibliography is not)."""
    f = []
    for d, ss in seiten.items():
        lit = next((i for i, s in enumerate(ss, 1) if re.search(r'(?m)^References\s*$', s)), None)
        if lit is None: continue
        for i, s in enumerate(ss, 1):
            if i <= lit: continue
            for m in re.finditer(r'(?m)^(Figure|Table) (\d+):', s): f.append(f'{d}: {m.group(1)} {m.group(2)} auf S. {i}, hinter References (S. {lit})')
    return f
def ohne_kommentare(text):   # `ohne_kommentare` = without comments: the LaTeX text with all % comments removed
    return '\n'.join(re.sub(r'(?<!\\)%.*', '', z) for z in text.splitlines())

# ------------------------------------------------------------------ checks (each returns a list of findings, empty = clean)
def pruef_Z(q, zahlen_tex):
    definiert = set(re.findall(r'\\@namedef\{Z@([^}]+)\}', zahlen_tex))
    return sorted({f'{d}: \\Z{{{k}}} ist nicht definiert' for d, (txt, _) in q.items()
                   for k in re.findall(r'\\Z\{([^}]+)\}', ohne_kommentare(txt)) if k not in definiert})
def pruef_V(q, logs=None):
    labels = {}
    for d, (txt, t) in q.items():
        if t: labels.setdefault(t, set()).update(re.findall(r'\\tlabel\{([^}]+)\}', ohne_kommentare(txt)))
    f = []
    for d, (txt, t) in q.items():
        s = ohne_kommentare(txt)
        if t:
            f += [f'{d}: \\tref{{{x}}} ohne \\tlabel in Teil {t}' for x in re.findall(r'\\tref\{([^}]+)\}', s) if x not in labels.get(t, set())]
        f += [f'{d}: \\partref{{{p}}}{{{x}}} ohne \\tlabel in Teil {p}' for p, x in re.findall(r'\\partref\{([IV]+)\}\{([^}]+)\}', s)
              if x not in labels.get(p, set())]
    for o, log in (logs or {}).items():
        if re.search(r'undefined references|Reference `[^\']*\' on page \d+ undefined', log): f.append(f'{o}: LaTeX meldet undefinierte Verweise')
        # doubly defined labels: the second \label silently overrides the first, so every reference to it points to the
        # wrong place. Bibliography keys report the same warning harmlessly (xr reads the .aux files of the other parts,
        # the own one wins); they carry no ':'.
        doppelt = sorted({k for k in re.findall(r"Label `([^']+)' multiply defined", log) if ':' in k})
        if doppelt: f.append(f'{o}: LaTeX meldet doppelt definierte Labels: {doppelt[:5]}')
    zaehl = {}
    for d, (txt, t) in q.items():
        if t:
            for x in re.findall(r'\\tlabel\{([^}]+)\}', ohne_kommentare(txt)): zaehl.setdefault((t, x), []).append(d)
    f += [f'\\tlabel{{{x}}} in Teil {t} {len(ds)}-mal definiert ({", ".join(sorted(set(ds)))})' for (t, x), ds in sorted(zaehl.items()) if len(ds) > 1]
    return f
def pruef_G(q):
    return [f'{d}:{i}: Gedankenstrich' for d, (txt, _) in q.items()
            for i, z in enumerate(ohne_kommentare(txt).splitlines(), 1) if '\u2014' in z or '---' in z]
def pruef_P(texte):
    return [f'{d}: Adresse {a}' for d, txt in texte.items()
            for a in re.findall(r'[\w.+-]+@[\w-]+(?:\.[\w-]+)+', txt) if a.lower() not in ERLAUBT]
def pruef_L(teil_tab, band_tab):
    f = []
    for t, tab in teil_tab.items():
        for lab, nr in tab.items():
            if not lab.startswith(t + ':'): continue
            if lab not in band_tab: f.append(f'{lab}: fehlt im Band')
            elif band_tab[lab] != nr: f.append(f'{lab}: Teil {nr} ≠ Band {band_tab[lab]}')
    return f

def pruef_N(q):
    # Notation: primitive part = \Prim{d}, cyclotomic polynomial = \Kreis{n}; a direct \Phi (or Φ) in the text would silently
    # undo the separation. The preamble (where the macros are defined) is not part of q.
    return [f'{d}:{i}: \\Phi direkt — \\Prim{{d}} oder \\Kreis{{n}} verwenden' for d, (txt, _) in q.items()
            for i, z in enumerate(ohne_kommentare(txt).splitlines(), 1) if re.search(r'\\Phi(?![A-Za-z])', z) or 'Φ' in z]

def pruef_B(bib):
    # Bibliography source check: every entry in literatur.bib names its first-hand source in the comment lines directly above it
    # (`% quelle: …`, quelle = source). An entry without a source may have been written from memory or copied from another list.
    f, kommentar = [], []
    for z in bib.splitlines():
        s = z.strip()
        if s.startswith('%'):
            kommentar.append(s); continue
        m = re.match(r'@(\w+)\s*\{\s*([^,\s]+)\s*,', s)
        if m and not any(k.startswith('% quelle:') for k in kommentar):
            f.append(f'literatur.bib: {m.group(2)} ohne `% quelle:`-Zeile')
        kommentar = []
    return f

BRE = re.compile(r'\b(\w+is(?:e|es|ed|ing|ation|ations)|\w+our(?:s|ed|ing|ite|ites|able|ably|al|ally|hood|hoods|ly|less|ful)?|'
                 r'analys(?:e|ed|ing)|paralys(?:e|ed|ing)|(?:modell|labell|travell|cancell|totall|signall|levell|equall|fuell|channell|counsell|tunnell)\w*|'
                 r'centre\w*|fibre\w*|metre|metres|litre\w*|whilst|amongst|programme\w*)\b', re.I)
# `BRE` = patterns for British spellings; `BRE_AUSNAHMEN` = exceptions: words that are correct American English but touched by a
# pattern; each exception is a genuine American English word, not a suppression
BRE_AUSNAHMEN = {
    'equally', 'totally',   # the patterns equall*/totall* target equalled/totalled; the adverbs are American English
    'otherwise', 'precise', 'precisely', 'raise', 'raised', 'raises', 'raising', 'rise', 'rises', 'arise', 'arises', 'arising', 'promise', 'promised',
    'promises', 'promising', 'exercise', 'exercised', 'exercises', 'exercising', 'surprise', 'surprised', 'surprises', 'surprising', 'noise',
    'wise', 'likewise', 'concise', 'premise', 'premises', 'praise', 'advise', 'advised', 'advises', 'advising', 'revise', 'revised', 'revises',
    'revising', 'devise', 'devised', 'devises', 'comprise', 'comprises', 'comprised', 'comprising', 'compromise', 'compromised', 'enterprise',
    'expertise', 'supervise', 'supervised', 'televise', 'improvise', 'improvisation', 'disguise', 'disguised', 'franchise', 'merchandise', 'poise',
    'cruise', 'vise', 'paradise', 'treatise', 'demise', 'reprise', 'crises', 'clockwise', 'pairwise', 'termwise', 'stepwise', 'pointwise',
    'coordinatewise', 'elementwise', 'bitwise', 'componentwise', 'entrywise', 'rowwise', 'digitwise', 'blockwise',
    'our', 'ours', 'your', 'yours', 'four', 'fours', 'hour', 'hours', 'hourly', 'tour', 'tours', 'pour', 'poured', 'pouring', 'flour', 'sour',
    'devour', 'contour', 'contours', 'detour', 'detours', 'scour', 'dour', 'amour', 'paramour', 'velour', 'troubadour'}
def pruef_A(q):
    # drop the arguments of reference/number/citation commands (labels can carry German words, for example `…weise`)
    def text(s): return re.sub(r'\\(?:tlabel|tref|label|ref|eqref|cite|citep|citet|nocite|Z|input|include|usepackage|partref\{[IV]+\})\{[^}]*\}', ' ', ohne_kommentare(s))
    return [f'{d}:{i}: britisch „{w}"' for d, (txt, _) in q.items()
            for i, z in enumerate(text(txt).splitlines(), 1) for w in BRE.findall(z) if w.lower() not in BRE_AUSNAHMEN]

def alle_pruefungen(q, zahlen_tex, logs, teil_tab, band_tab, texte, seiten=None):   # run all checks, {check letter: findings}
    return {'Z': pruef_Z(q, zahlen_tex), 'V': pruef_V(q, logs), 'G': pruef_G(q),
            'P': pruef_P({**{d: txt for d, (txt, _) in q.items()}, **texte}), 'L': pruef_L(teil_tab, band_tab), 'N': pruef_N(q),
            'B': pruef_B((S/'gemeinsam'/'literatur.bib').read_text(encoding='utf-8')), 'A': pruef_A(q), 'F': pruef_F(seiten or {})}

# ------------------------------------------------------------------ positive controls
def pk(q, zahlen_tex, teil_tab, band_tab):   # `pk` = positive controls: plant one error per check, each must be reported
    ok = True
    def erwarte(name, fehler, muss):   # `erwarte` = expect: the findings `fehler` must contain a finding with the text `muss`
        nonlocal ok
        treffer = [x for x in fehler if muss in x]
        sag(f'   PK {name}: {"✅ gemeldet" if treffer else "🔴 NICHT gemeldet"} — {treffer[:2]}'); ok &= bool(treffer)
    q1 = dict(q); q1['PK/z.tex'] = ('\\Z{gibtsnicht}', 'I'); erwarte('Z unbekannte Zahl', pruef_Z(q1, zahlen_tex), 'gibtsnicht')
    import zahlen_bauen as zb_pk
    schl = ", ".join(f"'k{i}': {i}" for i in range(60))
    erwarte('Z doppelter Schluessel', zb_pk.pruefe_doppelte_schluessel("d = {" + schl + ", 'k7': 0}"), "'k7'")
    nk_dup = zb_pk.pruefe_doppelte_schluessel("d = {" + schl + "}"); sag(f'   NK Z ohne Doppelung: {"✅ still" if not nk_dup else "🔴 " + str(nk_dup)}'); ok &= not nk_dup
    q2 = dict(q); q2['PK/v.tex'] = ('\\tref{weg} \\partref{II}{weg2}', 'I')
    f2 = pruef_V(q2); erwarte('V toter \\tref', f2, '\\tref{weg}'); erwarte('V toter \\partref', f2, '{weg2}')
    lab_ex = next((x for (t, x) in sorted({(t, x) for d, (txt, t) in q.items() if t for x in re.findall(r'\\tlabel\{([^}]+)\}', ohne_kommentare(txt))})), None)
    q2b = dict(q); q2b['PK/dup.tex'] = ('\\tlabel{' + lab_ex + '}', next(t for d, (txt, t) in q.items() if t and ('\\tlabel{' + lab_ex + '}') in txt))
    erwarte('V doppeltes \\tlabel', pruef_V(q2b), '-mal definiert')
    erwarte('V LaTeX: doppelt definiertes Label', pruef_V(q, {'PK': "LaTeX Warning: Label `I:thm:x' multiply defined."}), 'doppelt definierte Labels')
    nkv = pruef_V(q, {'NK': "LaTeX Warning: Label `GMPECM706' multiply defined."}); sag(f'   NK V Literaturschluessel-Warnung: {"✅ still" if not nkv else "🔴 " + str(nkv)}'); ok &= not nkv
    q3 = dict(q); q3['PK/g.tex'] = ('an em dash \u2014 here', None); erwarte('G Gedankenstrich', pruef_G(q3), 'PK/g.tex')
    erwarte('P fremde Adresse', pruef_P({'PK': 'write to x@beispiel.invalid'}), 'x@beispiel.invalid')
    nk = pruef_P({'NK': 'contact berksa@tutamail.com'}); sag(f'   NK P erlaubte Adresse: {"✅ still" if not nk else "🔴 " + str(nk)}'); ok &= not nk
    stau = {'PK/stau.pdf': ['see Figure 1 here', 'text', 'References\n[1] x', 'Figure 1: a caption at the end']}; erwarte('F Float-Stau', pruef_F(stau), 'Figure 1 auf S. 4')
    nk = pruef_F({'NK/ok.pdf': ['see Figure 1 and Table 2', 'Figure 1: caption', 'Table 2: caption', 'References\n[1] x']}); sag(f'   NK F Floats vor References: {"✅ still" if not nk else "🔴 " + str(nk)}'); ok &= not nk
    b2 = dict(band_tab); lab = next((l for l in b2 if l.startswith('II:')), None)
    if lab: b2[lab] = '9.9.9'
    erwarte('L Nummer Teil ≠ Band', pruef_L(teil_tab, b2), lab or '—')
    q4 = dict(q); q4['PK/n.tex'] = ('the primitive part $\\Phi_{d}$ and $Φ_5$', 'I'); f4 = pruef_N(q4)
    erwarte('N direktes \\Phi', f4, 'PK/n.tex'); ok &= sum('PK/n.tex' in x for x in f4) == 1
    nk4 = pruef_N({'NK': ('$\\Prim{d}$, $\\Kreis{12}$, $\\varphi(4d)$, \\PhiX', None)}); sag(f'   NK N Makros: {"✅ still" if not nk4 else "🔴 " + str(nk4)}'); ok &= not nk4
    bib_pk = '% quelle: https://beispiel.invalid\n@article{MitQuelle,\n  title = {x}\n}\n% nur ein Kommentar\n@misc{OhneQuelle,\n  title = {y}\n}\n'
    f5 = pruef_B(bib_pk); erwarte('B Eintrag ohne Quelle', f5, 'OhneQuelle'); ok &= not any('MitQuelle' in x for x in f5)
    sag(f'   NK B Eintrag mit Quelle: {"✅ still" if not any("MitQuelle" in x for x in f5) else "🔴 gemeldet"}')
    q6 = {'PK/a.tex': ('the behaviour of the labelled factorisation, whilst in the centre; its rigour', 'I')}; f6 = pruef_A(q6)
    erwarte('A britische Form', f6, 'behaviour'); ok &= len(f6) == 6
    sag(f'   PK A alle sechs eingepflanzten Formen: {"✅" if len(f6) == 6 else "🔴 " + str(len(f6))}')
    p1 = S.parent/'p2_a1_familie'/'release_v13'/'arxiv_upload'/'main.tex'
    if p1.exists():
        nk6 = pruef_A({'Paper1/main.tex': (p1.read_text(encoding='utf-8'), None)})
        sag(f'   NK A Paper 1 (amerikanisch): {"✅ still" if not nk6 else "🔴 " + str(nk6[:4])}'); ok &= not nk6
    else:
        sag(f'   NK A Paper 1: 🔴 Datei fehlt ({p1}) — Negativkontrolle NICHT gelaufen'); ok = False
    # real LaTeX test: an unknown key must abort the build
    probe = S/'teil1'/'pk_zahl.tex'
    probe.write_text('\\documentclass{article}\\newcommand{\\teilnr}{I}\\input{../gemeinsam/praeambel}\\begin{document}\\Z{gibtsnicht}\\end{document}\n', encoding='utf-8')
    rc, log = pdflatex('teil1', 'pk_zahl')
    for suffix in ('.tex', '.log', '.aux', '.pdf', '.out'):
        (S/'teil1'/f'pk_zahl{suffix}').unlink(missing_ok=True)
    flach = log.replace('\n', '')          # TeX wraps log lines after 79 characters: join them first, then search
    i = flach.find('Package serie Error')
    abbruch = rc != 0 and 'nicht in zahlen.tex definiert' in flach
    sag(f'   PK LaTeX bricht bei unbekannter Zahl ab: {"✅" if abbruch else "🔴"} (Exit {rc}) — Log: {flach[i:i + 90] if i >= 0 else "keine serie-Fehlermeldung gefunden"}')
    ok &= abbruch
    return ok

if __name__ == '__main__':
    t0 = time.time(); sag('=' * 90); sag('bauen.py — vier PDFs aus einer Quelle   ' + time.strftime('%Y-%m-%d %H:%M')); sag('=' * 90)
    werte, logs = bauen_alles()
    if werte is None:
        (AUS/'BUILD_PROTOKOLL.txt').parent.mkdir(exist_ok=True); (AUS/'BUILD_PROTOKOLL.txt').write_text('\n'.join(prot) + '\n', encoding='utf-8'); sys.exit(1)
    q = quellen(); zahlen_tex = (S/'gemeinsam'/'zahlen.tex').read_text(encoding='utf-8')
    teil_tab = {t: aux_tabelle(S/o/f'main_{o}.aux') for t, o in TEILE}; band_tab = aux_tabelle(S/'band'/'main_band.aux')
    seiten = pdf_seiten(); texte = pdf_texte(seiten)
    erg = alle_pruefungen(q, zahlen_tex, logs, teil_tab, band_tab, texte, seiten)
    sag('Pruefungen:')
    for k, name in (('Z', 'Zahlen aus den Daten'), ('V', 'Verweise'), ('G', 'Gedankenstriche'), ('P', 'Adressen (PII)'), ('L', 'Nummer Teil = Band'),
                    ('N', 'Notation \\Prim/\\Kreis'), ('B', 'Bib-Eintraege mit Quelle erster Hand'), ('A', 'American English (Muster)'),
                    ('F', 'kein Float hinter References')):
        sag(f'   {k} {name}: {"✅" if not erg[k] else "🔴 " + str(erg[k][:4])}')
    sag('   Labels je Teil: ' + ' · '.join(f'{t}: {sum(l.startswith(t + ":") for l in teil_tab[t])}' for t, _ in TEILE)
        + f' · Band: {len(band_tab)}')
    sag('   Beispiele: ' + ' · '.join(f'{l} = {band_tab.get(l)}' for l in ('I:thm:hoehe', 'II:prop:csechs', 'III:rem:probe')))
    gruen = not any(erg.values())   # `gruen` = green: no check reported anything
    if len(sys.argv) > 1 and sys.argv[1] == 'pk':
        sag('Positiv-Kontrollen (je ein eingepflanzter Fehler):'); gruen &= pk(q, zahlen_tex, teil_tab, band_tab)
        import probleme_bauen; sag('   Liste der offenen Probleme:'); gruen &= probleme_bauen.pk() == 0
    sag(f'{"BUILD GRUEN" if gruen else "BUILD ROT"} — PDFs in ausgabe/ · {time.time() - t0:.0f} s')
    if gruen and (W/'konkordanz_bauen.py').exists():
        k = subprocess.run([sys.executable, str(W/'konkordanz_bauen.py')], capture_output=True, text=True, encoding='utf-8')
        sag('   Konkordanz: ' + (k.stdout.splitlines()[0] if k.returncode == 0 and k.stdout else '🔴 ' + (k.stderr or '')[-200:]))
    (AUS/'BUILD_PROTOKOLL.txt').write_text('\n'.join(prot) + '\n', encoding='utf-8')
    sys.exit(0 if gruen else 1)
