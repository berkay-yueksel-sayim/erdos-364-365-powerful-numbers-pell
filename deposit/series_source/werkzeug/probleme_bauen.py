# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# probleme_bauen.py
# Purpose: builds the LIST OF OPEN PROBLEMS (Problem A, B, C, ...) of the volume from ONE source into two outputs, so that the
#   volume and the deposit file cannot drift apart (nobody writes them by hand). Source: the `question` environments of the three
#   parts (wording, title, label; the number comes from the .aux file of the part).
# Reads: <part>/inhalt.tex and the .aux file of each part, gemeinsam/zahlen.tex, the files checked by `problemvon_pruefen`.
# Writes: gemeinsam/probleme.tex, gemeinsam/probleme_def.tex, ausgabe/OPEN_PROBLEMS.md (relative to the folder above werkzeug/).
# Usage: python probleme_bauen.py (build) or python probleme_bauen.py pk (run the controls).
# Controls (pk): a question without a label, a \Z without a value and an unresolved \tref must make the build abort; the number of
#   problems must equal the number of question environments of the evidence-type route; an unknown \problemvon must be reported.

# Output 1: gemeinsam/probleme.tex — read by the volume in its framework introduction (\tref becomes \ref{<part>:...}, otherwise
#   the wording).
# Output 2: ausgabe/OPEN_PROBLEMS.md — the deposit file (LaTeX math stays $...$; \tref/\partref become numbers, \cite becomes
#   [key]).
# The letters run in the order of the parts: A, B, ... The discussion (what is known, what an answer would bring) stands ONLY
# in the part behind the question; both outputs point there. So there is exactly one place for every sentence.
import re, sys, pathlib, datetime
sys.stdout.reconfigure(encoding='utf-8')
W = pathlib.Path(__file__).resolve().parent; S = W.parent   # `W` = this folder (werkzeug), `S` = the folder above it
sys.path.insert(0, str(W)); import belege_bauen as bb   # `belege_bauen` = evidence builder (shared tool module of the build)
TEILE = bb.TEILE   # `TEILE` = parts: list of (part name, folder)
# `TEILNAME` = display names of the parts
TEILNAME = {'I': 'Part I (triples)', 'II': 'Part II (pairs)', 'III': 'Part III (further observations)'}

# `fragen` = questions
def fragen(S=S):
    """[(`teil`, `buchstabe`, `nummer`, `titel`, `label`, `koerper`, `abschnitt`)] = (part, letter, number, title, label, body, section)
    in document order; and a list of errors."""
    aus, fehler, i = [], [], 0   # `aus` = result list, `fehler` = errors, `i` = running index of the problem letter
    for teil, ordner in TEILE:   # `ordner` = folder of the part
        # `nummern` = numbers by label from the .aux file
        t = bb.ohne_kommentare((S/ordner/'inhalt.tex').read_text(encoding='utf-8')); nummern = bb.aux_nummern(ordner)
        for m in re.finditer(r'\\begin\{question\}', t):
            j = m.end(); titel, j = bb.optionales(t, j)
            ende = t.find('\\end{question}', j); koerper = t[j:ende]
            lab = re.search(r'\\tlabel\{([^}]*)\}', koerper); lab = lab.group(1) if lab else None
            koerper = re.sub(r'\\tlabel\{[^}]*\}|\\belegart\{[^}]*\}', '', koerper)
            koerper = re.sub(r'\s+', ' ', koerper).strip()
            nr = nummern.get(f'{teil}:{lab}') if lab else None
            if not lab: fehler.append(f'{ordner}: question „{titel[:40]}" ohne \\tlabel')
            if not titel.strip(): fehler.append(f'{ordner}: question {lab} ohne Titel')
            if lab and not nr: fehler.append(f'{ordner}: question {lab} hat keine Nummer in der .aux (erst bauen)')
            abschnitt = nr.split('.')[1] if nr else '?'
            aus.append((teil, chr(ord('A') + i), nr, titel.strip(), lab, koerper, abschnitt)); i += 1
    return aus, fehler

# `zahlen` = numbers: the \Z key -> value table read from gemeinsam/zahlen.tex (empty if the file does not exist)
def zahlen(S=S):
    p = S/'gemeinsam'/'zahlen.tex'
    if not p.exists(): return {}
    return dict(re.findall(r'\\@namedef\{Z@([^}]*)\}\{(.*)\}$', p.read_text(encoding='utf-8'), re.M))

# `tex_fassung` = TeX version: the content of gemeinsam/probleme.tex for the list `liste`, optionally preceded by a notation block
def tex_fassung(liste, notation=''):
    z = ['% probleme.tex — ERZEUGT von werkzeug/probleme_bauen.py. NICHT von Hand aendern; Quelle sind die question-Umgebungen der Teile.']
    if notation: z.append(notation)
    for teil, b, nr, titel, lab, koerper, abschnitt in liste:
        k = re.sub(r'\\tref\{([^}]*)\}', lambda m: f'\\ref{{{teil}:{m.group(1)}}}', koerper)
        z.append(f'\\par\\medskip\\noindent\\textbf{{Problem {b}}} ({titel}; Question~\\ref{{{teil}:{lab}}} in Part~{teil}).\\enspace {k}')
    return '\n'.join(z) + '\n'

# `UMLAUT` = LaTeX accent commands -> Unicode characters (for the Markdown output)
UMLAUT = {'{\\"o}': 'ö', '{\\"a}': 'ä', '{\\"u}': 'ü', '\\"o': 'ö', '\\"a': 'ä', '\\"u': 'ü', '{\\H{o}}': 'ő', '\\H{o}': 'ő', '\\ss{}': 'ß'}

# `md_text` = Markdown text of one question body: resolves \tref, \partref, \Z and \cite; returns (text, list of errors)
def md_text(koerper, teil, alle_nummern, Z):
    k = koerper
    fehler = []
    def tref(m):
        nr = alle_nummern.get(f'{teil}:{m.group(1)}')
        if not nr: fehler.append(f'\\tref{{{m.group(1)}}} in Teil {teil} unaufgeloest')
        return nr or '??'
    def partref(m):
        nr = alle_nummern.get(f'{m.group(1)}:{m.group(2)}')
        if not nr: fehler.append(f'\\partref{{{m.group(1)}}}{{{m.group(2)}}} unaufgeloest')
        return nr or '??'
    def zz(m):
        if m.group(1) not in Z: fehler.append(f'\\Z{{{m.group(1)}}} ohne Wert'); return '??'
        return Z[m.group(1)]
    k = re.sub(r'\\tref\{([^}]*)\}', tref, k)
    k = re.sub(r'\\partref\{([^}]*)\}\{([^}]*)\}', partref, k)
    k = re.sub(r'\\Z\{([^}]*)\}', zz, k)
    k = re.sub(r'\\cite\[([^\]]*)\]\{([^}]*)\}', lambda m: f'[{m.group(2)}, {m.group(1)}]', k)
    k = re.sub(r'\\cite\{([^}]*)\}', lambda m: f'[{m.group(1)}]', k)
    k = re.sub(r'\\emph\{([^}]*)\}', r'*\1*', k)
    k = k.replace('\\Prim{', 'M_{')
    for a, b in UMLAUT.items(): k = k.replace(a, b)
    k = k.replace('~', ' ').replace('\\S', '§').replace('\\%', '%').replace('\\&', '&')
    return k, fehler

# `md_fassung` = Markdown version: the content of ausgabe/OPEN_PROBLEMS.md (header, notation block, one section per problem);
#   returns (text, errors)
def md_fassung(liste, alle_nummern, Z):
    heute = datetime.date.today().isoformat()   # `heute` = today's date (written into the generated file)
    # `teile_mit` / `teile_ohne` = parts with / without questions
    teile_mit = [t for t, _ in TEILE if any(x[0] == t for x in liste)]; teile_ohne = [t for t, _ in TEILE if t not in teile_mit]
    # `woher` = where the questions come from
    woher = 'Parts ' + ' and '.join(teile_mit) + (f' (Part {" and ".join(teile_ohne)} states none)' if teile_ohne else '')
    z = ['# Open problems of the series', '',
         f'Generated on {heute} by `werkzeug/probleme_bauen.py` from the question environments of {woher}; the letters are those of the',
         'framework introduction of the collected volume, and the wording is the wording of the parts. The discussion of each problem, that is,',
         'what is known and what an answer would give, follows its statement in the part named. Mathematics is written in LaTeX between `$`;',
         'a bracket such as `[Reinhart2024, p. 6]` is a key of the bibliography that all parts share. The questions of Erdős themselves,',
         'whether three consecutive powerful numbers exist and whether the number of pairs below $x$ is bounded by a power of $\\log x$, are',
         'the questions the parts start from and are not repeated here.', '',
         '## Notation used below', '',
         '- A positive integer is *powerful* if every prime that divides it divides it at least twice.',
         '- For a squarefree $m \\ge 2$, $\\varepsilon_m = T_1 + U_1 \\sqrt{m}$ is the fundamental solution of $x^2 - m y^2 = 1$, and',
         '  $\\varepsilon_m^k = T_k + U_k \\sqrt{m}$; $T_k(m)$, $U_k(m)$ name the kernel $m$. $m\'$ is the product of the primes of $m$ that do not divide $U_1$.',
         '- The *rank* $\\alpha_m(p)$ of a prime $p \\nmid m$ is the smallest $k$ with $p \\mid T_k(m)$; $M_d$ is the part of $T_d(m)$ supported on the',
         '  primes of rank exactly $d$ (the *primitive part*).',
         '- *Wing (ii)* of Part I: the kernels $m \\equiv 7 \\pmod 8$ with $m \\mid U_1(m)$; they are the weak point of the lower bound for a triple.',
         '- In Part II, the *families* are Walker\'s Pell families of pairs $n, n + 1$, $F_f$ is the growth factor from one pair of the family $f$',
         '  to the next,',
         '  $w(m) = 1/(2m\' \\log_{10} \\varepsilon_m)$ is the rate of the tower of the kernel $m$, $P_2(N)$ counts the pairs $n - 1, n + 1$ of powerful',
         '  numbers with $n \\le N$, and $E_3(N)$ counts the odd squarefree $m$ with $m \\mid U_1(m)$, $T_1(m)$ even and $T_1(m) \\le N$.',
         '- Reinhart\'s condition (C) is Definition 1.4 of [Reinhart2024].', '']
    fehler = []
    for teil, b, nr, titel, lab, koerper, abschnitt in liste:
        k, f = md_text(koerper, teil, alle_nummern, Z); fehler += f
        z += [f'## Problem {b}: {titel}', '', f'*Stated as Question {nr} in {TEILNAME[teil]}, Section {abschnitt}.*', '', k, '']
    return '\n'.join(z), fehler

# `notation_aus_md` = notation from the Markdown: the same notation block as in OPEN_PROBLEMS.md, converted to LaTeX from THIS
#   source
def notation_aus_md(md):
    # (no second hand-written version, which would drift). Italic *x* becomes \\emph{x}, [key] becomes \\cite{key}.
    t = md.split('## Notation used below', 1)[1].split('\n## ', 1)[0]
    items = []
    for zeile in t.split('\n'):   # `zeile` = line
        if zeile.startswith('- '): items.append(zeile[2:])
        elif zeile.startswith('  ') and items: items[-1] += ' ' + zeile.strip()
    assert len(items) >= 5, ('Notationsblock nicht gefunden', len(items))
    def tex(s):
        s = re.sub(r'\*([^*$]+)\*', lambda m: '\\emph{' + m.group(1) + '}', s)
        return re.sub(r'\[([A-Za-z]+[0-9]{4}[a-z]?)\]', lambda m: '\\cite{' + m.group(1) + '}', s)
    return '\\paragraph{Notation used in the problems.}\n\\begin{itemize}\n' + '\n'.join('\\item ' + tex(i) for i in items) + '\n\\end{itemize}'

# `bauen` = build: collect the questions, write the three output files, and return a one-line summary (the build aborts on errors
#   if `pruefen`)
def bauen(pruefen=True):
    liste, fehler = fragen()
    alle_nummern = {}   # all numbers by "<part>:<label>"
    for teil, ordner in TEILE: alle_nummern.update(bb.aux_nummern(ordner))
    md, f2 = md_fassung(liste, alle_nummern, zahlen()); fehler += f2
    fehler += problemvon_pruefen(liste)
    if fehler and pruefen:
        for f in fehler: print('❌', f)
        raise SystemExit(f'probleme_bauen: {len(fehler)} Fehler')
    (S/'gemeinsam'/'probleme.tex').write_text(tex_fassung(liste, notation_aus_md(md)), encoding='utf-8')
    # in addition, each problem's letter via its label: \problemvon{II:q:ezwei} gives H. A text that means a SPECIFIC problem
    # names it this way and never via \problemeletzter (the letter of the last problem shifts as soon as a part gains a question).
    # `zeilen` = lines of probleme_def.tex
    zeilen = ['% ERZEUGT von werkzeug/probleme_bauen.py: letzter Buchstabe, Anzahl, und der Buchstabe je Label',
              f'\\providecommand{{\\problemeletzter}}{{{liste[-1][1]}}}\\providecommand{{\\problemeanzahl}}{{{len(liste)}}}',
              '\\providecommand{\\problemvon}[1]{\\ifcsname prob@#1\\endcsname\\csname prob@#1\\endcsname\\else\\textbf{??}\\fi}']
    zeilen += [f'\\expandafter\\def\\csname prob@{teil}:{lab}\\endcsname{{{b}}}' for teil, b, nr, titel, lab, koerper, abschnitt in liste]
    (S/'gemeinsam'/'probleme_def.tex').write_text('\n'.join(zeilen) + '\n', encoding='utf-8')
    (S/'ausgabe').mkdir(exist_ok=True); (S/'ausgabe'/'OPEN_PROBLEMS.md').write_text(md, encoding='utf-8')
    return f'{len(liste)} Probleme ' + ''.join(b for _, b, *_ in liste) + ' (' + ', '.join(f'{t}: {sum(1 for x in liste if x[0] == t)}' for t, _ in TEILE) + ')'

# `problemvon_pruefen` = check the \problemvon references: returns the list of errors (empty if all are fine)
def problemvon_pruefen(liste, texte=None):
    # every \problemvon{<part>:<label>} in the sources must point to a question of the list (otherwise LaTeX prints a bold "??")
    bekannt = {f'{t}:{lab}' for t, b, nr, ti, lab, k, a in liste}   # `bekannt` = known labels
    if texte is None:
        texte = {str(q.relative_to(S)): q.read_text(encoding='utf-8') for q in [S/'band'/'rahmen.tex', S/'band'/'main_band.tex', S/'gemeinsam'/'serie.tex']
                 + [S/o/'inhalt.tex' for _, o in TEILE] + [S/o/f'main_{o}.tex' for _, o in TEILE]}
    return [f'{d}: \\problemvon{{{x}}} zeigt auf keine question' for d, t in texte.items() for x in re.findall(r'\\problemvon\{([^}]*)\}', t) if x not in bekannt]

def pk():
    import tempfile, shutil
    liste, fehler = fragen(); print(f'  PK Bestand: {len(liste)} Fragen, {len(fehler)} Fehler'); ok = not fehler
    alle = {}
    for teil, ordner in TEILE: alle.update(bb.aux_nummern(ordner))
    # (1) a \tref to a foreign label must be reported
    _, f = md_text('see Corollary~\\tref{cor:gibtesnicht}', 'I', alle, zahlen()); print(f'  PK unaufgeloester \\tref gemeldet: {"✅" if f else "❌"}'); ok &= bool(f)
    # (2) a \Z without a value must be reported
    _, f = md_text('\\Z{gibtesnicht}', 'I', alle, zahlen()); print(f'  PK \\Z ohne Wert gemeldet: {"✅" if f else "❌"}'); ok &= bool(f)
    # (3) a question without a label, in a copy, must be reported
    tmp = pathlib.Path(tempfile.mkdtemp())
    for teil, ordner in TEILE:
        (tmp/ordner).mkdir(); src = S/ordner
        for name in ('inhalt.tex', f'main_{ordner}.aux'):
            if (src/name).exists(): shutil.copy(src/name, tmp/ordner/name)
    t = (tmp/'teil1'/'inhalt.tex').read_text(encoding='utf-8'); t = t.replace('\\tlabel{q:walsh}', '', 1)
    (tmp/'teil1'/'inhalt.tex').write_text(t, encoding='utf-8')
    _, f = fragen(tmp); print(f'  PK question ohne Label gemeldet: {"✅" if any("ohne" in x for x in f) else "❌"}'); ok &= any('ohne' in x for x in f)
    shutil.rmtree(tmp, ignore_errors=True)
    # (4) cross-check: as many problems as question environments with evidence type open (route of the evidence-type table)
    n_open = sum(1 for teil, ordner in TEILE for env, *_ in bb.route_b(teil, ordner)[0] if env == 'question')
    print(f'  PK Anzahl = questions der Belegart-Route: {"✅" if n_open == len(liste) else "❌"} ({n_open} = {len(liste)})'); ok &= n_open == len(liste)
    # (5) a \problemvon with an unknown label must be reported; the actual files must be silent
    f5 = problemvon_pruefen(liste, {'PK': 'see Problem~\\problemvon{II:q:gibtsnicht}'}); f5n = problemvon_pruefen(liste)
    print(f'  PK \\problemvon unbekannt gemeldet: {"✅" if f5 else "❌"} · Bestand still: {"✅" if not f5n else "❌ " + str(f5n[:2])}'); ok &= bool(f5) and not f5n
    print('PK', '✅ alle' if ok else '❌'); return 0 if ok else 1

if __name__ == '__main__':
    if sys.argv[1:] == ['pk']: sys.exit(pk())
    print('Offene Probleme:', bauen())
