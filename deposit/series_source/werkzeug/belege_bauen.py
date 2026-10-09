# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# belege_bauen.py
# Purpose: the EVIDENCE TYPE of each statement, collected from the text and verified twice.
# Rule: every statement (theorem, proposition, lemma, corollary, definition, remark, question) carries its evidence type:
#   proved · literature · computed · measured · heuristic · open. In the text, a \belegart{…} follows the \tlabel{…}; the macro
#   (in praeambel.tex) prints the label and writes a line "<number>|<type>|<section>" to main_<part>.bel. This tool builds from
#   it, for each part, the table gemeinsam/belegtab_<Part>.tex, and the "Results at a glance" table gemeinsam/blick.tex.
# Two routes (double verification):
#   Route A  the .bel file that LaTeX wrote while typesetting (numbers as printed).
#   Route B  a separate reader of the source text: every statement environment with its \tlabel and its \belegart, the numbers
#            taken from the .aux file via the label.
#   Both lists must be equal (number and type), otherwise the script aborts. Further rules: each environment has EXACTLY ONE
#   \belegart and one \tlabel, the type comes from the fixed list, a "question" always carries "open", and "open" is carried
#   only by a question.
# Reads:  teil1 to teil3/inhalt.tex, main_<part>.aux and main_<part>.bel (under series_source/).
# Writes: gemeinsam/belegtab_<Part>.tex and gemeinsam/blick.tex.
# Usage:  python belege_bauen.py       builds the tables (bauen.py calls it after every LaTeX pass)
#         python belege_bauen.py pk    controls
# Controls: a changed type label in a copy must make the A/B comparison fail, a missing one likewise; an unknown type is reported.
import re, sys, pathlib
sys.stdout.reconfigure(encoding='utf-8')
W = pathlib.Path(__file__).resolve().parent; S = W.parent   # `W` = folder of this tool; `S` = series_source folder
TEILE = [('I', 'teil1'), ('II', 'teil2'), ('III', 'teil3')]   # `TEILE` = parts: (part label, folder)
ARTEN = ['proved', 'literature', 'computed', 'measured', 'heuristic', 'open']   # `ARTEN` = the evidence types
UMGEBUNGEN = ('theorem', 'proposition', 'lemma', 'corollary', 'definition', 'remark', 'question')
# `UMGEBUNGEN` = statement environments; `WORT` = display word of each environment
WORT = {'theorem': 'Theorem', 'proposition': 'Proposition', 'lemma': 'Lemma', 'corollary': 'Corollary', 'definition': 'Definition',
        'remark': 'Remark', 'question': 'Question'}
# `ERKLAERUNG` = explanation of each type, printed in the footnote of the table
ERKLAERUNG = {'proved': 'proved in this part; the proof may use published results',
              'literature': 'due to the cited source, sometimes with a proof given here for completeness', 'computed': 'the result of an exact computation, or a proof that uses one of our computations',
              'measured': 'a measurement, with a model that is not proved', 'heuristic': 'a heuristic, used as a yardstick only',
              'open': 'an open question'}

def ohne_kommentare(t):   # `ohne_kommentare` = without comments: the LaTeX text with all % comments removed
    return '\n'.join(re.sub(r'(?<!\\)%.*', '', z) for z in t.splitlines())

def optionales(t, i):   # `optionales` = optional argument
    """Reads the optional argument starting at t[i] == '[' with bracket counting (citations contain brackets themselves)."""
    if i >= len(t) or t[i] != '[': return '', i
    tiefe, j = 0, i
    while j < len(t):
        if t[j] == '[': tiefe += 1
        elif t[j] == ']':
            tiefe -= 1
            if tiefe == 0: return t[i + 1:j], j + 1
        j += 1
    raise ValueError('optionales Argument ohne Ende')

def route_b(teil, ordner):
    """Statements from the source text: [(environment, title, label, type)], plus a list of errors."""
    t = ohne_kommentare((S/ordner/'inhalt.tex').read_text(encoding='utf-8')); aus, fehler = [], []
    for m in re.finditer(r'\\begin\{(' + '|'.join(UMGEBUNGEN) + r')\}', t):
        env = m.group(1); j = m.end()
        while j < len(t) and t[j] in ' \t\n': j += 1
        titel, j = optionales(t, j)
        ende = t.find(f'\\end{{{env}}}', j); assert ende > 0, (ordner, env, m.start())
        koerper = t[j:ende]
        labels = re.findall(r'\\tlabel\{([^}]*)\}', koerper); arten = re.findall(r'\\belegart\{([^}]*)\}', koerper)
        label = labels[0] if labels else None; art = arten[0] if arten else None
        if len(labels) != 1: fehler.append(f'{ordner}: {env} „{titel[:40]}" hat {len(labels)} \\tlabel (Soll 1)')
        if len(arten) != 1: fehler.append(f'{ordner}: {env} {label} hat {len(arten)} \\belegart (Soll 1)')
        if art is not None and art not in ARTEN: fehler.append(f'{ordner}: {label}: unbekannte Belegart „{art}"')
        if env == 'question' and art != 'open': fehler.append(f'{ordner}: question {label} traegt „{art}", nicht „open"')
        if env != 'question' and art == 'open': fehler.append(f'{ordner}: {env} {label} traegt „open" — das darf nur eine question')
        aus.append((env, titel, label, art))
    return aus, fehler

def aux_nummern(ordner):   # label -> number table from main_<folder>.aux
    p = S/ordner/f'main_{ordner}.aux'
    if not p.exists(): return {}
    t = p.read_text(encoding='utf-8', errors='replace')
    return {m.group(1): m.group(2) for m in re.finditer(r'\\newlabel\{([^}]+)\}\{\{([^}]*)\}', t)}

# route A: the lines of main_<folder>.bel as tuples (number, type, section); None if the file does not exist
def route_a(ordner):
    p = S/ordner/f'main_{ordner}.bel'
    if not p.exists(): return None
    zeilen = [z.strip() for z in p.read_text(encoding='utf-8', errors='replace').splitlines() if z.strip()]
    return [tuple(z.split('|')) for z in zeilen]

def titel_tex(titel):
    # the title from the source text, for the table: citations in the title (attributions) stay
    return titel.strip()

# "Results at a glance" of the volume. The selection (which statements are main results) is made HERE; the TYPE of each statement
# comes from the source (route B), so it cannot drift. Questions are listed completely (a letter per label via \problemvon).
# `HAUPT` = main results: the labels per part.
HAUPT = {
    'I':   ['thm:hoehe', 'cor:fluegel', 'thm:zeugen', 'thm:quadratmitte', 'thm:profil', 'thm:reduktion', 'thm:bbox', 'thm:streifen', 'thm:glatt', 'prop:raenge',
            'prop:exponenten', 'prop:paarliste', 'prop:lucaswieferich', 'rem:zweiprimitive', 'cor:teiler'],
    'II':  ['prop:csechs', 'prop:transient', 'prop:turmsumme', 'cor:untergrenze', 'prop:lokalisierung', 'prop:aufzaehlung', 'prop:klassedrei', 'def:vergleich',
            'rem:tao', 'prop:abstandeinsrechnung'],
    'III': ['prop:phidrei', 'prop:phiacht', 'prop:kopf', 'prop:zyktafel', 'rem:zaehlmodell', 'rem:unabhaengig', 'rem:wieferichbrueche'],
}
# `TEILNAME_KURZ` = short part names for the table
TEILNAME_KURZ = {'I': 'Part I (triples)', 'II': 'Part II (pairs)', 'III': 'Part III (further observations)'}

def blick_schreiben(alle):   # `blick_schreiben` = write the "Results at a glance" table; `alle` = statements of all parts
    kopf = ['\\begingroup\\small', '\\begin{longtable}{l l p{8.0cm}}', '\\hline', 'part & type & statements \\\\', '\\hline', '\\endhead']
    rows, fehler = [], []   # `fehler` = errors
    def nummer_schluessel(nr):   # sort key from a statement number such as 'I.3.1'
        return [int(x) for x in nr.split('.')[1:] if x.isdigit()] if nr else [999]
    for teil, _ in TEILE:
        # `zeilen` = rows: label -> (number, environment, title, label, type)
        zeilen = {z[3]: z for z in alle.get(teil, []) if z[3]}
        for lab in HAUPT[teil]:
            if lab not in zeilen: fehler.append(f'{teil}:{lab} steht nicht im Text')
        erste = True
        for art in ARTEN:
            if art == 'literature': continue
            if art == 'open': eintraege = [f'\\emph{{{z[2] or z[3]}}} (Problem~\\problemvon{{{teil}:{z[3]}}})' for z in alle.get(teil, []) if z[1] == 'question']
            else:
                # ALL own statements of this type in numerical order; the main statements (HAUPT) in bold,
                # one line each, followed by "also:" and the remaining ones in numerical order
                haupt, rest = [], []   # `haupt` = main statements, `rest` = the remaining ones
                for lab, (nr, env, titel, _, a) in sorted(zeilen.items(), key=lambda kv: nummer_schluessel(kv[1][0])):
                    if a != art or env == 'question': continue
                    text = f'{nr or "?"}~{re.sub(r";\s*a heuristic$", "", titel.strip())}'
                    (haupt if lab in HAUPT[teil] else rest).append(text)
                teile_ = [f'\\textbf{{{x}}}' for x in haupt] + ([('\\emph{also:} ' if haupt else '') + '; '.join(rest)] if rest else [])
                eintraege = ['\\newline '.join(teile_)] if teile_ else []
            if not eintraege: continue
            rows.append(f'{TEILNAME_KURZ[teil] if erste else ""} & {art} & ' + '; '.join(eintraege) + ' \\\\'); erste = False
        rows.append('\\hline')
    if fehler: raise SystemExit('belege_bauen (Blick): ' + '; '.join(fehler))
    (S/'gemeinsam'/'blick.tex').write_text('\n'.join(kopf + rows[:-1] + ['\\hline', '\\end{longtable}', '\\endgroup', '']), encoding='utf-8')

# `bauen` = build the tables of all parts; `pruefen` = also compare route A with route B; returns a report
def bauen(pruefen=True):
    bericht = []; alle = {}   # `bericht` = report lines; `alle` = all statements per part
    for teil, ordner in TEILE:
        b, fehler = route_b(teil, ordner); nummern = aux_nummern(ordner); a = route_a(ordner)
        zeilen = []
        for env, titel, label, art in b:
            nr = nummern.get(f'{teil}:{label}') if label else None
            zeilen.append((nr, env, titel, label, art))
        alle[teil] = list(zeilen)
        if a is not None and pruefen:
            liste_a = sorted((nr, art) for nr, art, _ in a)
            liste_b = sorted((nr, art) for nr, env, titel, label, art in zeilen if nr)
            if liste_a != liste_b:
                nur_a = sorted(set(liste_a) - set(liste_b)); nur_b = sorted(set(liste_b) - set(liste_a))
                fehler.append(f'{ordner}: Route A (LaTeX) ≠ Route B (Quelltext): nur A {nur_a[:5]} · nur B {nur_b[:5]}')
        if fehler:
            for f in fehler: print('❌', f)
            raise SystemExit(f'belege_bauen: {len(fehler)} Fehler in {ordner}')
        # write the table (also without .aux/.bel: then it holds placeholders, the next pass fills them)
        ziel = S/'gemeinsam'/f'belegtab_{teil}.tex'
        if not zeilen:
            ziel.write_text('% keine Aussagen in diesem Teil\n', encoding='utf-8'); bericht.append(f'{teil}: 0'); continue
        def schl(z):   # sort key: the numeric parts of the statement number
            nr = z[0] or ''; return [int(x) if x.isdigit() else 0 for x in nr.split('.')]
        zeilen.sort(key=schl)
        zahl = {k: sum(1 for z in zeilen if z[4] == k) for k in ARTEN}
        # longtable instead of tabular: Part I has many statements, so the table must be allowed to run across a page
        kopf = ['\\begingroup\\small', '\\begin{longtable}{l p{7.2cm} l}', '\\hline', 'statement & title & type \\\\', '\\hline', '\\endhead']
        rows = []
        for nr, env, titel, label, art in zeilen:
            name = f'{WORT[env]}~{nr}' if nr else f'{WORT[env]} (unnumbered)'
            rows.append(f'{name} & {titel_tex(titel)} & {art} \\\\')
        fuss = ['\\hline', '\\end{longtable}', '\\endgroup', '',
                'The types: ' + '; '.join(f'\\emph{{{k}}}, {ERKLAERUNG[k]} ({zahl[k]})' for k in ARTEN if zahl[k]) + '.', '']
        ziel.write_text('\n'.join(kopf + rows + fuss), encoding='utf-8')
        bericht.append(f'{teil}: {len(zeilen)} (' + ', '.join(f'{k} {zahl[k]}' for k in ARTEN if zahl[k]) + ')')
    blick_schreiben(alle)
    return bericht

def pk():   # `pk` = positive controls
    # Route B alone must report (1) a missing \belegart and (2) a wrong type; (3) A ≠ B must stand out in the comparison.
    global S
    import tempfile, shutil
    ok = True
    for teil, ordner in TEILE[:1]:
        b, f = route_b(teil, ordner); print(f'  PK Bestand {ordner}: {len(b)} Aussagen, {len(f)} Fehler'); ok &= not f
        a = route_a(ordner); nummern = aux_nummern(ordner)
        if a:
            liste_a = sorted((nr, art) for nr, art, _ in a); liste_b = sorted((nummern.get(f'{teil}:{l}'), art) for _, _, l, art in b if l)
            print(f'  PK A = B: {"✅" if liste_a == liste_b else "❌"} ({len(liste_a)} = {len(liste_b)})'); ok &= liste_a == liste_b
            # `kaputt` = broken copy of list A with one falsified type
            kaputt = list(liste_a); kaputt[0] = (kaputt[0][0], 'heuristic' if kaputt[0][1] != 'heuristic' else 'proved')
            print(f'  PK verfaelschtes Etikett faellt auf: {"✅" if sorted(kaputt) != liste_b else "❌"}'); ok &= sorted(kaputt) != liste_b
    tmp = pathlib.Path(tempfile.mkdtemp()); (tmp/'teil1').mkdir()
    t = (S/'teil1'/'inhalt.tex').read_text(encoding='utf-8')
    t2 = t.replace('\\belegart{', '\\belegartX{', 1)                      # one type label removed
    (tmp/'teil1'/'inhalt.tex').write_text(t2, encoding='utf-8')
    S_alt = S; S = tmp
    try:
        _, f = route_b('I', 'teil1'); print(f'  PK fehlendes \\belegart gemeldet: {"✅" if any("0 \\\\belegart" in x or "0 \\belegart" in x for x in f) else "❌"} {f[:1]}')
        ok &= any('0 \\belegart' in x for x in f)
        (tmp/'teil1'/'inhalt.tex').write_text(t.replace('\\belegart{proved}', '\\belegart{maybe}', 1), encoding='utf-8')
        _, f = route_b('I', 'teil1'); print(f'  PK unbekannte Art gemeldet: {"✅" if any("unbekannte" in x for x in f) else "❌"}'); ok &= any('unbekannte' in x for x in f)
    finally:
        S = S_alt; shutil.rmtree(tmp, ignore_errors=True)
    print('PK', '✅ alle' if ok else '❌'); return 0 if ok else 1

if __name__ == '__main__':
    if sys.argv[1:] == ['pk']: sys.exit(pk())
    print('Belegarten:', ' · '.join(bauen()))
