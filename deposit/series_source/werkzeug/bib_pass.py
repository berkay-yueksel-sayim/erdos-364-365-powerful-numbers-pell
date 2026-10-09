# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# bib_pass.py
# Purpose: first-hand bibliographic data for literatur.bib. For every entry of the list `EINTRAEGE` (entries) below it fetches the
#   metadata from the primary source: DOI via Crossref (otherwise DataCite, for Zenodo and institutional DOIs); arXiv via the
#   abstract page (citation_* metadata, journal reference, primary category). It builds a BibTeX candidate with a source line and
#   compares every number of our own record (volume, issue, pages, year) with the source; one verdict line per entry.
# It does NOT write to literatur.bib: a candidate is moved there only after review (title capitalization, special characters and
#   name forms are judgment calls). Works without a DOI or arXiv number are reported and have to be verified by hand.
# Network: read-only requests for public metadata; no full text, no code; neutral user agent.
# Usage: python bib_pass.py [Key1 Key2 ...]   (without keys: all entries)
#   writes ausgabe/BIB_PASS.txt and gemeinsam/literatur_kandidaten.bib (relative to the folder above this script)
import sys, re, json, time, html, pathlib, unicodedata, urllib.request, urllib.error, urllib.parse
sys.stdout.reconfigure(encoding='utf-8')
S = pathlib.Path(__file__).resolve().parent.parent
BERICHT = S/'ausgabe'/'BIB_PASS.txt'; KAND = S/'gemeinsam'/'literatur_kandidaten.bib'
UA = {'User-Agent': 'bib-pass/1.0 (metadata check)'}

# (key, source, our record as in the bibliography - only the numbers in it are compared);
# source = 'doi:...', 'arxiv:...' or 'suche:...' (`suche` = search: Crossref query string, for works without a DOI)
EINTRAEGE = [
 ('Sentance1981',          'doi:10.1080/00029890.1981.11995246', 'Amer. Math. Monthly 88(4), 272–274, 1981'),
 ('MollinWalsh1987',       'doi:10.1080/00150517.1987.12429724', 'Fibonacci Quart. 25(1), 34–37, 1987'),
 ('Walsh1988',             'doi:10.11575/PRISM/13488',           'MSc thesis, 1988'),
 ('ChimShoreySinha2019',   'doi:10.5486/PMD.2019.8397',          '2019'),
 ('ChimNairShorey2018',    'doi:10.46298/hrj.2019.5117',         '2018'),
 ('McIntoshRoettger2007',  'doi:10.1090/S0025-5718-07-01955-2',  '2007'),
 ('RossShenCai2025',       'doi:10.1080/00150517.2026.2656703',  '2025'),
 ('Fathi2026',             'doi:10.5281/zenodo.22178088',        '2026'),
 ('Ballot2019',            'doi:10.1080/00150517.2019.12427650', '2019'),
 ('Chan2012',              'doi:10.1017/S144678871200016X',      'J. Aust. Math. Soc. 93, 43–51, 2012'),
 ('Rizzi2026',             'doi:10.5281/zenodo.22030563',        '2026'),
 ('Sayim2026EMW',          'doi:10.5281/zenodo.22699364',        'v1.3, 2026'),
 ('Ma2026',                'arxiv:2608.23418',                   '2026'),
 ('MianSiddique2026',      'arxiv:2609.25011',                   '2026'),
 ('Reinhart2024',          'arxiv:2402.09827',                   '2024'),
 ('AnithaEtAl2021',        'arxiv:2101.04901',                   '2021'),
 ('Kym2026',               'arxiv:2605.24909',                   '2026'),
 ('She2025',               'arxiv:2507.16828',                   'Integers 25, A103, 2025'),
 ('Tao2026',               'arxiv:2603.27990',                   '2026'),
 ('Reuss2014',             'arxiv:1212.3150',                    '2014'),
 ('vanDoorn2026',          'arxiv:2605.06697',                   '2026'),
 ('Pandey2023',            'arxiv:2303.06610',                   '2023'),
 # works without a DOI in our record: Crossref search with strict acceptance (all record numbers + first author)
 ('MollinWalsh1986b',      'suche:Mollin Walsh On powerful numbers International Journal of Mathematics and Mathematical Sciences 1986', 'IJMMS 9, 801–806, 1986'),
 ('BennettWalsh1999',      'suche:Bennett Walsh The Diophantine equation b^2X^4-dY^2=1 Proceedings of the American Mathematical Society 1999', 'Proc. AMS 127(12), 3481–3491, 1999'),
 ('Granville2012',         'suche:Granville Primitive prime factors in second-order linear recurrence sequences Acta Arithmetica', 'Acta Arith. 155, 431–452, 2012'),
 ('Silverman1988',         'suche:Silverman Wieferich criterion and the abc-conjecture Journal of Number Theory 1988', 'JNT 30, 1988'),
 ('LaishramShorey2012',    'suche:Laishram Shorey Baker explicit abc-conjecture and applications Acta Arithmetica 2012', 'Acta Arith. 155, 2012'),
 ('StephensWilliams1988',  'suche:Stephens Williams Computation of real quadratic fields with class number one Mathematics of Computation 1988', 'Math. Comp. 50, 619–632, 1988'),
 ('Carmichael1913',        'suche:Carmichael On the numerical factors of the arithmetic forms alpha^n pm beta^n Annals of Mathematics 1913', 'Ann. of Math. 15, 30–70, 1913'),
 ('Yabuta2001',            'suche:Yabuta A simple proof of Carmichael theorem on primitive divisors Fibonacci Quarterly 2001', 'Fibonacci Quart. 39(5), 439–443, 2001'),
 ('Yabuta2007',            'suche:Yabuta The ABC-conjecture and the powerful numbers in Lucas sequences Fibonacci Quarterly 2007', 'Fibonacci Quart. 45(4), 362–365, 2007'),
 ('Snyder2000',            'suche:Snyder An alternate proof of Mason theorem Elemente der Mathematik 2000', 'Elem. Math. 55, 93–94, 2000'),
 ('Golomb1970',            'suche:Golomb Powerful numbers American Mathematical Monthly 1970', 'Amer. Math. Monthly 77(8), 848–852, 1970'),
 ('Lehmer1930',            'suche:Lehmer An extended theory of Lucas functions Annals of Mathematics 1930', 'Ann. of Math. 31, 419–448, 1930'),
 ('Lucas1878',             'suche:Lucas Theorie des fonctions numeriques simplement periodiques American Journal of Mathematics 1878', 'Amer. J. Math. 1, 184–240, 1878'),
 ('AktasMurty2017',        'suche:Aktas Murty Fundamental units and consecutive squarefull numbers International Journal of Number Theory 2017', 'IJNT 13(1) 243–252, 2017'),
 ('Bateman1954',           'suche:Bateman Squarefull integers American Mathematical Monthly 1954', 'Amer. Math. Monthly 61, 477–479, 1954'),
 ('Ribenboim2001',         'suche:Ribenboim On square factors of terms of binary recurring sequences and the ABC conjecture Publicationes Mathematicae Debrecen', 'Publ. Math. Debrecen 59, 459–469, 2001'),
 ('Rout2019',              'suche:Rout Lucas non-Wieferich primes in arithmetic progressions Functiones et Approximatio', 'Funct. Approx. Comment. Math. 60(2), 167–175, 2019'),
 ('ChenDing2017',          'suche:Chen Ding non-Wieferich primes in arithmetic progressions Proceedings of the American Mathematical Society 2017', 'Proc. AMS 145, 1833–1836, 2017'),
 ('Ding2019',              'suche:Ding non-Wieferich primes Comptes Rendus Mathematique 2019', 'C. R. Math. 357, 483–486, 2019'),
 ('BajpaiBennettChan2024', 'suche:Bajpai Bennett Chan Arithmetic progressions in squarefull numbers International Journal of Number Theory', 'IJNT 20, 19–45, 2024'),
 ('Walker1967',            'suche:Walker On the Diophantine equation mX^2-nY^2=pm1 American Mathematical Monthly 1967', 'Amer. Math. Monthly 74(5), 504–513, 1967'),
 ('Ballot2008',            'suche:Ballot On the 1/3 density of odd ranked primes in Lucas sequences Uniform Distribution Theory', 'Unif. Distrib. Theory 3(2), 129–145, 2008'),
 ('SorensonWebster2017',   'doi:10.1090/mcom/3134',              'Math. Comp. 86(304), 985–1003, 2017'),
 ('Stewart1977',           'suche:Stewart On divisors of Fermat Fibonacci Lucas and Lehmer numbers Proceedings of the London Mathematical Society 1977', 'Proc. London Math. Soc. (3) 35, 425–447, 1977'),
 ('Stewart2013',           'doi:10.1007/s11511-013-0105-y',      'Acta Math. 211, 291–314, 2013'),
 ('MurtySeguin2019',       'doi:10.1016/j.jnt.2019.02.016',      'J. Number Theory 201, 1–22, 2019'),
 ('HarringtonJones2023',   'doi:10.1017/S0004972723000138',      'Bull. Aust. Math. Soc. 108(3), 373–378, 2023'),
 ('CrandallDilcherPomerance1997', 'suche:Crandall Dilcher Pomerance A search for Wieferich and Wilson primes Mathematics of Computation 1997', 'Math. Comp. 66, 433–449, 1997'),
 ('Cohn1997',              'doi:10.4064/aa-78-4-401-403',        'Acta Arith. 78, 401–403, 1997'),
 ('Poonen2003',            'suche:Poonen Squarefree values of multivariable polynomials Duke Mathematical Journal 2003', 'Duke Math. J. 118, 353–373, 2003'),
 ('Granville1998',         'doi:10.1155/S1073792898000592',      'IMRN 1998, 19, 991–1009'),
 ('Hooley1967',            'suche:Hooley On the power free values of polynomials Mathematika 1967', 'Mathematika 14, 21–26, 1967'),
 ('Wieferich1909',         'suche:Wieferich Zum letzten Fermatschen Theorem Journal fur die reine und angewandte Mathematik 1909', 'J. reine angew. Math. 136, 293–302, 1909'),
 ('FelliniMurty2026',      'arxiv:2508.08472',                   '2025'),
 ('Chan2025',              'arxiv:2503.21485',                   '2025'),
 ('Chan2024',              'doi:10.5281/zenodo.11352677',        '2024'),
 ('BennettWalsh2024',      'doi:10.5281/zenodo.11352598',        '2024'),
 ('BiluHanrotVoutier2001', 'doi:10.1515/crll.2001.080',          'J. reine angew. Math. 539, 75–122, 2001'),
 ('Park2012',              'arxiv:1208.5353',                    '2012'),
 ('Luca1999',              'suche:Luca Arithmetic functions of Fibonacci numbers Fibonacci Quarterly 1999', 'Fibonacci Quart. 37(3), 265–268, 1999'),
 ('BugeaudLucaMignotteSiksek2005', 'suche:Bugeaud Luca Mignotte Siksek On Fibonacci numbers with few prime divisors Proceedings of the Japan Academy 2005', 'Proc. Japan Acad. Ser. A 81, 17–20, 2005'),
 ('Pongsriiam2019',        'suche:Pongsriiam Fibonacci Quarterly 2019',  'Fibonacci Quart. 57(5), 130–144, 2019'),
 ('SomerKrizek2015',       'suche:Somer Krizek On primes in Lucas sequences Fibonacci Quarterly 2015', 'Fibonacci Quart. 53(1), 2015'),
 ('Sanna2016',             'suche:Sanna The p-adic valuation of Lucas sequences Fibonacci Quarterly 2016', 'Fibonacci Quart. 54, 118–124, 2016'),
 ('Ballot2019Errata',      'doi:10.1080/00150517.2019.12427639', 'Fibonacci Quart. 57(4), 366, 2019'),
 ('Bouchard2026',          'doi:10.5281/zenodo.21854972',        'Zenodo, 2026 (Versions-DOI)'),
 ('KotyadaMuthukrishnan2018', 'suche:Kotyada Muthukrishnan Non-Wieferich primes in number fields and abc-conjecture Czechoslovak Mathematical Journal 2018', 'Czechoslovak Math. J. 68(2), 445–453, 2018'),
 ('RibenboimWalsh1999',    'suche:Ribenboim Walsh The ABC conjecture and the powerful part of terms in binary recurring sequences Journal of Number Theory 1999', 'J. Number Theory 74, 134–147, 1999'),
 ('StewartYu2001',         'suche:Stewart Yu On the abc conjecture II Duke Mathematical Journal 2001', 'Duke Math. J. 108, 169–181, 2001'),
 ('Yokoi1970',             'suche:Yokoi fundamental unit of real quadratic fields with norm 1 Journal of Number Theory 1970', 'J. Number Theory 2, 106–115, 1970'),
]
# Key filter: `python bib_pass.py Key1 Key2 ...` fetches only these entries; report and candidate file then contain
# only them (the candidate file is just a scratch area anyway; accepted entries live in literatur.bib).
if len(sys.argv) > 1:
    EINTRAEGE = [e for e in EINTRAEGE if e[0] in sys.argv[1:]]
    assert EINTRAEGE, 'kein bekannter Schluessel'
# `EIGENNAMEN` = proper names that must keep their capitalization in titles
EIGENNAMEN = {'Pell', 'Pellian', 'Diophantine', 'Lucas', 'Wieferich', 'Fibonacci', 'Erdős', 'Erdos', 'Mollin', 'Walsh', 'Möbius', 'Baker',
              'Lehmer', 'Mordell', 'Mason', 'Stothers', 'Golomb', 'Cramér', 'Riemann', 'Lean', 'Wall', 'Sun', 'Ankeny', 'Artin', 'Chowla',
              'Wolstenholme', 'Carmichael', 'Sylvester', 'Zsigmondy', 'Bang', 'Beckon', 'Sentance', 'Walker', 'Tao', 'Chan'}
# AMS abbreviations (MathSciNet) as in literatur.bib - amsplain expects them; for a new journal, add it here
# (`ABK` = abbreviations)
ABK = {'The American Mathematical Monthly': 'Amer. Math. Monthly', 'The Fibonacci Quarterly': 'Fibonacci Quart.', 'Publicationes Mathematicae Debrecen': 'Publ. Math. Debrecen', 'Mathematics of Computation': 'Math. Comp.', 'Journal of the Australian Mathematical Society': 'J. Aust. Math. Soc.', 'International Journal of Mathematics and Mathematical Sciences': 'Internat. J. Math. Math. Sci.', 'Proceedings of the American Mathematical Society': 'Proc. Amer. Math. Soc.', 'Acta Arithmetica': 'Acta Arith.', 'Journal of Number Theory': 'J. Number Theory', 'Elemente der Mathematik': 'Elem. Math.', 'Comptes Rendus. Mathématique': 'C. R. Math. Acad. Sci. Paris', 'International Journal of Number Theory': 'Int. J. Number Theory', 'The Annals of Mathematics': 'Ann. of Math. (2)', 'American Journal of Mathematics': 'Amer. J. Math.', 'Functiones et Approximatio Commentarii Mathematici': 'Funct. Approx. Comment. Math.', 'Hardy-Ramanujan Journal': 'Hardy-Ramanujan J.', 'Michigan Mathematical Journal': 'Michigan Math. J.'}

# `hole` = fetch: HTTP GET with up to three attempts; returns (status code, body), or (0, b'') after repeated failure
def hole(url, accept=None):
    h = dict(UA)
    if accept: h['Accept'] = accept
    for versuch in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=30) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            return e.code, b''
        except Exception:
            time.sleep(2 + 3 * versuch)
    return 0, b''

def schuetze(titel):
    """Brace proper names and words with several capital letters (acronyms) so that no bibliography style lowercases them."""
    def w(m):
        x = m.group(0)
        stamm = re.split(r"['’]", x)[0]  # possessive: Baker's -> {Baker}'s (otherwise amsplain would lowercase it)
        if stamm in EIGENNAMEN and stamm != x: return '{' + stamm + '}' + x[len(stamm):]
        return '{' + x + '}' if (x in EIGENNAMEN or sum(c.isupper() for c in x) >= 2) else x
    # TeX in the title (arXiv delivers e.g. Erd\H{o}s) stays untouched: no word directly after \ { } and none that touches
    # { or \ (otherwise Erd\H{o}s-Mollin-Walsh would become Erd\H{o}{s-Mollin-Walsh}). Multi-line Crossref titles (MathML)
    # are flattened to single spaces first.
    titel = ' '.join(titel.split())
    # Math with capital letters goes into {...}: BibTeX also lowercases inside $...$
    # (otherwise $b^2X^4-dY^2=1$ would become $b^2x^4-dy^2=1$)
    titel = re.sub(r'(?<!\{)(\$[^$]*[A-Z][^$]*\$)', r'{\1}', titel)
    return re.sub(r"(?<![\\{}A-Za-z])[A-Za-zÀ-ÿĀ-žő'’-]+(?![{\\A-Za-z])", w, titel)

# `tex` = TeX-escape & and %, convert en and em dashes to -- and ---
def tex(s):
    return (s.replace('&', r'\&').replace('%', r'\%').replace('–', '--').replace('—', '---')) if s else s

# `aus_crossref` / `aus_datacite` / `aus_arxiv` = from Crossref / DataCite / arXiv: turn the retrieved record into
# a BibTeX type and a field dictionary `f` (`aus_arxiv` also returns the journal reference text)
def aus_crossref(m):
    typ = {'journal-article': 'article', 'proceedings-article': 'inproceedings', 'book-chapter': 'incollection',
           'book': 'book', 'dissertation': 'phdthesis', 'posted-content': 'misc'}.get(m.get('type'), 'misc')
    f = {'author': ' and '.join(f"{a.get('family', '')}, {a.get('given', '')}".strip(', ') for a in m.get('author', [])),
         'title': schuetze(html.unescape(re.sub(r'<[^>]+>', '', (m.get('title') or [''])[0]))),
         'journal': ABK.get((m.get('container-title') or [''])[0], (m.get('container-title') or [''])[0]), 'volume': m.get('volume', ''), 'number': m.get('issue', ''),
         'pages': (m.get('page') or '').replace('-', '--'), 'year': str(((m.get('issued') or {}).get('date-parts') or [['']])[0][0]),
         'doi': m.get('DOI', '')}
    if typ != 'article': f['booktitle'] = f.pop('journal')
    return typ, f

def aus_datacite(a):
    typ = {'Dissertation': 'phdthesis', 'Preprint': 'misc', 'Text': 'misc', 'Software': 'misc', 'Dataset': 'misc'}.get(
        (a.get('types') or {}).get('resourceTypeGeneral', ''), 'misc')
    f = {'author': ' and '.join(c.get('name', '') for c in a.get('creators', [])),
         'title': schuetze((a.get('titles') or [{}])[0].get('title', '')), 'howpublished': a.get('publisher', ''),
         'year': str(a.get('publicationYear', '')), 'doi': a.get('doi', ''), 'note': 'Version ' + str(a.get('version')) if a.get('version') else ''}
    return typ, f

def aus_arxiv(seite, aid):
    t = seite.decode('utf-8', 'replace')
    meta = lambda k: [html.unescape(x) for x in re.findall(rf'<meta name="{k}" content="([^"]*)"', t)]
    jref = re.findall(r'class="tablecell jref">(.*?)</td>', t, re.S)
    prim = re.findall(r'class="primary-subject">[^(]*\(([^)]+)\)', t)
    vers = re.findall(r'\[v(\d+)\]', t)
    f = {'author': ' and '.join(meta('citation_author')), 'title': schuetze(' '.join((meta('citation_title') or [''])[0].split())),
         'year': (meta('citation_date') or [''])[0][:4], 'eprint': aid, 'archiveprefix': 'arXiv', 'primaryclass': prim[0] if prim else '',
         'note': f'arXiv:{aid}' + (f' (v{max(map(int, vers))})' if vers else '')}
    return 'misc', f, (' '.join(re.sub(r'<[^>]+>', '', jref[0]).split()) if jref else '')

# `zahlen` = numbers: the set of digit strings in s;  `ohne_akzent` = without accents
def zahlen(s):
    return set(re.findall(r'\d+', s or ''))

def ohne_akzent(s):
    return ''.join(c for c in unicodedata.normalize('NFKD', s) if not unicodedata.combining(c))

# `zeilen` = report lines, `bib` = BibTeX candidate entries; `sag` = say: print a line and keep it for the report
zeilen, bib = [], []
def sag(x=''): print(x, flush=True); zeilen.append(x)
sag(f'BIB-PASS · {time.strftime("%Y-%m-%d %H:%M")} · {len(EINTRAEGE)} Eintraege · Quelle je Eintrag in der Zeile darueber')
# `quelle` = source, `unser` = our record
for key, quelle, unser in EINTRAEGE:
    art, kennung = quelle.split(':', 1)
    jref = ''
    # `herkunft` = origin: URL of the retrieved data; `jref` = journal reference found on the arXiv page
    if art == 'doi':
        st, b = hole('https://api.crossref.org/works/' + kennung)
        if st == 200:
            typ, f = aus_crossref(json.loads(b)['message']); herkunft = 'https://api.crossref.org/works/' + kennung
        else:
            st, b = hole('https://api.datacite.org/dois/' + kennung, 'application/vnd.api+json')
            if st != 200:
                sag(f'🔴 {key}: {kennung} weder bei Crossref noch bei DataCite (HTTP {st})'); continue
            typ, f = aus_datacite(json.loads(b)['data']['attributes']); herkunft = 'https://api.datacite.org/dois/' + kennung
    elif art == 'suche':
        # Work without a DOI in our record: Crossref search. A hit is accepted only if ALL numbers of our record (volume, issue,
        # pages, year) occur in its metadata AND its first author occurs in the record or the key. Otherwise the three best
        # hits are listed for a manual decision. (`treffer` = hits, `gefunden` = found, `erst` = first author)
        st, b = hole('https://api.crossref.org/works?' + urllib.parse.urlencode({'query.bibliographic': kennung, 'rows': 3}))
        treffer = json.loads(b)['message']['items'] if st == 200 else []
        gefunden = None
        for m in treffer:
            typ, f = aus_crossref(m)
            qz = set().union(*(zahlen(f.get(k, '')) for k in ('volume', 'number', 'pages', 'year')))
            erst = ohne_akzent(((m.get('author') or [{}])[0].get('family') or '').lower())
            if erst and erst in ohne_akzent((unser + ' ' + key).lower()) and not (zahlen(unser) - qz):
                gefunden = (typ, f); break
        if not gefunden:
            sag(f'⚠️ {key}: Suche ohne bestaetigten Treffer — von Hand entscheiden:')
            for m in treffer:
                t_, f_ = aus_crossref(m)
                sag(f'     ? {f_["doi"]} · {f_.get("author", "")[:40]} · {f_.get("year")} · {f_.get("journal") or f_.get("booktitle", "")} '
                    f'{f_.get("volume", "")}({f_.get("number", "")}) {f_.get("pages", "")} · {f_.get("title", "")[:70]}')
            time.sleep(1); continue
        typ, f = gefunden; herkunft = f'https://api.crossref.org/works/{f["doi"]} (gefunden ueber die Suche „{kennung}")'
    else:
        st, b = hole('https://arxiv.org/abs/' + kennung)
        if st != 200:
            sag(f'🔴 {key}: arXiv {kennung} HTTP {st}'); continue
        typ, f, jref = aus_arxiv(b, kennung); herkunft = 'https://arxiv.org/abs/' + kennung
    # `quellzahlen` = numbers found in the source; `fehlt` = numbers of our record missing there; `urteil` = verdict
    quellzahlen = set()
    for k in ('volume', 'number', 'pages', 'year'):
        quellzahlen |= zahlen(f.get(k, ''))
    fehlt = sorted(zahlen(unser) - quellzahlen - zahlen(jref), key=int)
    urteil = '✅' if not fehlt else f'⚠️ unsere Zahlen {fehlt} nicht in der Quelle'
    sag(f'{urteil} {key} · {typ} · {f.get("author", "")[:60]} · {f.get("year")} · {f.get("journal") or f.get("howpublished") or f.get("note", "")}'
        f' {f.get("volume", "")}({f.get("number", "")}) {f.get("pages", "")}' + (f' · Journal-Referenz: {jref}' if jref else ''))
    sag(f'     Titel: {f.get("title", "")[:150]}')
    felder = ',\n'.join(f'  {k:<13}= {{{tex(v)}}}' for k, v in f.items() if v)
    bib.append(f'% quelle: {herkunft} · abgerufen {time.strftime("%Y-%m-%d")} · unser Record: {unser}'
               + (f'\n% arXiv-Journal-Referenz: {jref}' if jref else '') + f'\n@{typ}{{{key},\n{felder}\n}}\n')
    time.sleep(1)
BERICHT.write_text('\n'.join(zeilen) + '\n', encoding='utf-8')
KAND.write_text('% literatur_kandidaten.bib — ERZEUGT von werkzeug/bib_pass.py. NICHT eingebunden. Ein Eintrag wandert erst nach Durchsicht\n'
                '% nach literatur.bib (Titel-Schreibweise, Namen, Sonderzeichen sind Urteil).\n\n' + '\n'.join(bib), encoding='utf-8')
sag(f'→ {BERICHT.name} · {KAND.name} ({len(bib)} Kandidaten)')
