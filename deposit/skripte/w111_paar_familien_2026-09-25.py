# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w111_paar_familien_2026-09-25.py
# -*- coding: utf-8 -*-
# Purpose: from which Pell families do the 39 pairs (n, n+1) of consecutive powerful numbers come? Each powerful number is
#   n = a^2 b^3 = (ab)^2 * b with b squarefree, so b is the SQUAREFREE PART s(n). A pair n, n+1 therefore solves
#   d*X^2 - b*Y^2 = 1 with b = s(n), d = s(n+1), Y = sqrt(n/b), X = sqrt((n+1)/d) (side conditions b | Y, d | X).
#   The "family" of a pair is (b, d). Square families: b = 1 or d = 1 (one of the two is a square, Erdos #365 question 1).
# Reads: the module w114_a060355_quelle_2026-09-26.py (terms of OEIS A060355 read from its data file; third-party data,
#   terms 27-39 are not stored in the archive); factorization with sympy.factorint (library as a black box).
# Writes: ergebnisse/w111_paar_familien_result.json and ergebnisse/w111_paar_familien_output.txt. No arguments.
# Control (PK, aborts on failure): for all 39 pairs, n and n+1 are powerful (fully factorized) and
#   d*X^2 - b*Y^2 = 1 with b | Y, d | X.
#
# Question. Are the 39 pairs up to 3.9*10^21 a few repeating families (then the number per family grows like log x / log(factor)),
#   or do NEW families keep appearing? Both are compatible with "linear in log x"; the split shows WHERE the slope comes from.
# Expectation, stated in advance: (8, 9), (288, 289), (9800, 9801), ... from x^2 - 8y^2 = 1 form a family (b, d) = (2, 1) with the
#   factor (3 + sqrt 8)^2 ~ 34 per term; w103 counted 9 of 24 from x^2 - 8y^2 = 1. [assumption] About half of the pairs in a few
#   square families, the rest in single families (Walker type).


import sys, json, math, time, pathlib, re
from collections import defaultdict
import sympy
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'
aus = []
def sag(s=''):
    print(s, flush=True); aus.append(s)
sag('='*100); sag('w111 — PELL-FAMILIEN DER 39 PAARE   ' + time.strftime('%Y-%m-%d %H:%M')); sag('='*100)
import importlib.util  # the terms come from w114 (read from the b-file), not hard-coded in this script
_s = importlib.util.spec_from_file_location('w114', pathlib.Path(__file__).resolve().parent/'w114_a060355_quelle_2026-09-26.py')
w114 = importlib.util.module_from_spec(_s); _s.loader.exec_module(w114)
A = w114.glieder()
assert len(A) == 39 and A[0] == 8
# `s_und_pruef` = squarefree part s(n) and the factorization of n; asserts that n is powerful
def s_und_pruef(n):
    f = sympy.factorint(n)
    assert all(e >= 2 for e in f.values()), ('nicht powerful', n, f)
    return math.prod(p for p, e in f.items() if e % 2 == 1), f
# `fam` = families: (b, d) -> list of the n; `zeilen` = rows of the result file
fam = defaultdict(list); zeilen = []
for i, n in enumerate(A, 1):
    b, fn = s_und_pruef(n); d, fn1 = s_und_pruef(n + 1)
    Y2, X2 = n // b, (n + 1) // d; Y, X = math.isqrt(Y2), math.isqrt(X2)
    assert Y*Y == Y2 and X*X == X2 and d*X*X - b*Y*Y == 1 and Y % b == 0 and X % d == 0
    fam[(b, d)].append(n)
    # all 39 pairs are our own computation (w166), so every row carries the value; `lg_n1` stays in every row for
    # the build tool zahlen_bauen (key paar_grenze) and for w101.
    zeilen.append(dict(i=i, n=str(n), lg_n1=round(math.log10(n + 1), 6), b=b, d=d, X=str(X), Y=str(Y)))
sag(f'PK ✅  alle 39 Paare: n und n+1 powerful (vollstaendig faktorisiert), d·X² − b·Y² = 1 mit b | Y, d | X.')
sag()
sag(f'{len(fam)} Familien (b, d) = (quadratfreier Anteil von n, von n+1):')
ordn = sorted(fam.items(), key=lambda kv: kv[1][0])
fam_res = []
for (b, d), ns in ordn:
    art = 'Quadrat-Familie' if 1 in (b, d) else 'KEIN Quadrat (Walker-Typ)'
    fak = [ns[j+1]/ns[j] for j in range(len(ns) - 1)]
    sag(f'  (b, d) = ({b}, {d})  Feld ℚ(√{b*d}) · {len(ns):>2} Paar(e) · erstes n = {w114.etikett(A.index(ns[0]) + 1, ns[0])} · {art}' +
        (f' · Faktor je Glied ≈ {", ".join(f"{x:.4g}" for x in fak[:4])}{" …" if len(fak) > 4 else ""}' if fak else ''))
    fam_res.append(dict(b=b, d=d, anzahl=len(ns), n=[w114.etikett(A.index(x) + 1, x) for x in ns], quadrat=1 in (b, d)))
sag()
# new families per decade
# `dek` = per decade e (number of digits minus 1): [number of pairs, number of pairs from a NEW family]; `gesehen` = seen families
dek = defaultdict(lambda: [0, 0]); gesehen = set()
for n in A:
    key = next(k for k, v in fam.items() if n in v); e = len(str(n)) - 1
    dek[e][0] += 1
    if key not in gesehen: dek[e][1] += 1; gesehen.add(key)
sag('Je Groessenordnung (Dekade von n): Paare / davon aus einer NEUEN Familie')
sag('  ' + ' · '.join(f'10^{e}: {a}/{nf}' for e, (a, nf) in sorted(dek.items())))
# `q` = pairs in square families, `nq` = pairs in non-square families
q = sum(len(v) for k, v in fam.items() if 1 in k); nq = 39 - q
sag(f'Quadrat-Familien liefern {q} von 39 Paaren, Nicht-Quadrat-Familien {nq}. Familien mit ≥ 2 Paaren: {sum(1 for v in fam.values() if len(v) >= 2)}; '
    f'Einzelfamilien: {sum(1 for v in fam.values() if len(v) == 1)}.')
# result record: `paare` = the 39 pairs (rows), `familien` = families, `je_dekade` = per decade, `quadrat_paare` /
# `nicht_quadrat_paare` = number of pairs in square / non-square families
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), quelle=w114.QUELLE,
           paare=zeilen, familien=fam_res, je_dekade={str(e): v for e, v in sorted(dek.items())}, quadrat_paare=q, nicht_quadrat_paare=nq)
(ERG/'w111_paar_familien_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w111_paar_familien_output.txt').write_text('\n'.join(aus) + '\n', encoding='utf-8')
