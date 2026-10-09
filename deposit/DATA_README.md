# DATA_README: data and code of the series on consecutive powerful numbers

**Author:** Berkay Yüksel Sayim (ORCID 0009-0004-4993-7352, contact: berksa@tutamail.com).
**Licenses:** data, sources and documentation CC BY 4.0 (see `LICENSE-paper`), code Apache 2.0 (see `LICENSE` and `NOTICE`).
**Integrity:** `MANIFEST_sha256.txt` lists the SHA-256 hash and the size of every file in this folder.
**Release:** data release 1 (October 2026), belonging to the v1.0 papers. The archive of further scripts and results of the investigation, and English comments, follow with v1.1.

## What is here

- `ergebnisse/` (results): the result files (`*_result.json`) and raw text outputs (`*_output.txt`) of the runs. `CORE_FILES.txt` lists the 89 files that the papers cite. The folder also holds the primality certificates (`w164_zertifikate/`, `w165_zertifikate/`, `w95_primzertifikate.json`, `w136_zertifikate.json`) with the runs that check them, the decimal expansions of the open building blocks (`w147_offene_bloecke_dezimal.txt`, `w155_algebraische_teile_dezimal.txt`), and the inputs and outputs of the included scripts.
- `skripte/` (scripts): the Python scripts that wrote these files. The name of a script starts with its run label (for example `w152`) and usually ends with the date on which it was written; its result files carry the same label.
- `series_source/`: the LaTeX sources of the volume (`band/`) and of Parts I to III (`teil1/` to `teil3/`), the shared files (`gemeinsam/`, among them `zahlen.tex` with every computed number, `literatur.bib`, `probleme.tex`), the figures (`abbildungen/`) and the build tool (`werkzeug/`).
- `werkzeug/` (tools): container recipes for GMP-ECM and PARI/GP. The programs themselves are not part of the deposit.
- `ecpp/`: our own prover for elliptic curve primality certificates (method of Atkin and Morain, written from the published description). `e1` to `e3` are the building blocks and a single proof step, `e4` the search strategy, `e5_teil1.py` and `e5_teil2.py` the runs on the numbers of Parts I and II. Its certificates are in `ecpp/ergebnisse/`; each one is accepted by our checker `skripte/w164_ecpp_pruefer_2026-10-02.py` and by PARI/GP's `primecertisvalid` (`ecpp/pari_pruefen.gp`, output `ecpp/ergebnisse/e5_pari_pruefung_output.txt`). It produced certificates for the 6 numbers of Part I and 17 of the 19 numbers of Part II that PARI/GP certified; the two largest (1191 and 1687 digits) are beyond its present reach. The papers rest on the PARI/GP certificates; these are a second, independent route.

## How a number gets into the papers

The LaTeX sources refer to each computed number by a key, `\Z{key}`. The build tool `series_source/werkzeug/zahlen_bauen.py` knows for each key the result file(s) it comes from and how the printed value is computed. It checks every file it reads against `MANIFEST_sha256.txt` (hash and size), stops if a file is missing or differs, and writes `series_source/gemeinsam/zahlen.tex`. Run it from any directory:

```
python series_source/werkzeug/zahlen_bauen.py
```

## Controls

The scripts run positive controls (a known case must come out right, for example the factorization T₇(7) = 2³ · 29 · 197 · 2857 of Mollin and Walsh) and, where possible, negative controls (a deliberately wrong input must be rejected), and abort if a control fails. In the outputs these appear as `PK` or `PK+` (positive control) and `PK-`, `PK−` or `NK` (negative control).

## Conventions

- Kernel m: the squarefree part of n² − 1. T₁(m) + U₁(m)√m is the smallest solution of x² − my² = 1 with y > 0, and T_k + U_k√m its k-th power.
- m′: the product of the primes p | m with p ∤ U₁(m). A lattice point is a pair (m, k) with k an odd multiple of m′, so that m | U_k.
- Status of a lattice point: `parity` (T_k odd, no candidate), `witness` (a prime p with v_p(T_k) = 1 is given, so T_k is not powerful), `open` (no witness found in the searched range).
- Residue classes: kernels m ≡ 7 (mod 8) carry the candidates for triples (middles n ≡ 0 mod 4); kernels m ≡ 3 (mod 8) carry pairs at distance two (middles n ≡ 2 mod 4).

## Glossary of German names

| German | English | German | English |
|---|---|---|---|
| ergebnisse | results | skripte | scripts |
| werkzeug | tool(s) | gemeinsam | shared |
| band | volume | teil | part |
| inhalt | content | rahmen | framework (introduction of the volume) |
| zahlen | numbers | belege, belegtab | evidence, table of evidence types |
| probleme | problems | literatur | references |
| abbildungen | figures | bauen | build |
| Kern | kernel | Mitte | middle (of a triple) |
| Zeuge | witness | Leiter, Turm | ladder, tower |
| Paar, Paare | pair, pairs | Familie | family |
| Abstand, Luecke | distance, gap | Zertifikat | certificate |
| Primteil, primitiver Teil | primitive part | Rang | rank (of apparition) |
| Pruefung, pruefen | check | Stelle(n) | digit(s), or position |
| gerade, ungerade | even, odd | offen | open |

## Environment

Python 3.12 with gmpy2 2.3, SymPy 1.14 and NumPy 2.4. The factor searches use GMP-ECM 7.0.6 and the primality certificates PARI/GP 2.17.2, both run as installed programs inside the container images described in `werkzeug/`; a separate script of ours checks every certificate. All code in this deposit was written for this project.

## Edits for this release

The code comments and docstrings of the scripts were translated into English, and internal working notes were removed; German names of variables and functions were kept and are explained in the comments. The LaTeX sources carry no comments. The code itself is unchanged (checked by comparing the syntax trees), except for a few lines that set paths or names: the data path of `zahlen_bauen.py` and the paths in `ecpp/` point to this folder, `bib_pass.py` sends a neutral user agent, and `bauen.py` skips a concordance step whose tool is not included. In a few printed messages, internal labels were replaced by neutral wording, in the script and in its stored output: `w20`, `w25`, `w35`, `w48`, `w55`, `w58`, `w61`, `w103`, `w112`, `w134`, `w155`, `w165`. No computed value changed; the build tool reproduces every number of the papers from this folder. Some outputs record a hash of the script that wrote them; that hash refers to the script as it was run, before its comments were translated, and therefore differs from the hash of the file published here.

## Not yet included

The archive of the investigation (further scripts, results and negative results that the papers do not print) and English versions of the comments and outputs will follow with v1.1. Until then the archive is part of the v1.0 deposit description in the papers but not of this repository.
