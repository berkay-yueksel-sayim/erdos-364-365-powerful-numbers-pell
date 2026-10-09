# Consecutive powerful numbers, Pell equations, and Erdős problems #364 and #365

A positive integer is powerful if every prime that divides it divides it at least twice. Pairs of consecutive powerful numbers exist in abundance, but no three consecutive powerful numbers are known; Erdős asked whether they exist. This repository holds a collected volume, its three parts, and the data and code behind the computed numbers.

> **Status: preprint (v1.0) with data release 1 (October 2026).** The papers are the v1.0 PDFs of the Zenodo records. The data deposit holds every result file that the papers cite, the inputs these result files depend on, the scripts that wrote them, and the LaTeX sources with the build tool. The remaining archive of the investigation and an English translation of the code comments will follow with v1.1.

## Papers

| | Title | Zenodo (all versions) | File |
|---|---|---|---|
| Volume | Consecutive Powerful Numbers and Pell Equations: Triples, Pairs, and Further Observations | [10.5281/zenodo.23127546](https://doi.org/10.5281/zenodo.23127546) | `papers/volume_consecutive-powerful-numbers-and-pell-equations.pdf` |
| Part I | Triples of Consecutive Powerful Numbers: Pell Addresses, a Height Inequality, and Witnesses | [10.5281/zenodo.23127702](https://doi.org/10.5281/zenodo.23127702) | `papers/part-I_triples-of-consecutive-powerful-numbers.pdf` |
| Part II | Pairs of Consecutive and Near-Consecutive Powerful Numbers: Pell Families, Towers, and Localization | [10.5281/zenodo.23127704](https://doi.org/10.5281/zenodo.23127704) | `papers/part-II_pairs-of-consecutive-and-near-consecutive-powerful-numbers.pdf` |
| Part III | Further Observations on Powerful Numbers: Gaps, Cyclotomic Values, and Independence Tests | [10.5281/zenodo.23127706](https://doi.org/10.5281/zenodo.23127706) | `papers/part-III_further-observations-on-powerful-numbers.pdf` |

The PDFs are identical to those of the Zenodo records; `SHA256SUMS.txt` lists their SHA-256 hashes.

Related earlier preprint: *Nonexistence of consecutive powerful triplets around cubes with mixed prime factorizations*, [10.5281/zenodo.20654529](https://doi.org/10.5281/zenodo.20654529).

## Data and code

| Folder | Content |
|---|---|
| `deposit/ergebnisse/` | result files (JSON) and the raw text output of the runs, including the primality certificates and the decimal expansions of the open building blocks; the 89 files that the papers cite are listed in `deposit/CORE_FILES.txt` |
| `deposit/skripte/` | the Python scripts that wrote these files; they state what they expect and run positive controls, and where possible negative controls, before they compute |
| `deposit/series_source/` | LaTeX sources of the four PDFs and the build tool (`werkzeug/`) |
| `deposit/werkzeug/` | recipes of the container images for GMP-ECM and PARI/GP (the programs themselves are not included) |
| `deposit/ecpp/` | our own prover for elliptic curve primality certificates and its certificates, a second route next to the PARI/GP certificates |
| `deposit/MANIFEST_sha256.txt` | SHA-256 hash and size of every file in `deposit/` |
| `deposit/DATA_README.md` | conventions, glossary, reproduction |
| `deposit/OPEN_PROBLEMS.md` | the open problems of the volume |

**Reproducing the printed numbers.** `python deposit/series_source/werkzeug/zahlen_bauen.py` checks every result file it reads against the manifest, recomputes each number in its printed form, and writes `deposit/series_source/gemeinsam/zahlen.tex`. It stops if a file is missing or differs from the manifest.

**Language and edits.** The code comments are in English; German names of folders, files, variables and the raw text outputs are kept, and `deposit/DATA_README.md` has a glossary. Internal working notes were removed; the code itself is unchanged. In a few printed messages, internal labels were replaced by neutral wording, in the script and in its stored output; no computed value changed (list in `deposit/DATA_README.md`).

## Corrections since v1.0

These will be incorporated in v1.1. None of them changes a result or a number of the series.

1. **OEIS A135735 and the condition m | U₁(m).** Part I, Remark 2.5, and the passages that repeat it say that the squarefree m with m | U₁(m) are exactly the members of OEIS A135735. This is too general. The members of A135735 are squarefree d that divide the coefficient v of the fundamental unit u + vω of the maximal order of ℚ(√d) (A. Reinhart, Remark 5.5 in arXiv:2305.09267, published in Acta Arith. 211 (2023); the search goes back to Stephens and Williams, Math. Comp. 50 (1988)). Outside the residue classes used in the series the two notions can differ: m = 2 (with T₁ + U₁√2 = 3 + 2√2) and m = 1221562 are squarefree and satisfy m | U₁(m), but are not in A135735. For squarefree m ≡ 3 (mod 4) the two notions agree, because then ℤ[√m] is the maximal order and x² − my² = −1 has no solution, so T₁ + U₁√m is the fundamental unit. This covers every kernel m ≡ 7 (mod 8) of a triple and every m with T₁(m) even, since T₁(m) even forces m ≡ 3 (mod 4).
2. **Citation.** A. Reinhart, *A counterexample to the Pellian equation conjecture of Mordell*, is published in Acta Arith. 215 (2024), 85–95, [doi:10.4064/aa240214-3-4](https://doi.org/10.4064/aa240214-3-4); the papers cite the arXiv version.
3. **Range of the extended search.** A later paper, A. Reinhart, *A counterexample to the conjecture of Ankeny, Artin and Chowla* ([arXiv:2410.21864](https://arxiv.org/abs/2410.21864), Section 3; to appear in Algebra & Number Theory), reports the search for squarefree d with d | y up to 2.3 · 10¹⁴. The papers of the series use the earlier range up to 5.325 · 10¹³, with its caveat that no independent double check of the search was done; the bounds that rest on this range will be reviewed in v1.1.
4. **Reinhart's condition (C).** Problem D of the volume (Part I, Question 9.4; `deposit/OPEN_PROBLEMS.md`) restates condition (C) of Reinhart (Definition 1.4 of the paper in item 2) as "is there an index k with T_k(d) powerful and d | U_k(d)". The definition also requires d ≡ 7 (mod 8), k odd and U_k odd (for d ≡ 7 (mod 8) the last condition is the same as T_k even). The restatement is therefore not equivalent; the open problem is Reinhart's condition as he defines it.
5. **Our own check of all 241 citations against the sources (October 2026)** found these further imprecisions; none of them changes a result:
   - Corollary 5.1 of Ross, Shen and Cai is about the Möbius duals $\prod_{d \mid n} F_d^{\mu(n/d)}$, not about the primitive parts of the Fibonacci numbers; Part I states it correctly in one place and incorrectly in another.
   - Corollary 6.6 of Kym holds for indices above a bound and outside an exceptional set, and for primitive divisors that do not divide A; Part I omits these hypotheses and calls the arXiv preprint published.
   - That the family of Golomb's pair 12167, 12168 is infinite follows from Walker's Theorem 3.5 (p. 115), not from p. 111.
   - Luca's bound $\omega(F_n) \ge \tau(n) - 2$ (p. 267) is stated for $6 \nmid n$.
   - The bound $10^{51075}$ for a triple follows from Theorem 1.8 of Chim, Shorey and Sinha together with the bound they give for a triple just before it; it is not stated there as a theorem.
   - The heuristic of van Doorn is a product of factors over squarefree m, not over primes.
   - Locators: Reinhart's Definition 1.4 is on p. 2 of the arXiv version (his remark on 4099215, 39028039587479 and 7 is on p. 6); She's conjecture is on p. 8 of the Integers version; Tao's paper is now at arXiv version 3 (same numbering); Ryan 2026 should be cited by its version DOI 10.5281/zenodo.21456694.

## Author

Berkay Yüksel Sayim, Independent Researcher, Germany · ORCID [0009-0004-4993-7352](https://orcid.org/0009-0004-4993-7352) · berksa@tutamail.com

These are partial results and do not settle Erdős problem #364. Generative AI tools were used in this work; the use is disclosed in the preprints.

## License

The papers, the data, the LaTeX sources and the documentation are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), see `LICENSE.md` and `deposit/LICENSE-paper`. The code (the scripts, the build tool and the container recipes) is licensed under the [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0), see `deposit/LICENSE` and `deposit/NOTICE`. Both licenses allow use for any purpose, including commercial use; both require that the source is credited. Programs of others that the scripts call as tools, such as GMP-ECM and PARI/GP, are not part of this repository and keep their own licenses.
