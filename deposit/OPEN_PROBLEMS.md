# Open problems of the series

Generated on 2026-10-04 by `werkzeug/probleme_bauen.py` from the question environments of Parts I and II and III; the letters are those of the
framework introduction of the collected volume, and the wording is the wording of the parts. The discussion of each problem, that is,
what is known and what an answer would give, follows its statement in the part named. Mathematics is written in LaTeX between `$`;
a bracket such as `[Reinhart2024, p. 6]` is a key of the bibliography that all parts share. The questions of Erdős themselves,
whether three consecutive powerful numbers exist and whether the number of pairs below $x$ is bounded by a power of $\log x$, are
the questions the parts start from and are not repeated here.

## Notation used below

- A positive integer is *powerful* if every prime that divides it divides it at least twice.
- For a squarefree $m \ge 2$, $\varepsilon_m = T_1 + U_1 \sqrt{m}$ is the fundamental solution of $x^2 - m y^2 = 1$, and
  $\varepsilon_m^k = T_k + U_k \sqrt{m}$; $T_k(m)$, $U_k(m)$ name the kernel $m$. $m'$ is the product of the primes of $m$ that do not divide $U_1$.
- The *rank* $\alpha_m(p)$ of a prime $p \nmid m$ is the smallest $k$ with $p \mid T_k(m)$; $M_d$ is the part of $T_d(m)$ supported on the
  primes of rank exactly $d$ (the *primitive part*).
- *Wing (ii)* of Part I: the kernels $m \equiv 7 \pmod 8$ with $m \mid U_1(m)$; they are the weak point of the lower bound for a triple.
- In Part II, the *families* are Walker's Pell families of pairs $n, n + 1$, $F_f$ is the growth factor from one pair of the family $f$
  to the next,
  $w(m) = 1/(2m' \log_{10} \varepsilon_m)$ is the rate of the tower of the kernel $m$, $P_2(N)$ counts the pairs $n - 1, n + 1$ of powerful
  numbers with $n \le N$, and $E_3(N)$ counts the odd squarefree $m$ with $m \mid U_1(m)$, $T_1(m)$ even and $T_1(m) \le N$.
- Reinhart's condition (C) is Definition 1.4 of [Reinhart2024].

## Problem A: finiteness of the kernels of wing (ii)

*Stated as Question I.9.1 in Part I (triples), Section 9.*

Are there only finitely many squarefree $m \equiv 7 \pmod 8$ with $m \mid U_1(m)$ and $T_1(m)$ even?

## Problem B: the primitive part of $T_d(m)$

*Stated as Question I.9.2 in Part I (triples), Section 9.*

Let $m \equiv 7 \pmod 8$ be squarefree with $T_1(m)$ even, let $d \ge 5$ be odd, and let $M_{d}$ be the part of $T_d(m)$ supported on the primes of rank exactly $d$. Can $M_{d}$ be powerful?

## Problem C: a middle that is a power of $2$

*Stated as Question I.9.3 in Part I (triples), Section 9.*

Can $2^a - 1$ and $2^a + 1$ both be powerful for some $a \ge 2$?

## Problem D: Reinhart's condition (C)

*Stated as Question I.9.4 in Part I (triples), Section 9.*

Does Reinhart's condition (C) [Reinhart2024, Definition 1.4] hold for some squarefree $d$, that is, is there an index $k$ with $T_k(d)$ powerful and $d \mid U_k(d)$? For $d \mid U_1(d)$ the second condition holds for every $k$.

## Problem E: convergence of the rates

*Stated as Question II.5.3 in Part II (pairs), Section 5.*

Does $\sum_f 1/\log F_f$, taken over all families, converge?

## Problem F: growth of the tower rate

*Stated as Question II.6.7 in Part II (pairs), Section 6.*

Does $c_M = \sum_{m \le M} w(m)$, the sum over odd squarefree $m$ with $T_1(m)$ even of $w(m) = 1/(2m' \log_{10} \varepsilon_m)$, stay bounded as $M \to \infty$?

## Problem G: the kernels that divide $U_1$

*Stated as Question II.8.1 in Part II (pairs), Section 8.*

Is there a $\theta < 61/180$ with $E_3(N) \ll N^{\theta}$? Since $E_3(N) \le P_2(N)$, Blomer and Schöbel give $E_3(N) \ll N^{\gamma}$ for every $\gamma > 61/180$ [BlomerSchoebel2013, Theorem 3]; the question is whether the squarefree odd $m$ with $m \mid U_1(m)$, $T_1(m)$ even and $T_1(m) \le N$ are rarer than that.

## Problem H: the most frequent gap

*Stated as Question III.5.1 in Part III (further observations), Section 5.*

Let $g(x)$ be the most frequent gap between consecutive powerful numbers up to $x$. Is $g(x)$ powerful for all large $x$, and does the share of powerful values among the $N$ most frequent gaps up to $x$ tend to $1$ when $N$ is fixed and $x \to \infty$?

## Problem I: powerful values of $x^4 + 1$

*Stated as Question III.5.2 in Part III (further observations), Section 5.*

Is $x^4 + 1$ powerful for some integer $x \ge 1$?
