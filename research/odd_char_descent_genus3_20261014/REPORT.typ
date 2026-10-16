#set page(paper: "a4", margin: (x: 1.5cm, y: 1.4cm), numbering: "1 / 1")
#set text(size: 8.6pt)
#set par(justify: true, leading: 0.5em)
#set heading(numbering: none)
#show heading.where(level: 1): set text(size: 12pt)
#show heading.where(level: 2): set text(size: 9.6pt)
#show heading: set block(above: 0.8em, below: 0.4em)
#let T = text(fill: rgb("#2a78d6"), weight: "bold")[T]
#let E = text(fill: rgb("#1a8a5f"), weight: "bold")[E]
#let H = text(fill: rgb("#c26b00"), weight: "bold")[H]
#let O = text(fill: rgb("#52514e"), weight: "bold")[O]

= Odd-characteristic Weil descent onto genus 3

2026-10-14. This is proposal 5 of `problems/ghs-to-smith-pipeline-characteristic-mismatch.md`.
It combines a literature reading with an exact native enumeration. It
solves no discrete logarithm and claims no cost. Tool:
`examples/odd_char_ghs_genus.rs` (std-only Rust). Runs `20261014T1758Z` and
`20261014T1800Z` are under `research/odd_char_descent_genus3_20261014/results/`.
Protocol `PROTOCOL.md` was frozen at commit `4ad7dd5db`, with one amendment.
Tags: [#T] theorem (source or derived here), [#E] exact computation, [#H]
heuristic, [#O] out of scope or not established. Conductor is not
installed on this host. The full report with every quotation is
`REPORT.md`.

== Answer
- [#T, source] *Yes, for $n = 3$.*
  - The curves are $E\/FF_(q^3): y^2 = e g(x)(x-alpha)(x-alpha^q)$, with
    $alpha in FF_(q^3) without FF_q$ and $g in FF_q [x]$ of degree 1 or 2.
  - They have a GHS cover of degree 8 over $PP^1(x)$ by a genus-3
    *hyperelliptic* curve over $FF_q$, for every odd $q$.
  - There are $Theta(q^2)$ isomorphism classes of such $E$.
  - Sources: Momose–Chao, ECC 2007, slide 24; Joux–Vitse, Eurocrypt 2012,
    §4.1, eq. (3), which calls them "the only" such curves, citing Diem
    2003, Momose–Chao 2005 and Thériault 2003.
- [#E] On all 81 instances ($p <= 31$, three $alpha$ per $p$, three choices
  of $h$), Joux–Vitse's explicit $H: y^2 = 4N^2 h(F\/(4N))$ (App. A.1) has
  $L_H(T) = 1 - t T^3 + p^3 T^6 = L_(op("Res") E)$. This also matches the
  repository's own $FF_7$ instance ($328$; $N = 8, 50, 296$).
- [#T, derived] $g = 1 + 2^(m-2)(r-4)$ for every degree-2 $x$, and the
  genus is independent of $q$. Genus 3 needs $(m,r) in {(2,6),(3,5)}$ and
  a branch-set period of 2 or 3. If $FF_q (j) = FF_(q^n)$ with $n >= 4$,
  then every $x$ gives $g >= 5$.
- [#T] The problem record is *not* closed as a proved obstruction, because
  a family exists. The characteristic-2 pipeline as stated stays blocked:
  Smith's construction "does not apply in characteristics 2 and 3" (Smith
  2009, p. 2).

== Sources (exact locations)
- *Chao–Momose, ECC 2007 slides.*
  - sl. 17: Condition (C), "$g = d g_0$, the smallest possible genus".
  - sl. 20: for char $!= 2$ and even $d$, covers exist only at
    $(d,n) = (2,2)$; for odd $d$, at $(5,4)$, at $d = 2^n - 1$, and
    (decomposable) at $(3,3)$.
  - sl. 24: the hyperelliptic family.
  - sl. 25: the non-hyperelliptic Types 1 and 2, covering degree 4.
  - sl. 26: when the cover is hyperelliptic.
  - sl. 33: $\#{lambda} approx q^3\/2$.
  - Proofs not read. [#T]
- *Joux–Vitse, Eurocrypt 2012 (author PDF).*
  - p. 8, eqs. (3) and (4): $Theta(q^2)$ curves, all with $4 | \#E$.
    Isogeny-walk conjecture [#H].
  - pp. 14–15: the explicit $H$.
  - p. 16 (App. C): for $n = 6$ and $j in.not FF_(q^3) union FF_(q^2)$,
    "the minimum genus of the Fq-cover obtained by the GHS method is 9". [#T]
- *Diem, Compositio 147 (2011).*
  - p. 94: "$g(C) = sum_L g(L)$, where L runs over all subextensions of
    k(C)|k(x) of degree two" (Kani–Rosen).
  - Prop. 4.10(a): $g <= (2^n-1)(2n-1)$.
  - p. 95: for odd $q$, $k(E)|k(x)$ is tamely ramified at 4 places.
  - p. 78: Gaudry's $tilde(O)(q^(2-2\/n))$ is heuristic [#H].
- *Hess, LMS JCM 7 (2004).*
  - p. 168: the odd-characteristic generalisation is due to Diem.
  - p. 171: "analogous results can be obtained for Kummer extensions".
  - p. 174, Thm 8: the genus is exponential in $m$ (Artin–Schreier). [#T]
- *Smith, arXiv:0806.2995.*
  - p. 2: the restriction to characteristic at least 5.
  - Prop. 1 and Def. 1 (pp. 4–5): tractable subgroups.
  - Lemma 1 (p. 6): $\#S(H)$ by factorisation pattern.
  - Props. 2 and 4 (pp. 10, 15): two square conditions.
  - Hypotheses 2–3 [#H].
- *Tian, arXiv:2012.07173.*
  - p. 2: prime $n in.not {3,5,7}$: "the attack fails"; genus 3 at $n = 3$
    "necessary condition … three rational points of order 2".
  - p. 18: $(ell,ell,ell)$-isogenies $op("Res")(E) -> op("Jac")(C)$ for
    prime order. [#T]
- *Abstracts only:* ePrint 2009/236, 2015/805, 2004/073, 2005/119 and
  2001/054.
- [#O] *Not read*, with no copy reachable from this host:
  - Diem, "The GHS attack in odd characteristic", JRMS 18 (2003);
  - Thériault, JRMS 18 (2003);
  - Diem–Thomé, J. Cryptology (2008);
  - Gaudry, J. Symb. Comput. (2009);
  - the ePrint PDFs (HTTP 403).

== Theorems (derived here; proofs in REPORT.md) and enumeration
Set $S = union.big_(i<n) sigma^i B$, $r = |S|$, and
$m = dim op("span")_(FF_2) {1_(sigma^i B)}$; let $k$ be the period of $B$.
- *Theorem 1.* $g = 1 + 2^(m-2)(r-4)$, by Riemann–Hurwitz for a tame
  $(ZZ\/2)^m$ cover. It equals Diem's subcover sum; [#E] checked on every
  type for $n <= 14$.
- *Theorem 2.*
  - $g = 3 <==> (m,r) in {(2,6),(3,5)}$.
  - $m = 2 ==> k in {2,3}$, and $(m,r) = (3,5) ==> k = 3$.
  - $k >= 4 ==> g >= 5$.
- *Theorem 3.* For even $n >= 6$ no type of any period has $g = n$: the
  genus would force $m = 2$ and $r = n+3$, but $r <= 8$.
- [#E] 601 types for $n <= 8$. Exact branch-set counts for 144 pairs
  $(q,n)$ ($q = p, p^2$, $p <= 31$) add up to $q^(4n) - q^(2n)$; the
  completeness check also passes at $q = 5$ for $n <= 14$. 11 unit tests
  reproduce:
  - Momose–Chao's minimal-cover list;
  - the four $n = 3$ genus-3 shapes;
  - Joux–Vitse's two genus-9 shapes;
  - Diem's bound;
  - Smith's Lemma 1 table.

#table(
  columns: 15, inset: 2.5pt, align: center, stroke: 0.4pt + rgb("#bbbbbb"),
  [$n$], [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12], [13], [14],
  [min $g$, period $n$], [1], [2], [3], [5], [5], [9], [7], [33], [81], [49], [1793], [25], [9217], [41],
)

#figure(image("figures/genus_vs_n.svg", width: 100%),
  caption: [(a) Genera of period-$n$ branch types, the minimum, the generic
  $(n-1)2^n+1$, and $g = n$; genus 3 is circled and occurs only at $n = 2, 3$.
  (b) $log_10$ share of degree-4 branch sets with genus 3. Run `20261014T1800Z`.])

== The hyperelliptic family and Smith (third arm)
- [#E] The Weierstrass cycle type and Smith's count:
  - $(3,3,1,1)$ for $h = x$ and $h = x(x-1)$, so $\#S(H) = 3$ rational
    tractable subgroups;
  - $(6,2)$ for irreducible $h$, so $\#S(H) = 1$;
  - all 81 instances also have $4 | \#E(FF_(p^3))$.
- [#T, derived] Why: $sigma$ normalises the bielliptic group
  ${1, phi, phi^sigma, phi^(sigma^2)}$, so it acts on each free orbit with
  exactly one fixed point.
- [#O] Smith's Props. 2 and 4 are not checked. It is not checked that
  $J_H\/S$ is a smooth plane-quartic Jacobian. No isogeny is built.
- [#T] Types 1 and 2 give plane quartics directly, so Smith matters only
  for curves whose sole genus-3 cover is hyperelliptic. The census finds 20
  such curves at $q = 5$ and 42 at $q = 7$.
- [#O] Conflict with Tian's "three rational points of order 2": with
  irreducible $h$, $\#E(FF_(p^3))[2] = 2$ on all 27 instances, yet the
  L-check passes. Not resolved, because the cited sources were not read.

== Weil restriction directly ($n = 3$)
- [#T, standard] With its own polarization, $op("Res")(E)$ is geometrically
  $E times E^sigma times E^(sigma^2)$ with the product polarization, so it
  is never a smooth-curve Jacobian.
- [#T, source + E] It is $FF_q$-isogenous to a genus-3 Jacobian
  - via $(2,2,2)$ GHS covers (hyperelliptic family; Types 1 and 2);
  - or via Tian's $(ell,ell,ell)$-isogenies for prime order.
- [#T] The target can be hyperelliptic when $4 | \#E$. For prime order,
  Tian (citing [19]) reports non-hyperelliptic only.
- [#O] A full characterisation is not established here.

== Genus-$k$ covers, $k = 6, 10, 14$ (`docs/curves/jacobian-transfers.json`)
- [#T, source] Joux–Vitse p. 16 give minimum GHS genus 9 at $k = 6$.
- [#T, source] Momose–Chao sl. 20: no $g = d$ cover for even
  $d != 2$.
- [#T + #E] By Theorem 3 there is no genus-$k$ type at all. For
  $FF_p (j) = FF_(p^k)$ the minima are $9, 49, 41$, for every odd $p$.
- *Plainly:* the 846 `cover_descent` rows priced at $g = k$ use a genus the
  GHS construction does not deliver.
- [#O] Non-GHS covers are not ruled out. The file was not edited.

== Does not show
- [#O] Costs, logarithms or rho comparisons.
- [#O] $FF_q$-models for even $n$.
- [#O] Smith's Props. 2 and 4.
- [#O] Characteristic 2 or 3.
- [#O] Lean: not formalized, because there is no toolchain on this host.
