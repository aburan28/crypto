// Results report. Every number is read from data/*.json; nothing here
// computes a result (rounding is presentation only).
#set page(paper: "a4", margin: (x: 1.6cm, y: 1.5cm), numbering: "1")
#set text(size: 9pt)
#set par(justify: true, leading: 0.55em)
#set heading(numbering: "1.")
#show heading: set text(size: 10.5pt)
#let s7 = json("/data/summary_p7.json")
#let s11 = json("/data/summary_p11.json")
#let vd = json("/data/verdict.json")
#let f(x, n: 4) = str(calc.round(x, digits: n))
#let tbl(..args) = text(size: 7.6pt, table(inset: 3pt, stroke: 0.4pt + luma(170), ..args))
#let names = (
  ("X1_cycle_type", [X1: Frobenius cycle type on the Weierstrass points (finest $J[2]$ invariant)]),
  ("X2_sylow2", [X2: 2-Sylow subgroup of $J(FF_p)$ (class group)]),
  ("X3_j4_profile", [X3: $\#J(FF_(p^k))[4]$, $k | 2m$ (frozen descent)]),
  ("X3plus_j4_sets", [X3+: halvable sets, canonical (frozen descent)]),
  ("X3c_j4_profile_corrected", [X3c: X3 corrected (Amendment 1)]),
  ("X3cplus_j4_descent_sets_corrected", [X3c+: X3+ corrected (Amendment 1)]),
)
#align(center)[#text(size: 13pt, weight: "bold")[Which invariant finer than $L$ carries Smith's `full_success`?] \
`research/smith_class_invariant_20261015` · 2026-10-15 · protocol `PROTOCOL.md` frozen at `1461f9445`,
Amendment 1 at `f6917f05b`, both before any null was run]

= Question
At $p = 7$ and $p = 11$ Smith's genus-3 `full_success` clusters within $L$-polynomial isogeny classes beyond what
independent coins per $FF_p$-isomorphism class give (`research/smith_iso_granularity_20261015`, observed class share
$C$ = #f(s7.C) / #f(s11.C), below all 999 iso-level permutations). This study asks whether a Galois-module invariant
finer than $L$ (the cycle type of Frobenius on $J[2]$, the 2-Sylow subgroup of $J(FF_p)$, or Frobenius on $J[4]$)
explains the clustering: conditioning the iso-level permutation null on it should bring $C$ inside the null's central
95%, while an $L$-class placebo of equal cardinality should not.

= Verdicts (frozen decision rule)
#tbl(columns: (31%, 23%, 23%, 23%),
  [*candidate*], [*$p = 7$*], [*$p = 11$*], [*overall*],
  ..names.map(((k, lab)) => (lab, ..(s7, s11).map(s => {
    let c = s.candidates.at(k)
    [#f(c.null.median) [#f(c.null.lo), #f(c.null.hi)]; gap share #f(c.gap_share, n: 3); *#c.verdict*]
  }), [*#vd.at(k).overall*])).flatten(),
)
Cells: null median [2.5%, 97.5%] of $C$ over 999 permutations within (stratum, invariant) cells; gap share
$g_X = ("med" N_0 - "med" N_X) / ("med" N_0 - C)$ with $N_0$ the unconditioned iso-level null
(median #f(s7.null_unconditioned.median) / #f(s11.null_unconditioned.median)). "Invalid" is the frozen rule's
outcome for X3/X3+, whose descent instrument failed control K4 (Section 4).

#figure(image("/figures/class_invariant_nulls.svg", width: 100%),
  caption: [Observed $C$ (dashed) against each conditioned null (blue: central 95%, dot median) and its $L$-class
  placebo (orange), 999 permutations each, seed 20261015. Source `figures/plot.py`, data `data/summary_p*.json`.])

= Measured
#tbl(columns: (40%, 30%, 30%),
  [*quantity*], [*p = 7*], [*p = 11*],
  [tractable isomorphism classes (invariants computed on all)], [#s7.controls.tractable_classes], [#s11.controls.tractable_classes],
  [observed $C$; unconditioned null 2.5–97.5% (K1 reproduction)],
  [#f(s7.C); #f(s7.null_unconditioned.lo)–#f(s7.null_unconditioned.hi) (#s7.k1_reproduces_frozen_null)],
  [#f(s11.C); #f(s11.null_unconditioned.lo)–#f(s11.null_unconditioned.hi) (#s11.k1_reproduces_frozen_null)],
  ..names.map(((k, lab)) => ([cells; share of classes in single-class cells: #k],
    ..(s7, s11).map(s => { let c = s.candidates.at(k); [#c.cells; #f(c.share_in_singleton_cells, n: 4)] }))).flatten(),
  [X3c classes with an ambiguous $k > 1$ entry], [#s7.controls.counts.x3c_classes_with_ambiguous_k],
  [#s11.controls.counts.x3c_classes_with_ambiguous_k],
)

#include "findings.typ"

= Validation
#let ok(c, k) = if c.at(k) [pass] else [*fail*]
#tbl(columns: (52%, 24%, 24%),
  [*control*], [*p = 7*], [*p = 11*],
  [K1: unconditioned null reproduces `summary_p*.json` of PR \#1928 exactly], [#s7.k1_reproduces_frozen_null], [#s11.k1_reproduces_frozen_null],
  [K2: X1/X3/X3c invariants equal on a second orbit member (every class)], [#ok(s7.controls, "k2_pass_j4")], [#ok(s11.controls, "k2_pass_j4")],
  [K2: X2 equal on a second orbit member (1,000-class sample)], [#ok(s7.controls, "k2_pass_sylow")], [#ok(s11.controls, "k2_pass_sylow")],
  [K3: every class certified; $sum v_i = v_2(\#J(FF_p))$; rank $= dim J(FF_p)[2]$], [#ok(s7.controls, "k3_pass")], [#ok(s11.controls, "k3_pass")],
  [K4: frozen descent count at $k = 1$ equals class group (mismatching classes)],
  [#ok(s7.controls, "k4_pass") (#s7.controls.counts.k4_mismatch)], [#ok(s11.controls, "k4_pass") (#s11.controls.counts.k4_mismatch)],
  [K5: frozen profile within the $L$-polynomial bounds (failing entries)],
  [#ok(s7.controls, "k5_pass") (#s7.controls.counts.k5_failures)], [#ok(s11.controls, "k5_pass") (#s11.controls.counts.k5_failures)],
  [K4′: corrected $k = 1$ count brackets the class group], [#ok(s7.controls, "k4c_pass")], [#ok(s11.controls, "k4c_pass")],
  [K5 for X3c], [#ok(s7.controls, "k5c_pass")], [#ok(s11.controls, "k5c_pass")],
  [instrument errors (field membership, characters, relation labels)],
  [#(s7.controls.counts.j4_instrument_errors + s7.controls.counts.sylow_instrument_errors)],
  [#(s11.controls.counts.j4_instrument_errors + s11.controls.counts.sylow_instrument_errors)],
)
K4 failed by design of the mathematics, not of the code: when every Frobenius orbit on the roots has even length,
the even-degree descent map $J(K)\/2J(K) -> L^*\/K^* L^(*2)$ has a kernel of order exactly 2 (Amendment 1 gives the
cohomological argument). The failures sit exactly on the all-even cycle types (8, 4.4, 2.6, 2.2.4, 4.2.2, 2.2.2.2),
and the corrected count never disagrees with the class group (K4′).

= Method
All computation is native Rust (`examples/smith_class_invariant.rs`, `compute` then `summarize`), on the representative
model `y^2 = t f(x)` of every tractable isomorphism class of the frozen rows `iso_p{7,11}.csv.gz`.
- *X2.* Class group of the function field: generators all places of degree $<= 3$, relations $"div"(a)$ for monic
  irreducible $a$ and $"div"(y - v)$ for random $v$ of degree $<= 4$ whose norm $v^2 - t f$ is 3-smooth and
  squarefree; Smith form over $ZZ\/2^(e+1)$, accepted only when its order is $2^(2e+1)$, which with
  $"Pic" = ZZ ⊕ J(FF_p)$ certifies the 2-Sylow subgroup.
- *X3, X3c.* Roots of $f$ in its splitting field $FF_(p^m)$ ($m <= 15$); descent images
  $delta_j(T_S) = product_(i in S, i != j)(e_i - e_j) dot f'(e_j)^([j in S])$, and for swap-type points a rational
  representative from $y - beta g_S(x)$ with $beta^(q+1) = t$; quadratic characters by norms. X3c resolves the order-2
  kernel with the Kummer pairing on the stable part and takes the $k = 1$ count from X2.
- *Nulls.* 999 permutations of `full` within (stratum, invariant) cells, seed 20261015; placebo: each $L$-class gets a
  hashed bucket among as many values as the invariant takes in the stratum.

= Scope
Exact enumeration at $p = 7, 11$, degree-8 models only. X3c+ is an intrinsic invariant built from 2-descent data and
is not claimed to be a complete invariant of Frobenius on $J[4]$; X3c leaves #s7.controls.counts.x3c_ambiguous_entries /
#s11.controls.counts.x3c_ambiguous_entries profile entries as two-valued. A negative verdict excludes these
invariants as explanations at these primes; it does not show that no 2-adic invariant matters. Nothing here builds
or prices an isogeny walk or says anything at cryptographic sizes.

= Reproduce
```
cargo build --release --example smith_class_invariant
B=target/release/examples/smith_class_invariant; D=research/smith_class_invariant_20261015/data
$B compute --p 7 --out $D; $B compute --p 11 --out $D     # 4 min and ~1 h on 2 vCPU
$B summarize --dir $D 7 11
python research/smith_class_invariant_20261015/figures/plot.py   # layout only
```
Run report: `reports/run-20261015T2332Z/report.pdf`. Prior: PR \#1928, `research/smith_iso_granularity_20261015/`.
Source: B. Smith, J. Cryptology 22 (2009), arXiv:0806.2995.
