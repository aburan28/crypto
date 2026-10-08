# ISO-1 weak isogeny classes: trace census and invariant audit (2026-10-07)

## Requested question and scope

The requested base primes are **p = 37 through about 200**, with
q = p² and curves over F_(q³) = F_(p⁶).  The requested result is a label for
every ordinary full-2-torsion isogeny class at every such p, followed by a
formula that can replace both the cover-branch reach census and random
seed grinding.  This round completed **p = 37** in that range and used
p = 11, 13, 17 as development and held-out checks.  The remaining primes
and a necessary-and-sufficient formula remain open.

| Requested item | Status | Evidence or gap |
|:--|:--|:--|
| Every trace at p = 37 | Measured, probabilistic point-count assignment | [p37 trace rows](p37_twist_derived.csv), 50,654 Hasse candidates and 49,284 ordinary rows |
| Every trace at p = 41 through about 200 | Not completed | Direct normalized-representative enumeration is too costly at that range with the present point counter |
| Fit 2-splitting, Frobenius-order 2-depth, trace residues, class-number parity | Completed on p = 11, 13; tested on p = 17 and p = 37 | [held-out p37 fit](fit_p11_p13_to_p37.txt), [held-out p17 fit](fit_p11_p13_to_p17.txt) |
| Actual endomorphism-ring 2-volcano levels | Not measured | A curve's level needs its endomorphism order; `v₂(f_pi)` is only a class-wide upper bound |
| Exact criterion replacing reach census and seed sieve | Not established | The strongest simple condition has 290 false positives at p = 37; see counterexample below |

## Field, representatives, and trace assignment

The field implementation is `F_p[u]/(u²-w)` followed by
`F_(p²)[theta]/(theta³-s)`, with the first nonsquare `w` in F_p and the first
noncube `s` in F_(p²) chosen by the code in `jv_cover.rs`.  A weak curve is
`y² = x(x-alpha)(x-sigma(alpha))`, where `sigma` is the p²-Frobenius and
`alpha` is outside F_(p²).  Translation and F_(p²) square scaling leave
`2q²+2q` normalized representatives.  For p = 37 this is **3,751,060**.
The run evaluated one representative from each quadratic-twist pair and
derived the opposite trace, so it made 1,875,530 point-count calls.  A
full four-branch run and the twist-derived run agreed byte for byte at
p = 13; at p = 7 two individual point-count assignments failed independent
square-table audits: provisional [346 became -686](p7_trace346_audit.txt)
and [554 became -430](p7_trace554_audit.txt).  The trace labels are consequently
**probabilistic computational labels**, not independently certified
cardinalities of every representative.  The exact normalized representative
count and the twist trace identity are algebraic.
One positive p = 11 witness at provisional trace 38 had
[exact square-table trace 38](p11_trace38_audit.txt); its fresh
twist-derived run agreed with the full run on all weak versus zero labels.

The pre-existing census used `w` from F_p as its supposed nonsquare in
F_(p²).  Every nonzero F_p element is a square in F_(p²), so that run
duplicated one branch and omitted the other.  The source now selects an
actual F_(p²) nonsquare.  Prior p = 7–31 weak-class counts and reach
fractions are historical and require correction; p = 23 and 31 have not
yet been rerun with the corrected representative set.

The trace `t = p⁶+1-#E(F_(p⁶))` identifies the F_(p⁶) isogeny class.
The CSV includes every `t ≡ 2 (mod 4)` in the Hasse interval, including
zero-witness rows.  The `trace_status` column separates ordinary traces
from boundary and p-divisible traces.  These are class diagnostics, not
individual curve or DLP subgroup manifests, so no `IC1` candidate or
speedup measurement is asserted.

## Measured results

| p | q | normalized weak representatives | ordinary trace rows | weak rows | weak rows with v₂(f_pi)=1 | zero rows with v₂(f_pi)>=2 | random full-2 curves in weak class (n=4,000, Wilson 95%) |
|---:|---:|---:|---:|---:|---:|---:|:--|
| 11 | 121 | 29,524 | 1,210 | 542 | 0 / 606 | 62 / 604 | 2,363 / 4,000 = 0.5908 [0.5754, 0.6059] |
| 13 | 169 | 57,460 | 2,028 | 928 | 0 / 1,014 | 86 / 1,014 | 2,324 / 4,000 = 0.5810 [0.5656, 0.5962] |
| 17 | 289 | 167,620 | 4,624 | 2,198 | 0 / 2,312 | 114 / 2,312 | 2,480 / 4,000 = 0.6200 [0.6048, 0.6349] |
| **37** | **1,369** | **3,751,060** | **49,284** | **24,352** | **0 / 24,642** | **290 / 24,642** | **2,433 / 4,000 = 0.6083 [0.5930, 0.6233]** |

The [p = 37 run receipt](p37_run.txt) and CSV sum to 3,751,060 representative counts.  Its 50,654
candidate rows include 1,370 nonordinary or boundary rows, none with an
observed weak representative.  The p = 37 census took **1,106.807 s**
with `RAYON_NUM_THREADS=4`, and recorded **563,856,133,786** charged F_p
multiplications.  This is a class-label diagnostic on a contended host;
it is not an isolated CPU speed comparison.

The two denominators differ: **24,932/49,284 ordinary trace rows (50.59%)**
have no observed weak representative, while **1,567/4,000 sampled random
full-2-torsion curves (39.18%)** land in those zero rows.  Larger
isogeny classes carry more sampled curves, so the class fraction and the
curve-weighted reach fraction should not be interchanged.

![Observed class strata and invariant flow](class_strata.svg)

The [editable visual source](../../examples/iso1_class_census.rs) reads
the four census CSVs.  The bars show complete observed trace counts, so
sampling intervals do not apply to them.  The 4,000-curve reach fractions
above do have binomial sampling intervals.

## Invariant fit and counterexamples

Write `D = t²-4p⁶ = f_pi² D_K`, with `D_K` the fundamental CM
discriminant.  The field and class trace determine `f_pi`; it is the
**Frobenius-order** conductor, not the conductor of a particular curve's
endomorphism ring or its 2-volcano level. It bounds the possible 2-adic
conductor level within the class, but the census has no per-curve
endomorphism-order proof. For ordinary `t ≡ 2 (mod 4)`,
the following two tests are arithmetically equivalent:

`v₂(f_pi) >= 2`  iff  `(t/2)² ≡ p⁶ (mod 16)`  iff
`t/2 ≡ ±p³ (mod 8)`.

All observed weak classes at p = 11, 13, 17, 37 pass this condition.
The first equivalence is an arithmetic identity for these trace
discriminants; **necessity for all weak curves at arbitrary p is an
empirical conjecture, not a proof here**.  The condition is not
sufficient: the p = 37 held-out table gives 24,352 true positives,
290 false positives, zero false negatives, and 24,642 true negatives,
or 48,994/49,284 correct labels (99.41%).  Splitting of 2 and
maximal-order class-number parity fit substantially worse; combining
those variables and a finer 2-adic residue did not repair the false
positives.

Two p = 37 classes exhibit the obstruction to a formula from just the
proposed small invariants.  At `t = -92,218` there are **zero** observed
weak representatives, whereas `t = 38,854` has **24**.  Both have
`v₂(f_pi)=3`, ramified 2, even maximal-order class number, and the same
trace modulo **2¹⁷**; the traces differ by exactly 131,072.  PARI/GP
[`qfbclassno` output](cm_counterexample_gp.txt) gives class numbers 912
and 3,168 respectively, both even.
At p = 17, `t = 9,714` has zero weak representatives while `t = 9,586`
has 24, with the same depth 3, inert 2, even class-number parity and
trace residue modulo 64.  Thus no rule using only those feature values
can classify these rows.  A trace value itself remains a class invariant;
the missing result is a compact predictive criterion.

The residual p = 37 zero rows concentrate near the Hasse edge:
[234 of 290](p37_edge_bins.csv) lie in the outer fifth of `|t|/(2p³)`, where
234/4,930 high-depth rows have no observed weak representative (4.75%).
The inner fifth has 6/4,926 (0.12%).  This association is measured;
it does not prove that small isogeny-class size causes each absence.
There are central counterexamples, and the more extreme `t=-100,410`
has 12 weak representatives while `t=-92,218` has zero.

As a control, 445/2,000 sampled p = 11 curves with a square Legendre
parameter still had Frobenius depth 1.  The observed condition is
therefore more specific than the simple fact that the norm-one
parameter is a square.

## Reproduction and limits

The executable is `examples/iso1_class_census.rs` built with
`cargo build --release --example iso1_class_census` on arm64, rustc
1.93.1, from base commit `9fbdf13c17184ff0866e61fa326adcfd477b041c`
plus the changes in this worktree.  The p = 37 run command was:

```sh
RAYON_NUM_THREADS=4 target/release/examples/iso1_class_census 37 research/iso1_weak_classes_20261007/p37_twist_derived.csv --derive-twists
```

The smaller CSVs used the full four-branch mode without
`--derive-twists`.  The fit used p = 11 and 13 for training and p = 17
or 37 as holdout.  Reach samples used `sample P CSV 4000`, with
`StdRng` seed `1 ^ 0xE8AC7`.  The p = 37 executable's pre-report SHA-256
was `923d687011fd804adb2fae9c105509f4b2b337b0187a7b0d43d9cdba67864b55`.
The program records failed and zero-witness trace rows in each CSV.

At p = 199, the same normalization would require 3,136,557,604
representatives, 1,568,278,802 direct point counts with twist symmetry.
Scaling the p = 37 run by the representative count and the point
counter's approximate `p^(3/2)` baby-step cost gives roughly **130 days**
for that prime alone on comparable resources.  This is an extrapolation,
not a timed p = 199 run.  The full p = 37–200 request therefore needs
a different exact class criterion or a substantially faster trace
algorithm.  Neither was established in this round.

The original Joux–Vitse paper gives the weak model and explicitly treats
independence of weak form and isogeny class as an assumption, rather than
proving universal reach: [Joux and Vitse, §4.1](https://eprint.iacr.org/2011/020.pdf).
The broader Legendre-family reach theorem concerns all order-divisible-by-four
classes and does not identify this norm-one weak subfamily:
[Auer and Top](https://arxiv.org/abs/math/0106273).
