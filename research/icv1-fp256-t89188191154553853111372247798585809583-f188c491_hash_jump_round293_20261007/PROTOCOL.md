# P-256 decorrelated hash-jump factor-base screen, round 293: protocol

Status: preregistered before implementation or execution on 7 October 2026.

## Question and hypotheses

Round 189 is the best of the 289 executed two-delta screens, but its complete
17-sign orbit covers only `0.44017512787806307` of its emitted states and has a
factor-base-corrected selector-stage ratio of `2.320918327943466` to Pollard
rho.  This round asks two separate questions:

1. Can a different choice of all 17 columns repair the two-delta base?
2. If not, can sparse deterministic large jumps between long `+G` runs remove
   the arithmetic-band collision while retaining a sub-rho selector stage?

The preregistered hypotheses are:

- **H1 (two-delta closure):** an exact equal-Hamming-weight support ceiling is
  below the selector gate for every 17-column tuple on the frozen Round 189
  base.  If true, tuple search is closed without pretending that an
  unexecuted combinatorial search found a better tuple.
- **H2 (decorrelation):** at least one hash-jump candidate has exact distinct
  coverage at least `0.9686171464465044` and both its ideal-local and
  delta-frequency-corrected selector ratios below one.
- **H3 (promotion):** even if H2 passes, the candidate is not promoted unless
  it also has a same-family structured residual degree at most five and a
  demonstrably non-generic complete path.  A randomized translation walk is
  rho under another name, not parity for index calculus.

This is a bounded stage experiment.  It does not authorize an unplanted
full-depth P-256 relation unless every promotion gate below passes.

## Frozen identity and inputs

- curve: `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- comparison factor base: `FB1h2f8621cda105`;
- columns `B = 131458`, relation arity and sign depth `17`;
- frozen Round 189 receipt SHA-256:
  `8630dd6ba930c5df7b3b1dd1307bf2fd2289be1d872f0562b5a408c9c4305d0d`;
- frozen executed-screen manifest SHA-256:
  `8723757f965eec5cfb1310de1091c8668815212c31823fd30b976de48ebadf45`;
- Round 29 exact-memoryless closure SHA-256:
  `9be6e5fc38644c34ad746dec9a6542da7e35ef2d8e11eb39a83de7af0f54f443`;
- Round 30 hidden-state receipt SHA-256:
  `8927791b3b149c60a0d4cfb43b77aca21ede9ff21c9c9451a534247ec1f56ad1`;
- Round 31 two-delta receipt SHA-256:
  `9b9ee16c5a24e868d1b9304f08aa80ca39ac23b1b4586ef730534e739172c435`;
- local ideal-oracle ratio `0.964336477130181` and selector coverage gate
  `0.9686171464465044`.

All dependencies are byte-hash checked before a result is written.

## H1: tuple-independent two-delta certificate

Let `R=7021`, let `C=226` be Round 189's emitted interval length minus one,
and let `a_j` be a two-delta coefficient.  With the mechanical rare-edge word,

```text
R a_j = R alpha + rho_j (mod n),  0 <= rho_j < B.
```

For a selected column with residual common capacity `c_j`, the sign-boundary
weight is `w_j = 2a_j+c_j`.  Its anchor-free scaled value is

```text
v_j = 2 rho_j + R c_j.
```

The mechanical-word endpoint condition confines every `v_j` to an interval
of width at most `W=B+R`.  In Hamming layer `k`, complement symmetry therefore
confines the scaled subset starts to width at most
`min(k,17-k) W`.  Splitting by residue modulo `R`, an optimistic union ceiling
for that layer is

```text
min(binomial(17,k)(C+1),
    R(C + ceil(min(k,17-k)W/R) + 2)).
```

The native certificate must derive every term with integer arithmetic, sum all
18 layers, and independently confirm the coefficient scaling and endpoint
identities on all 131,458 frozen coefficients.  It may ignore overlap between
different layers, so it is an upper bound.  The corresponding ideal and
corrected ratio floors charge Round 189's frozen `29,884,448` additions and
must be reported even if they exceed one.

## H2: sparse hash-jump candidates

Evaluate exactly three registered rare-edge counts:

```text
R in {1021, 2047, 4093}.
```

Rare positions use the same mechanical word as the two-delta family.  Common
edges have delta one.  For a given `R` and jump-attempt number, the first
`R-1` rare deltas are independent SHA-256 reductions under the domain

```text
<curve>/hash-jump-round293/<R>/<jump-attempt>/<rare-ordinal>
```

and the final rare delta is the unique value that closes the coefficient
cycle modulo the P-256 subgroup order.  Reject the whole attempt if any rare
delta is zero, one, duplicated up to orientation where that changes the exact
delta inventory, if the coefficient cycle is not distinct, or if no anchor in
the registered 4096-attempt domain makes all columns nonidentity and unique up
to sign.  Record every rejected attempt and reason.  Candidate IDs use
`LDHJR{R}h<first-12-hex-of-coefficient-SHA256>` and the same digest/byte-count
storage fields as the existing `LD2R...` screening receipts; no `FB1` identity
is assigned unless all promotion gates pass.

For each accepted base, hash-draw distinct 17-column tuples under the frozen
Round 33 column domain and accept the first whose common-run capacity reaches

```text
ceil(219 * 6935 / R).
```

This preserves the previous high-capacity selection pressure across rare-edge
densities.  Record draws, rejected tuples, columns, capacities and digest.

Run the complete `2^17` sign orbit with the sign-oriented monotone traversal.
Use exact arbitrary-precision modular interval union, not a hash estimate.
Every factor-base edge and tuple endpoint is replayed in subgroup scalar
arithmetic; direct enumeration on two toy groups and native sign depths 8, 10
and 12 must agree with interval union with zero false positives and negatives.

## Boundary and accounting

The common delta occurs `B-R` times.  Every accepted rare delta must have its
exact multiplicity recorded.  The optimistic exact-delta collision
probability is

```text
p_delta = sum_d (count(d)/B)^2,
hidden_ratio = 0.964336477130181 / sqrt(p_delta).
```

For capacity `C`, charge at least

```text
2^17 C + (2^17-1) + 16 + 17 + 2^17
```

P-256 addition equivalents: path steps, sign boundaries, initial sum,
boundary precomputation and target corrections.  Report:

- emitted and exact distinct states;
- ideal-local and delta-frequency-corrected additions per exact distinct
  state, normalized to rho;
- factor-base construction, delta generation, tuple draws, sorting, merging,
  validation and replay operations separately;
- wall/CPU time, peak RSS, materialized bytes and artifact hashes;
- the cold one-target ratio with setup left **unset** unless every foreign
  operation has a measured addition-equivalent conversion.

No discarded branch is credited as exhaustive coverage.

## Algebraic and non-generic gates

For each candidate, prove or reject all of the following independently of its
coverage:

1. the 131,458 points are unique up to sign, hence have 131,458 distinct
   short-Weierstrass x-coordinates;
2. an exact univariate membership polynomial for the explicit x-set has degree
   at least 131,458 unless a separately verified lower-degree defining
   structure is supplied;
3. structured residual degree of regularity is at most five on a same-family
   complete reference and does not merely omit the membership constraint;
4. the transition does not factor only through a represented group sum as the
   partitioned translation classified by Round 29;
5. relation collection, rank allowance, sparse linear algebra and final
   recovery are all priced in the same unit as rho.

A candidate may pass the selector gate while failing every promotion gate.
That outcome is reported as a stage diagnostic, never as rho parity.

## Success and stop conditions

The execution gate requires exact construction, zero replay failures, complete
full-depth interval unions for all three registered candidates, exact controls,
deterministic JSON, and peak materialized storage below `2^50` bytes.

The selector gate requires exact distinct fraction at least
`0.9686171464465044` and both selector-stage ratios below one.  Promotion also
requires structured residual degree at most five, projected usable relation
cost below `2^103`, collection of 138,031 independent rows below `2^120`
including duplicates/rank and sparse linear algebra, peak storage below
`2^50`, and a complete non-generic DLP below matched one-target rho.

Stop without an unplanted full-depth relation if no candidate passes every
promotion gate.  Publish the exact boundary table and dominant obstruction,
including a selector-stage success if it only recreates generic rho or destroys
the low-degree factor-base representation.

## Deliverables

Implement and test in Rust.  Preserve the protocol commit before executing the
experiment.  Emit a compact JSON result, isolation receipt and result report;
update every affected P-256 dashboard panel; validate the generated site; and
publish the work in a stacked pull request.
