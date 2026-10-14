# Protocol: why does the (7,4) weak curve of #1757 have two weak lines?

Frozen 2026-10-13, before any new code was run on any input below. Follow-up
to `research/weak_subfield_hyperelliptic_decomposition_20261012/` (merged in
#1757), whose README leaves one item open: the E3 weak curve at (q,n) = (7,4)
has two weak affine lines, one of which is not a scaling line through the known
`β`; every other field has exactly one. That study's frozen records are not
edited. Conductor: no control plane was reachable from this session
(`coord_start_work` returned 404), so there is no task id; the intended scope
is `dir:research/weak_line_7_4_followup_20261013` only. Structural diagnostic
on toy fields; no speedup or sub-rho claim.

## Setting

`K = F_{q^n}`, `f(x) = β⁵ g(x/β)` with `g ∈ F_q[x]` monic quintic, `β ∈ K \ F_q`.
A weak line is `L = {β't + b : t ∈ F_q}` with `c·f(β't+b) ∈ F_q[t]` for some
`c ∈ K*`; `census` in `weak_curve_lab.rs` counts these lines as sets.
Write `ψ1 = g`, the F_q-polynomial of the known line. Define
`Aut(ψ) = {ν = (ζ,c) ∈ Aff(K) : ψ(ζs + c) = ζ⁵ψ(s)}`, the affine symmetries of
the monic quintic (the hyperelliptic involution is trivial on `x`). `φ` is the
`q`-Frobenius.

## Competing hypotheses and tests

- **H1 (affine symmetry).** `g` has a nontrivial `ν ∈ Aut(g)`, and a second
  line follows from Lemma A below. Test: enumerate `Aut(g)` over `K` exactly.
  Falsifier: `Aut(g)` trivial while the census finds two lines, or the
  Lemma-A count differs from the census count.
- **H2 (quadratic subfield).** The second line comes from `F_{q²} ⊂ F_{q⁴}`: a
  twist trivialised over `F_{q²}`. Test: take the K-affine map `μ` between the
  two lines; H2 predicts its scalar lies in `F_{q²}`, and that the two
  `F_q`-forms `y² = ψ_i(t)` have Frobenius data related by a twist. Falsifier:
  the scalar of `μ` lies outside `F_{q²}`. H2 is not exclusive with H1: the
  data may show H1 acting through `F_{q²}`. H2 without any symmetry of `g` has
  no mechanism if Lemma A holds; it is then falsified together with H3.
- **H3 (coincidence).** The second line is an accident not tied to symmetry:
  `Aut(g)` is trivial and the line is nonetheless weak. Falsifier: `Aut(g)` is
  nontrivial. If H3 survived, `census` would be suspect, and an independent
  check would be added (the line is verified directly: `c·f(β't+b) ∈ F_q[t]`).

**Lemma A (derived here, tested by the controls).** Let `ψ1 ∈ F_q[t]` be a
monic quintic and `f = ψ1∘σ1⁻¹` for an affine `σ1` (the known line). Lines
`σ1(μF_q)`, `μ ∈ Aff(K)` modulo right multiplication by `Aff(F_q)`, are weak
iff `ν := μ^φ μ⁻¹ ∈ Aut(ψ1)` (apply `φ` to `c·ψ1(μt) ∈ F_q[t]`), and by
Hilbert 90 for `Aff = G_a ⋊ G_m` a given `ν` occurs iff its norm
`N(ν) = φ^{n-1}(ν) ∘ … ∘ φ(ν) ∘ ν` is the identity. Hence

    #weak lines of f  =  #{ ν ∈ Aut(g) : N(ν) = 1 }.

The composition order is fixed by the code, and the count is checked against
the census (so a wrong orientation shows up as a mismatch). Consequences to
test: (i) the count is at least 1 (take `ν = 1`); (ii) more than one line
forces `Aut(g) ≠ 1`; (iii) any `ν ≠ 1` fixes `∞` and permutes the 5 roots, so
has order dividing 4 or 5 (char 5 excluded: translations); (iv) for `ν = (−1,
2x0)` over `F_q` (`g` odd about `x0 ∈ F_q`), `N(ν) = ν^n`, so the extra line
exists iff `n` is even, and then it is `x0 + δF_q` with `δ ∈ F_{q²}`,
`δ^q = −δ`.

## Experiments

- **F1 reproduction.** Rebuild the E3 (7,4) weak curve from the same stream
  (tag 300+74) and compare its JSON record with the recorded line byte for byte.
- **F2 the two lines.** For each weak line: `β` index, canonical `b`, `c`, the
  monic `ψ_i`, Frobenius data `(a1,a2)`, `#J` over `F_q`, `F_{q²}`, `K`
  (Frobenius-polynomial route and, for `K`, point counts). The K-affine map
  `μ` between the lines; tests of `ψ2 = λ ψ1(αt+γ)` over `F_q`, `F_{q²}` and
  `K`; `Aut(g)` by a full scan over `ζ ∈ K*` (`c` is forced by the `x⁴`
  coefficient, `c = g4(ζ−1)/5`), with each element's field of definition; the
  Lemma-A prediction.
- **F3 control census.** Fresh weak curves from stream tags `7000+…` with the
  E3 sampler (monic `g ∈ F_q[x]`, nonzero lower coefficients, squarefree,
  `β ∉ F_q`). Even n: (3,2) (7,2) (11,2) (13,2) (3,4) (7,4) (11,4) (13,4)
  (3,6); odd n: (3,3) (7,3) (11,3) (3,5). 300 curves per field (2000 at
  (7,4)). Per curve: census line count, Lemma-A count, `|Aut(g)|`. The census
  is the shortcut `b = −f4/5` (exact: the `x⁴` coefficient of `f(β't+b)/β'⁵`
  lies in `F_q` iff `b ∈ −f4/5 + β'F_q`), cross-checked against the exhaustive
  `census` on every curve when `|K| ≤ 2401` and on the first 5 curves
  elsewhere. Char 5 is excluded (Lemma A and translations).
- **F4 population by `g`.** The line count depends on `g` and `(q,n)` only. For
  all monic squarefree `g ∈ F_q[x]` with nonzero lower coefficients (q = 3, 7,
  11, 13) and n = 2, 3, 4 (and 6 for q = 3), the Lemma-A count, with `ζ`
  restricted to the roots of unity of order 2, 4, 5 in `K` (justified by
  consequence (iii), cross-checked against the full scan on F2/F3 curves).
  Output: fraction of `g` with more than one line, by symmetry type.
- **F5 positive control.** Constructed `g` odd about a random `x0 ∈ F_q`
  (`(x−x0)⁵ + a(x−x0)³ + b(x−x0)`, `a,b ≠ 0`; squarefree) at (3,2) (7,2) (3,4)
  (7,4) (3,6): Lemma A predicts 2 lines; at (3,3) (7,3) (3,5): 1 line. And
  `x⁵+c` over `F_{11}` at n = 2, 3, 4 (`ζ₅ ∈ F_11`, `N(ν) = ν^n`): 1 line, in
  contrast with E3b's (3,4) result (5 lines). The detector must see the effect
  the explanation predicts, or the explanation is untested. Census counts are
  by the exhaustive `census` for `|K| ≤ 2401`.

## Success criteria

H1 is supported if the Lemma-A count equals the census count on every F2, F3
and F5 curve with the (7,4) curve having `Aut(g) ≠ 1`. It is falsified by any
mismatch (reported, and the lemma abandoned) or by trivial `Aut(g)` on the
(7,4) curve. Counts are reported as exact tallies; there is no significance
claim, only agreement or disagreement of an exact prediction.

## Not measured

Timing; fields above `|K| = 117649`; relation spans of the second line beyond
the 1.2·10⁷ cap of the parent lab.

## Language

Computation in Rust (`weak_line_lab.rs`, dependency-free, `rustc --edition 2021
-O`; a fork of `weak_curve_lab.rs` whose field, polynomial, Cantor and census
code is kept verbatim). Python only plots recorded values.
