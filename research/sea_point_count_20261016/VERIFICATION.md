# What is verified, and what remains open

This verification concerns the combinatorial change in PR #1923. It is not a
claim that the elliptic-curve implementation has been fully formally verified.
The TLA+ files are specifications, not generated Rust or extracted code. Lean
checks mathematical statements about an explicitly lawful additive group;
there is no automatic refinement proof connecting that structure to padd.

## Model and code map

| Rust behavior | TLA+ action/invariant | Lean statement | Coverage and premise |
| --- | --- | --- | --- |
| `progression_bsgs` retains all baby indices | `BabyStep`, `BabyComplete` | `bsgs_complete` | Complete multimap; lawful group; width positive |
| Giant steps enumerate bounded indices | `GiantStep`, `SearchPrefix` | `bsgs_complete` | Quotient/remainder decomposition, including partial last blocks |
| Binary-search intersection requires sorted hits | `BabyOrdered`, `HitsOrdered` | `giant_blocks_ordered`, `signed_progression_ordered` | Increasing baby insertion, strictly separated giant blocks, positive trace stride |
| Old one-index table drops collisions | `Legacy=TRUE` | Not asserted correct | Negative control loses the true candidate for a nonzero point in a group of order two |
| Atkin sets filter search hits | `allowed` in `FinishSearch` | Candidate-preservation premise | The true residue must be in each Atkin set; computation of those sets is not proved |
| Full search or scalar refinement on the next point | `NextSample`, `TrueSurvives` | `refinement_preserves`, `direct_equals_full_intersection` | Every sample test accepts the true order; the heuristic only selects a strategy |
| Successive sampled tests retain truth | `TrueSurvives` across rounds | `refinements_preserve` | Induction over any number of tests, under the stated sample premise |
| Return when one candidate survives | `ReturnUnique`, `UniqueSound` | `unique_sound` | True candidate retained, singleton candidate set |
| Bounded ambiguous exit | `ReturnAmbiguous`, `Terminates` | No unbounded termination claim | Finite round bound and weak fairness of `Next` |
| Order certificate forces a unique Hasse-interval value | `Certification.Sound` | `interval_unique` | Genuine common divisor exceeds interval width; both orders lie in interval |
| Sampled annihilation alone | `Certification.NontrivialSampleSound` | Deliberately not a theorem | Negative control: claimed four annihilates a nonzero point in a group of order two |

`GroupOps` states the ordinary group laws as structure fields; it does not
assume the search result. `smul_add` and `subtract_left` are proved from those
laws. `bsgs_complete` then proves the full equivalence between enumerated hits
and the scalar equation, for arbitrary interval length and positive width.
It is not merely a theorem about one tested parameter set.

The library's exhaustive point tests cover the Rust search on actual small
prime-field curves, including repeated infinity. These executable comparisons
support the implementation/model correspondence, but do not constitute an
unbounded refinement proof of the Rust source or compiler output.

The model abstracts performance cutoffs and the cost heuristic by allowing
either valid refinement strategy. It exercises BSGS at small interval counts
to explore its loop structure; Rust normally chooses the linear walk below
its 1024-candidate cutoff. Thus the count-two legacy trace is a loop-level
counterexample, not a claim that the public SEA entry point chooses BSGS for
two candidates. The Rust exhaustive regression also tests length 1025, and
the Lean search theorem has no finite count bound.

## Checked configurations

- `Fixed.cfg`: cyclic groups of order 1 through 6; interval counts 1 through 8;
  strides 1, 2, 3; every true index and first point; three forms of
  truth-preserving filter; up to two follow-up samples; both strategies.
  165,750 distinct states, 209,680 generated; invariants and termination passed.
- `Wide.cfg`: groups 1 through 9; counts 1 through 12; the same strides,
  points, indices and filters; first search and bounded return. 310,500 distinct
  states; invariants and termination passed.
- `Certification.cfg`: interval endpoints, actual and claimed orders, and point
  representatives through order 16. 58,956 distinct states; soundness passed.
- `Legacy.cfg`: expected invariant failure, exit 12. The retained trace has
  group order 2, nonzero point 1, stride 2, true index 1, and count 2. Index 1
  is discarded when both baby keys are infinity; the filter then loses truth.
- `SampleOnly.cfg`: expected invariant failure, exit 12, using a nonzero point.
  This explains why the previous CM constructor's sampled checks were
  insufficient. The corrected constructor requires exhaustive counting or
  successful curve/twist Hasse certificates, and returns its evidence.

These are exhaustive finite-state checks within the specified bounds, using
TLC's state fingerprinting. The tool's collision-risk estimates are retained in
the logs. They are not an unbounded theorem about elliptic curves.

## Proof tools and trust boundary

Lean 4.19.0, commit `6caaee842e94`, with no package dependencies beyond its
standard library. `lake-manifest.json` records an empty dependency list.
The final proof suite contains ten theorem declarations. Axiom inspection
reports only standard Lean `propext` and `Quot.sound` where used; no reported
theorem depends on `sorryAx`, admitted proofs, or a custom research axiom.
The initial unsuccessful induction proof and its failed build log are retained
separately; they are superseded by the successful proof, not reported as checked.
An independent `lean4checker --fresh` executable is not installed here.

TLC comes from the v1.8.0 release URL with SHA-256
`15d103b7e5617780bcc904c0af12f1ad242e2b345c0301136c5be4588f17edac`.
Its reported build is 2026.10.15.145915, revision `932d35e`. CI verifies that
exact artifact hash, reruns positive configurations, and requires both negative
controls to fail with the expected invariant and exit status.

The proofs concern group search, set refinement, and interval uniqueness. No
new polynomial identity is asserted, so the polynomial-evaluation
`identity.certificate/v1` format is not substituted for these Lean proofs.

## Open implementation obligations

1. Prove that the actual bigint/field representation and padd implement a
   lawful elliptic-curve group, including inverses, infinity and exceptional
   cases. Link the proof to concrete Rust operations and integer bounds.
2. Establish the correctness of trace parity, Elkies cycles and Atkin
   candidate generation. The current proof assumes those constraints retain
   the true trace; it does not prove that assumption for the Rust code.
3. Connect point-order divisor accumulation and factor validation to genuine
   divisibility of the actual group order. Probable-prime tests are a separate
   trust assumption; Hasse uniqueness alone cannot validate a false divisor.
4. Prove the CM class-polynomial/rounding method or replace it with certified
   interval or exact CRT construction, and prove the curve/twist selection.
   The existing small exhaustive regression does not prove this in general.
5. Establish machine-level refinement, memory safety assumptions, RNG contract,
   overflow bounds, and any timing-resistance requirements. This research code
   is not claimed constant-time.

These obligations remain open. Passing TLC, Lean, unit tests, or CI does not
discharge them. The [problem record](../../problems/sea-rust-group-refinement-proof.md)
preserves the full implementation verification requirement and its concrete
next proof step. Final tool outputs are in [formal/logs/final](formal/logs/final/).

## Replay

Run `bash formal/verify.sh` from any directory, with Java 21, curl, sha256sum,
timeout and the pinned Lean toolchain on PATH. `TLA_JAR` may point to the
already downloaded checked JAR. `FORMAL_LOG_DIR` selects a fresh output path;
committed evidence under `formal/logs` and `lean/logs` is preserved.
