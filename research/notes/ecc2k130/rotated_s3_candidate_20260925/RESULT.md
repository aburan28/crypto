# Exact affine recursive-S3 candidates on the frozen rotated PDP corpus

**Decision: the bounded semantic diagnostic passes, with one concrete
affine-chain incompleteness case.** The producer and independent verifier
enumerated all rational factor-x masks, all affine recursive-S3 root paths,
all signed factor-point tuples, and all 64 exact frozen Q+T branches. All
18 affine candidate paths on the two arms lift to signed witnesses of their
specific full target; no all-affine signed witness was missed. One additional
n13 witness reaches O at a prefix and is absent from the affine-only chain.
Thus a solver reporting no affine-chain model cannot by itself establish PDP
UNSAT. This is a solver-free toy semantic result, not an S6/S7 exporter or an
ECC2K-130 ECDLP timing.

The inputs are the SHA-pinned [#767 corpus](../rotated_pdp_corpus_20260925/RESULT.md)
and [#770 semantic gate](../rotated_m56_export_gate_20260925/RESULT.md).
The pre-outcome protocol, complete source, 600-second/512-MiB per-child caps,
and inherited archive hashes were placed in draft [PR #774](https://github.com/aburan28/crypto/pull/774)
before the first n13 result. An initial independent replay failed solely
because JSON parsed parity-counter keys as strings while its in-memory
comparison used integers. Its original freeze, producer output and failure
receipt remain in [`evidence/failure_0`](evidence/failure_0/README.md).
The corrected verifier and new freeze were committed before a fresh, separate
full run; no producer mathematics or target changed. Final `FROZEN.json`
SHA-256 is `1edd91441f0d17e882c8b5905bf15abacfef1daf2cc435aab210a13e1aa818de`.

| Toy arm | Rational x masks | Signed point tuples | Affine S3 candidate paths / distinct masks | Exact all-affine tuples / masks | Exact O-prefix tuples / exceptional-only masks | Candidate-only masks / nonrational-prefix paths / rational-prefix paths without signed target | Producer / independent verifier wall | Full-DLP S / rho |
|:--|--:|--:|--:|--:|--:|--:|--:|:--|
| n13, m5 | 243 | 3,125 | 8 / 8 | 8 / 8 | 1 / 1 | 0 / 0 / 0 | 0.212 / 0.160 s | unset / unset |
| n19, m6 | 4,096 | 117,649 | 10 / 10 | 10 / 10 | 0 / 0 | 0 / 0 / 0 | 14.024 / 9.934 s | unset / unset |

The correctness boundary was fixed before measurement: **all** exact
all-affine witnesses must be recovered, **zero** candidate-only masks or
unliftable paths may be called relations, and O-prefix witnesses must be
reported separately. The observed affine-path completeness ratio is 8/8
and 10/10; path-level and mask-level false-positive counts are both zero on
these fixed targets. The one exceptional n13 mask is `[0,0,0,2,1]` for
target index 12, `Q+O=(7256,3272)`. The complete #770 point witness is
`[(0,1),(0,1),(0,1),(6433,7897),(217,374)]`; its first two points sum to O.
The fourth and fifth points make the exact full target. It has no affine c2
coordinate, so an affine-only chain cannot express this witness. The negative
targets have zero affine paths **and** were independently certified absent
by the complete point oracle; their certification does not come from S3
UNSAT.

The explicit trace condition accepts 121/243 or 122/243 n13 masks according
to the fixed target's parity, and exactly 2048/4096 n19 masks per parity.
Across the 32 branches per arm, it accepts half the mask-target pairs.
All 18 actual affine candidate paths pass the condition, so it removed no
candidate path on this corpus. This is a rational-PDP input prefilter count,
not measured SAT construction/solve savings; all four torsion branches
together require both parities.

The n13 producer solved 2,007 S3 equations and expanded 1,764 prefix paths;
the n19 producer solved 72,736 equations and expanded 68,640 paths. The root
case counts were respectively 27/1,284/696 and 256/53,784/18,696 for
degenerate-no-root/two-root/unique-square-root; neither arm encountered a
linear root or trace-one quadratic on these rational factor domains. Every
returned root was substituted into S3. The independent verifier used a
GF(2)-pivot right inverse for `z²+z=h` instead of the producer's half trace,
and it separately rebuilt the complete signed-point oracle with bit-serial
Fermat field arithmetic. It matched each target's candidate path list,
point multiplicity, mask set, exception class, branch counts and root-case
counts. Its n3 exhaustive S3 root self-test also passed.

The complete native work ledger, including the reference oracle, is:

| Toy arm | Producer field inversions / mul / squares | Verifier field inversions / mul / squares | Each side's point additions | Peak producer / verifier RSS | Candidate stage cost relative to rho |
|:--|:--|:--|--:|:--|:--|
| n13, m5 | 15,556 / 43,549 / 90,727 | 2,130 / 114,220 / 52,990 | 12,532 | 26.5 / 27.2 MiB | unset |
| n19, m6 | 709,684 / 1,906,068 / 4,693,955 | 81,445 / 5,603,556 / 2,558,866 | 588,277 | 26.5 / 51.9 MiB | unset |

The field counters use different native scopes: the verifier's `square`
also enters its `mul` counter, while the producer's square is separate.
They are deliberately **not** combined into a common operation count or
an index-calculus S. The [raw evidence](evidence/README.md) preserves the
four full outputs, stdout/stderr, UTC child intervals, source/input hashes,
exit codes, and the first replay failure. Both full-run children per arm
finished within the fixed caps; the full sequential runner took 24.5 s on
the local host. Wall time is a practicality receipt, not an attack speedup.

This stage is classified as **accounting/semantic gate**: it tightens the
conditions a later exporter must satisfy but has no measured full-DLP gain.
The next PR should encode O/inverse prefix branches alongside affine S3,
produce direct S6/S7 and chain equations from the same fixed point targets,
and prove model-to-point equivalence before solver timing. At n131 it must
use the [literal public P/Q import](../challenge_point_import_20260925/RESULT.md),
charge export, solve, failed queries, rank and log recovery, and retain a
matched rho comparator. These toy counts are not a transfer estimate for
the public challenge point.
