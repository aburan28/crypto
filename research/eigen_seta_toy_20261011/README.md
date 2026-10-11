# Eigen-Séta toy implementation status

The source is `isogeny_algos/src/eigen_seta_toy.rs`. Run
`cargo run --release --example eigen_seta_toy` from `isogeny_algos/`.

This is a bounded correctness experiment on the public CM curve of
discriminant -7. It is **not** an implementation of the full proposal and
must not protect data. Its public sender-index search takes at most 81
isogenies, so anyone can derive the toy shared key. A separate
secret-dependent single-chain decoder now exercises the proposed lollipop
mechanism at `p=21503`, `N=128`, `d=3`. A bounded two-halves oracle at
`p=16127`, `N=64`, `d=9` does the same by exhaustively tabulating odd-torsion
images. Neither key hides the CM endomorphism.

| Requested part | Status | Evidence and gap |
| --- | --- | --- |
| Key generation | Toy complete | Computes a genuine norm-2 CM endomorphism and its 2-power eigenbasis; the CM curve and its endomorphism ring are public. Random-curve key generation is absent. |
| Encapsulation | Toy complete | Builds the 3-power Vélu isogeny, applies deterministic odd diagonal masks, and hashes a canonical ciphertext transcript. |
| Decapsulation | Secret-dependent toy complete on two fixtures | `toy_trapdoor_decaps` recovers the sender index from Kani component x-coordinates and a 3-adic DLP in both the single-chain and bounded two-halves oracle regimes. `toy_decaps` remains the public exhaustive reference path. The toy secret is publicly reconstructible. |
| Single Kani chain and component x-maps | Toy complete on the tested special-key fixture | Product-isomorphism steps, one terminal theta quotient, bounded split-chart search, and odd-torsion image recovery. Independent direct `φ θ' φ̂` evaluation agrees with the recovered `α` x-coordinates for all three sender indices. |
| Two-halves first quotients | Toy complete on the tested `p=16127,N=64,d=9` fixture | Both original and Rosati-dual N-halves reach Jacobians. Their final Igusa–Clebsch keys match; their theta nulls and swapped-source payload images agree projectively for all nine sender indices. The 4-division table for graph lifts has been replaced by two-adic digit lifting. |
| Dual of the first half | Partial | `toy_first_half_dual_prefix` starts from the Jacobian theta null, evaluates two quotient steps, and stops at an intermediate product with four steps remaining, retaining the lifted kernel images. Continuing needs a product-aware point representation; this does not yet yield the Kani matrix `F`. |
| Second Kani half and final component maps | Exponential toy oracle only; scalable evaluator not implemented | A complete `E_C[d]²` first-half image table inverts the quotient on odd torsion, and an 11,520-element theta-coordinate search identifies candidate middle automorphisms. This yields x-actions and trapdoor decapsulation for `d=9`; it is not a dual theta/Richelot chain, uses O(d⁴) state, and cannot support full-size parameters. |
| Full-size parameters and security | Not established | The toy uses 16,127- and 21,503-element prime fields and tiny sender degrees. No constant-time arithmetic, CCA transform, random-curve keygen, or security reduction exists here. |

The design note derives independent sender index `s` and masks `a,b` from a
seed, then says to recover the seed after recovering `s`. That inversion is
not specified. This toy instead uses `s` itself as the seed and derives the
masks from it. Re-encryption is consequently well-defined, but public
enumeration of `s` is trivial at these parameters.

## Single-chain decoder experiment (2026-10-11, preregistered)

Hypothesis: for the CM(-7) toy key at `p=21503`, `N=2^7`, `d=3`, the
single-chain Kani relation `N=d²·7+c²+c'²` holds with `(c,c')=(8,1)`, and
the existing Rust theta evaluator, after its explicit product steps, reaches
a split quotient. Freeze the sender indices `s=0,1,2`, keygen RNG seed
`0x44c0`, and the source revision recorded in the resulting commit. This is
a diagnostic of the special-key regime only; it does not discharge the
two-halves random-curve requirement.

Success gate: all three ciphertexts produce a split final codomain, and any
recovered component action agrees with direct `φ θ' φ̂` evaluation on the
tested odd-torsion basis. Stop the experiment if the chain returns an
inconsistent quotient or cannot produce point images; retain that failure.
Count keygen, chain, component recovery, and verification separately if timing
is reported. No asymptotic or security conclusion follows from this fixture.

Outcome: the initial Rust theta chain failed after an intermediate split,
because it assumed every nonterminal quotient stayed a Jacobian. An explicit
product-isomorphism step resolved this fixture's split, leaving six product
steps and one theta step. All three sender indices split to two factors with
the ciphertext curve's j-invariant. A bounded search over level-2 theta
coordinate changes factored every tested payload; the recovered `α` x-action
matched an independently constructed dual 3-isogeny and the direct lollipop.
`toy_trapdoor_decaps` recovered all three sender indices across four toy keys.
Tampered ciphertexts were rejected. These are correctness observations, not
security or performance measurements.

The Sage reference cited in the supplied design note is tracked on repository
branch `research/endring-pke-strategy-20261010` at
`research/eigen_seta_prototype_20261011/`. It implements special-key
single-chain decoding but explicitly leaves the random-curve two-halves
method open. The Rust implementation uses its mathematical behavior as a
reference and does not execute or extend its Python code.

## Two-halves odd-torsion oracle experiment (preregistered 2026-10-11)

Fixture: `p=16127`, `N=64`, `d=9`, sender index `s=0`, keygen seed
`0x44c0`. Hypothesis: the missing middle-Jacobian automorphism can be
identified among bounded level-2 theta-coordinate changes that fix the
first-half null. As a toy diagnostic, enumerate `E_C[9]²` through the first
half to invert it on odd torsion, then test whether a null-preserving change
composed with the dual map has the prescribed off-diagonal `[c']` action.
Acceptance requires the resulting `α` x-coordinates to agree with an
independent direct `φ θ' φ̂` computation on an odd-torsion basis. Stop if
the chart search exceeds 50,000 projective states, the Kummer image table is
inconsistent, or no candidate satisfies both component constraints; retain
the exact failure. This exponential odd-torsion enumeration is **not** an
implementation of the scalable dual theta chain and cannot establish the
full scheme even if it succeeds. No performance comparison is planned.

Outcome: at sender index `0`, the first-half Kummer table contained 3,281
projective images from 6,561 `E_C[9]²` inputs. The coordinate group search
closed after 11,520 matrices. Two candidate middle actions satisfied the
public off-diagonal component constraints: the identity gave zero for both
tested `α` images, while a nontrivial length-7 word gave x-coordinates
`(11351,6234)` and `(5971,11365)`. Independent direct evaluation of
`φ θ' φ̂ + [c]` agreed with the nontrivial candidate on both basis points.
The test command was
`cargo test --release --lib two_halves_middle_automorphism_oracle -- --ignored --nocapture`
in `isogeny_algos/` (toolchain 1.99.0). The decoder subsequently recovered
all nine sender indices for each of three toy keys; each successful recovery
was checked by deterministic re-encapsulation. All nine first-half action
tests also included a candidate matching direct lollipop evaluation. This
oracle cannot replace the missing dual-chain algorithm: its `d⁴` table would
be infeasible when `d≈2^256`, and the bounded chart search is not constant
time. No timing or security claim is made.

For the tested `p=16127`, `N=2^6`, `d=3^2` fixture, each first Kani half
needs an initial diagonal quotient followed by factorwise 2-isogenies before
the existing theta chain applies. The final quotient is a Jacobian. Its
matching dual first-half quotient is not a complete two-halves evaluation;
the bounded odd-torsion oracle fills it only at the tiny test size. Passing
either KEM round trip does not fill the scalable requirement.

## Jacobian-start dual-prefix diagnostic (2026-10-11; recorded after run)

On the existing `p=16127, N=64, d=9`, keygen seed `0x44c0`, sender index
`0` fixture, the first attempt to run the entire dual quotient from the
middle Jacobian returned no result. The failure was not at the initial
isogeny: a second-step theta null had a zero coordinate, requiring a change
of theta chart. After adding bounded chart changes, the dual quotient ran
two `(2,2)` steps and then reached an intermediate product with four steps
remaining. `toy_first_half_dual_prefix` now returns that product null and
the images of the two lifted kernel generators instead of silently treating
the product as a construction failure. The experiment was not preregistered
as a separate protocol; it is an implementation diagnostic, not evidence of
a scalable decoder.

The existing split theta coordinates retain x-coordinates of each factor
but lose relative y-signs. Continuing through the product needs explicit
point lifts or translated Kummer coordinates; selecting arbitrary signs
would not define the same kernel. The desired second half additionally
requires the middle-Jacobian automorphism and final component maps. These
are implementation gaps, not a proof of mathematical impossibility. The
first-half graph lifting itself now uses binary digit extraction rather
than an `N²`-entry 4-division table; the `d⁴` odd-torsion oracle remains.

Validation: `cargo test --release --lib` in `isogeny_algos/` passes the
Jacobian-prefix fixture and the pre-existing toy decapsulation tests. The
fixture checks the split boundary, not the dual composition `[N]` or the
intended Kani matrix `F`, neither of which is fully evaluated here.
