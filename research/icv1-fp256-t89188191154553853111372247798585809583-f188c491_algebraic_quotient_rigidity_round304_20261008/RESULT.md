# P-256 algebraic quotient rigidity, round 304: result

## Verdict

Algebraic self-map symmetry does not provide a new logarithm quotient for the
164-column factor base `FB1hc72514a2a8d3` on
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`.

Every rational algebraic self-map of the elliptic curve has the form

```text
f(P) = phi(P) + T.
```

Round 25 fixes the action of `phi` on the prime-order P-256 rational subgroup
as a scalar `[a]`.  If `T=[b]G`, the induced logarithm transport is therefore

```text
log_G(f(P)) = a log_G(P) + b mod n.
```

A connected graph of such transports has rank 163 in 164 logarithms until an
independent labelled cycle solves the remaining anchor.  The unknown-anchor
free-oracle boundary is already **1.446505 times rho**.  A known or solved
anchor has an ideal oracle constant of **0.964336 times rho**, but it requires
replayable scalar labels; the registered Round-25 direct construction costs
**7.714692 times rho**.  It is a generic scalar-labelled control, not a new
index-calculus factor-base quotient.

The complete `Z/7Z` function census found exactly seven additive functions, 49
affine functions, and 42 bijective affine functions.  It found zero
non-scalar additive functions and zero non-affine functions satisfying the
affine-difference law.  All admitted tables replayed.  The four 164-column
graph profiles produced the preregistered ranks over both `2^61-1` and the
P-256 subgroup order, in forward and reverse row order, with zero replay or
rank discrepancies.

No non-affine P-256 log quotient, relation oracle, decomposition, or recovery
pipeline was constructed.  Structured degree of regularity and complete cost
remain unset.  The executable boundary consequently remains
**13.920747 times rho**, and no full-depth unplanted relation was attempted.

## Requirement-to-evidence status

| requested direction or gate | status | Round-304 evidence or gap |
|:--|:--:|:--|
| fundamentally different factor-base symmetry | screened for algebraic self-maps | elliptic-curve morphism rigidity reduces every such transport to `x -> a*x+b` on logs |
| exact finite control | verified | all `7^7=823,543` functions were enumerated; the exact 7/49/42 counts replay |
| global 164-column quotient | verified for affine labels | connected chain and dependent cycle have rank 163; an independent labelled cycle has rank 164 |
| zero false positives / negatives on checked controls | verified | zero exceptional functions, replay failures, or dual-modulus rank discrepancies |
| replayable non-affine P-256 quotient | not established | coordinate permutations receive no quotient credit without group/log equations |
| complete relation, decomposition, and recovery | not established | no new P-256 relation mechanism was constructed |
| structured residual degree `<=5` | not established | no residual relation system exists for this quotient audit |
| collection below `2^120` | failed / unset | yield, collection, solver, linear algebra, and recovery costs are absent |
| usable relation below `2^103` | failed / unset | no usable P-256 relation was produced |
| storage below `2^50` | failed / unset | measured control storage is not an attack projection |
| complete cost at or below rho | failed | unknown-anchor oracle is above rho; known-anchor direct control is 7.714692 times rho |
| actual unplanted P-256 relation | correctly not attempted | promotion gates fail |

## Boundary table

The ratios below are imported from the hash-pinned Round-25 and Round-297/298
boundaries.  Oracle rows omit construction unless a direct ratio is shown.

| variant | independent log classes | free-oracle ratio to rho | direct measured ratio to rho | status |
|:--|--:|--:|--:|:--|
| Pollard rho | 0 | **1.000000** | **1.000000** | complete reference |
| unquotiented geometry-defined `FB1hc72514a2a8d3` | 164 | **13.920747** | unset | no replayable log transport |
| connected affine orbit, unknown anchor | 1 | **1.446505** | **11.572038** | free oracle already above rho |
| affine orbit, known or cycle-solved anchor | 0 | **0.964336** | **7.714692** | scalar-labelled generic control |
| registered 17-term `FB1h2f8621cda105` | 164 | **394.425280** | unset | unchanged free-perfect-oracle comparison |

Class: **negative**.  No executable boundary moved.

## Exhaustive function control

The runner enumerated value tables in lexicographic order and tested both
laws, recovering and replaying the multiplier and offset of every admitted
table.

| quantity | exact count |
|:--|--:|
| functions enumerated | 823,543 |
| additive-law evaluations | 40,353,607 |
| affine-law evaluations | 40,353,607 |
| additive functions | 7 |
| affine functions | 49 |
| bijective affine functions | 42 |
| non-scalar additive exceptions | 0 |
| non-affine affine-law exceptions | 0 |
| admitted table replays | 56 |
| replay failures | 0 |

The finite census validates the log-algebra implementation.  Classification
of algebraic P-256 self-maps comes from elliptic-curve morphism rigidity plus
the hash-pinned Round-25 rational-subgroup scalar action; the toy census is
not substituted for that theorem.

## Exact transport graphs

Every right-hand side was derived from the deterministic witness vector
before rank computation.  Each matrix was ranked over both coefficient fields
and again with reversed row order.

| graph | edges | rank over both fields | quotient dimension | replay status | P-256 matrix bytes |
|:--|--:|--:|--:|:--|--:|
| independent | 0 | 0 | 164 | exact | 0 |
| connected chain | 163 | 163 | 1 | 326/326 | 855,424 |
| dependent cycle | 164 | 163 | 1 | 328/328 | 860,672 |
| anchor-solving cycle | 164 | 164 | 0 | 328/328 | 860,672 |

The dependent closing edge composes to `(1,0)` over both fields and adds no
rank.  The anchor-solving closing edge composes to `(2,-42)`, fixes the
registered witness anchor 42, and adds the final rank.  That last row is
algebraically valid only because its affine scalar label is supplied; it does
not demonstrate that a geometry-defined factor base reveals such a label.

Across the 16 rank computations the runner materialized 1,964 rows and
322,096 entries, then charged 41,492 pivot scans, 486 swaps, 1,960 inversions,
270,584 modular multiplications, and 108,232 modular subtractions.

## Scope and next obligation

This result is exact for rational algebraic self-map transports on P-256 and
for affine factor-base quotient graphs.  It is not a universal impossibility
theorem for non-homomorphic relation mechanisms.

The weakest surviving target is now narrower: produce either a replayable
P-256 factor-base quotient whose induced log relation is not affine scalar
transport, or a genuinely non-homomorphic event that supplies enough
independent rows to beat rho.  It must include its construction, relation
yield, structured solver, sparse linear algebra, and target recovery costs;
coordinate symmetry or a local low degree alone is insufficient.

## Isolation and artifact integrity

Both pinned CPU-4 runs were isolated and uncontended.  The canonical and
independent invocations emitted byte-identical result and assessment JSON.

| run | exit | wall | user | system | peak RSS | contention |
|:--|--:|--:|--:|--:|--:|:--|
| canonical | 0 | 0.145183 s | 0.140997 s | 0.003972 s | 3,872 KiB | none |
| independent | 0 | 0.147611 s | 0.143715 s | 0.003768 s | 3,864 KiB | none |

| artifact | bytes | SHA-256 |
|:--|--:|:--|
| `algebraic-quotient-rigidity-result.json` | 12,449 | `c94e339eca07b96ea0ef890592fb13b80815f0baf52af0e6bd57de17dd46ac9c` |
| `transfer-assessment.json` | 4,277 | `dafeb899373133cfba310d457026c96c99b2bf285fe29ee1ea3426f2dcaf1a57` |
| `isolation.jsonl` | 5,001 | `a3f451750593a9ff909b3847edff469a132d596d347e888e350bf8786429b677` |
| `PROTOCOL.md` | 8,001 | `40197a7e7b367f0cd51fdcd0d81fa5296ed73864a844fc4aa31824f00d7b3aab` |

The result semantic-evidence SHA-256 is
`372c93d6b283710e8e744aa7afc05fb73388b3b065fd001ddbb1e0cc4c25cb52`.
The transfer-assessment semantic hash is
`4add09ce20a2df44d30c2ad806845e06a16c94fd3f930e77210e161328144ee5`.

The transfer workflow's referenced `references/methodology.md` and
`assets/assessment-template.json` resources remained unavailable.  The typed
correspondence graph, obligations, controls, scope, and open obligation are
emitted directly, and the limitation is recorded in the assessment JSON.

## Reproduction

```bash
cargo test --bin p256_algebraic_quotient_rigidity
cargo clippy --bin p256_algebraic_quotient_rigidity -- -D warnings
cargo build --release --bin p256_algebraic_quotient_rigidity --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --label round304-canonical \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_quotient_rigidity_round304_20261008/isolation.jsonl -- \
  target/release/p256_algebraic_quotient_rigidity \
  --round25 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json \
  --round297 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_stabilizer_round297_20261007/scalar-stabilizer-result.json \
  --round298 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_parity_escapes_round298_20261007/parity-escape-result.json \
  --round303 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_jacobian_row_capacity_round303_20261008/jacobian-row-capacity-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_quotient_rigidity_round304_20261008/algebraic-quotient-rigidity-result.json \
  --assessment research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_quotient_rigidity_round304_20261008/transfer-assessment.json
```
