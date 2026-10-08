# Protocol E-3: rank from self-collisions in the root index

Frozen 2026-10-07, before any instrument is built.  **Stage diagnostic**
until priced; the registered class depends on the decision rule below.
`S` is reported on both axes, operations and memory.  Status: **PENDING**.

## Derivation (stated before measuring)

Two states of the root index with the same canonical root are a four-term
relation among factor-base points, `±π^{i}P_a ± π^{j}P_b ± π^{k}P_c ±
π^{l}P_d = O`, with no target drawn and no probe spent.  With `S = K²·n`
states and about `r/(2n)` canonical root classes, the expected number of
unordered colliding pairs is

    S² / (2 · r/(2n))  =  K⁴ n³ / r,

less the degenerate ones (the same pair twice, or a pair and its mate).
The rank stage needs about `K` independent rows, so collisions alone pay
for the rank when `K ≳ (r/n³)^{1/3}`.  At the ladder sizes:

| curve | log₂ r | `(r/n³)^{1/3}` | states at that `K` |
|:--|--:|--:|--:|
| `icv1-f2m41-tm2308219-7f48b14a` | 39.0 | about 200 | 1.6·10⁶ |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | about 520 | 1.4·10⁷ |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | about 890 | 4.8·10⁷ |
| `icv1-f2m83-t6151469093347-cdcc5432` | 52.9 | about 2,500 | 5.2·10⁸ |

The §8a gate curve `icv1-f2m83-tm6151469093347-debefd74` has `r` near
`2^{81}`, where the required `K` is about `1.6·10⁶` and the index `2^{44}`
states; it is out of reach for this protocol and is listed to say so.

What the trade is: the guided rank at `K = 600` on the `a = 1` curve spent
about `1.5·10⁹` probes; an index of `5.2·10⁸` states costs about that many
`S₃` solves to build, needs no rank probes at all, and at the larger `K`
cuts online probes per relation by `(K'/K)²`, about 17×.  Against that,
memory grows by `(K'/K)²` too, about 17× before E-1's packing, and the
linear algebra has 2,500 columns instead of 600.  The zero-collision probe
on branch `cursor/ic-boundary-experiments-d111`
(`examples/koblitz_zero_collision_probe.rs`, unmerged) already counts
collisions and their rank; it has no cost accounting and no decision rule.

None of this changes the product law: the self-collision rank is a point
inside the pair-table family, priced by memory.  It can only move the
frontier along the memory axis.

## Instrument (Rust, to build in the follow-on PR)

1. Extend the compact-orbit producer with `KIC_RANK=collisions`: build the
   index, sort or hash canonical roots, emit every non-degenerate
   collision as a relation row, run the incremental solver, and stop at
   full rank or at index exhaustion.
2. Run `K ∈ {0.8, 1.0, 1.2, 1.5} × (r/n³)^{1/3}` on the four reachable
   curves above, three seeds each.  Record collisions found, the
   degenerate fraction, the rank reached, and whole-process instructions
   and resident bytes under Callgrind as in the strong-rho ladder.
3. Run the guided rank at the ledger's `K` on the same curves as the
   control, same unit.
4. Where collisions reach full rank, solve the frozen single targets at
   the new `K` and record online probes.

## Predictions (pass/fail)

- **C1 (law).**  Measured non-degenerate collisions over `K⁴n³/r` lie in
  `[0.5, 2.0]` at every `K` and curve.
- **C2 (degeneracy).**  At least 90% of collisions are non-degenerate.
- **C3 (rank).**  Full rank is reached at `K ≤ 1.2 × (r/n³)^{1/3}` on
  every curve, and never at `0.8×`.
- **C4 (online).**  Online probes to the frozen targets fall by at least
  `0.7 × (K'/K)²` against the ledger's `K`.
- **C5 (verification).**  Every column log and relation verifies in the
  group; every frozen target's scalar is reproduced.

## Decision rule (registered)

If C1, C2, C3 and C5 pass, price the row: whole-process instructions and
resident bytes against the guided rank at the ledger's `K` and against the
strong rho.  Adopt as a producer option if instructions fall by at least
1.4× at no more than 4× the memory after E-1's packing, classed
**engineering**; otherwise record the measured point on the frontier's
memory axis, classed **boundary**.  C1 failing by more than 2× in either
direction halts the protocol for review of the class count `r/(2n)`.

## Stop condition and inadmissible moves

Bounded: four curves, four `K` values, three seeds, one control.

Inadmissible: counting degenerate collisions as relations; reporting
probes saved without the memory spent; comparing against the guided rank
at a different `K` without saying so; using the `a = 0` gate curve's
numbers from the `a = 1` curve; reading a memory-bought constant as a
change to the ratio-to-rho exponent.
