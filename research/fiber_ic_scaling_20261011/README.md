# Does fiber-aware index calculus ever cross rho? No: the gap grows as `#E^{1/6}`

**Date.** 2026-10-11.  **Status.**
- Measured: three audited `ecbench` sessions, 816 records, all verified,
  every audit's 12 replays reproduced.
- The protocol ([`PROTOCOL.md`](PROTOCOL.md)) was frozen before the first
  session.  Amendment 1 was added after reading the second; it is labelled
  wherever it is used.

**Result class.**
- Whole-method cold cost against matched rho, one target per run
  (AGENTS.md §2, §8).  This is a negative result.
- Synthetic public curves.  No deployed-curve claim, no key recovery.

PDF: [`report.pdf`](report.pdf).

## The question

The three conductor-gap rounds found levers that cut index-calculus cost on
cofactor curves:
- the fiber-closed pair table
  ([round 1](../isogeny_conductor_gap_ic_20261010/README.md));
- partial fibers and a corrected closure
  ([round 2](../isogeny_conductor_gap_round2_20261011/README.md));
- the 2-torsion quotient
  ([round 3](../isogeny_conductor_gap_round3_20261011/README.md)).

Round 2's best cell was 12.8× rho.  Asked whether this speeds up anything
end to end, this study measures the best fiber-aware IC against rho across
field sizes at a fixed cofactor.  The question is whether the gap closes.

## Answer

![Envelope IC ÷ rho against size](figures/envelope_vs_size.svg)

For each curve, the envelope is the cheapest fiber-closed MITM arm over
factor-base sizes `B` = 4 to 256 (`prime-abscissa:size=B`,
`mitm-fiber:fiber=closed,negation_folded=1,fiber_gen=coset`), divided by
`rho.negation` on the same curve, both as mean `S`.

| cofactor | curves | envelope ÷ rho at 18 → 22 → 26 → 30 bits (per curve) | slope of log₂ ratio on log₂ #E | 95 % interval |
|--:|--:|:--|--:|:--|
| 4 | 8 | 8.4, 8.8 → 11.0, 11.3 → 20.6, 17.9 → 35.1, 37.0 | **0.177** | 0.145 to 0.209 |
| 72 (with amendment 1) | 8 | 20.1, 13.4 → 17.7, 21.2 → 42.5, 37.6 → 60.4, 63.1 | **0.168** | 0.115 to 0.221 |

Both slopes are the `1/6` the protocol derived.  IC with a fiber-closed
table costs about `N^{2/3} h^{−1/3}` per target, and rho costs
`(N/h)^{1/2}`.  So the ratio grows as `(N·h)^{1/6}`, about 1.12× per bit of
`#E`.
- From 18 to 30 bits it went from 8.4× to 37× at `h = 4`, and from 13× to
  63× at `h = 72`.
- The fitted lines reach a ratio of 1 only at `log₂ #E` ≈ 0.8 and −5.5,
  which are not curves.
- Extrapolated upward, the `h = 4` fit gives about 2^45× at a 256-bit
  group.  That figure is derived and was not measured.

**Fiber-aware relation generation does not speed up the end-to-end
discrete-log solve on these curves at any size, and it falls further behind
rho as the curves grow.**  The levers of rounds 1 to 3 lower the IC
constant.  They cannot change the exponent, because the pair-table method's
`N^{2/3}` is worse than rho's `N^{1/2}`.

## Against the predictions

- **P-S1, slope ≈ 1/6 and > 0: holds** for both cofactors.  It holds for
  `h = 72` only with amendment 1, as follows.
  - On six of the eight cofactor-72 curves, the pre-registered sweep
    (`B ≥ 16`) put the envelope at its floor, `B = 16`.  The derivation's
    optimum there is `B` ≈ 4 to 10, so those curves were over-priced.
  - The pre-registered `h = 72` fit alone gives slope −0.073 (interval
    −0.183 to 0.038).  That is an artifact of the unresolved envelope, not
    a closing gap.
  - Amendment 1 reran the same curves and targets at `B ∈ {4, 8}`.  Over
    the combined sweep the envelope is resolved on every curve: `B = 4` or
    8 below 26 bits, 16 or 8 at 26, and 32 at 30.
- **P-S2, cofactor 72 above cofactor 4 by about `18^{1/6}` ≈ 1.6: holds
  roughly.**  The measured ratio at matched size is 1.94, 1.74, 2.09 and
  1.71 at 18, 22, 26 and 30 bits.
- **P-S3, the fiber beats the blind table by at least `h^{1/3}`: fails as
  stated.**
  - The blind table at `B = 128` over the closed envelope is:
    - at `h = 4`: 4.8 to 5.0 at 18 bits, 1.7 to 1.9 at 22, 1.4 at 26,
      2.3 to 2.8 at 30;
    - at `h = 72`: 5.8 to 6.9, 2.3 to 2.5, 1.2 to 1.4, 2.2 to 3.0.
  - The blind arm was swept at one size only, so these over-state the
    fiber's gain at most sizes.  Even so, they fall below 1.6 (`h = 4`) and
    4.2 (`h = 72`) from 22 or 26 bits on.
  - The fiber's constant-factor gain over the blind table is smaller than
    the crude model's.  It does not change the exponent either way.

The `h = 4` envelope at 30 bits sits at the top of the sweep (`B = 256`).
Its neighbour `B = 128` is within 2 to 8 percent, so the cost curve is flat
there and the slope is barely affected.

## Sessions

| session | spec | records | audit | isolation earned | run report |
|:--|:--|--:|:--|:--|:--|
| `scaling_h4` (`ECBS1h4a5009a8f716`) | [`specs/scaling_h4.json`](specs/scaling_h4.json) | 336 | ok, 12 replays reproduced | L2 × 231, L1 × 105 | [`run-20261011T2037Z`](reports/run-20261011T2037Z/README.md) |
| `scaling_h72` (`ECBS1hbf4321e2bec8`) | [`specs/scaling_h72.json`](specs/scaling_h72.json) | 336 | ok, 12 replays reproduced | L2 × 197, L1 × 139 | [`run-20261011T2039Z`](reports/run-20261011T2039Z/README.md) |
| `scaling_h72_small_b` (`ECBS1h54835df66b11`), amendment 1 | [`specs/scaling_h72_small_b.json`](specs/scaling_h72_small_b.json) | 144 | ok, 12 replays reproduced | L2 × 88, L1 × 56 | [`run-20261011T2044Z`](reports/run-20261011T2044Z/README.md) |

**Curves.**
- Random prime-field curves `y² = x³ + ax + b` with `#E = h·r`, `r` prime:
  two per size at 18, 22, 26 and 30 bits, for `h = 4` and 72.
- Found by [`examples/fiber_scaling_curves.rs`](../../examples/fiber_scaling_curves.rs)
  (order by baby-step giant-step over the Hasse interval, confirmed on a
  second point).  `ecbench` re-checked each on load.
- Slugs are in `TABLES.md`.

**Analysis.**  [`examples/fiber_scaling_tables.rs`](../../examples/fiber_scaling_tables.rs)
writes [`TABLES.md`](TABLES.md) and [`tables.json`](tables.json):
- per-curve means;
- the envelope;
- least-squares slopes with Student-t intervals;
- the cofactor comparison.

`S` is `ecbench`'s counted cost over the floor.  It charges every group
operation, including oracle set-up and the pair table.  Both sides are
lower bounds:
- the IC arms leave table lookups, square roots, row operations, Legendre
  symbols and target-guard probes unpriced;
- rho leaves its canonicalisations unpriced.

Pricing the IC side's unpriced work would raise the ratios, and pricing
rho's canonicalisations would lower them.  No wall
times are reported.

**Graphs.**
- [`figures/envelope_vs_size.svg`](figures/envelope_vs_size.svg), from
  `figures/plot.py`, which plots `tables.json`.
- Each run report's `s_by_arm.svg`: every arm's `S` against `log₂ #E`.

## What this leaves open

- **The amortised regime.**  Fiber-aware IC is not measured here when a
  table is reused across many targets.  There the set-up `π·s` leaves the
  per-target cost, and the relation phase `∝ N/(B·s)` can fall below rho's
  `√(N/h)` once memory allows `B·s > √(N·h)`.
  - That is the multi-target regime the scoreboard keeps separate from
    one-target results.  This study does not admit or measure it.
  - The derivation says the fiber multiplies the online gain by `s` for the
    same table memory.
- **Other decompositions.**  Only the `m = 2` pair table was tested.
  Algebraic decompositions (Semaev and Weil descent) were shown in rounds 1
  and 3 to be dominated by solver work that `S` leaves unpriced.

## Reproduce

```bash
cargo build --release --bin ecbench --example fiber_scaling_curves --example fiber_scaling_tables
./target/release/examples/fiber_scaling_curves --cofactors 4,72 --bits 18,22,26,30 --per-size 2 \
  --seed 20261011 --out research/fiber_ic_scaling_20261011
research/fiber_ic_scaling_20261011/run_sessions.sh
./target/release/examples/fiber_scaling_tables --dir research/fiber_ic_scaling_20261011 \
  --out research/fiber_ic_scaling_20261011/TABLES.md --json research/fiber_ic_scaling_20261011/tables.json \
  --combine scaling_h72_combined=scaling_h72+scaling_h72_small_b
```
