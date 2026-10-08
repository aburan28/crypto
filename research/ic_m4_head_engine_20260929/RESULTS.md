# The m = 4 exponent audit on the head engine: results

This was written after the run. [PREREGISTRATION.md](PREREGISTRATION.md) is unchanged, and
the raw readout is [runs/readout.txt](runs/readout.txt).

## Run

| | |
|:--|:--|
| engine | `m4_exponent_audit-4ff512f2`, sha256 `307bcd96…6545` (`runs/binary.sha256`) |
| reproduction control | exact, before registration (`baseline/REPRODUCTION-4ff512f2.txt`) |
| wall | 2026-09-29 16:06:34Z to 16:28:01Z |
| isolation | the whole run under the benchmark lock; parallel cells pinned one per CPU (1, 2, 3), degree cells on CPU 3 |
| targets | identical to the first audit's: same seeds, same `k` on every cell (checked) |

## Registered verdict

**Closed for this engine at these sizes.** `ĉ = 1.060` bits per unit `n`, bootstrap band
[1.058, 1.093], with `B = 10,000` and no invalid replicates. No Semaev target censored.

## All four predictions hold

1. **The verdict is closed.** This was predicted.
2. **The head slope is not lower.** It is 1.060, above the first audit's 0.985 and inside
   the predicted [0.90, 1.40]. The head engine is cheaper at every cell, by a factor that
   shrinks with `n`, so the slope is steeper:

   | cell | frozen `2809b498`, lower median | head `4ff512f2` | frozen / head |
   |:--|--:|--:|--:|
   | `K0n9l2` | 68,235 | 30,805 | 2.22 |
   | `K1n9l2` | 39,941 | 12,624 | 3.16 |
   | `K1n11l3` | 719,041 | 454,647 | 1.58 |
   | `K0n13l3` | 1,493,937 | 985,448 | 1.52 |
   | `K0n15l4` | 5,259,157 | 3,345,342 | 1.57 |
   | `K1n15l4` | 5,079,215 | 3,299,274 | 1.54 |
   | `K1n17l4` | 17,475,739 | 12,787,163 | 1.37 |
   | `K0n19l5` | 50,076,470 | 33,575,448 | 1.49 |
   | `K1n19l5` | 80,892,472 | 54,785,935 | 1.48 |

3. **The controls pass.**
   - enumeration null 0.828, inside the registered [0.60, 1.00];
   - random-system null censored at every cell, a pass;
   - 0 oracle disagreements;
   - 0 group-check failures.
4. **Cheaper at `n = 9` on both curves.** The head engine is cheaper there by 2.22× and
   3.16×.

## Reported, not decisive

- **Refuted-only fit.** `ĉ = 1.060`, band [1.056, 1.062].
- **Prime `n` only.** `ĉ = 0.845`, band [0.843, 0.965]. The reading is closed, as in the
  primary fit. The same prime-versus-all gap appeared in the first audit (0.833 against
  0.985), on the same targets.
- **Degree arm.** It is identical in outcome to the first audit's:
  - `D = 6` on the first unsatisfiable target at `ℓ = 2`, on both curves;
  - both `ℓ = 3` cells were killed by the 300 s per-cell CPU limit before finishing a
    target. That is machine protection, so the result is censored, not negative;
  - 0 count faults.

## What this says

The post-freeze engine changes are constant-factor gains, about 1.4–1.6× from `n = 11` up,
with more at `n = 9`. They do not move the exponent down; on these nine cells they move it
slightly up.

`m = 4` stays closed for both engines at `n ≤ 19`. Per the first audit's §10, that is no
asymptotic statement, and it does not close the algebraic route. The route's next test is
the symmetry-lever slope test (`research/ic_symmetry_lever_slope_20260929/`). It decides
whether an `m = 5` audit with that lever is worth building.
