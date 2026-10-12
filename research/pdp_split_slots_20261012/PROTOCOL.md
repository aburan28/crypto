# Protocol: orbit-union τ-slots for the split and hybrid PDP, measured in ecbench

**DRAFT, not yet registered: the grid is fixed after the sizing runs. No measured run has started.** Only the smoke
runs listed at the end preceded it. Conductor task `T-2`. This round
follows `research/pdp_split_dreg_20261011` (round 1, PR #1709).

## Question

Round 1 used the standard subspace `V = span(1, z, …, z^{ℓ−1})`. This
round asks whether the Frobenius closure `∪ⱼ σʲ(V)` lowers the cost of
the split and hybrid point-decomposition systems, measured end to end in
ecbench rather than projected. The closure is folded one column per
signed Frobenius orbit, and each query solves one system per τ-slot tuple.

### What the closure can and cannot change (stated before running)

- **Columns.** The closure folds about `n·|V|` points into about the same
  column count as `V`. So the number of relations needed does not shrink
  by `n`.
- **Systems per relation.** A τ-slot system with summand `i` over
  `σ^{jᵢ}(V)` meets more decompositions per solve than the plain `V^m`
  system, up to the ordering factor `m!` (= 6 at `m = 3`). Round 1's lead
  L5 measured 1.34× and 1.79× at `m = 2`.
- **Cost per system.** Each system has the same unknowns and degree as on
  `V`, so its cost should be unchanged.

## Code

- `pdp_split.rs`: per-summand bases (`build_split_system_bases`,
  `run_hybrid_bases`, `exhaustive_decompositions_bases`), `slot_tuples`,
  `closure_census`, with tests.
- `src/bin/ecbench/pdp_split_oracle.rs`: the oracle `pdp-split:m,plan,guess,slots,max_free`.
  It is registered beside `pdp3-koblitz`, together with the existing
  `koblitz-closure` base.
  - Charged: the sign lift's group additions.
  - Counted, not charged: F4's word XORs, under the `f4-f2` unit "word XORs
    (elimination only)". So `S` is a lower bound, as for every algebraic
    oracle in ecbench.
  - Test: `tests/ecbench.rs::pdp_split_verifies_on_standard_and_closure_bases`.

## Part A — ecbench sessions (primary)

**Curves:** the round-1 Koblitz curves.

| curve | dimension `ℓ` |
|---|---|
| `icv1-f2m11-t67-05f5aa36` (`K_1`, n = 11) | 4 |
| `icv1-f2m13-t181-515ee569` (`K_0`, n = 13) | 4 |
| `icv1-f2m17-tm101-00378d4e` (`K_1`, n = 17) | 5 |
| `icv1-f2m19-tm797-9c54981b` (`K_1`, n = 19) | 4 |

Round 1 found no decomposition in 256 targets at n = 17, `ℓ = 4`. That is
why n = 17 uses `ℓ = 5`.

**Arms per curve,** all `m = 3`, guess 4, `max_free` 24, `max_trials` 200000:

| arm | method |
|---|---|
| `rho-strong` (reference) | `rho.signed_frobenius_strong` (lanes 32, dp_bits 4, step_cap_factor 2000) |
| `rho-frob` | `rho.signed_frobenius` |
| `std-unsplit-g4` | `koblitz-standard-subspace:dimension=ℓ` × `pdp-split:plan=3,slots=one` |
| `clo-unsplit-g4` | `koblitz-closure:dimension=ℓ` × `pdp-split:plan=3,slots=all` |
| `std-chain-g4` | `koblitz-standard-subspace:dimension=ℓ` × `pdp-split:plan=2+1,slots=one` |
| `clo-chain-g4` | `koblitz-closure:dimension=ℓ` × `pdp-split:plan=2+1,slots=all` |

The chain arms run only where the smoke timings below fit the budget;
`SPECS.md` lists the arms frozen per curve.

**Sizes and isolation.** Eight public one-target workloads per curve, one
round, no warm-up. The figures are operation counts and replay
identically, so no wall time is reported. Isolation level L0 is
accepted.

**Seeds.** Target seed `202610121101`; measurement seed `202610121102`.

**Audit.** `ecbench verify --replay-all` for every session.

**Figures.**

- **Recorded S.** A lower bound.
- **Priced S (sensitivity).** `S + F4 word XORs · ns_per_word_xor / √r`
  per run, at `docs/ic/calibration.json`'s pinned ratio for the curve,
  which was calibrated on another host.
- **Ratios on the same workload:**
  - priced S ÷ strong-rho S, and ÷ signed-Frobenius-rho S;
  - closure ÷ standard twin.

  Medians over the 8 workloads, with min–max.
- **Analysis:** native `pdp_split_analyze --ecbench …`.

## Part B — stage diagnostic

**Grid.** Round 1's `m = 3` cells: `3:1:11:4`, `3:0:13:4`, `3:1:17:4`,
`3:1:19:4`, `3:0:13:5`, `3:1:17:5`. Each cell is run with
`pdp_split_dreg --base closure`:

- `--slots-per-target 4 --targets 4`;
- `--guesses 0,2,4`, plans `3` and `2+1`;
- seed `20261012`, `--budget-secs 300 --max-free 24`;
- `--lambda-systems 512` (yield per system, exhaustively counted over
  random (target, slot-tuple) pairs).

**Checks.** Every run is checked against exhaustive enumeration for its
slot tuple.

**Projection.** `(columns + 1) / λ_sys · X̄ · ratio / √r`, compared with
round 1's standard-base projection on the same cell. This part gives the
per-system view (solving degree, yield, cost) beside Part A's end-to-end
measurement. Its projection stays labelled a projection.

`m = 4` is not run on the closure: about 2.7 h per cell at four slot
tuples.

## Predictions (written first)

- **R1.** On every curve, `clo-unsplit-g4` solves fewer systems per run
  than `std-unsplit-g4`. Median over workloads, at least 1.5× fewer.
- **R2.** On every curve, the priced S of `clo-unsplit-g4` is below that of
  `std-unsplit-g4`: the median per-workload ratio is below 1.
- **R3.** The yield per system λ_sys on the closure is at least 1.5× the
  standard per-target λ in every Part B cell where both are nonzero.
- **R4.** Solving degree per system is the same on the closure as on the
  standard base, for the same plan and g, in every Part B cell (up to ±1).
- **R5.** No IC arm's median priced S is below the strong-rho median S on
  its curve. **Falsification target:** any curve where one is.
- **R6.** Every ecbench run verifies, and every Part B run is complete and
  exact.

## Inadmissible

The following are inadmissible after any measured run:

- changing curves, dimensions, arms, seeds, guess, `max_free` or pricing;
- dropping failed, exhausted or regressing runs;
- quoting priced S as recorded S, or either as a speedup claim;
- any claim about deployed curves, asymptotics or ECC2K-130.

## Smoke runs (sizing only; scratch, not committed)

- `K_1`/n = 11, ℓ = 4, 2 targets: all 8 runs verified after fixing the
  oracle's closure check, which had wrongly required `x = 0` in the base.
  Replay OK.
  - unsplit: 58–59 systems per run on the standard base, 25–29 on the
    closure.
- n = 13, ℓ = 4, 1 target, systems per run (standard → closure):
  - unsplit: 119 → 41, in 2.1 s → 0.6 s;
  - chain: 119 → 41, in 21.6 s → 7.7 s.
- n = 17 (ℓ = 5) and n = 19 (ℓ = 4): see `SPECS.md`.
