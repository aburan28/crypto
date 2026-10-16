# Frozen n41/n53 shared-rank K8/K16 counted-unit panel

Status: preregistered. This protocol, `SPEC.json`, `PLAN.json` and
`FREEZE.json` are committed and pushed before the first measured execution.
No seed, cap, arm, workload or decision rule may change after an outcome is
seen; an amendment is a new commit that keeps every original row.

## Question

The n37 thread ended with
[`carry_both_host_noise_exceeds_gate`](../ecbench_n37_native_online_wall_20261004/RESULT.md):
the K16 one-target online lead cannot be promoted without an L2 host, and any
n41/n53 transfer "still requires separate matched full-rank workloads and
charged cold costs". This round carries the **same four arms with identical
method parameters** to the next two registered sizes of the same family and
asks, in the harness's counted unit, whether the n37 figures survive:

- `icv1-f2m41-tm2308219-7f48b14a` (Koblitz `a=0,n=41`, subgroup order
  `r = 2^39.00`, 82 automorphisms, EC1 `EC1N41Ce0he09550ab560a`);
- `icv1-f2m53-tm56619371-dac20a85` (Koblitz `a=0,n=53`, `r = 2^44.26`, 106
  automorphisms, EC1 `EC1N53Ce0hb097de99be9a`).

Both slugs are in `docs/curves/registry.json`. The field degrees 41 and 53
are prime, so neither field has a proper intermediate subfield over GF(2)
(AGENTS.md §8b); `ord_41(2) = 20` and `ord_53(2) = 52`, so the nontrivial
cyclotomic block splits at 41 and is irreducible at 53, as at 83 and 131.

## Frozen inputs

`SPEC.json` (`ecbench.spec/v1`, spec id `ECS1hb5e99ec6937e`) fixes:

- **Workloads.** 16 public one-target workloads per curve from
  `hash_to_subgroup_v1` at target seed `202610055101`; 32 workloads, each a
  separate one-target run. No solver receives a logarithm; `[k]G = Q` is the
  only check.
- **Arms**, byte for byte the n37 native-wall spec's arms:
  `rho-strong` = `rho.signed_frobenius_strong` (lanes 32, dp_bits 8,
  step_cap_factor 2000), the reference; `ic-k8` and `ic-k16` =
  `ic.shared_rank` with `compact-orbit-scan:columns=8,raw_x_cap=1000000`
  and `columns=16`, rank seed `202610042032`, `rank_max_trials` 100000,
  `target_max_attempts` 512; `ic-k16-control`, the identical K16 method as
  the A/A control. The method ids are `ECM1h1d9961ee4601`,
  `ECM1h90d7515e0e7a` and `ECM1h670d3d9b99a0`.
- **Measurement.** Five measured rounds and one warmup, arm order
  alternating at seed `202610055102`, `isolation_required: L2`, 600 s per
  child. 768 executions, 640 of them measured.
- **Statistics.** Ratios are ratios of sums over matched (workload, round)
  pairs; intervals are 20,000-resample two-stage (workload-block) percentile
  bootstraps at seed `202610055103`, produced by
  `ecbench compare --resamples 20000 --seed 202610055103 --save`.
- **Primary row per curve.** Target index 0 of each curve, designated
  before any run: `W15dcf1a74cb8` (n41) and `W32013b060ceb` (n53). The other
  15 workloads per curve are secondary same-size replications, reported
  individually and in the ratio of sums; never a batch.

`PLAN.json` is the harness's expansion of the spec. `FREEZE.json` records
the file hashes, the primary workloads and the overlap checks below.

### Target overlap

Every planned target point and workload id was checked against:

1. every committed `ecbench` session whose plan names an `n=41` or `n=53`
   Koblitz curve (`research/ecbench_calibration_20261002/sessions/koblitz`,
   `.../followup-koblitz`, `research/ecbench_pair_claw_20261003/sessions/koblitz`;
   88 prior workloads, all planted-scalar targets): 0 point and 0 id overlaps;
2. the four frozen public points of PR #1353
   (`experiments/koblitz-base-size-cold-panel-20261004/TARGETS.md`, not on
   `main`): 0 overlaps;
3. a text search of `research/`, `experiments/`, `docs/`, `examples/`,
   `src/` and `tools/` on `main` for each abscissa as `0x`-hex, bare hex and
   decimal: 0 files outside this directory. This covers the ledger's own
   n41/n53 one-target rows (`research/ic_single_target_20260930/runs/k0n41`,
   `.../k0n53`), the exponent runs and the SAT autolab n53 panels, whose
   points are stored as decimal or hex coordinates.

`experiments/koblitz-n41-n53-single-target-20261004` does not exist on
`main` or on any open branch (GitHub code search returned no match), so no
check against it was possible.

### Pilot, disclosed

To set the per-child cap, one throwaway `ecbench` session ran on this Mac
with target seed `999001` (one public point per curve, one round, no
warmup, the same three distinct methods). It is kept under [`pilot/`](pilot/)
and is **not** part of the panel. Its outcomes, on the busy host:

| Cell | Outcome | Process wall | Note |
|:--|:--|--:|:--|
| n41 rho | verified | 0.04 s | `S = 0.135` |
| n41 K8 | error: "target-blind rank setup did not verify every base column" | 58.7 s | the rank search is target-blind and run at the fixed rank seed |
| n41 K16 | verified | 25.9 s | 1,312 usable points, 21,226 rank trials, 225 target attempts, `S = 42.57` |
| n53 rho | verified | 0.19 s | `S = 0.176` |
| n53 K8 | same error | 87.6 s | |
| n53 K16 | same error | 185.8 s | |

The cap is 600 s, 3.2× the largest pilot process wall and 23× the largest
verified solve. The rank setup is target-blind with a fixed seed, so a cell
that fails in the pilot is **predicted** to fail identically on all 96
executions of that cell; those rows are run anyway and kept. The pilot's
`S` values are not panel results and are not cited anywhere else. Rough
wall budget from the pilot: about 15 hours on this busy Mac, dominated by
the predicted-failing n53 K16 cells; the cap is not shrunk to fit.

## Hypotheses, stated before the run

The n37 counted figures are: cold counted IC/rho lower-bound quotient K8
**4.5501 [4.110, 5.040]** (16 untouched points,
[`ecbench_n37_k8_k16_20261004`](../ecbench_n37_k8_k16_20261004/RESULT.md))
and K16 **5.0344 [4.346, 5.866]** (eight points,
[`ecbench_n37_rank_columns_20261004`](../ecbench_n37_rank_columns_20261004/RESULT.md)),
with K16/K8 1.0175 [0.973, 1.059] on the same 16 points; and mean counted
target-only (online) GAE K16 **477.33** against K8 **3,384.60**, i.e. a K16/K8
target-only quotient of about 0.141.

- **H1.** The counted cold IC/rho lower-bound quotient for K8 and for K16
  stays above 1 at n41 and at n53. H1 is *retained* for an arm at a curve if
  the arm's secondary interval lies wholly above 1.0, *falsified* if it lies
  wholly below 1.0, and *undecided* otherwise. A cell in which the arm has no
  verified solve on any workload yields no quotient; H1 is then *not
  evaluable* there and the cell is reported as unknown with its failure
  counts, which is itself the finding that the frozen n37 parameters do not
  produce a verified solve at that size.
- **H2.** K16's target-only counted work relative to K8 keeps the n37
  ordering (K16 below K8). Per curve, the statistic is the ratio of sums of
  the five `target_*` phase GAE over paired verified K16 and K8 runs with
  the same bootstrap. H2 is *retained* if the interval lies wholly below 1.0,
  *falsified* if wholly above 1.0, *undecided* otherwise, and *not evaluable*
  if either arm has no verified run at that curve.

## Decision rules

Per curve, exactly the n37 native-wall protocol's two rows, applied to the
counted unit first and to L0 wall only as an exploratory echo:

1. **Primary target-zero row.** For the designated workload, the mean over
   the five paired measured rounds of each arm's cold counted GAE and the
   resulting K8/rho, K16/rho and K16/K8 quotients. These are lower bounds
   while native work (field arithmetic, hashing, allocation, table probes,
   modular combination) is counted but unpriced.
2. **Secondary ratio of sums.** Over all 16 workloads × 5 rounds, the
   per-curve `ecbench compare` ratio of sums of `S` with its 20,000-resample
   workload-block 95% percentile interval at seed `202610055103`, for
   `rho-strong`→`ic-k8`, `rho-strong`→`ic-k16`, `ic-k8`→`ic-k16` and
   `ic-k16`→`ic-k16-control`. A comparison involving a failed run is
   `incomplete` and is reported as such; a quotient is quoted only from
   verified pairs.
3. **Counted-unit decision**, after `ecbench verify --replay-all` reproduces
   every verified measured run:
   - `counted_quotients_above_one`: every IC arm verified on every
     workload of both curves and every IC/rho interval lies wholly above 1;
   - `frozen_parameters_do_not_transfer`: at least one IC arm has zero
     verified solves at some curve, and every IC/rho interval that exists
     lies wholly above 1;
   - `counted_lead_candidate`: some IC/rho interval lies wholly below 1
     (a lower-bound quotient; it would still need native pricing,
     independent replay and an L2 host before meaning anything);
   - `undecided` otherwise.
4. **A/A.** The counted K16/control ratio must be exactly 1.0 over the
   verified pairs (determinism check, ecbench README §8); its L0 wall
   deviation is reported descriptively only.
5. **Wall.** Every run on this host earns L0. The target-zero online-wall
   ratios and the secondary online ratio of sums with the same bootstrap
   are computed and labelled exploratory; `online_speedup` and
   `fully_priced_cold_speedup` stay `null`. The L2 gate of
   `AGENTS.md` and the n37 protocol remains pending and is addressed by
   [`L2_RUNBOOK.md`](L2_RUNBOOK.md), which was not executed in this round.

Classification by AGENTS.md §3 is expected to be **accounting**: no
admitted speedup, counted lower bounds, a new size priced for the first time
with the n37 method unchanged.

## Inadmissible

Changing any method parameter or `K`; changing the cap after a result;
replacing or dropping a target; dropping a failed, timed-out or
out-of-memory row; dividing counted lower bounds into a speed claim;
quoting an L0 wall ratio as a speedup; pooling the two curves into one
quotient; inferring n83 or ECC2K-130 behaviour from n41/n53.

## Cost accounting

768 child executions on one Apple M4 Pro (macOS, L0, host shared with
other long-running experiments), estimated ~15 wall hours from the pilot;
`ecbench verify --replay-all` replays every verified measured run once more
(the failing cells are not replayed by construction); the analyzer and
comparisons are seconds. No cloud or paid resource is used in this round.

## Reproduction

```sh
cargo build --release --bin ecbench
target/release/ecbench plan --spec research/ecbench_n41_n53_shared_rank_20261005/SPEC.json --json
target/release/ecbench run --spec research/ecbench_n41_n53_shared_rank_20261005/SPEC.json \
  --out research/ecbench_n41_n53_shared_rank_20261005/sessions/mac-l0 --cpus none --allow-busy --wait
target/release/ecbench verify --dir research/ecbench_n41_n53_shared_rank_20261005/sessions/mac-l0 \
  --replay-all --exit-code --out research/ecbench_n41_n53_shared_rank_20261005/AUDIT.json
for pair in "rho-strong ic-k8" "rho-strong ic-k16" "ic-k8 ic-k16" "ic-k16 ic-k16-control"; do
  set -- $pair
  target/release/ecbench compare --dir research/ecbench_n41_n53_shared_rank_20261005/sessions/mac-l0 \
    --a "$1" --b "$2" --resamples 20000 --seed 202610055103 --save
done
cargo build --release --example ecbench_n41_n53_shared_rank_analyze
target/release/examples/ecbench_n41_n53_shared_rank_analyze \
  research/ecbench_n41_n53_shared_rank_20261005 \
  research/ecbench_n41_n53_shared_rank_20261005/sessions/mac-l0 \
  research/ecbench_n41_n53_shared_rank_20261005/AUDIT.json /tmp/n41-n53-decision.json
cmp /tmp/n41-n53-decision.json research/ecbench_n41_n53_shared_rank_20261005/DECISION.json
```
