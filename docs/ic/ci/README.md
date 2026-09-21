# ic end-to-end benchmark (CI gate)

`.github/workflows/ic-e2e-benchmark.yml` runs the **whole** index-calculus
method — select → collect → logs → solve, then the signed-Frobenius ρ baseline
on the same known-answer targets in the same process — on three frozen ledger
rungs, and gates the result against the reference in this directory. It is the
CI form of the rule in `AGENTS.md` §2 that end-to-end speed is the measure of
speed: nothing on this page is a stage number.

| rung | base | `n` | log₂ r | targets | ρ/IC whole-process, dev host → CI runner | crosses ρ e2e |
|:--|:--|--:|--:|--:|:--|:--|
| `docs/ic/params/k0n31.json` | pruned divisor, 35 columns | 31 | 20.5 | 32 | 0.03 → 0.04 | no |
| `docs/ic/params/k0n41-subgroup.json` | subgroup, 5248 points | 41 | 39.0 | 32 | 2.70 → 2.26 | yes |
| `docs/ic/params/k0n53-subgroup.json` | subgroup, 15264 points | 53 | 44.3 | 32 | 2.50 → 2.03 | yes |
| `docs/ic/params/k0n61-subgroup-wide.json` | subgroup, 36112 points, compact table, 2-summand descent | 61 | 47.2 | 32 | 1.87 → 1.77 | yes |

The workflow gates against **`ic-e2e-reference-v4.json`**, frozen from the
artifact of the workflow's own run on a GitHub `ubuntu-latest` runner so that
the wall ratios are matched to the host they are compared on. The earlier
references stay in the tree as "before" marks: `v1` (three rungs, dev host),
`v2` (three rungs, runner), `v3` (four rungs, dev host). Every pinned counter is
identical across all of them and across the differently compiled `ic` binaries
— including the 652,056,328-entry compact pair table at `n = 61` — while the
whole-process ratios are 6–20 % lower on the runner, which is the size of the
host effect this gate must not mistake for a regression. Each reference's
`frozen_from` block records its commit, `ic` binary hash and host. The four
rungs take about 3.5 min on the runner's four cores; the `n = 61` rung peaks
at 4.5 GB of resident memory.

Read the ladder down the rows: at a fixed batch of 32 targets the whole-process
margin over ρ **shrinks** as the subgroup grows, 2.26 → 2.03 → 1.77 from 39 to
47 bits, while the charged (descent-only) ratio rises 10.8 → 26 → 76. That is
the §5 lesson in one table: the stage number improves as the method's
end-to-end standing worsens, because precomputation grows faster than the
descent saves.

## What the gate checks

Per rung, fail closed (`scripts/ic_e2e_benchmark.py check`):

1. **Correctness.** Workflow status `complete`; every target recovered,
   verified as `[d]G = Q`, and equal to the planted scalar; every ρ walk
   verified on the same target. A row without a verified answer is not a
   result, on either side.
2. **Pinned counters, exact.** Factor-base points and columns, pair-table
   pairs, collection trials, summands scanned, relations, descent trials, ρ
   iterations and ρ group additions must equal the reference bit for bit.
   Every one of them is seeded and deterministic (checked across repeated
   runs). A drift is not a failure of the code, it is a change of algorithm,
   and the gate exists to make that change deliberate: freeze a new reference
   beside the old one (below) rather than let it pass silently.
3. **End-to-end wall ratio.** `ρ seconds / (select + collect + logs + pair
   table + descent seconds)`, both sides on the same host in the same process,
   may not fall below the frozen ratio by more than the tolerance (default 50%,
   `--wall-ratio-tolerance`), and a rung that crossed ρ end to end when frozen
   must still cross. Wall time is a practicality note under `AGENTS.md` §6,
   which is why it is gated only as a paired same-host ratio with a wide
   tolerance and never as an absolute number.
4. **ρ is a real opponent.** ρ's group additions per target per √r must lie
   within a factor of 3 of the expected `√(π/4n)` for the signed-Frobenius
   walk on `A = 2n` classes. A broken baseline makes every ratio meaningless,
   which is how an earlier revision produced a spurious crossover.

The summary table is written to the job summary and, with every report,
stderr, the resumable run state and the manifest, uploaded as an artifact.

## What passing means, and does not

Passing says: the method still recovers every logarithm, still does exactly
the work it did when frozen, and its end-to-end standing against ρ on these
rungs has not regressed on this hardware. That is **engineering at most** in
the `AGENTS.md` §3 sense, never an advance.

It does not say the method is faster than ρ in the repository unit. The IC
side's `S = ops / √r` is reported as **null** on every row, on purpose: the
pinned IC counters are in mixed units — group additions for the pair table and
for ρ, table lookups for the summand scan, scalar multiplications for the
probes — and no measured conversion between them is recorded here. `AGENTS.md`
§8 says to leave such a quantity null rather than infer it. The `whole_process`
verdicts at degrees 41 and 53 are same-host wall statements about a batch of 32
targets and carry exactly the caveats of
`docs/ic/runs/koblitz-subgroup-bases-20260912.json` (`what_this_is_not`).
The charged ρ/IC column is advisory and is not gated: on its own it is the
stage number this gate was built not to be fooled by.

## Run history: the paired interval a runtime claim needs

`AGENTS.md` §8 lets a *runtime* claim stand only on paired baseline/candidate
reruns on matched hardware with a 95 % interval that excludes no improvement.
Each passing CI run is one such paired rerun (rho and IC in the same process,
same runner image), and `scripts/ic_e2e_history.py` keeps them:

```bash
gh run download <run-id> --dir /tmp/art            # the workflow's artifact
python3 scripts/ic_e2e_history.py add --run-dir /tmp/art/ic-e2e-benchmark-<run-id> \
    --run-id <run-id> --url https://github.com/aburan28/crypto/actions/runs/<run-id>
python3 scripts/ic_e2e_history.py report --runs docs/ic/ci/runs/*.json \
    --reference docs/ic/ci/ic-e2e-reference-v4.json
```

`add` writes a compact record (host, `ic` hash, commit, per-rung counters and
whole-process seconds) to `docs/ic/ci/runs/<run-id>.json` and never overwrites;
`report` pools only runs whose counters are identical to each other and to the
frozen reference — different counters are a different algorithm — and gives,
per rung, the per-run rho/IC ratios, a t-based 95 % interval, and the same for
rho − IC in seconds. The rung's claim "faster than rho end to end" is set only
with at least three runs and both intervals excluding no improvement. It is a
statement about wall time on that runner; `S` stays null.

Current history (5 runs on the `ubuntu-latest` image, see `runs/`):

| rung | log₂ r | runs | rho/IC per run | 95 % CI | runtime claim |
|:--|--:|--:|:--|:--|:--|
| `k0n31` | 20.5 | 5 | 0.036, 0.036, 0.037, 0.037, 0.038 | [0.036, 0.037] | no (IC slower) |
| `k0n41-subgroup` | 39.0 | 5 | 2.228, 2.234, 2.256, 2.206, 2.251 | **[2.210, 2.260]** | faster than rho, end to end |
| `k0n53-subgroup` | 44.3 | 5 | 1.996, 2.064, 2.034, 2.004, 1.957 | **[1.961, 2.062]** | faster than rho, end to end |
| `k0n61-subgroup-wide` | 47.2 | 3 | 1.774, 1.696, 1.784 | **[1.632, 1.871]** | faster than rho, end to end |

Add each new passing run's artifact; the table above is regenerated from the
`report` output, not edited by hand.

## The `S` column: instruction counts, the unit the rest of the ledger uses

The gate leaves the index-calculus `S` null because its counters are in mixed
units. `scripts/ic_e2e_instructions.py` fills that column the way the
beats-rho tournament thread does: **valgrind callgrind instructions (Ir)**,
`S = Ir / (targets · √r)`, one single-threaded pass per rung
(`RAYON_NUM_THREADS=1`, so rayon workers do not spin under callgrind's
serialised threading and the count is work, not scheduling). `Ir_rho` is the
inclusive cost of `koblitz_signed_frobenius_rho_with_progress` (all 32 calls,
setup included) read from `callgrind_annotate --inclusive=yes`; `Ir_ic` is the
rest of the process, so start-up, curve construction, target resolution and
JSON output are charged to the index-calculus side. Single-threaded totals
reproduce to about 2 × 10⁻⁵. Callgrind is ~50× slower than native, so this is
a frozen measurement re-run when the `ic` binary changes, not a CI step.

```bash
python3 scripts/ic_e2e_instructions.py measure --ic target/release/ic --output /tmp/ir \
    --params docs/ic/params/k0n31.json docs/ic/params/k0n41-subgroup.json docs/ic/params/k0n53-subgroup.json
python3 scripts/ic_e2e_instructions.py report --instructions /tmp/ir/instructions.json
```

Frozen: [`instructions/ladder-20260920-n31-41-53-61.json`](instructions/ladder-20260920-n31-41-53-61.json)
(`ic` built from `a6ad6236`, valgrind 3.22.0, x86_64):

| rung | log₂ r | `S`, index calculus | `S`, rho | rho / IC in Ir | rho / IC in wall (runner, 95 % CI) | IC Ir in the pair table / in the scan |
|:--|--:|--:|--:|--:|:--|--:|
| `k0n31` | 20.5 | 167,091 | 5,912 | 0.035 | [0.036, 0.037] | 38 % / 2 % |
| `k0n41-subgroup` | 39.0 | 1,155 | 1,542 | **1.335** | [2.210, 2.260] | 57 % / 41 % |
| `k0n53-subgroup` | 44.3 | 2,292 | 2,761 | **1.204** | [1.961, 2.062] | 47 % / 52 % |
| `k0n61-subgroup-wide` | 47.2 | 2,723 | 2,724 | **1.001** | [1.632, 1.871] | 86 % / 11 % |

**Read the two ratio columns together.** The wall ratio compares a pipeline
whose pair table and collection run on every core against a single-threaded
rho; the instruction count removes that parallelism, and the end-to-end
margin is 1.34, 1.20 and **1.00** in work where the runner clock said 2.26,
2.03 and 1.77. At 47 bits the whole-process advantage in the ledger's unit
is gone: what the clock shows there is the pipeline's four cores against
rho's one, and the pair table alone is 86 % of the pipeline's instructions.
That is the AGENTS.md §6 rule — wall time is never the headline — with the
number attached, and it is why the gate pins counters and shows the charged
column as advisory only. In the repository unit the pipeline beats matched
rho end to end by 20–34 % at 39–44 bits on a 32-target batch, reaches parity
at 47, and loses by 28× at 20. Class: engineering (no ratio to the counting
floor moved; nothing here does what a generic algorithm cannot). `S_rho` at
n = 31 (5,912) sits beside the tournament thread's 5,388 on comparable
cells, which is the check that the two threads are on one axis. Callgrind
totals sum every thread; `RAYON_NUM_THREADS=1` still leaves one rayon worker,
to which the table build is handed, so the totals are whole-process work.

## Running it locally

```bash
cargo build --release --bin ic
python3 scripts/ic_e2e_benchmark.py run --ic target/release/ic --output /tmp/ic-e2e \
    --params docs/ic/params/k0n31.json docs/ic/params/k0n41-subgroup.json docs/ic/params/k0n53-subgroup.json docs/ic/params/k0n61-subgroup-wide.json
python3 scripts/ic_e2e_benchmark.py check --output /tmp/ic-e2e --reference docs/ic/ci/ic-e2e-reference-v4.json
python3 -m unittest discover -s scripts -p 'test_ic_e2e_benchmark.py'
```

On a machine unlike the runner the counters will still be identical and the
wall ratios may sit outside the tolerance in either direction; that is the
host effect, not a verdict. Freeze a local reference for local work.

## Re-freezing after a deliberate change

A change that moves a counter — a different probe schedule, a different
descent, a retuned ρ — is allowed; what is not allowed is landing it without
recording the before and after. In the same pull request:

```bash
python3 scripts/ic_e2e_benchmark.py freeze --output /tmp/ic-e2e \
    --reference-out docs/ic/ci/ic-e2e-reference-v5.json --note "why v5 supersedes v4"
```

then point `IC_E2E_REFERENCE` in the workflow at the new file. `freeze`
refuses to overwrite, so the superseded reference stays in the tree as the
"before" mark (`AGENTS.md` §7: a superseded figure moves, it does not
vanish). Classify the change in the PR by the §3 test — a lower counter with
a lower total is engineering; a lower counter with a higher total is
relabelling — and update the scoreboard if the figures it cites moved.
Freeze from CI hardware, as v2 and v4 were: let the workflow run once against the old
reference (it will fail on the drift, which is the point), download the run
artifact, `freeze` from it, and land the new reference in the same PR. The
manifest inside the artifact records the host and binary hash.

Adding a rung is the same operation: add its parameter file (with
`baseline.rho` on) to the `run` step and to a new reference. The gate refuses a
rung that was run but not frozen and a frozen rung that was not run.
