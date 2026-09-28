# ic end-to-end benchmark (CI gate)

`.github/workflows/ic-e2e-benchmark.yml` runs the **whole** index-calculus
method — select → collect → logs → solve, then the signed-Frobenius ρ baseline
on the same known-answer targets in the same process — on three frozen ledger
rungs, and gates the result against the reference in this directory. It is the
CI form of the rule in `AGENTS.md` §2 that end-to-end speed is the measure of
speed: nothing on this page is a stage number.

| rung (v6) | base | `n` | log₂ r | targets | relations | summands scanned | descent trials | IC whole process, same host: v5 rung → v6 rung | ρ/IC whole-process, v5 → v6 (same host) | crosses ρ e2e |
|:--|:--|--:|--:|--:|--:|--:|--:|:--|:--|:--|
| `docs/ic/params/k0n31.json` | pruned divisor, 35 columns | 31 | 20.5 | 32 | 78 → 78 | 277,760 → 277,760 | unchanged | unchanged | 0.08 → 0.06 | no |
| `docs/ic/params/k0n41-subgroup-wide.json` | subgroup, 6560 points, 80 columns | 41 | 39.0 | 32 | 67 → 83 | 2,770,944 → 2,287,616 | 726 → 346 | 0.316 s → 0.258 s | 8.34 → 10.28 | yes |
| `docs/ic/params/k0n53-subgroup-wide.json` | subgroup, 30528 points, 288 columns | 53 | 44.3 | 32 | 147 → 334 | 31,694,848 → 16,064,512 | 7,803,621 → 2,057,777 | 1.775 s → 1.013 s | 15.60 → 27.28 | yes |

**`ic-e2e-reference-v6.json`** is what the workflow gates against. It
supersedes v5 because the degree-41 and degree-53 rungs now use **wider
factor bases** — 6000 and 30000 requested points (6560 and 30528
materialised) against 4759 and 15000 — with the collection window scaled
with the base (207 and 954 summands), on the same curves, targets, seeds,
solver, descent summand count and aimed, windowed collection as the v5
rungs. The v5 rung files are unchanged.

The width was chosen by a sweep of the base size on each curve (ρ off,
same host; `docs/ic/perf/OPTIMIZATION_PLAN.md` log), which only paid off
once the folded pair-table build stopped dominating: keyed once with the
scan's own kernels and placed by partitioning, the build at degree 53 fell
from 0.09 s to 0.04 s at 15,000 points and from 1.03 s to 0.33 s at
45,000. A wider base needs more relations (columns grow like `|F|`) but
finds them in fewer summands (pair sums grow like `|F|²`), and it cuts the
descent's trials by 3.8× at degree 53 and 2.1× at degree 41. The sweep's
minimum was flat between about 20,000 and 30,000 points at degree 53 and
around 6,000 at degree 41; unit size and window scale moved it by less
than run-to-run spread. By `AGENTS.md` §3 this is **engineering**: a
parameter of the same method, re-chosen because a stage got cheaper; the
generic-group floor does not move.

The IC figures in the table are from one host and one `ic` binary, v5 rung
beside v6 rung; v6 itself is frozen on a dev host, as v3–v5 were, and its
`frozen_from` records the commit, the binary hash and the host.

| rung (v5) | base | `n` | log₂ r | targets | relations | summands scanned | IC whole process, same host: v4 rung → v5 rung | ρ/IC whole-process, v4 → v5 (dev) | crosses ρ e2e |
|:--|:--|--:|--:|--:|--:|--:|:--|:--|:--|
| `docs/ic/params/k0n31.json` | pruned divisor, 35 columns | 31 | 20.5 | 32 | 78 → 78 | 277,760 → 277,760 | 0.074 s → 0.073 s | 0.04 → 0.08 | no |
| `docs/ic/params/k0n41-subgroup-aimed.json` | subgroup, 5248 points, 64 columns | 41 | 39.0 | 32 | 124 → 67 | 15,072,256 → 2,770,944 | 0.421 s → 0.297 s | 3.39 → 7.75 | yes |
| `docs/ic/params/k0n53-subgroup-aimed.json` | subgroup, 15264 points, 144 columns | 53 | 44.3 | 32 | 458 → 147 | 242,514,432 → 31,694,848 | 4.41 s → 1.65 s | 3.00 → 14.84 | yes |

The v4 rungs were `docs/ic/params/k0n41-subgroup.json` and
`docs/ic/params/k0n53-subgroup.json`; their files are unchanged, and the
ledger records that cite them still mean what they meant. Earlier columns
of the ρ/IC ratio, for the v4 rungs: v1 (dev) 2.70 / 2.50, v2 (runner)
2.23 / 2.00, v3 (dev) 1.79 / 1.51, v4 (dev) 3.39 / 3.00.

The v5 rungs, above as the "before" column, are the table below. v5
(`ic-e2e-reference-v5.json`) superseded v4 for two reasons, one of them a
gate fix.

*The rungs collect the way the measurements say to.* The degree-41 and
degree-53 rungs keep their curve, base, targets, pair-table tier and
descent, and collect with a window of about `|F|/32` summands aimed at the
least-mentioned columns: one planned unit, then further units until every
column is determined (`k0n41-subgroup-aimed.json`, window 164 and units of
1024 probes; `k0n53-subgroup-aimed.json`, window 477 and units of 2048).
The window and the aim were measured and classified in
`docs/ic/runs/koblitz-collection-window-20260921.json` (a full `m = 3`
scan meets each triple three times and keeps one; a window keeps all
three chances) and `docs/ic/runs/koblitz-collection-aim-20260922.json`
(aimed at the least-mentioned columns, relations fall to within 3 % of
the counting floor). Both left the ledger rungs sweeping, as an open
question; v5 answers it for the gate. At degree 53 the rung scans 7.65×
fewer summands for 147 relations instead of 458, and the whole IC process
on one host falls from 4.41 s to 1.65 s. Descent trials, ρ and the stored
pairs are identical to v4's. By `AGENTS.md` §3 this is **engineering**:
the counters and the total both fell, the generic-group floor did not
move, and no decomposition is found that the sweep could not find.

*The collection counters now see extension units.* They used to be read
from the collect stage, which reports only the planned units; units
collected inside the logs stage, because the planned ones left a column
undetermined, were work the gate could not see. A v5 rung does almost all
of its collecting there, so reading the old way would have pinned 2 of 67
relations at degree 41. `scripts/ic_e2e_benchmark.py` now reads the logs
stage's totals, which equal the collect stage's whenever nothing is
extended: the v4 rungs check identically against v4 under the new read.
A reused logs stage, or a logs total below the collect stage's, is
refused.

v5 is frozen on a dev host, as v3 and v4 were; `frozen_from` records the
commit, the `ic` binary hash and the host. The IC whole-process figures
in the table are both from that host and the same `ic` binary, v4 rung
beside v5 rung.

v4 superseded v3 because the workflow now computes its own probing volume
and passes it to `PairSumTable::build_within_for`
(`docs/ic/runs/koblitz-probe-volume-20260921.json`), so the tier is
priced for the run being built rather than for the volume the cost model
was calibrated at. All three rungs move to the **folded** tier: they
probe between 1.45× and 1,266× less than that calibration volume, and
less probing to amortise the build over is when the fold wins. Priced in
group additions at each rung's own degree, that is 1.371× fewer
operations at degree 41 and 1.229× at degree 53.

v3 superseded v2 because the tier came to be chosen by a measured cost
model rather than a fixed ladder
(`docs/ic/runs/koblitz-tier-crossover-20260921.json`), and all three
rungs built the **compact** tier where v2 built the full one.

That change went through this gate undetected the first time, which is
worth recording. Every pinned counter was identical, because
`pair_table_stored_pairs` is `|F|(|F|+1)/2` for the full tier *and* the
compact one — they hold the same pairs and differ only in how a pair is
stored. The gate reported "counters identical" while the algorithm had
changed, which is the one thing it exists to refuse. The tier is
therefore pinned beside the counters now, a report without one is
refused, and `TierPinTests` fails if that check is removed. A reference
frozen before tier reporting (v1, v2) carries no tier and skips that
comparison, which is why v3 had to be frozen rather than v2 amended.

v3's ρ/IC ratios are lower than v2's for two reasons that should not be
confused: it is frozen on a dev host rather than the runner (the v1→v2
column shows that effect on its own), and the compact tier really is a
different algorithm. The counters and the tier are host-independent; the
wall ratios are not, and are gated only as a paired same-host ratio.

The previous references stay in the tree as the "before" marks.
`ic-e2e-reference-v2.json` was frozen from the
artifact of the workflow's own first run on a GitHub `ubuntu-latest` runner
so that the wall ratios are matched to the host they are compared on.
`ic-e2e-reference-v1.json` is the same three rungs frozen on a 4-vCPU dev
host: every pinned counter is
identical between the two hosts (and between the two differently compiled
`ic` binaries), while the whole-process ratios are 17–20 % lower on the
runner, which is the size of the host effect this gate must not mistake for
a regression. Each reference's `frozen_from` block records its commit, `ic`
binary hash and host. The three rungs take about 60 s on the runner's four
cores.

## What the gate checks

Per rung, fail closed (`scripts/ic_e2e_benchmark.py check`):

1. **Correctness.** Workflow status `complete`; every target recovered,
   verified as `[d]G = Q`, and equal to the planted scalar; every ρ walk
   verified on the same target. A row without a verified answer is not a
   result, on either side.
2. **Pinned counters, exact, and the pair-table tier.** Factor-base points
   and columns, pair-table pairs, collection trials, summands scanned,
   relations, descent trials, ρ iterations and ρ group additions must equal
   the reference bit for bit — and so must the tier, which none of those
   counters can see: `pair_table_stored_pairs` is `|F|(|F|+1)/2` for the
   full and the compact tier alike.
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
    --reference docs/ic/ci/ic-e2e-reference-n61-v4.json
```

`add` writes a compact record (host, `ic` hash, commit, per-rung counters and
whole-process seconds) to `docs/ic/ci/runs/<run-id>.json` and never overwrites;
`report` pools only runs whose counters are identical to each other and to the
frozen reference — different counters are a different algorithm — and gives,
per rung, the per-run rho/IC ratios, a t-based 95 % interval, and the same for
rho − IC in seconds. The rung's claim "faster than rho end to end" is set only
with at least three runs and both intervals excluding no improvement. It is a
statement about wall time on that runner; `S` stays null.

Current history (5 runs on the `ubuntu-latest` image, see `runs/`), all on the
four-rung `n = 61` ladder archived below, not on the v5 aimed rungs:

| rung | log₂ r | runs | rho/IC per run | 95 % CI | runtime claim |
|:--|--:|--:|:--|:--|:--|
| `k0n31` | 20.5 | 5 | 0.036, 0.036, 0.037, 0.037, 0.038 | [0.036, 0.037] | no (IC slower) |
| `k0n41-subgroup` | 39.0 | 5 | 2.228, 2.234, 2.256, 2.206, 2.251 | **[2.210, 2.260]** | faster than rho, end to end |
| `k0n53-subgroup` | 44.3 | 5 | 1.996, 2.064, 2.034, 2.004, 1.957 | **[1.961, 2.062]** | faster than rho, end to end |
| `k0n61-subgroup-wide` | 47.2 | 3 | 1.774, 1.696, 1.784 | **[1.632, 1.871]** | faster than rho, end to end |

Add each new passing run's artifact; the table above is regenerated from the
`report` output, not edited by hand.

## Running it locally

```bash
cargo build --release --bin ic
python3 scripts/ic_e2e_benchmark.py run --ic target/release/ic --output /tmp/ic-e2e \
    --params docs/ic/params/k0n31.json docs/ic/params/k0n41-subgroup-aimed.json docs/ic/params/k0n53-subgroup-aimed.json
python3 scripts/ic_e2e_benchmark.py check --output /tmp/ic-e2e --reference docs/ic/ci/ic-e2e-reference-v5.json
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
Freeze from CI hardware, as v2 and n61-v4 were: let the workflow run once against the old
reference (it will fail on the drift, which is the point), download the run
artifact, `freeze` from it, and land the new reference in the same PR. The
manifest inside the artifact records the host and binary hash.

Adding a rung is the same operation: add its parameter file (with
`baseline.rho` on) to the `run` step and to a new reference. The gate refuses a
rung that was run but not frozen and a frozen rung that was not run.

## Archived `n = 61` references (from #506)

#506 froze a four-rung ladder that adds `docs/ic/params/k0n61-subgroup-wide.json`
to the three v4-era rungs. Its two references were named `v3` and `v4` on that
branch, which collided with main's `v3` and `v4`; they are kept under
**`ic-e2e-reference-n61-v3.json`** (dev host) and
**`ic-e2e-reference-n61-v4.json`** (GitHub `ubuntu-latest` runner). They are
archived inputs, not the gated reference: the workflow gates against
`ic-e2e-reference-v5.json` above. The run records in `runs/` that this ladder
produced keep the `reference` path they were recorded with, so their `v3`/`v4`
mean these `n61` files. The records and a report over them come from
`scripts/ic_e2e_history.py`.

| rung | base | `n` | log₂ r | targets | ρ/IC whole-process, dev host → CI runner | crosses ρ e2e |
|:--|:--|--:|--:|--:|:--|:--|
| `docs/ic/params/k0n31.json` | pruned divisor, 35 columns | 31 | 20.5 | 32 | 0.03 → 0.04 | no |
| `docs/ic/params/k0n41-subgroup.json` | subgroup, 5248 points | 41 | 39.0 | 32 | 2.70 → 2.26 | yes |
| `docs/ic/params/k0n53-subgroup.json` | subgroup, 15264 points | 53 | 44.3 | 32 | 2.50 → 2.03 | yes |
| `docs/ic/params/k0n61-subgroup-wide.json` | subgroup, 36112 points, compact table, 2-summand descent | 61 | 47.2 | 32 | 1.87 → 1.77 | yes |

Every pinned counter is identical across the two `n61` references, including
the 652,056,328-entry compact pair table at `n = 61`. The four rungs take about
3.5 min on the runner's four cores; the `n = 61` rung peaks at 4.5 GB of
resident memory.

Read the ladder down the rows: at a fixed batch of 32 targets the whole-process
margin over ρ **shrinks** as the subgroup grows, 2.26 → 2.03 → 1.77 from 39 to
47 bits, while the charged (descent-only) ratio rises 10.8 → 26 → 76. That is
the §5 lesson in one table: the stage number improves as the method's
end-to-end standing worsens, because precomputation grows faster than the
descent saves.
