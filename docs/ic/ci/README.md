# ic end-to-end benchmark (CI gate)

`.github/workflows/ic-e2e-benchmark.yml` runs the **whole** index-calculus
method — select → collect → logs → solve, then the signed-Frobenius ρ baseline
on the same known-answer targets in the same process — on three frozen ledger
rungs, and gates the result against the reference in this directory. It is the
CI form of the rule in `AGENTS.md` §2 that end-to-end speed is the measure of
speed: nothing on this page is a stage number.

| rung | base | `n` | log₂ r | targets | frozen ρ/IC whole-process | crosses ρ e2e |
|:--|:--|--:|--:|--:|--:|:--|
| `docs/ic/params/k0n31.json` | pruned divisor, 35 columns | 31 | 20.5 | 32 | 0.03 | no |
| `docs/ic/params/k0n41-subgroup.json` | subgroup, 5248 points | 41 | 39.0 | 32 | 2.70 | yes |
| `docs/ic/params/k0n53-subgroup.json` | subgroup, 15264 points | 53 | 44.3 | 32 | 2.50 | yes |

The frozen figures are `ic-e2e-reference-v1.json`, whose `frozen_from` block
records the commit, the `ic` binary hash and the host they came from. The
three rungs take about 50 s on four cores.

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

## Running it locally

```bash
cargo build --release --bin ic
python3 scripts/ic_e2e_benchmark.py run --ic target/release/ic --output /tmp/ic-e2e \
    --params docs/ic/params/k0n31.json docs/ic/params/k0n41-subgroup.json docs/ic/params/k0n53-subgroup.json
python3 scripts/ic_e2e_benchmark.py check --output /tmp/ic-e2e --reference docs/ic/ci/ic-e2e-reference-v1.json
python3 -m unittest discover -s scripts -p 'test_ic_e2e_benchmark.py'
```

## Re-freezing after a deliberate change

A change that moves a counter — a different probe schedule, a different
descent, a retuned ρ — is allowed; what is not allowed is landing it without
recording the before and after. In the same pull request:

```bash
python3 scripts/ic_e2e_benchmark.py freeze --output /tmp/ic-e2e \
    --reference-out docs/ic/ci/ic-e2e-reference-v2.json --note "why v2 supersedes v1"
```

then point `IC_E2E_REFERENCE` in the workflow at the new file. `freeze`
refuses to overwrite, so `v1` stays in the tree as the "before" mark
(`AGENTS.md` §7: a superseded figure moves, it does not vanish). Classify the
change in the PR by the §3 test — a lower counter with a lower total is
engineering; a lower counter with a higher total is relabelling — and update
the scoreboard if the figures it cites moved. For a reference frozen on CI
hardware rather than a dev host, download the run artifact and `freeze` from
it; the manifest inside records the host.

Adding a rung is the same operation: add its parameter file (with
`baseline.rho` on) to the `run` step and to a new reference. The gate refuses a
rung that was run but not frozen and a frozen rung that was not run.
