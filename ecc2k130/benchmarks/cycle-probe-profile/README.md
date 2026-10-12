# Exact table-v3 cold-probe diagnostic

This default-off profile counts why every exact v3 cold probe returns on the
selected B16 fast2 kernel. It records total history hints, fast2 exact-cycle
hits, general exact cycles by length 1 through 8, probes still open after eight
steps, distinguished-point aborts, exceptional-denominator aborts, affine
`next()` calls, and selected anchor exits.

The terminal counters partition hints:

```
hints = fast2_hits + sum(general_cycles_by_length_1_8)
      + open_eight_step_probes + dp_aborts
      + exceptional_denominator_aborts
```

`anchor_exits` is a subset of exact cycles. `affine_next_calls` counts actual
full affine `next()` calls; fast2's second-denominator check is not counted as
one. The device combines increments by active warp, distinct outcome, and
distinct next-call count before issuing global atomics.

The frozen job builds the current RTX PRO 6000 B16 preset with
`CYCLE_PROFILE=1`, runs the synthetic host reconciliation controls, then starts
three fresh benchmark-mode runs with identical run id, geometry, effective
DP weight 0, and 64 steps per
launch. The only changed input is the launch budget: 1 (early), 8 (medium), and
32 (long). `counts.tsv` retains exact update counts, log hashes, and the full
JSON counter line for comparison. The medium and long runs contain the same
deterministic prefix as the early run.

Launch on Modal only when GPU execution is authorized:

```sh
modal run --detach modal_job.py \
  --job benchmarks/cycle-probe-profile/gpujob.sh \
  --out /tmp/ecc2k-cycle-probe-profile --gpu RTX-PRO-6000
```

This is an instrumented diagnostic. Its elapsed time and displayed iteration
rate are not performance evidence and must not be compared with uninstrumented
throughput results.

## Measured counts

The long prefix completed 3,154,116,608 updates and recorded 7,834,496 hints
(0.248390% of updates). Of those hints, 2,876,286 (36.7131%) were exact fast2
cycles, 4,955,914 (63.2576%) remained open after all eight steps, 2,282 closed
at length four and 14 at length six. No other general cycle length, DP abort,
or exceptional denominator occurred. The probe executed 42,532,810 affine
`next()` calls, 5.4289 per hint and 0.013485 per complete update.

The early and medium deterministic prefixes reconcile independently and show
the hint rate rising from 0.210672% to 0.224326% and then 0.248390%. The exact
rows are retained in [`counts.tsv`](counts.tsv). The 6,024-byte raw archive has
SHA-256 `e08ffba7345e1713b4610ddba498e4d01881ade97b45e347715acadd2e2591bf`.
[`independent-audit.json`](independent-audit.json) checks source hashes, host
controls, log hashes, update budgets and counter reconciliation. The executable
bytes were not retained, so its recorded binary hash is provenance only.

These counts describe `--bench`, which unconditionally uses DP weight 0. The
first frozen job declared 48 in its wrapper metadata, but the client reset it;
all raw backend markers say 0. The source and this note now state the effective
value. No collection-mode inference is made.
