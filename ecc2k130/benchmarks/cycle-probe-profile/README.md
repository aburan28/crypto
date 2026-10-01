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
three fresh runs with identical run id, geometry, threshold, and 64 steps per
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
