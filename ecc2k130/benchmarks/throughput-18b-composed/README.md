# Composed 18 B/s table-walk result

Target declared before confirmation: beat the current **17.9 B complete scalar
updates/s** table-walk baseline on one RTX PRO 6000, with the same walk,
product count and geometry and with device reports replayed by the scalar
reference.

## Boundary and result

Unit: billions of completed scalar updates per second (`finished:` only).
The derived one-affine-addition floor remains 22--25 B/s. The measured
reference is the inline-both table walk. Every row is engineering: the
algorithm and its ratio to the generic-group bound are unchanged.

| variant | median B/s | / 17.9 target | / paired baseline | correct | class |
|---|---:|---:|---:|---|---|
| inline-both table walk | 17.949369 | 1.0028 | 1.0000 | 300/300, 0 dropped | reference |
| + ALU polynomial square + paired CLMAD schedule | 17.975480 | 1.0042 | 1.0015 | scout replay passed | engineering |
| **+ reduced-input polynomial-to-ONB conversion** | **18.109759** | **1.0117** | **1.0089** | **300/300, 0 dropped** | **engineering** |
| one-addition floor | 22--25 | 1.23--1.40 | 1.23--1.39 | derived | floor |

The winning candidate measured 18.045--18.200 B/s in six balanced-order
repetitions; every paired ratio exceeded one. Its paired geometric mean is
1.00978 with exploratory 95% paired log-ratio t interval
`[1.00666, 1.01291]`. The interval is a runtime diagnostic, not an ECDLP
operation-count claim.

## What composed

- `PACKED_ALU_SQUARE=1` moves the per-update polynomial square from the
  `CLMAD`/FP64 pipe to logic. Nsight put that pipe at 71.6% in the baseline.
- `PACKED_PAIR_CLMUL=1` issues the twelve independent carryless products of a
  paired multiplication before either reduction. With ALU square enabled,
  the carryless pipe has room for this scheduling.
- `PACKED_FROM_REDUCED=1` uses the five-word inverse transform for the reduced
  polynomial state consumed by the selector, instead of zero-extending it
  through the nine-word unreduced-product transform.

Each lever existed independently and had previously been neutral or negative.
Their composition is the measured result. Products remain 5.3125 per update;
no work moved outside the timer.

## Reproduce

Build both arms with the common flags from
[`../throughput-25b-inline/build-config.json`](../throughput-25b-inline/build-config.json).
The candidate additionally sets:

```text
PACKED_ALU_SQUARE=1 PACKED_PAIR_CLMUL=1 PACKED_FROM_REDUCED=1
```

Run each binary with:

```sh
./ecc2k130 --curve 131 --packed --bench --steps 1024 --launches 32 --verify 0
./ecc2k130 --curve 131 --packed --threads 1024 --steps 16 --launches 16 \
  --dp-weight 50 --verify 300
```

[`result.json`](result.json) freezes the samples, order, hashes, ratios and
configuration. `balanced-comparison.log` is the raw six-repetition comparison;
the two verification logs retain the correctness rows. The initial build and
scout logs are retained alongside them.

## Rejected profiler follow-up

Nsight reported 56% excessive shared-memory wavefronts from random table
lookups, but moving selection LUTs to the read-only global/L2 path regressed to
16.001 B/s. Compacting shared addends to 33,536 bytes enabled three blocks per
SM at 80 registers and regressed further to 15.432 B/s. Those opt-in layouts
remain source-visible as negative engineering evidence; neither is part of the
winning configuration.

After the winning composition, `PACKED_TOP_HOIST=1` was also retested because
inlining had lowered the old build's register pressure. It measured 18.132773
B/s against an alternating 18.141865 B/s control, ratio 0.99950, and won only
one of six pairs. Its 300-report replay passed, but it remains rejected.
`rejected-top-hoist.log` preserves that follow-up.
