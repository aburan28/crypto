# 20 B/s on this RTX PRO 6000

Written before the timings. Walk-rate engineering, not an ECDLP exponent claim.

Eighteen already timed every *priced static* leftover on a Modal RTX PRO 6000
and stopped at 16.667 B/s
([benchmarks/eighteen](../eighteen/README.md)).
This round is a different question: can **this** instance, with clocks locked,
CUDA 13.3.73, Nsight Compute, and one unmeasured scheduling lever, clear
**20 B complete scalar updates/s**?

## Boundaries

- **Floor:** one affine addition per step, ≈ 22–25 B/s
  ([THROUGHPUT-30B.md](../../THROUGHPUT-30B.md),
  [ITERATION-FUNCTION.md](../../ITERATION-FUNCTION.md)). 20 is under the
  floor, so the hardware does not forbid it.
- **Reference:** table walk + byte pivot, 16.667 B/s median on a Modal
  RTX PRO 6000 (driver 580.95.05, automatic 96,256 workers, unlocked clocks
  2385–2407 MHz). Same SKU, not this allocation.
- **Unit:** billions of completed scalar updates per second (`finished:`
  line only, via `codegen.benchreport.parseRate`).
- **Class:** engineering. The table walk is a different iteration and
  cannot join the live campaign. Pair-ILP does not change the walk.

## What this instance adds

| lever | priced effect | why eighteen did not settle it |
|---|---|---|
| SM clocks locked at 2430 MHz | 2430/2400 ≈ +1.25% → ~16.88 B/s | Modal samples ran 2385–2407 MHz |
| Driver 595.91.07 / CUDA 13.3.73 | unknown | eighteen used 580.95.05 |
| `PACKED_PAIR_ILP=1` | ~9% *if* ptxas dual-issues the second product's `clmad` with the first reduction | static SASS is unchanged; a count is not a rate |
| Nsight Compute | diagnose ALU vs CLMAD vs stall | previous G7 counter collection failed |

None of those is a 20% instruction cut. The honest prior is **~18.4 B/s**
if pair-ILP pays its full priced overlap *and* the clock lock converts;
20 still needs the pipes to issue more than the 16.67 run did.

## Acceptance, written before the run

Success is a verified median **> 20.0 B/s** on this RTX PRO 6000 Blackwell
Server Edition, identity matching the requested walk, every sample valid,
same 5.3125 field products per update.

Inadmissible: changing the GPU SKU, counting two GPUs as one, dropping
identity, changing the product count, quoting a model as a rate, or
scraping boost-clock progress lines instead of `finished:`.

A row without a verified answer is not a result. Report replay of 300
device distinguished points is required before timing each binary.

If every arm stays ≤ 20.0, this thread stops. The remaining distance is
the product and the reduction, which are the floor.

## Recipe

From `ecc2k130/`, GPU idle, CUDA 13.3.73 on `PATH`:

```sh
bash benchmarks/throughput-20b-gpu/run.sh
```

Automatic workers: that is the geometry that produced 16.56 / 16.67 B/s.
Three binaries, three alternating repetitions, 50,465,865,728 updates per
sample. Clocks locked at 2430 MHz when the driver allows it.
