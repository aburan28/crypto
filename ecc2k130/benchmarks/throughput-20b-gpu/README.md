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
| Nsight Compute | diagnose ALU vs CLMAD vs DRAM stall | previous G7 counter collection failed |
| `PACKED_L2_PERSIST=1` | hide compact-field traffic in the 80 MiB persisting L2 window | Nsight on this card showed DRAM ~65% / ALU ~35% at a 16-step diagnostic; eighteen priced the ALU and never asked L2 |

Pair-ILP is not a 20% instruction cut. The honest ALU prior is still **~18.4 B/s**
if dual-issue pays *and* the clock lock converts. Persist is a different
question: if the walk is DRAM-bound, 20 B/s is under the DRAM ceiling this
card can post (~26 B/s if the DRAM share of the 16.7 B/s run went to zero).
If Nsight's ALU pipe is not the logic pipe the static count uses, persist
is a no-op. The receipt decides.

Slot prefetch (`PACKED_SLOT_PREFETCH`) and slot unroll (`UNROLL_SLOTS=2`)
are knobs, default off. They do not change the product count.

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

## Result (2026-09-17)

This RTX PRO 6000 Blackwell Server Edition, driver 595.91.07, CUDA 13.3.73,
automatic 96,256 workers, 50,465,865,728 updates per sample, clocks locked
at 2430 MHz (SM still 2340–2422 under load). Every binary replayed 300
device reports with 0 dropped. Product count 5.3125. Walk-rate engineering,
not an ECDLP exponent claim.

| variant | median B/s | / 20 | / this shipping 14.436 | / 22 B floor | / eighteen 16.667 | class |
|---|---:|---:|---:|---:|---:|---|
| shipping product | 14.436 | 0.722 | 1.000 | 0.656 | 0.866 | reference |
| table walk + byte pivot | 16.474 | 0.824 | 1.141 | 0.749 | 0.988 | engineering |
| + `PACKED_PAIR_ILP=1` | 16.617 | 0.831 | 1.151 | 0.755 | 0.997 | engineering |
| + `PACKED_L2_PERSIST=1` | 17.081 | 0.854 | 1.183 | 0.776 | 1.025 | engineering |
| + `UNROLL_SLOTS=2` | 17.298 | 0.865 | 1.198 | 0.786 | 1.038 | engineering |

No arm clears 20. Pair-ILP was priced at ~9% if ptxas dual-issued and
measured **+0.9%**. Persist was the DRAM lever Nsight reopened and
measured **+2.8%** on the ILP control (SM held 2422 MHz at ~530 W; the
ILP control fell to 2340–2370 MHz at ~555 W). Slot unroll 2 added
**+1.3%** at the same 2422 MHz. The best median is 17.298 B/s, 0.865 of
the target.

This instance's shipping walk is 14.436 against eighteen's 15.116 on a
Modal card of the same SKU; the table+pivot control here is 16.474
against eighteen's 16.667. The clock lock does not hold 2430 under load.

Rejected scouts, not in the table (no 300-report replay): slot prefetch
regressed; a global table at 80 registers resident 3 blocks/SM at ~15.5
B/s (1.5× working set, SM ~2280 MHz); the same global table at 104
registers / 2 blocks was 16.08 B/s, so the LDS copy still wins; unroll 4
matched unroll 2's register count and was slightly slower.

The campaign default stays the shipping walk. This thread stops: every
arm is ≤ 20.0, and the leftover distance is the product and the
reduction.

Frozen source: [summary.json](summary.json), [result.json](result.json).

## Recipe

From `ecc2k130/`, GPU idle, CUDA 13.3.73 on `PATH`:

```sh
bash benchmarks/throughput-20b-gpu/run.sh
```

Automatic workers: that is the geometry that produced 16.56 / 16.67 B/s.
Four binaries plus the two DRAM/scheduling arms that Nsight reopened:
shipping, table+pivot, pair-ILP, L2 persist, and slot-unroll-2 on persist.
Three alternating repetitions, 50,465,865,728 updates per sample. Clocks
locked at 2430 MHz when the driver allows it. Verify uses `--dp-cap 262144`
so automatic occupancy cannot overflow the default 65,536-report buffer.
