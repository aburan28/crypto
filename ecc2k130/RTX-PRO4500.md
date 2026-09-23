# RTX PRO 4500 (EC2 g7) — the rate, and the decision it settles

> **Update, [POWER-BOUND.md](POWER-BOUND.md):** both receipts below are for
> builds the tree has since superseded on the RTX PRO 6000. The 4500 sits on
> its 165 W limit, so its rate ranks builds by updates per joule; measured on
> a 6000, `make gpu-rtx-pro4500` (the 20 B/s build plus one fused pass per
> step) is 1.072× per joule over the 20 B/s build, which is itself 1.25× the
> table-walk receipt's. **Predicted 7.2 – 7.4 B/s on a 4500, against 5.552
> and 5.107 here — a prediction, not a receipt** (EC2 launches are blocked in
> the session that made it; POWER-BOUND.md §5 has the command that measures
> it). The break-evens below were computed against the 6000's shipping rate
> and have to be redone with both cards on their best builds.

Every audited throughput figure in this tree is for the **RTX PRO 6000**:
14.637530 B/s benchmarking and 14.1 B/s collecting
([THROUGHPUT-30B.md](THROUGHPUT-30B.md), [aws/README.md](aws/README.md)). The
fleet currently running the campaign is on **g7.2xlarge**, which carries the
**RTX PRO 4500**. One `--bench` figure for that part now exists, taken as the
control arm of the iteration-function comparison
([ITERATION-FUNCTION.md](ITERATION-FUNCTION.md) §6,
[`benchmarks/table-walk/comparison.json`](benchmarks/table-walk/comparison.json)):
**5.107 B/s** median of three, shipping walk, automatic workers, and it is
below every break-even in the table below. The `run.sh` receipt described next
is still owed: the figure above was not taken through it, and the collecting
row is still empty.

That gap was not cosmetic. It decides which instance type the campaign should
buy, and the recollected figure sat within a few percent of the break-even,
on the wrong side of it; the measured one sits 20–25% below it.

## Measure it

On an **idle** g7 instance, from the `ecc2k130` directory:

```bash
bash benchmarks/rtx-pro4500/run.sh
```

Writes `benchmarks/rtx-pro4500/result.json` with each command's verbatim output,
following [benchmarks/hardware-limits/result.json](benchmarks/hardware-limits/result.json).

It reads the rate with the tree's own `codegen.benchreport` — `parseRate`, which
accepts only the single synchronised `finished:` line and rejects a run without
one, and `summarizeSamples`, which takes the median across repeats. That is the
method behind the RTX PRO 6000 figures, so the two are comparable by
construction rather than by intention. The periodic progress lines also carry
`M it/s` and read high during boost-clock warmup; a maximum over those would
have sat the 4500 above a same-method comparison by enough to cross both
break-evens below.

The script **refuses to run while another process holds the GPU**. A benchmark
sharing a device with a collecting worker steals SMs from it and has them stolen
back; both figures come out low and neither measures anything. Stop the worker
first — it checkpoints (`--checkpoint`) — or accept that `ALLOW_CONTENTION=1`
produces a number that is *not* comparable to the RTX PRO 6000 figure.

It reuses the arithmetic and layout options of the RTX PRO 6000 preset
(`RTX_PRO6000_ENV` in the [Makefile](Makefile)), which are field-arithmetic and
storage choices rather than GPU-model choices — passed in the names a local
`make` consumes (`BATCH`, `THREADS`, `MINBLOCKS`, `PACKED_*`), not the
`ECC_PACKED_*` spelling, which is the Modal image interface and would leave
every option at its Makefile default. `result.json` records the resulting
`-DECC_*` defines so a reader can confirm the preset reached the build instead
of trusting that it did. It deliberately does **not**
reuse that preset's 385,024 workers: that count is four times automatic on a
188-SM part and has no claim to be right on a smaller one. Automatic scales with
`multiProcessorCount` (`src/main.cu`), so the script starts there.

## The decision this feeds

Spot prices observed **2026-09-16, us-east-1 and us-west-2**, via
`describe-spot-price-history` (cheapest AZ per type):

| region | type | GPU | $/hr |
|---|---|---|---:|
| us-west-2 | g7e.2xlarge | RTX PRO 6000 | 1.2128 |
| us-west-2 | g7.2xlarge | RTX PRO 4500 | 0.5705 |
| us-east-1 | g7e.2xlarge | RTX PRO 6000 | 1.5671 |
| us-east-1 | g7.2xlarge | RTX PRO 4500 | 0.6858 |
| us-east-2 | g7e.2xlarge | RTX PRO 6000 | 1.8920 |
| us-east-2 | g7.2xlarge | RTX PRO 4500 | 1.1847 |

Cost is a property of the work, so the only figure that ranks instance types is
**iterations per dollar**. Against the 6000's 14.1 B/s collecting rate, g7 wins
its region only above a threshold:

| region | g7 must exceed | for g7 to beat g7e per dollar |
|---|---:|---|
| us-west-2 | **6.63 B/s** | 14.1 × 0.5705 / 1.2128 |
| us-east-1 | **6.17 B/s** | 14.1 × 0.6858 / 1.5671 |
| us-east-2 | **8.83 B/s** | 14.1 × 1.1847 / 1.8920 |

A figure of roughly 6 B/s had been stated from recollection; it appears in no
receipt in this tree and is not treated as evidence here. The measured
`--bench` rate is **5.107 B/s** (5.151 / 5.107 / 5.069 over three alternating
repetitions, `--bench --steps 1024 --launches 32`, automatic worker count,
`RTX_PRO6000_ENV` arithmetic and layout, CUDA 13.3.1 in the
`nvidia/cuda:13.3.1-devel` container on an otherwise idle `g7.2xlarge`;
[raw log](benchmarks/table-walk/raw/rtx4500-lut-build-verify-bench.log)).
The card sat at its 165 W limit at 1.92–1.97 GHz throughout. That is below
the break-even in every region: per dollar, g7e beats g7 by 6.63 / 5.107 =
**1.30×** in us-west-2 and 6.17 / 5.107 = **1.21×** in us-east-1. The
collecting figure and the `run.sh` receipt are still to be taken; neither can
move the verdict unless collecting on the 4500 loses less to DP handling than
it does on the 6000 (14.64 → 14.1, 3.7%).

| | value |
|---|---|
| RTX PRO 4500, `--bench` | **5.107 B/s** (median of 3; `benchmarks/table-walk/comparison.json`, not yet re-taken through `run.sh`) |
| RTX PRO 4500, collecting | *pending* |
| RTX PRO 6000, `--bench` | 14.637530 B/s |
| RTX PRO 6000, collecting | 14.1 B/s |

## What is already settled

**The Ada types are a separate, open question.** EC2 g6 (L4) and g6e (L40S)
are `sm_89` and have no rate here either; [ADA-L4-L40S.md](ADA-L4-L40S.md)
carries that decision, the harness that measures it, and the one preset knob
(`PACKED_CLMAD`) that does not transfer off Blackwell. Take g6/g6e and g7/g7e
quotes in the same `benchmarks/ada/prices.sh` run: quotes from different days
are not a comparison at a break-even this narrow.

**Region dominates instance type.** us-west-2 is about 25% cheaper than
us-east-1 for both parts, and that gap is larger than the g7-vs-g7e gap in
either region. At the recollected rate the fleet's present configuration
(4 × g7 in us-east-1, ~$2.74/hr) buys 3.15e13 iterations per dollar, against
4.19e13 for g7e in us-west-2 — about **1.33x**, or ~$490/month at equal
throughput. That conclusion holds whichever side of the break-even the 4500
lands on, because it is a statement about price, not about the GPU.

Two costs to check before moving region: the DP bucket
(`ecc2k130-<account>`) becomes cross-region, so price egress on upload; and
`aws/README.md` records g7e as offered in us-west-2, us-east-1 and us-east-2.

## Expected shape of the answer

The kernel issues 74.1 lane-instructions per SM-clock on the 6000, 95% of the
best rate that part's hardware probe ever measured
([THROUGHPUT-30B.md](THROUGHPUT-30B.md)). An issue-bound kernel scales with
SMs × clock, so the 4500's rate should be close to the 6000's scaled by their
SM counts and sustained clocks. A measurement far *below* that scaling is the
interesting outcome — it would mean something other than issue rate binds on
the smaller part, and the tuning that reached 74.1 lanes/SM-clock does not
transfer. A measurement at or above it means there is no missing performance to
chase on g7, and the choice is purely the price arithmetic above.

It landed on the scaling: 14.64 × (82 / 188) × (1.95 / 2.44 GHz) = 5.10 B/s
against 5.107 measured, so the kernel is issue-bound on the 4500 exactly as on
the 6000 and the part's only handicap is its 165 W power limit, which holds
the clock at 1.95 GHz where the 6000 sustains 2.42. There is no g7-specific
performance to chase; the decision is the price arithmetic, and it favours
g7e.
