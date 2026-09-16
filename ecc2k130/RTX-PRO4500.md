# RTX PRO 4500 (EC2 g7) — rate pending, and the decision that waits on it

Every throughput figure in this tree is for the **RTX PRO 6000**: 14.637530 B/s
benchmarking and 14.1 B/s collecting ([THROUGHPUT-30B.md](THROUGHPUT-30B.md),
[aws/README.md](aws/README.md)). The fleet currently running the campaign is on
**g7.2xlarge**, which carries the **RTX PRO 4500** — and this tree records no
rate for that part at all.

That gap is not cosmetic. It decides which instance type the campaign should
buy, and the recollected figure sits within a few percent of the break-even,
on the wrong side of it.

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

**The measured rate is not yet recorded.** A figure of roughly 6 B/s has been
stated from recollection; it appears in no receipt in this tree and is not
treated as evidence here. If it is right, g7e is marginally the better buy in
every region. If the part reaches 6.7 B/s, g7 wins in us-west-2. The margin is
small enough in both directions that only a receipt settles it.

| | value |
|---|---|
| RTX PRO 4500, `--bench` | *pending — fill from `benchmarks/rtx-pro4500/result.json`* |
| RTX PRO 4500, collecting | *pending* |
| RTX PRO 6000, `--bench` | 14.637530 B/s |
| RTX PRO 6000, collecting | 14.1 B/s |

## What is already settled

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
