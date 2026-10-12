# The RTX PRO 4500: a power-bound part, priced from the RTX PRO 6000

Question: can the table-walk kernel be made faster on the **RTX PRO 4500
Blackwell** (EC2 `g7.2xlarge`, the fleet's card), which the tree has measured
twice (shipping walk 5.107 B/s, table walk 5.552 B/s, both at the old
256 × 2 geometry; [RTX-PRO4500.md](RTX-PRO4500.md),
[ITERATION-FUNCTION.md](ITERATION-FUNCTION.md) §6.3) and never with any of
the work since?

**What could not be done in this session:** run anything on a 4500. EC2
`RunInstances` returns `Blocked` in every region while AWS Health carries
two open risk events on the account (`AWS_RISK_CREDENTIALS_EXPOSURE_SUSPECTED`,
`AWS_RISK_ACCOUNT_CONSOLE_COMPROMISE`); Modal rents no RTX PRO 4500; there is
no RunPod key. Every 4500 figure below that is not one of the two receipts
above is a **prediction**, marked as such, from measurements on an RTX PRO
6000 and the model of §2. The command that turns each prediction into a
measurement is in §5.

## 1. Why the 6000 can price the 4500

The two parts have the same SM (`sm_120`): ITERATION-FUNCTION.md §4.5
measured every instruction stream within 1% per SM-clock, and the per-clock
table/shipping ratio was 1.164 on the 4500 against 1.161 on the 6000. What
differs is the envelope: 82 SMs at a **165 W** limit (2.01 W per SM) against
188 at 600 W (3.19 W per SM). Both 4500 receipts sat on that limit at
1.82 – 1.97 GHz, and the kernel that did more per clock was clocked 6.6%
lower by the card, so the fixed-power gain (+8.7%) was half the per-clock
gain (+16.4%). **On the 4500 the unit that ranks builds is updates per
joule, not updates per SM-clock.** The 6000 does not reach its limit on
most of these kernels (560 W of 600, 2422 MHz held), so on it the two
units separate and both can be read.

The container cannot cap the 6000's power to the 4500's per-SM budget
(`nvidia-smi -pl 378`: "Insufficient Permissions"), so the prediction is
energy-based rather than an emulation.

## 2. The model and its calibration

Energy per update is measured on the 6000 with the board power sampled
every 100 ms *during* each bench (`--steps 1024 --launches 64`, automatic
workers; the loaded window is the samples at ≥ half the trace's peak, less
the first five), divided by the `finished:` rate.
[benchmarks/power-bound/gpujob.sh](benchmarks/power-bound/gpujob.sh),
[summarize.py](benchmarks/power-bound/summarize.py).

Prediction: at a fixed power cap two builds' rates stand in the ratio of
their updates per joule. It is first-order — it ignores that the 4500 runs
the lower-voltage end of the same V/f curve, and that static power is a
larger share of 165 W than of 600 W — so it is calibrated against the one
pair the tree has on both cards:

| pair | 6000 updates/J ratio (this note) | 4500 measured rate ratio at 165 W | model error |
|---|---:|---:|---:|
| table walk / shipping walk, both 256 × 2 | 28.44 / 25.50 = **1.115** | 5.552 / 5.107 = **1.087** | **+2.6%** (optimistic) |

The error bar every prediction below carries is that 2.6%, and the two
anchors (scaling from the 4500's shipping receipt or from its table
receipt) give the ends of each predicted range.

## 3. Survey on the 6000

One RTX PRO 6000 (Modal, driver 580.95.05, nvcc 13.3.73, 2026-09-22), three
alternating repetitions per binary, medians;
[rtx-pro-6000/summary.json](benchmarks/power-bound/rtx-pro-6000/summary.json).
Every binary re-walked 300 of 300 device reports with 0 dropped except
640 × 1, whose automatic 1.93 M walks overflowed the 262,144-record report
buffer on the first chunk (the guard of ITERATION-FUNCTION.md §6.2 firing;
the same build verified 300/300 with the reference's distinguished points
at a forced walk count in AUTOSWEEP.md §4).

| build | B/s on the 6000 | board W | SM MHz | **M updates / J** | / 20 B/s build | 4500 predicted, B/s (from shipping – table anchor) |
|---|---:|---:|---:|---:|---:|---|
| shipping walk, 256 × 2 (the 4500's 5.107 receipt) | 14.467 | 567 | 2422 | 25.50 | 0.717 | *5.107 measured* |
| table walk, 256 × 2 (the 4500's 5.552 receipt) | 16.498 | 582 | 2383 | 28.44 | 0.800 | *5.552 measured* |
| **20 B/s build**, `gpu-rtx-pro6000-20b`, 512 × 1 | 20.095 | 565 | 2422 | **35.57** | 1 | **6.94 – 7.12** |
| 20 B/s build at 384 × 1 | 18.533 | 509 | 2422 | 36.42 | 1.024 | 7.11 – 7.29 |
| 20 B/s build at 640 × 1 | 20.127 | 592 | 2398 | 34.11 | 0.959 | 6.66 – 6.83 |
| **20 B/s build + `TABLE_FUSED`** | 20.151 | **529** | 2422 | **38.08** | **1.071** | **7.43 – 7.63** |
| 20 B/s build + `PACKED_ALU_SQR` | 20.031 | 583 | 2422 | 34.35 | 0.966 | 6.70 – 6.88 |
| two chains, 256 × 32 | 13.895 | 429 | 2422 | 32.44 | 0.912 | 6.33 – 6.50 |

Reading it:

- **The 20 B/s build should be worth +25 – 28% on a 4500** over its best
  receipt (6.94 – 7.12 against 5.552; +21 – 24% from §4's session, where
  the same binary drew 3% more power), because the one-block geometry and
  the scheduling of ONE-BLOCK-GEOMETRY.md buy energy per update as well as
  rate: 35.6 against 28.4 M updates/J. That is a prediction for a card the
  fleet runs, from work that was only ever measured on the other card. It
  moves the 4500 above RTX-PRO4500.md's per-dollar break-even in us-east-1
  (6.17 B/s) and us-west-2 (6.63).
- **Rate and energy rank the variants differently.** 640 × 1 is the fastest
  build on a cool 6000 and the worst of the one-chain builds per joule — which
  is why the warm 6000 took it back in clock in AUTOSWEEP.md §4. 384 × 1 is
  8% slower on the 6000 and 2.4% better per joule. And `TABLE_FUSED`, neutral
  on the 6000's rate (+0.3%), draws 36 W less for it — x and y are read once
  per update instead of twice (190 → 122 bytes; ONE-BLOCK-GEOMETRY §5), and
  state traffic is energy even when it is not time — so it is **+7.1% per
  joule**, the one lever here that the 6000's rate could not see and a
  power-capped card would.
- The two-chain kernel draws the least power and does the least per joule:
  idle datapaths still cost static power.

## 4. The fused pass, measured for the 4500

**Target, declared before the run:** `TABLE_FUSED` is the 4500 build iff,
alternating with the 20 B/s build on one 6000, its updates per joule are
**≥ 1.03×** in the median with every paired repetition above 1.0, its rate
is not lower than 0.99×, 300/300 re-walked, and its distinguished-point set
on the forced 1,540,096 walks is byte-identical to the reference's. Three
more fused forms are measured beside it (its load pipeline, 384 × 1, the
squaring on the carry-less unit) because the pass changes the energy mix
and the geometry optimum may move with it.

One RTX PRO 6000 (Modal, 2026-09-22), five alternating repetitions,
[rtx-pro-6000-fused/summary.json](benchmarks/power-bound/rtx-pro-6000-fused/summary.json).
All five binaries re-walked 300 of 300 with 0 dropped, and all five produced
the same 1,480,482 distinguished points on the forced 1,540,096 walks (hash
`6cb064cd…`, the same as every RTX PRO 6000 binary since TWO-CHAINS.md).

| build (20 B/s build + …) | B/s on the 6000 | rate paired / ref | board W | M updates / J | **per joule, paired / ref (min – max)** | 4500 predicted, B/s | class |
|---|---:|---:|---:|---:|---|---|---|
| — (reference, `gpu-rtx-pro6000-20b`) | 20.045 | 1 | 583 | 34.37 | 1 | 6.71 – 6.88 | reference |
| **`TABLE_FUSED=1`** = **`make gpu-rtx-pro4500`** | 20.124 | 1.004 (≥ 1.003) | 544 | 37.01 | **1.072 (1.063 – 1.111)** | **7.23 – 7.41** | **engineering; meets the target; the 4500 build** |
| fused + load pipeline | 19.693 | 0.982 | 542 | 36.29 | 1.067 (1.013 – 1.076) | 7.08 – 7.27 | rate below 0.99: fails the target |
| fused at 384 × 1 | 18.789 | 0.938 | 498 | 37.72 | 1.095 (1.083 – 1.123) | 7.36 – 7.55 | rate below 0.99: fails the target |
| fused + squaring on the carry-less unit (`PACKED_ALU_SQUARE=0`) | 19.097 | 0.953 | 497 | **38.42** | **1.119 (1.105 – 1.135)** | **7.50 – 7.70** | fails the rate clause; **the second arm for the 4500 itself** |

The reference drew 583 W in this session against 565 W in §3's, on another
6000 on another hour, so its updates per joule moved 3.4% between sessions;
the paired column, taken within one session and one repetition, is what
the target reads, and the predictions use this session's figures.

**Against the target:** `TABLE_FUSED` is 1.072× per joule in the median,
every repetition above 1.06, at 1.004× the rate, 300/300, identical points —
it is the 4500 build, `make gpu-rtx-pro4500`, and on the 6000 itself it is a
free +0.4% (within the noise of that card's rate, as ONE-BLOCK-GEOMETRY §5
found). Predicted on a 4500: **7.2 – 7.4 B/s, +30 – 33% over the 5.552
receipt and +42 – 45% over the fleet's shipping-walk 5.107**, with the
model's +2.6% optimism already inside the range's low end.

The rate clause was there so the 4500 build would not be slower on a card
that is not power-bound, and it is what rules the last row out as the
*declared* pick. But that row is the most efficient build measured, and why
is worth recording: the λ² squaring as five `CLMAD`s instead of ~60 logic
ops is 5% slower on the 6000's saturated carry-less unit and 4% cheaper in
energy — the unit is slow on this die, not expensive. On a part that is
power-bound rather than unit-bound the model puts it 4% above the fused
build, which is inside the model's error bar, so it is the second arm of the
4500 run (`RTX4500_EXTRA=PACKED_ALU_SQUARE=0`), not a claim.

## 5. The measurement owed

On one idle `g7.2xlarge` (RTX PRO 4500), once EC2 launches work again:

```sh
cd ecc2k130
python3 aws/bench_job.py --job benchmarks/power-bound/gpujob.sh --out /tmp/rtx4500 \
    --region us-west-2 --instance-type g7.2xlarge --env SET=rtx4500,REPS=5
python3 benchmarks/power-bound/summarize.py /tmp/rtx4500 benchmarks/power-bound/rtx-pro-4500
```

It builds the two receipts' builds as controls (shipping and table walk at
256 × 2), the 20 B/s build, `gpu-rtx-pro4500` and its `PACKED_ALU_SQUARE=0`
arm; verifies each (300/300; the four table-walk builds held to one
distinguished-point set at the forced `82 × 512 × 16` walks); and benches
them alternating with power sampled during each run. What it decides:

- whether the prediction holds — the 4500's own rate ratio of
  `gpu-rtx-pro4500` to the 5.552 control against the predicted 1.30 – 1.33;
  a miss beyond the 2.6% calibration error is a failure of §2's model, to be
  recorded as such;
- which of the two fused arms is the 4500's build;
- whether the 4500 now beats the RTX PRO 6000 per dollar: at the
  predicted 7.3 B/s it clears RTX-PRO4500.md's break-evens in us-east-1
  (6.17) and us-west-2 (6.63), but those were computed against the 6000's
  *shipping* rate, and the 6000 has the same levers — the per-dollar
  comparison has to be redone with both cards on their best builds.

## 6. Classification

| change | class | evidence |
|---|---|---|
| energy per update as the ranking unit for a power-capped part | method; calibrated | model +2.6% on the one pair measured on both cards |
| `TABLE_FUSED` for the 4500 (`make gpu-rtx-pro4500`) | **engineering, predicted +30 – 33% on a 4500; +7.2% per joule measured on the 6000** | §4; the rate on a 4500 is a prediction until §5 runs |
| fused + `PACKED_ALU_SQUARE=0` | engineering, second arm | +11.9% per joule, −4.7% rate on the 6000 |
| 640 × 1, two chains, `ALU_SQR` | worse per joule on the 6000 | §3 |
