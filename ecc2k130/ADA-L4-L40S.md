# Ada (EC2 g6 / g6e) — what transfers off sm_120, and what has to be measured

Every throughput figure in this tree is for the **RTX PRO 6000 Blackwell Server
Edition**, `sm_120`: 14.637530 B/s benchmarking and 14.1 B/s collecting
([THROUGHPUT-30B.md](THROUGHPUT-30B.md), [aws/README.md](aws/README.md)). EC2
**g6** carries the **L4** and **g6e** the **L40S**, both `sm_89` (Ada). This
tree records no rate for either, and until this change it could not have
produced a correct client for them at all.

That is two separate gaps, and only the first is a bug.

## 1. The fleet could not build an Ada client

`aws/build.sh` baked `PACKED_CLMAD=1` into the knob string unconditionally,
while `aws/README.md` invites `ARCHES="89 90 120"` "if a non-Blackwell instance
type joins later". Those two instructions are in conflict, and the conflict is
silent: the knob is not a portability flag, it is a *tuning* decision bought
with a pipe balance nobody has measured outside Blackwell.

`clmad` costs about 37.7 ALU operations' worth of work per carryless operation
— that price is measured, in `benchmarks/clmad-price`. It is still the right
trade on the RTX PRO 6000 for one reason, recorded in the comment in
[include/packed131.h](include/packed131.h): Nsight Compute put the walk at
**87.3% on the ALU pipe and 51.4% on the pipe that carries `clmad`**. The walk
is short of ALU there and has carryless capacity to spend, so moving work onto
the dearer unit wins. Change the ratio of units and that argument does not
survive; it does not even point in a known direction.

So `CLMAD` no longer travels with `ARCHES`. `aws/build.sh` defaults it to 1
only for a Blackwell-only build and to 0 as soon as anything older is in
`ARCHES`, and `CLMAD=1` is available to override once a receipt justifies it.

The same change closed a second hazard it created. The published S3 key was
`bin/<source-sha>/`, and `ARCHES` and `CLMAD` now vary independently of the
source: two builds differing only in those would have overwritten each other,
leaving `campaign.json` pointing at one binary and its manifest describing the
other. The key is now bound to the source, the architectures **and** the knobs,
and the manifest carries that `buildSha256` alongside the source sha.

## 2. Nothing has measured the Ada rate

```bash
bash benchmarks/ada/run.sh          # on the g6/g6e instance, or: make bench-ada
```

Writes `benchmarks/ada/result.json` with each command's verbatim output,
following [benchmarks/rtx-pro4500/result.json](benchmarks/rtx-pro4500/run.sh).
It reads the rate with the tree's own `codegen.benchreport` — `parseRate`,
which accepts only the single synchronised `finished:` line and rejects a run
without one, and `summarizeSamples`, which takes the median across repeats.
That is the method behind the RTX PRO 6000 figures, so the two are comparable
by construction rather than by intention.

It **measures both arms of the CLMAD question** on the same device allocation
rather than taking either side of §1 on trust, and reports the ratio and which
preset the receipt selects. An arm whose build fails is recorded as a failed
arm — `clmad` is PTX 9.3 / CUDA 13.3, and a toolkit or target that rejects it
fails at compile time, which is a result — and the other arm still yields a
rate.

Three refusals rather than a quietly wrong number: it will not run while
another process holds the GPU (a benchmark and a collecting worker steal SMs
from each other and neither figure measures anything), it will not record an
L4/L40S result on a part that is not one, and it never turns an invalid run
into a number. It builds `sm_89` only; the Makefile default is a
five-architecture fat binary and ptxas takes minutes per architecture on this
kernel.

`make bench-g6-modal` and `make bench-g6e-modal` run the same preset through
Modal (`ECC_GPU=L4` / `L40S`). `make gpu-g6` / `make gpu-g6e` emit the thin
`sm_89` software-product client. The T4 / g4dn sibling is
[T4-G4DN.md](T4-G4DN.md). `make bench-ada-modal` is the same L40S entry
(`ADA_GPU=L4` selects the L4). Workers stay **automatic** in both:
the RTX PRO 6000 preset's 385,024 is four waves on a 188-SM part and has no
claim on a 142-SM (L40S) or 58-SM (L4) one, and `autoThreads`
([include/packedengine.cuh](include/packedengine.cuh)) scales with
`multiProcessorCount` and with the occupancy the build actually admits.

Everything else in the preset — batch 16, 256-thread blocks, two resident
blocks, 256-worker state tiles, compact storage, weighted prefixes, shared
Frobenius masks — is a field-arithmetic or storage choice rather than a GPU
choice, and is carried over unchanged, exactly as `benchmarks/rtx-pro4500`
carries it. The geometry survives the move on a register argument: the kernel
uses 122 registers per thread ([NATIVE-CARRYLESS.md](NATIVE-CARRYLESS.md)),
and 2 × 256 threads/SM against Ada's 64K registers/SM leaves 128 per thread, so
the two-block residency the preset asks for is admissible on `sm_89` as well.
The shared Frobenius masks are 448 words, 1,792 bytes per block, far inside
Ada's per-block shared limit.

## What the measurement is expected to show

Stated as a prediction so that the receipt can contradict it.

**H1.** The walk is issue-bound on integer logic — it is at 74.1
lane-instructions per SM-clock on the 6000, 95% of the best rate that part's
hardware probe ever measured ([THROUGHPUT-30B.md](THROUGHPUT-30B.md)).

**H2.** An Ada SM retires integer logic at half the per-clock rate of an
`sm_120` SM, Ada splitting its lanes between an FP32-only half and an
FP32/INT32 half where Blackwell's are unified.

Under H1 and H2 together, rate scales as SMs × clock × per-SM issue width, and
the 6000's 14.637530 B/s over 188 SMs — 77.9 M it/s per SM — gives, at equal
clock:

| part | SMs | scaled expectation |
|---|---:|---:|
| L40S (g6e) | 142 | ~5.5 B/s |
| L4 (g6) | 58 | ~2.3 B/s, **before** power |

Those are arithmetic from H1 and H2, not measurements, and neither hypothesis
is established in this tree — H2 in particular is a vendor-architecture claim
this repository has never tested. The L4 figure carries a further caveat that
the L40S figure does not: the L4 is a **72 W** part and the L40S a 350 W one,
so the L4's sustained clock is a much larger discount on its boost clock, and
the scaling above assumes equal clocks. Expect the L4 to land below 2.3 B/s;
how far below is the question.

A measurement far below the scaling on the **L40S** is the interesting outcome
— it would mean something other than issue rate binds on Ada, and the tuning
that reached 74.1 lanes/SM-clock does not transfer. A measurement at or above
it means there is no missing performance to chase on g6e and the choice is
purely the price arithmetic below.

## The decision this feeds

Cost is a property of the work, so the only figure that ranks instance types is
**iterations per dollar**:

```
iterations per dollar  =  rate (it/s) × 3600 / ($/GPU-hour)
```

g6/g6e beats the g7/g7e the campaign runs on today exactly when

```
ada_rate  >  blackwell_rate × ada_price / blackwell_price
```

using the *collecting* rate on both sides, not `--bench` — `--bench` suppresses
distinguished-point handling and the collecting rate is lower (14.1 vs
14.637530 B/s on the 6000). Fetch a matched set of spot prices with:

```bash
bash benchmarks/ada/prices.sh            # g6, g6e, g7, g7e; cheapest AZ per type
```

| | value |
|---|---|
| L4 (g6), `--bench` | *pending — fill from `benchmarks/ada/result.json`* |
| L4 (g6), collecting | *pending* |
| L40S (g6e), `--bench` | *pending* |
| L40S (g6e), collecting | *pending* |
| RTX PRO 4500 (g7) | *pending — see [RTX-PRO4500.md](RTX-PRO4500.md)* |
| RTX PRO 6000 (g7e), `--bench` | 14.637530 B/s |
| RTX PRO 6000 (g7e), collecting | 14.1 B/s |

No price table is printed here because none has been observed for g6/g6e in
this tree, and [RTX-PRO4500.md](RTX-PRO4500.md) already records what happens
when a rate is carried on recollection: the g7-vs-g7e break-even sits within a
few percent and only a receipt settles it. The same discipline applies here,
with one more reason — an Ada quote and a Blackwell quote taken on different
days are not a comparison, so take all four in one `prices.sh` run.

**Region still dominates instance type.** That conclusion from
[RTX-PRO4500.md](RTX-PRO4500.md) is a statement about price, not about the GPU,
and it holds whatever the Ada rates turn out to be: us-west-2 ran about 25%
cheaper than us-east-1 for both Blackwell parts, a wider gap than the one
between the parts in either region. Before moving region, price egress on the
DP bucket (`ecc2k130-<account>`), which becomes cross-region.

## If the receipt selects g6/g6e

Build and publish an Ada client with:

```bash
BUCKET=... ARCHES="89" ./aws/build.sh /path/to/ecc2k130      # CLMAD defaults to 0
BUCKET=... ARCHES="89" CLMAD=1 ./aws/build.sh /path/...      # only with a receipt
```

The binary key now includes the knobs, so an Ada client and a Blackwell client
coexist in the bucket instead of overwriting one another. `bootstrap.sh` reads
the expected carryless marker off that build's `manifest.json` instead of
requiring 1, so a CLMAD-free Ada client boots; a binary that disagrees with its
own manifest is still refused, and a missing manifest still demands 1. `campaign.json` still
carries one `binaryKey` and one geometry, so a **mixed-architecture** fleet is a
further change and not merely a further build — the batch/threads/minBlocks
assertion in `build.sh` is what would catch an attempt to serve both from one
entry, and it should keep catching it until the campaign record grows a
per-architecture binary map.
