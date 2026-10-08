# Frobenius steps through shared-memory nibble tables

`PACKED_SIGMA_TABLE=1` (`ECC_PACKED_SIGMA_TABLE`) computes the fused sigma
walk's `x + sigma^j(x)` and `y + sigma^j(y)` as one GF(2)-linear map on the
polynomial coordinates the state already uses, read out of a table in shared
memory, instead of converting `y` to the normal basis, running the paired
56-stage masked swap network on `x` and `y`, and converting both results
back. It is implemented for the fused schedule (`SIGMA_FUSED=1`) and changes
nothing about the iteration function, the distinguished-point rule, reports
or checkpoints: the map is the same field automorphism, only evaluated in
different coordinates. `make gpu-rtx-pro6000-sigma-table` builds it.

The headline below compares the table at 512 x 1 against the preset at
256 x 2 with four waves of workers; the population sweep further down
attributes +3.6% of that to the one-block geometry and +1.8% to the table,
and finds another +5.9% in the L2 persisting window at one wave; the
headline protocol on that recommended configuration gives **16.86 B/s**,
paired median **1.102** against the preset at its own best population.

Under the fused headline protocol on one RTX PRO 6000 (five A/A pairs, five
alternating A/B pairs, 64 launches per sample, replay and sorted-corpus gates),
the verdict is **promote as compatible engineering**: A/B paired median ratio
**1.023930** (minimum 1.023803, maximum 1.024086) against an A/A maximum drift
of 0.099%, session medians **15.884369** against 15.513653 B complete scalar
updates/s. An earlier 32-launch screen without an A/A panel measured a paired
median of 1.025148 over five pairs. Both clear the gate the fused tuning star
preregistered (paired geometric mean at least 1.015, every pair at least
1.005), which none of that star's ten arms met. The 26 B/s objective remains
unmet.

## Why this lever

[ROOFLINE.md](ROOFLINE.md) finds the integer logic pipe the busiest unit of
the packed kernel, with the carry-less unit second and the load/store pipe
far behind, so the opening is work that moves from logic to lookups. In the
fused select, the Frobenius step was the largest block of logic that was not
multiplication: the paired swap network (about 500 instructions for `x` and
`y` together, with 56 mask loads from the [shared copy](SHARED-SIGMA.md)),
one `fromPolynomial131` for `y` and two `toPolynomial131` after the network.
A nibble-table evaluation of the same map costs about 150 logic operations
and 37 shared-memory loads per coordinate, and needs no conversions at all.
The conversion of `x` to the normal basis stays, because the weight that
drives the distinguished-point rule and the jump index is a normal-basis
property; `y` is converted only when a report is written.

## Arithmetic

With `T` the normal-to-polynomial conversion and `sigma` the Frobenius
permutation of the normal basis, `L_j = T (I + sigma^j) T^-1` is a fixed
131 x 131 matrix over GF(2) for each `j` in 3..10. Writing the polynomial
coordinates as 33 nibbles, `L_j(x)` is the XOR of 33 table rows, one per
nibble position and nibble value. Rows hold output bits 0..127 as one
16-byte entry (8 maps x 33 positions x 16 values x 16 bytes = 67,584 bytes);
output bits 128..130 are the parities of the input against three 131-bit
masks per map (480 bytes). Each nibble costs one funnel shift, one `LOP3`
that masks the nibble into a row offset and ORs in the map's base address,
and one `LDS.128`; the rows XOR into four accumulators, and the three top
bits cost 15 `LOP3`, three `POPC` and a few combines.

The table is built on the host by `buildSigmaTable` in `packed131.h` from the
existing `fromPolynomial131`, `sigma131` and `toPolynomial131`, so it agrees
with the verified arithmetic by construction; the engine uploads it once and
every walk block stages it into dynamic shared memory before its first step
(the staging precedes the per-thread exit, so the barrier is uniform). At
68 KB per block the kernel runs one block per SM, so the measured
configuration uses 512-thread blocks with one resident block, which keeps the
preset's 128-register budget and 16 resident warps. The compiler drops the
now-unreferenced shared mask copy, so the candidate's static shared
footprint is zero and its dynamic footprint 68,352 bytes.

## Validation

* `make test-packed` and `make test-packed-network` check, in every host
  build, that the table reproduces `toPolynomial131(x + sigma^j(x))` for all
  131 polynomial basis vectors, dense and random inputs and all eight maps,
  and that the result converts back to the independent reference's
  `x + sigma^j(x)` (5,160 cases).
* `src/testpackedcuda.cu`, compiled with the knob, runs the device lookup
  through shared memory on 3,112 cases and compares it with both the host
  table and the network path; its `packed arithmetic sigma table: 1` marker
  identifies the build.
* `codegen/testpackedclient.py` runs the complete client with report replay,
  restarts, byte-identical resume and the overdue guard on the table binary.
* The job's replay gate re-walks 300 reports from each binary on the host
  (96,256 workers, 95 steps, seven launches, cutoff 48) and requires the two
  1,709,477-record corpora to be the same sorted multiset; they were,
  SHA-256 `7e0e8c90a7a1256c4f2c7f5ad4f23d5334de025a82590639f25516c0be25de92`.

## Measurement

[benchmarks/sigma-table/gpujob.sh](benchmarks/sigma-table/gpujob.sh) ran
through `modal_job.py` on one RTX PRO 6000 Blackwell Server Edition
(`GPU-d0d6311b-74da-213d-36ff-462e0b4e66d2`, driver 610.57.04, CUDA 13.3.73,
600 W limit). Both binaries use the fused preset's knobs
(`PRO6000_SIGMA_FUSED_KNOBS`: batch 16, native carry-less multiply, compact
state, weighted prefixes, shared masks, inline polynomial products); the
candidate adds `PACKED_SIGMA_TABLE=1` with `THREADS=512 MINBLOCKS=1`. Every
timed sample runs 385,024 workers, 6,160,384 walks, 1,024 steps and 32
launches: 201,863,462,912 complete scalar updates, `--verify 0` after the
gates above. Two warm-ups per binary are excluded.

| pair | order | control B/s | table B/s | ratio |
|---:|---|---:|---:|---:|
| 1 | control, table | 15.555107 | 15.946285 | 1.025148 |
| 2 | table, control | 15.511560 | 15.924801 | 1.026641 |
| 3 | control, table | 15.533519 | 15.904710 | 1.023896 |
| 4 | table, control | 15.497462 | 15.902108 | 1.026110 |
| 5 | control, table | 15.533439 | 15.897102 | 1.023412 |

Warm-ups were 15.640611 and 15.581936 (control) and 16.124771 and
16.051822 B/s (table). `nvidia-smi` sampled after each run shows the table
binary at 2,280–2,340 MHz against 2,362–2,407 MHz for the control, at
similar power (489–556 W against 509–564 W): the kernel is spending part of
its instruction saving on clock, which is consistent with the roughly 7.5%
the same change measured on the older software-multiplier preset (below).

| walk kernel, static | control | table |
|---|---:|---:|
| registers / local bytes / spills | 126 / 0 / 0 | 128 / 0 / 0 |
| shared bytes, static + dynamic | 1,792 + 0 | 0 + 68,352 |
| instructions | 5,264 | 5,104 |
| `LOP3` / `SHF` / `IMAD` | 2,410 / 1,174 / 584 | 2,359 / 1,043 / 494 |
| `LDS` / `CLMAD` | 56 / 79 | 140 / 79 |

### Headline protocol

[benchmarks/sigma-table/gpujob-headline.sh](benchmarks/sigma-table/gpujob-headline.sh)
(`make bench-rtx-pro6000-sigma-table`) is the fused headline job with the
table as candidate: the same frozen hardware and compiler checks, the GPU
arithmetic and client integration suites for both binaries, the 300-report
replay with sorted-corpus identity through the fused job's `corpus_identity`
tool, four excluded warm-ups, five alternating control/control pairs, five
alternating control/candidate pairs at `385024 x 16 x 1024 x 64` updates
each, and the native summarizer. It ran on the same GPU and driver as the
screen below.

| pair | A/A control_a | A/A control_b | A/B control | A/B table | A/B ratio |
|---:|---:|---:|---:|---:|---:|
| 1 | 15.498753 | 15.509355 | 15.511745 | 15.883911 | 1.023993 |
| 2 | 15.503252 | 15.518575 | 15.514757 | 15.884515 | 1.023833 |
| 3 | 15.507054 | 15.515976 | 15.513653 | 15.884899 | 1.023930 |
| 4 | 15.513678 | 15.512937 | 15.514261 | 15.883553 | 1.023803 |
| 5 | 15.516455 | 15.514619 | 15.510774 | 15.884369 | 1.024086 |

A/A maximum absolute drift 0.000988; A/B median 1.023930, minimum 1.023803;
decision `PROMOTE_ENGINEERING`, `goal_26b_met` false. Both replay corpora
held 1,709,477 sorted v1 records and the identity tool reported them
IDENTICAL. Receipts under
[benchmarks/sigma-table/headline/](benchmarks/sigma-table/headline/):
[result.json](benchmarks/sigma-table/headline/result.json),
[samples.tsv](benchmarks/sigma-table/headline/samples.tsv), the per-run logs,
build, arithmetic, integration, replay and corpus outputs, binary and source
hashes.

### Headline protocol, recommended configuration

The same job with `HEADLINE_PERSIST=1` (`make bench-rtx-pro6000-sigma-table`
style, candidate `PACKED_L2_PERSIST=1` at its one-wave population of 96,256
workers against the preset at its four-wave 385,024), one RTX PRO 6000
(driver 580.95.05): A/A maximum drift 0.053%, A/B paired ratios 1.102281,
1.102789, 1.101823, 1.103130 and 1.100996, **median 1.102281**, session
medians **16.864369** against 15.299614 B complete scalar updates/s, verdict
`PROMOTE_ENGINEERING`, identical 1,709,477-record corpora. Receipts under
[headline-persist/](benchmarks/sigma-table/headline-persist/). This is the
number of record for the table: **16.86 B/s, +10.2%** over the confirmed
fused preset measured in the same session.

### Screen

Receipts of the 32-launch screen: [result.json](benchmarks/sigma-table/result.json),
[samples.tsv](benchmarks/sigma-table/samples.tsv), the per-run logs, build,
arithmetic, integration, replay and corpus-identity outputs, binary and
source hashes, and [walk-static-mix.json](benchmarks/sigma-table/walk-static-mix.json)
(the 68 MB SASS dumps are not retained). The screen has no A/A panel and
uses 32 rather than 64 launches per sample; the headline protocol above is
the measurement of record.

### One-knob star on the table build

[benchmarks/sigma-table/gpujob-star.sh](benchmarks/sigma-table/gpujob-star.sh)
takes the table build (batch 16, 512 threads, one block per SM) as baseline
and changes one knob per arm. Every arm passes the arithmetic and
integration suites and reproduces the baseline's 1,709,477-record replay
corpus with the walk population held fixed. Three alternating pairs per arm
at 32 launches, one RTX PRO 6000 (driver 580.95.05, a different allocation
from the headline run):

| arm | baseline B/s | arm B/s | paired ratios | median | decision |
|---|---:|---:|---|---:|---|
| `BATCH=32` | 15.912–15.981 | 13.421–13.428 | 0.839832, 0.843342, 0.843873 | 0.843342 | reject |
| `PACKED_INV_POLY=2` | 15.903–15.912 | 15.979–15.984 | 1.004748, 1.004383, 1.004513 | 1.004513 | below gate |
| `PACKED_FROM_REDUCED=1` | 15.912–15.917 | 15.908–15.918 | 0.999410, 1.000312, 1.000275 | 1.000275 | reject |

Batch 32 doubles the per-thread state and loses a sixth of the rate, so the
preset's batch 16 stands. The polynomial-basis inversion is a consistent
small gain on top of the table, as it was in the fused star, but stays below
the 1.015 promotion gate; the five-word inverse conversion is flat. Receipts
under [benchmarks/sigma-table/star/](benchmarks/sigma-table/star/).

### Attribution, population and the persisting window

[benchmarks/sigma-table/gpujob-population.sh](benchmarks/sigma-table/gpujob-population.sh)
separates the three things the table build changed at once and adds the two
run-time levers ONE-BLOCK-GEOMETRY.md found for the table walk. Four builds,
all gates and identical corpora, three rotating rounds at 32 launches on one
RTX PRO 6000 (driver 580.95.05):

| build | workers | B/s (median of 3) |
|---|---:|---:|
| fused preset, 256 x 2 | 96,256 (one wave) | 15.109 |
| fused preset, 512 x 1 (no table) | 96,256 | 15.647 |
| table, 512 x 1 | 96,256 | 15.932 |
| table, 512 x 1, `PACKED_L2_PERSIST=1` | 96,256 | **16.877** |
| fused preset, 256 x 2 | 385,024 (four waves, 2 rounds) | 15.497 |
| table, 512 x 1 | 385,024 (2 rounds) | 15.887 |

So the one-block geometry alone is +3.6%, the table on top of it +1.8%, and
the persisting L2 window (80 MiB of the 104.7 MB one-wave blob, x/y/pchain
entirely inside it) a further +5.9%. The 256 x 2 preset is 2.6% faster at
four waves than at one, while the table build is indifferent to the
population, which is why the four-wave headline above credits the table with
only 2.4%. The recommended configuration is the table at 512 x 1 with the
persisting window at the automatic one-wave population: **16.88 B/s**, 8.9%
above the preset at its best population in the same session.

[benchmarks/sigma-table/gpujob-tagdenom.sh](benchmarks/sigma-table/gpujob-tagdenom.sh)
tried the table walk's other lever, `SIGMA_TAG_DENOM=1`: no denominator
field, the jump index tagged into the prefix's spare tail bits and
`x + sigma^j(x)` rebuilt from the table in the reverse pass. The walk stays
bit-identical (all gates, identical corpora) but the rate falls to **0.795**
of the table build, with or without the window: the rebuilt lookup lands on
the serial inverse chain, and the GPU sits at its maximum clock with 100 W
to spare, so this is latency, not work. Rejected; the knob stays for the
record. Receipts under [population/](benchmarks/sigma-table/population/) and
[tagdenom/](benchmarks/sigma-table/tagdenom/).

[benchmarks/sigma-table/gpujob-pipeslot.sh](benchmarks/sigma-table/gpujob-pipeslot.sh)
tried `SIGMA_PIPE_SLOT=1`: the fused reverse pass issues the next slot's
loads and its pair product (which depend only on the inverse chain) before
the current slot's square, second product, stores and selection, so the
carry-less unit has work during the ALU tail. Bit-identical walk, all gates,
128 registers with 8 bytes of local memory. Four rounds against the table
with the persisting window at one wave: **0.983893** (16.607 against 16.879
B/s, every round between 0.9836 and 0.9841). ptxas was evidently already
overlapping across the slot boundary as far as the registers allow, and the
explicit pipeline costs more in live state than it buys. Rejected; the knob
stays for the record. Receipts under [pipeslot/](benchmarks/sigma-table/pipeslot/).

### Where the update's time goes

`roofline.py` on both builds (receipts in `population/roofline-*.txt`), dynamic
lane-instructions per scalar update and the pipe ceilings at 2.415 GHz:

| pipe | fused preset | table build | ceiling (table) |
|---|---:|---:|---:|
| ALU (LOP3, SHF, IADD, ...) | 1,650 | 1,355 | 21.4 B/s |
| FMA-pipe integer | 236 | 143 | |
| LSU | 81 | 95 | 76 B/s |
| `CLMAD` | 38.1 | 38.1 | **19.3 B/s** at 1.62 lanes/SM-clk, 23.7 at 1.99 |

By function the table build spends 441 lanes in `reducePolynomial131`,
377 in `product131`, 335 in the two table lookups, 103 in the weight's
`fromPolynomial131`, 86 in the inversion's `toPolynomial131` and 38 in the
top-bit parities. The carry-less unit is the busiest pipe: 82% at the
four-wave 15.88 B/s and 87% at 16.88, against ONE-BLOCK-GEOMETRY.md's
1.62 lane-CLMADs per SM-clock. The sigma walk issues 38.1 `CLMAD`s per
update, five more than the table walk's 33.1, because its generated reducer
spends one `CLMAD` per product on the quotient; at that count the unit's
floor is 19.3 B/s, and 22 B/s needs the count back at 33 with the ALU no
heavier than today, which is what a low-weight pentanomial basis would buy
(THROUGHPUT-30B.md lever 1, now that the basis machinery is a table).

### Earlier preset

The same change was first measured on the CUDA 13.0 software-multiplier
preset of the `cursor/ic-boundary-experiments-d111` branch (batch 32, 256
threads, two blocks per SM, 6.9 B/s): three alternating pairs gave a median
**+7.532%** (6.742071 to 7.249909 B/s), batch 64 on top of the table added
0.17%, and batch 64 alone lost 0.86%. Receipt:
[legacy-preset/result.json](benchmarks/sigma-table/legacy-preset/result.json)
and the `compare.py` that produced it (it targets that branch's knobs, not
this tree's).
