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

Measured on one RTX PRO 6000 against the confirmed fused preset, five of
five alternating pairs favour the table, with a **paired median ratio of
1.025148** (minimum 1.023412, maximum 1.026641): session medians
**15.904710** against 15.533439 B complete scalar updates/s. That clears the
gate the fused tuning star preregistered (paired geometric mean at least
1.015, every pair at least 1.005), which none of that star's ten arms met.
The 26 B/s objective remains unmet.

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

Receipts: [result.json](benchmarks/sigma-table/result.json),
[samples.tsv](benchmarks/sigma-table/samples.tsv), the per-run logs, build,
arithmetic, integration, replay and corpus-identity outputs, binary and
source hashes, and [walk-static-mix.json](benchmarks/sigma-table/walk-static-mix.json)
(the 68 MB SASS dumps are not retained). This run has no A/A panel and uses
32 rather than 64 launches per sample; the headline protocol in
[benchmarks/sigma-fused/HEADLINE-PROTOCOL.md](benchmarks/sigma-fused/HEADLINE-PROTOCOL.md)
should be run before the number replaces the confirmed 15.436677 B/s.

### Earlier preset

The same change was first measured on the CUDA 13.0 software-multiplier
preset of the `cursor/ic-boundary-experiments-d111` branch (batch 32, 256
threads, two blocks per SM, 6.9 B/s): three alternating pairs gave a median
**+7.532%** (6.742071 to 7.249909 B/s), batch 64 on top of the table added
0.17%, and batch 64 alone lost 0.86%. Receipt:
[legacy-preset/result.json](benchmarks/sigma-table/legacy-preset/result.json)
and the `compare.py` that produced it (it targets that branch's knobs, not
this tree's).
