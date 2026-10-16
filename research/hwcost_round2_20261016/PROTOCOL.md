# Protocol: `hwcost` round 2: ISA, `fp64` branches, GPU pipes, IC in time and dollars

Frozen before any run or analysis below; the commit that adds this file is
the protocol commit. Round 1 (`research/hwcost_model_20261016/`) built and
validated the model; this round follows its four "next steps".

## Scope and boundary

Four parts, one study. Nothing here solves a discrete logarithm or changes an
index-calculus pipeline, so there is no index-calculus improvement and no new
gap-to-rho row (AGENTS.md success criterion). Part D compares IC and rho in
time and dollars from frozen records and ranks them only where wall time and
every §5a band weight agree (AGENTS.md §5a).

Host for A and B: this session's cloud container, as in round 1
(Intel Xeon @ 2.10 GHz, family 6 model 207, 4 vCPUs, KVM; isolation L0).
**No hardware performance counters on this host**: the guest kernel reports
`Performance Events: unsupported CPU family 6 model 207 no PMU driver,
software events only.` (dmesg), so branch misses cannot be counted; part B
tests the branch hypothesis causally instead.

## A. ISA: a `-C target-cpu=native` build

- Builds at the protocol commit: **baseline** `cargo build --release --bin
  hwcost` (x86-64 baseline target) and **native**
  `CARGO_TARGET_DIR=target-native RUSTFLAGS="-C target-cpu=native" cargo build
  --release --bin hwcost` (this CPU: BMI2, ADX, AVX2, AVX-512 incl. IFMA).
  Each binary records the target features it was compiled with.
- Runs, back to back in this order, each the full `hwcost measure` (round-1
  parameters: 7 reps, 100 ms, batches 2^4 … 2^22, fields below):
  **A1** baseline, **B1** native, **A2** baseline, **B2** native →
  `runs/isa-base-1`, `runs/isa-native-1`, `runs/isa-base-2`,
  `runs/isa-native-2`. Each is fitted with `hwcost fit --fit-max-batch 1024`.
- Fields: `fp64`, `fp64bf` (part B), `w1`–`w4` (as round 1).
- `Fp52` (the AVX-512 IFMA backend) is not an arm: its scalar arithmetic is
  `Fp64`'s element for element, and IFMA serves only its F4 row accumulator,
  which the rho iteration does not use.

Per arm (every op arm `field.op.mode`, kernel arm `field × B`, scaling arm,
memory and clock arm), with medians over each run's 7 reps: `r_k = B_k / A_k`
(k = 1, 2; inverted for bandwidth and clock arms, where higher is better)
and the same-build drift `d = |A2 / A1 − 1|` (`hwcost compare`).

| id | claim | rule |
| --- | --- | --- |
| I1 | classification of every arm | **native faster** if `r1, r2 < 1 − max(2d, 0.02)`; **native slower** if `r1, r2 > 1 + max(2d, 0.02)`; else **no detected effect**. Reported per arm and counted per class; no threshold. |
| I2 | model C on the native build | V1 (round 1's calibrated-model claim, median \|error\| on B ≥ 4096) ≤ 0.30 for every field, in B1 and in B2 |
| I3 | calibration transfer | native kernel predicted from native op costs + A1's fitted residuals (`hwcost fit --residual-from runs/isa-base-1`): median \|error\| over test batches ≤ 0.30 for every field, in B1 and B2 |

Round 1's claims (C0, V1–V3) are re-evaluated on A1 and A2 as a
replication and reported whatever they say.

## B. The `fp64` branch hypothesis

Round 1 found the `fp64` kernel 1.6× slower than `w1`'s although its
operations are cheaper in isolation, and named data-dependent branches in
`Fp64::sub` and `Fp64::redc` as an unverified lead. `fp64bf`
(`hwcost::kernel::Fp64Bf`) is the same field, the same Montgomery
representation and bit-for-bit the same results, with the reduction's final
conditional subtraction and the subtraction's wrap-around written as masks.
Arms `fp64` and `fp64bf` are interleaved inside every run of part A.

| id | claim | rule |
| --- | --- | --- |
| F0 | gate | every run checks that the `fp64bf` walk equals the `Fp64` walk point for point (`check_same_walk`), and the unit test compares 10,000 random operations; a mismatch aborts the run |
| F1 | the branches cost the time | `fp64bf` kernel median ≤ 0.85 × `fp64`'s at **every** B ≥ 4096, in **both** A1 and A2 |
| F2 | the a priori model then explains `fp64bf` | V2 (round 1's a priori claim) on `fp64bf` ≤ 0.30, in A1 and A2 |
| F3 | static evidence (descriptive) | conditional jumps counted in the disassembly (`objdump -d`) of the timed kernel closure for `fp64` and for `fp64bf` in the baseline binary |

F1 passing supports the lead (only the branches differ); F1 failing rejects
it. Either way the outcome is reported.

## C. A GPU profile

Data: `ecc2k130/benchmarks/roofline/calibration.json`, 46 builds (21 on the
RTX PRO 6000, 25 on the B200), each with dynamic lane-instruction counts per
update (ALU, FMA, CLMAD, issue, LSU) and a measured rate. Profiles
(`hwcost::gpu`, every constant with its source): SMs as cores, measured
pipe ceilings (lane-ops per SM-clock) from `hardware-limits` and the
fast-clmad probe, issue 128 (architectural, flagged), measured clocks, Modal
list prices. The update's pipes overlap: SM-clocks per update = the busiest
pipe's `count / ceiling` (`Phase::pipes_overlap`, `ecc2k130/roofline.py`'s
model). Two pipe models: **separate** (ALU and FMA apart; primary) and
**shared-int** (ALU+FMA on one measured 69.3-lane mixed datapath;
descriptive). Run: `hwcost gpu-check --fit-arm base --reference-gpu
rtx-pro-6000 --check-rate rtx-pro-6000:14.637530 --check-rate
rtx-pro-6000:15.115792` → `runs/gpu/`.

| id | claim (separate pipes) | rule |
| --- | --- | --- |
| G1 | a priori, no fitted parameter, all 46 builds | median \|error\| ≤ 0.30; also reported: builds with η = measured / predicted > 1 (the ceiling model refuted there) |
| G2 | η fitted on each GPU's `base` build, tested on that GPU's other builds | median \|error\| ≤ 0.30, and whether it is ≤ 0.10 |
| G3 | η of the RTX PRO 6000 `base` build carried to every B200 build | median \|error\| ≤ 0.30 |

Disclosure: the repository's `calibration.txt` already lists each RTX PRO
6000 build's η under similar models; the thresholds above are round 1's 30%,
not tuned to it. The expected ECC2K-130 cost at the shipping σ-walk's
measured rates (14.637530 and 15.115792 B/s, $3.031/h, 2^60.9 iterations) is
printed as a check of the cost arithmetic against the paper's $124k
(`paper/ecc2k130/sec-solvetime.tex`).

## D. IC priced in time and dollars

Data: `research/ic_frontier_campaign_20261015` Phase 1 (32/32 verified runs
per arm per curve, host Xeon 2.30 GHz 4 vCPU, L0; its §5a band files).
Pairs: **Koblitz** `ic-compact` vs `rho-frob` (8 curves, log2 r 11–39);
**prime** `ic-pair` vs `rho-neg` (6 curves, log2 r 20–40). Run:
`hwcost ic-price` → `runs/ic-koblitz/`, `runs/ic-prime/`.

Per curve: `P(w)` from the band (`mean_s[w] · √r`), the mean measured wall
time per solve (counted, verified runs), `t(w) = wall_rho / P_rho(w)` and the
band-predicted IC wall `P_IC(w) · t(w)`, the implied weight `w*`, and dollars
per solve at one vCPU of the illustrative c7i.xlarge ($0.1785/h ÷ 4).

| id | quantity (descriptive; no threshold) |
| --- | --- |
| D1 | per curve: is the measured IC/rho wall ratio inside the band's ratio range over w ∈ {0, 0.1, 1}? and `w*` |
| D2 | `log2 y = a + α log2 r` fitted per arm for wall time and each `P(w)` over the tier — Koblitz log2 r ∈ [27, 40] (4 curves), prime [23, 41] (5 curves) — and projected to log2 r ∈ {48, 64, 80, 96, 112, 128}: core-years, dollars, IC/rho ratio, and a rank only where wall time and all three weights agree |

Projections are extrapolations of fitted exponents, labelled as such, not
measurements.

## Commands

```sh
cargo build --release --bin hwcost
CARGO_TARGET_DIR=target-native RUSTFLAGS="-C target-cpu=native" cargo build --release --bin hwcost
D=research/hwcost_round2_20261016/runs
for r in base-1:target native-1:target-native base-2:target native-2:target-native; do
  n=${r%%:*}; t=${r##*:}; mkdir -p $D/isa-$n
  ./$t/release/hwcost measure --out $D/isa-$n 2> $D/isa-$n/stderr.txt
done
for n in base-1 native-1 base-2 native-2; do ./target/release/hwcost fit --run $D/isa-$n > $D/isa-$n/fit.txt; done
mkdir -p $D/compare; ./target/release/hwcost compare --base $D/isa-base-1,$D/isa-base-2 --cand $D/isa-native-1,$D/isa-native-2 --out $D/compare > $D/compare/summary.txt
for n in native-1 native-2; do mkdir -p $D/isa-$n/transfer; ./target/release/hwcost fit --run $D/isa-$n --residual-from $D/isa-base-1 --out $D/isa-$n/transfer > $D/isa-$n/transfer/fit.txt; done
./target/release/hwcost gpu-check --calibration ecc2k130/benchmarks/roofline/calibration.json --fit-arm base \
  --reference-gpu rtx-pro-6000 --check-rate rtx-pro-6000:14.637530 --check-rate rtx-pro-6000:15.115792 --out $D/gpu
C=research/ic_frontier_campaign_20261015
./target/release/hwcost ic-price --band $C/analysis/band-koblitz.json --sessions $C/phase1/sessions \
  --ic-arm ic-compact --rho-arm rho-frob --tier 27,40 --targets 48,64,80,96,112,128 --out $D/ic-koblitz
./target/release/hwcost ic-price --band $C/analysis/band-prime.json --sessions $C/phase1/sessions \
  --ic-arm ic-pair --rho-arm rho-neg --tier 23,41 --targets 48,64,80,96,112,128 --out $D/ic-prime
```

## Addendum 1 (frozen before its runs): branch-free code that survives the compiler

Written after parts A–D ran. F3 showed that `Fp64Bf`'s masks did not reach
the machine code: in the baseline build `step_fp64` and `step_fp64bf` are
the same instruction sequence up to register allocation (303 instructions,
32 conditional jumps, 4 `cmov`), so F1 (FAIL, worst ratio 1.05) did not test
branch removal. The native build, whose `fp64` kernel ran 2.3–2.4× faster,
emits 22 jumps and 14 `cmov` but also `mulx` and vector code, so it does not
isolate the branches either.

`fp64cm` (`hwcost::kernel::Fp64Cm`) writes the same two selections with
`std::hint::select_unpredictable`. **Precondition, checked before freezing**:
in the baseline build `step_fp64cm` compiles to 285 instructions, 20
conditional jumps and 12 `cmov` (`objdump -d`), against `step_fp64`'s 303,
32 and 4; its results are bit for bit `Fp64`'s (unit test and run gate).

Runs, baseline build at this addendum's commit, back to back:
`hwcost measure --fields fp64,fp64cm --out runs/cm-1`, then `runs/cm-2`
(all arm groups, round-1 parameters), each fitted with `hwcost fit`.

| id | claim | rule |
| --- | --- | --- |
| F1′ | branch-free code costs less on the baseline ISA | `fp64cm` kernel median ≤ 0.85 × `fp64`'s at every B ≥ 4096, in both cm-1 and cm-2 (`hwcost pair --a fp64 --b fp64cm --id F1prime`) |
| F2′ | the a priori model then explains `fp64cm` | V2 on `fp64cm` ≤ 0.30 in both runs |

F1′ passing supports the branch hypothesis on this ISA and attributes the
named share of the native build's `fp64` speedup to it; failing rejects it.

**Also post hoc (labelled as such in the report):** part D's fitted wall-time
exponents for rho are flattened by a fixed per-solve floor (≈3 ms; rho wall
α 0.26 Koblitz, 0.34 prime, against GAE α 0.46–0.47), so the pre-registered
wall-time projections understate rho at large r by orders of magnitude.
`ic-price` now also writes `projections_gae.csv`: the fitted `P(w)` exponents
converted to time at the matched rho's ns per priced GAE at the largest
measured size. Both are reported; the pre-registered one with this defect
named.

```sh
for n in 1 2; do mkdir -p $D/cm-$n; ./target/release/hwcost measure --fields fp64,fp64cm --out $D/cm-$n 2> $D/cm-$n/stderr.txt; done
for n in 1 2; do ./target/release/hwcost fit --run $D/cm-$n > $D/cm-$n/fit.txt; done
mkdir -p $D/cm-pair; ./target/release/hwcost pair --runs $D/cm-1,$D/cm-2 --a fp64 --b fp64cm --min-batch 4096 --threshold 0.85 --id F1prime --out $D/cm-pair > $D/cm-pair/summary.txt
```
