# Fused sigma B16 tuning star: frozen RTX PRO 6000 screen

Protocol source parent: `73556d4bcdee3e464f58608be5a8b2986496860d`.
The counter-free fused B16 path was frozen in `9bfe4ef05bf8f2b8a37ed8ae10941d92ededa41d`
and confirmed in `e0b0858b179fe28f753353e7a7313cfe23d71884` at
15.436677 billion complete scalar updates/s.  The later geometry experiment
retained B16/T256/minBlocks2.  This document freezes a one-knob tuning screen
before any star measurement is made.

Class: compatible kernel engineering.  The job runs a benchmark and bounded
correctness replay only.  It does not run a search, solver, collision recovery,
key recovery or multi-GPU aggregation.

## Fixed baseline

The hardware and toolchain are exactly one NVIDIA RTX PRO 6000 Blackwell
Server Edition, native `sm_120` code and CUDA 13.3.73.  A launch fails before
building unless the hardware, compiler and one clean 40-hex `SOURCE_REV`
match.  Every build uses batch 16, 256-thread blocks, minimum two resident
blocks, 96,256 correctness threads and 385,024 benchmark threads.

The common arithmetic and walk are:

```
PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
PACKED_GENERATED_PRODUCT=1 PACKED_INLINE_POLY=3 PACKED_CLMAD=1
PACKED_CLMAD_SQUARE=0 PACKED_KARAT3=0 PACKED_WEIGHTED_PREFIX=2
PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 PACKED_TOP_CLMAD=0
PACKED_STATE_TILE=256 PACKED_ADD_COMBINE=0 PACKED_TOP_HOIST=0
PACKED_ONB_INV=0 PACKED_ALU_SQR=0 PACKED_SQUARE_TABLE=0
PACKED_CHAINS=1 PACKED_SLOT_PREFETCH=0 PACKED_SLOT_PIPELINE=0
SIGMA_FUSED=1 WITNESS=0 WALK_TABLE=0 PHASE_PROFILE=0
```

The baseline additionally fixes the nine varied settings to
`PAIR_ILP=0`, `PAIR_CLMUL=0`, `L2_PERSIST=0`, `UNROLL_SLOTS=1`,
`FROM_REDUCED=0`, `INV_POLY=0`, `CLMUL_FLAT=0`, `ALU_SQUARE=0` and
`SIGMA_FUSED_LATE_Y=0`.

## Frozen star

Each candidate differs from that baseline in exactly the one setting shown.
The native `star_arms.hpp` table is the executable arm registry used by the
builder, log checker and summarizer.

| arm | only changed setting | reason it remains admissible |
|---|---|---|
| `pair-ilp` | `PACKED_PAIR_ILP=1` | changes the two-product buffer schedule |
| `pair-clmul` | `PACKED_PAIR_CLMUL=1` | interleaves the twelve carry-less issues before reduction |
| `l2-persist` | `PACKED_L2_PERSIST=1` | changes the access policy for the otherwise identical field blob |
| `unroll2` | `UNROLL_SLOTS=2` | changes only the fused reverse-loop unroll request |
| `from-reduced` | `PACKED_FROM_REDUCED=1` | uses the reduced-input inverse basis conversion |
| `inv-poly1` | `PACKED_INV_POLY=1` | inlines the polynomial-basis inversion links |
| `inv-poly2` | `PACKED_INV_POLY=2` | uses one out-of-line copy of the same link product |
| `clmul-flat` | `PACKED_CLMUL_FLAT=1` | groups ordinary product low halves before high halves |
| `alu-square` | `PACKED_ALU_SQUARE=1` | moves the polynomial square to the logic-pipe spelling |
| `late-y` | `SIGMA_FUSED_LATE_Y=1` | delays normal-Y materialization on the fused DP-zero hot path |

Prior table-walk and B200 measurements do not decide this fused RTX PRO 6000
screen because the walk, fusion and resource balance differ.  They do identify
combinations that are structurally redundant here, and those are pruned:

- pair-CLMUL is an alternative branch to pair-ILP in
  `mulPolynomialPair131`, so there is no pair-CLMUL plus pair-ILP arm;
- pair-CLMUL supplies its own interleaved product and therefore does not
  consume the ordinary flat-CLMUL schedule; there is no combination arm;
- `INV_POLY=1` and `=2` compute the same inversion chain with alternative
  inlining, so they are separate one-knob arms and are not composed;
- late-Y is valid only with `SIGMA_FUSED=1`, which is fixed in every arm;
- the square table is table-walk-only and is incompatible with this sigma
  star.

No candidate inherits a prior candidate's setting.  A later confirmation may
compose independently confirmed winners only under a separately frozen
protocol.

## Fail-closed correctness and identity

The job hashes the Makefile, all source/include/generated files and the native
protocol tools.  It records the exact source commit, full common and varied
Make settings, build command, binary SHA-256, ptxas resource excerpt, GPU UUID,
driver, clocks, power limit and memory.

Every execution log is reopened by `star_log_check.cpp`.  The checker requires
all common and varied runtime identity lines, including
`packed inline polynomial: 3`, and it requires the runtime marker
`packed launch bounds: 256 threads, 2 min blocks`.  It also records registers,
local bytes, static shared bytes and the L2 access-policy window.  The
L2-persist arm is rejected if a nonempty window was not actually installed.
The backend line must name the exact fixed workers, slots, DP weight and steps.

Before timing, every arm runs seven odd 95-step launches at DP weight 48 with
run id 37.  Each must replay 300/300 reports against the host reference with
zero drops and no mismatch or overflow.  The native corpus comparator sorts
and compares all eleven headerless-v1 corpora and retains a canonical stream
and SHA-256.  Any build, marker, resource, replay or corpus failure suppresses
all timing.

## Equal-work timing

Every warmup, A/A row and star row completes exactly
`385024 * 16 * 1024 * 32 = 201863462912` scalar updates.  The panel contains:

1. one excluded equal-work warmup for each of the eleven binaries;
2. five baseline/baseline pairs with alternating A/B position;
3. three baseline/candidate pairs per candidate, with candidate traversal
   forward, reverse and half-rotated and A/B position alternating by arm and
   round.

Each log is checked natively before its rate is admitted to the TSV.  The TSV
also binds the log SHA-256 and the observed SM clock, power draw and
temperature.  There are exactly 81 timing rows; a missing, duplicate or
wrong-work row invalidates the result.

## Frozen decision

The native summarizer first requires A/A maximum symmetric drift no greater
than 1%.  For each candidate it computes the three paired candidate/baseline
ratios.  Its geometric-mean threshold is

```
max(1.015, 1 + A/A max drift + 0.005)
```

and its minimum-pair threshold is

```
max(1.005, 1 + A/A max drift).
```

A candidate qualifies only if it clears both thresholds and its candidate
median exceeds its paired-baseline median.  The highest qualifying geometric
mean is selected for a separate five-pair confirmation; otherwise
`selected_arm` is `null`.  The screen does not promote a preset by itself.
The 26,000 M/s objective is marked met only when the selected candidate's
screen median reaches it, and still requires the separate confirmation before
an engineering promotion.

`test_star_native.sh` compiles the native tools with warnings as errors,
exercises select-none, winner, malformed-panel, excessive-noise and log-marker
cases, and runs the packed host arithmetic control for each applicable
arithmetic arm.  Neither the test nor the GPU job invokes Python.
