# Cheaper table-walk selection (in progress)

Question: [ONE-BLOCK-GEOMETRY.md](ONE-BLOCK-GEOMETRY.md) §6 leaves the
forward pass short of a fully fed carry-less unit because selection issues
~54 `LDS` per slot per warp (34 of them random byte lookups). Can a different
selection cut those lookups enough to move complete-update throughput, without
changing the walk?

This note declares the experiment before any GPU measurement. Class will be
**engineering** if `S` (here B complete updates/s) rises and the ratio to the
CLMAD floor moves only because the idle fraction shrinks — the floor itself is
unchanged (still 33.1 CLMAD/update → 22.3 B/s on the RTX PRO 6000).

## 1. Boundary

Unit: billions of complete scalar updates per second on one RTX PRO 6000, same
command as ONE-BLOCK-GEOMETRY (`--curve 131 --packed --bench --steps 1024
--launches 32 --verify 0` at automatic workers), with 300/300 re-walks and an
identical distinguished-point multiset against the control binary.

**Floor (derived):** the carry-less unit, unchanged: **22.3 B/s** at 2.42 GHz
for 33.1 CLMAD/update (ONE-BLOCK §1). No selection change alters that count.

**Reference:** `make gpu-rtx-pro6000-20b` as shipped: **20.078 B/s** median
(0.90 of the floor). Selection cost inside that kernel is the residual the
forward pass still pays after `TABLE_PIPE_SELECT`.

## 2. Falsification target (declared in advance)

Success on one RTX PRO 6000, paired against the reference rebuilt in the same
session:

* median complete-update rate ≥ **20.50 B/s** (about +2.1%, enough to clear
  noise and show the MIO cut), every pair favouring the candidate, and
* 300/300 device reports re-walked, byte-identical DP set on a forced common
  walk count, and
* shared-memory footprint still ≤ 48 KB per block (two-block geometry remains
  buildable; the 20 B/s build uses one block and has more headroom).

Abandon if the paired median is below **20.20 B/s**, if any correctness gate
fails, or if the layout forces spills / register growth that drops occupancy
below the reference's 16 warps per SM.

Inadmissible: changing the walk, the CLMAD count, the batch, or the DP rule;
quoting a phase-only cycle cut as a method speedup; measuring on a different
SKU without a matched reference rebuild.

## 3. Candidate under build: `TABLE_PHASE_POPC`

Replace the 17 byte-table lookups in `twPhase` with eight bit-plane popcounts
against precomputed masks of coordinates whose `L` bit is set — the form
`packedtablewalk.cuh` rejected for the *whole* selection when priced against
ALU, but which may win for *phase alone* now that the pivot is already on the
byte table (`TABLE_PIVOT_BYTES=1`) and the forward pass is MIO-bound rather
than ALU-bound (ONE-BLOCK §6).

Static price to fill before any GPU run (host-equivalent `TW_FN` path):

| piece | reference (byte table) | candidate (popc planes) |
|---|---|---|
| phase `LDS.U8` | 17 | 0 |
| phase ALU (priced in LOP3 slots) | ~72 (THROUGHPUT-20B §3) | 8×`POPC` + masks + reduce; measure with `kernel_cost.py` / host objdump |
| shared bytes for phase | 17×256 | 131×5 mask words (same as `maskLt` family) or less if folded into existing masks |

Host gate: `make test-table-walk-host` must pass with both
`TABLE_PHASE_POPC=0` and `1` against the scalar reference for every tested
point. The candidate stays **off by default**.

## 4. What was checked and rejected on the way here

| idea | result | class |
|---|---|---|
| Dense GFNI for bitsliced `multPrep`/`toOnb` (HOST-SCHEDULE §7) | Static: ≥17² GFNI + transpose ≈ 1.6k ops vs 1.2k XOR instr today; GFNI acts on bytes *inside* a register, bitsliced maps XOR whole words — needs a layout change first | accounting (priced miss) |
| Host `BATCH=48/64/128` vs 32 | Paired medians +1–2%, pairs not all >1 on this container; noise dominates (same failure mode as HOST-SCHEDULE §4) | inconclusive; default stays 32 |
| Host `WITNESS=1` vs 0 | Paired median on/off ≈ 0.996 | neutral |

## 5. Status

* Boundaries and target: this file.
* Implementation: `TABLE_PHASE_POPC` knob in `include/packedtablewalk.cuh`,
  host fill of the eight L-bit planes, `twPhase` popc path, Makefile flag,
  and `make test-table-walk-host` covering nibble/bytes × table/popc (4
  binaries). Host gate: **PASS**, 4096 points, zero mismatches on every arm.
  Shared bytes with `PIVOT_BYTES=1`: 48,732 → **44,540** (−4,192).
* GPU measurement: blocked here (no `nvcc`). When a card is available, run
  the 20 B/s reference paired with `TABLE_PHASE_POPC=1` and file the receipt
  under `benchmarks/cheaper-selection/`. Default stays off until that clears
  §2.

