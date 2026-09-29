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

## 3. Candidate: `TABLE_PHASE_POPC`

Replace the 17 byte-table lookups in `twPhase` with eight bit-plane popcounts
against precomputed masks of coordinates whose `L` bit is set — the form
`packedtablewalk.cuh` rejected for the *whole* selection when priced against
ALU, but which may win for *phase alone* now that the pivot is already on the
byte table (`TABLE_PIVOT_BYTES=1`) and the forward pass is MIO-bound rather
than ALU-bound (ONE-BLOCK §6).

### Static price (filled)

Slot weights from `benchmarks/clmad-price` on sm_120: LOP3/PRMT = 1,
POPC = 3.97. Host gate `make test-phase-cost` (same points, both arms):

| piece | byte table (`PHASE_POPC=0`) | popc planes (`PHASE_POPC=1`) |
|---|---|---|
| random `LDS.U8` | **17** | 0 |
| broadcast `LDS.U32` (plane masks) | 0 | 40 |
| ALU slots (static) | ≈ 33 source / ~72 with address math (THROUGHPUT-20B) | **AND 40 + POPC 40×3.97 + shift/add ≈ 245.8** |
| shared bytes (`PIVOT_BYTES=1`) | 48,732 | **44,540** (−4,192) |
| host ns/call (practicality) | 4.4 | 85.1 |

Reading the ALU column: the candidate is **~3–7× dearer** on the logic pipe
than the byte table, in exchange for deleting 17 random byte loads. It only
pays if those loads are what starve the carry-less unit in the forward pass
(ONE-BLOCK §6). That is a GPU question; the static table alone would abandon
the knob. The falsification target in §2 still stands — a card either clears
20.50 B/s or the candidate is abandoned below 20.20.

Host gate: `make test-table-walk-host` and `make test-phase-cost` both PASS.
Off by default.

### GPU job

```sh
make bench-cheaper-selection OUT=/tmp/cheaper-selection
# or: modal run --detach modal_job.py \
#        --job benchmarks/cheaper-selection/gpujob.sh \
#        --out DIR --gpu "RTX PRO 6000"
```

Paired `gpu-rtx-pro6000-20b` ± `TABLE_PHASE_POPC`, 300/300 re-walks, DP-set
identity, five alternating benches, profile dumps. Receipts under
`benchmarks/cheaper-selection/`.

## 4. What was checked and rejected on the way here

| idea | result | class |
|---|---|---|
| Dense GFNI for bitsliced `multPrep`/`toOnb` (HOST-SCHEDULE §7) | Static: ≥17² GFNI + transpose ≈ 1.6k ops vs 1.2k XOR instr today; GFNI acts on bytes *inside* a register, bitsliced maps XOR whole words — needs a layout change first | accounting (priced miss) |
| Host `BATCH=48/64/128` vs 32 | Paired medians +1–2%, pairs not all >1 on this container; noise dominates (same failure mode as HOST-SCHEDULE §4) | inconclusive; default stays 32 |
| Host `WITNESS=1` vs 0 | Paired median on/off ≈ 0.996 | neutral |

## 5. Status

* Boundaries and target: this file.
* Implementation: `TABLE_PHASE_POPC` in `include/packedtablewalk.cuh`, identity
  line `packed table phase popc:`, `make test-table-walk-host`,
  `make test-phase-cost`, GPU job `benchmarks/cheaper-selection/gpujob.sh`.
* Host gate: **PASS** (4096 points, zero mismatches; shared 48,732 → 44,540).
* Static price: **filled** — candidate is much heavier on ALU (~246 slots vs
  ~33–72); only a MIO-starvation win on the 6000 can save it. Default stays off.
* GPU measurement: blocked here (no `nvcc` / Modal token). Run
  `make bench-cheaper-selection OUT=…` on an RTX PRO 6000 to clear or abandon
  §2; file the receipt under `benchmarks/cheaper-selection/`.

## 6. If the card abandons popc

The next cheaper-selection lever to price (not built here) is cutting the
remaining **pivot** byte lookups or the selection's address math without
touching POPC — the byte phase table is already cheap on the ALU column, so
the leftover MIO is the 17 pivot `LDS.U8` plus sign/history, not the phase.

