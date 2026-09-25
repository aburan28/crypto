# Pricing the Gröbner stage and the whole logarithm in instructions

**Harness:**
- `examples/groebner_stage_bench.rs`, which gains `--rung i`;
- `examples/ir_calibration.rs` (new);
- `ic run --solver groebner --summands 3`;
- `valgrind --tool=callgrind`.

**Frozen evidence:** `research/groebner_instructions_20260925/`.

**Predecessors:**
- [`RESEARCH_CHAIN_SPLIT_ORDER.md`](RESEARCH_CHAIN_SPLIT_ORDER.md) (#690) and
  [`RESEARCH_SUPPORT_LOCAL_MULTIPLIERS.md`](RESEARCH_SUPPORT_LOCAL_MULTIPLIERS.md)
  (#712), whose results this round re-prices.
- §7 of the latter, which promised this round.

**Status: registered, measured (§5): G1 and G3 met, G2 not met as
registered (one rung, thread scheduling, §5.1); accounting.**
- §0–§4 were committed before any registered run (`31bedbdf`).
- §1 lists every run made before registration.
- Results are appended below §4 and do not edit it.
- **Class, stated in advance: accounting.**  No algorithm changes.  The code
  changes are harness only: a rung filter and a calibration example.

## The question

The Gröbner-stage rounds have priced the decomposition oracle in 64-bit word
operations:
- elimination XORs;
- specialisation reads and writes;
- linear elimination's row operations.

Everything else inside the stage is uncharged in both arms:
- building the decomposition system;
- building rows;
- sorting, hashing and allocation;
- scanning for leading terms;
- substituting generators.

Round 2 already found that an uncharged phase can grow and hide a
regression (its §6).  This round asks how much of the stage the unit sees,
and what the last two rounds and the whole logarithm cost when every
instruction is counted.  The whole logarithm has had a null `S` because
`ic run` counts only its oracle in operations.

## 0. The boundaries, and the unit

**The floor.**  Unchanged.
- For ECC2K-130 it is `m·2^131` candidate tuples
  ([`RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md) §5).
- For a whole logarithm on a ladder cell it is the generic floor
  `√(πr/2)` group operations, which is `S ≈ 1.25`.

**The reference.**  Pollard rho, counted (`ic_boundary::rho_reference`), run
on the same instance and priced in the same unit as the pipeline.

**The unit: instructions retired**, counted by callgrind (valgrind 3.22.0).
- **Stage instructions (`stage Ir`):** instructions executed inside
  `groebner_decompose` (`--toggle-collect='*groebner_decompose*'`).  This
  covers the whole oracle: system construction, the splitting solve with
  every phase in it, lifting and verification of every decomposition.
- **Process instructions (`process Ir`):** every instruction of an `ic run`
  process, from start to the verified logarithm.
- **Conversion to group-addition equivalents (GAE):**
  - `ir_calibration add` measures instructions per group addition on the
    fast arithmetic that the rho reference and the oracle ladder's GAE use.
    It runs `ic_boundary::calibrate_group`'s own loop, on the same binary and
    the same `(a, n)`.
  - `S = process Ir / (Ir per addition) / √r`.
  - Rho gets its `S` by the same formula, from `ir_calibration rho`.

**Why instructions.**
- They are **complete**: every cost the method incurs is in them, including
  the ones no counter was written for.
- They are **deterministic**: a repeated run of `K_0/2^23 m=3` differed by 28
  of 2,136,316,624, with the arms on one binary.
- They are **independent of the host**.

What they are not is independent of the build.  A different compiler or
binary moves them, so:
- every comparison is made on one binary;
- the binary's hash and the valgrind version are pinned in the manifest;
- word operations stay beside them as the build-independent column.

## 1. What was run before registration

All on `2026-09-24`/`25`, on `main` at `48229a70` plus the `--rung` filter,
disclosed so that the cells it touched are read as tuning data:

- **A profile of every rung of the five stage ladders**
  (`presweep/`, `presweep/table.md`) in the current default (support-local
  multipliers) and in round 2's reference (occurring-variable multipliers).
  - **Counted word operations are 0.2–2% of the stage's instructions**:
    51–586 instructions per counted operation, an 11× spread across rungs.
  - **In instructions, round 2 is `1.12×`** on its registered holdout against
    `1.50×` in word operations.  It is `1.01×` on the four largest cells
    (`K_0/2^31`, `K_1/2^29`) against `1.43–1.63×`, `1.17×` on the
    predecessor's T1′ against `1.65×`, and `1.00×` on the frozen ladder
    against `1.02×`.
  - One frozen rung, `K_1/2^15 m=2`, is `0.79×` in instructions
    (29.5 M → 37.6 M) and `1.04×` in word operations.
  - Where the instructions go, summed over every rung of the current default:
    specialisation `51%` (`ReducedBasis::insert` `29%`), sorting `22%`,
    root builds `18.5%`, linear elimination `11.6%`, system construction
    `8.7%`.  The counted elimination kernels are `8%`.
  - On `K_0/2^31 m=3`, the rewrite scatter inside specialisation is `36%` of
    the rung's instructions.  It is charged per live word but costs per set
    bit, about 186 instructions per counted operation.  The m4ri elimination
    is `9.6%`, at about 7.5 instructions per counted XOR.
- **Round 2's whole-log wall ratio** after its hash fix was `1.41×`
  (`K_0/2^13`, 15 seeds).  Instructions on one of those runs (seed 210,
  round 2's `profile/`) gave `1.46×` in the stage (332.9 M → 228.3 M) and
  `1.33×` for the whole process (417.3 M → 312.8 M).  On the R2 holdout, wall
  gave `1.08×` and instructions `1.12×`: instructions track wall, and word
  operations do not.
- **Calibration probes:**
  - `ir_calibration add`: 924 instructions per addition on `K_0/2^13` and
    680 on `K_0/2^9`, over 200,000 additions.
  - `ir_calibration rho` on `K_0/2^13`, seeds 201–203: 940 k, 706 k and
    555 k instructions for counted `S` of 22.4, 16.8 and 13.2.  Rho's
    instructions per counted GAE are 938, 937 and 940, so its overhead
    beyond the additions is about 1%.

## 2. The arms, and the suites

**Arms,** every control explicit and recorded, all on one binary:
- **A0**, the default before #690: `KIC_CHAIN_ORDER=layout KIC_LINEAR_ELIM=0
  KIC_F4_MULTIPLIERS=occurring KIC_F4_DROP=complete`.
- **A1**, the default between #690 and #712: `interleaved`, `1`,
  `occurring`, `complete`.
- **A2**, the current default: `interleaved`, `1`, `support`, `complete`.

**Stage ladders.**
- All five ladders: frozen, chain, chain-holdout, chain-holdout-2, R2 holdout.
- Every admissible rung, one callgrind run per rung per arm (`run.sh`).
- Three rungs per arm are run a second time for the determinism gate: the
  frozen `K_0/2^13 m=2`, chain-holdout-2 `K_0/2^23 m=3`, and R2
  `K_0/2^9 m=4`.

**Whole logarithms,** in process instructions: `ic run --degree n --curve-a 0
--summands 3 --solver groebner --random-target --batch 1`, one run per seed
per arm.
- `K_0/2^13`: seeds 201–210, plus holdout 301–305.
- `K_0/2^9`: seeds 201–205.

These are round 2's seeds.

**Reference,** in the same unit: `ir_calibration rho 0 n seed` on the same
twenty seeds, and `ir_calibration add 0 n 200000` for `n = 13` and `9`.

## 3. Gates, and the convention adopted

This round has no speedup target; it is accounting by construction.  Its
gates are:

- **G1, identity.**  Every counter of every rung equals the registered run of
  the same arm where one exists:
  - A1 and A2: the registered runs of round 2;
  - A0: round 1's registered reference on the four ladders it ran (it never
    ran the R2 holdout).
  - For the whole logarithms, A1's and A2's trials, relations and oracle word
    operations equal round 2's registered runs, seed by seed.
  - A failure voids that arm's row.
- **G2, determinism.**  Each rerun reproduces its stage instructions within
  one part in `10⁶`.
- **G3, correctness.**  Every whole logarithm recovers the planted `k` and
  verifies `[k]G = Q`, and every rho run recovers and verifies its `k`.

Reported for every rung and arm, not gated:
- stage instructions, word operations and their ratio;
- the phase shares;
- the ratios of round 1 (A0 / A1) and round 2 (A1 / A2) in both units.

For the whole logarithms, `S` for every run of every arm and of rho, and the
ratio to rho.

**The convention this round adopts,** stated before its numbers:
- From the next Gröbner-stage round on, every table carries stage
  instructions beside word operations.
- A claimed stage gain needs the stage-instruction total to fall on the
  registered holdout.
- Any per-rung no-regression clause applies in instructions as well as in
  word operations.
- A whole-logarithm claim is stated in `S` from process instructions, with
  the ratio to rho measured the same way.
- Word operations remain the build-independent diagnostic.

**Inadmissible:**
- choosing rungs or seeds after seeing their instructions;
- narrowing the collection toggle to leave a phase out;
- comparing arms across binaries;
- counting instructions of a run whose identity gate failed.

## 4. What the numbers will and will not change

- **Round 1's and round 2's registered results stand as registered,** in the
  unit they were registered in.  This round adds their instruction-priced
  ratios beside them and changes neither round's status.
- **The scoreboard shows both.**  If round 2's R2-holdout ratio in
  instructions is below its registered `1.3×` (§1 suggests `1.12×`), the page
  says so beside the `1.50×`.
- **The whole-logarithm `S` stops being null** for these two cells.  It
  remains a measurement at `n = 9` and `13` only.  Rho's `S` at these sizes
  is dominated by its jump-table setup (§1: 13–22) and is not the asymptotic
  `1.3`.  The ratio to rho is taken against rho measured here, not against
  `1.3`.
- **No phase gets a new counter.**  Whether to charge the specialisation
  scatter per set bit, or the system construction per term, is left to a
  later round.  That round would write the counter to land within a stated
  band of these instruction counts.

## 5. Results

*Appended after every registered run; §0–§4 are not edited.  Every number below
is read from `research/groebner_instructions_20260925/`: `tables.md` is
generated by `table.py` from the `summary.json` of every job, next to its
gzipped callgrind profile.  Every job ran on one binary built from
`31bedbdf`.*

**Bottom line.**
- **The word-operation unit sees 0.2–2% of the Gröbner stage's
  instructions**: 51–586 instructions per counted operation, 98 over the R2
  holdout and 258 over the frozen ladder.
- **Round 1 survives the re-pricing.**  In instructions it is `5.23×` on its
  holdout T1′, against its registered `4.91×` in word operations.  On the R2
  holdout, which it never ran, it is `10.6×` against `17.3×` in word
  operations.
- **Round 2 shrinks to `1.12×`** on its registered holdout, against `1.50×`.
  It is `1.01×` on the three `K_0/2^31` cells and `1.00×` on the frozen
  ladder.  One frozen rung regresses to `0.79×`.
- **The whole logarithm has a measured `S`.**
  - `K_0/2^13`: `S ≈ 19,100` against rho's `19.8` on the same seeds and in the
    same unit, **`964×` rho**.
  - `K_0/2^9`: `S ≈ 8,700` against `45.2`, `192×`.

The class is **accounting**: no algorithm moved.

### 5.1 The gates

- **G1, identity: met.**  Every field but timings, on every rung, equals the
  registered run of the same arm:
  - A1 and A2 against round 2, on all five ladders;
  - A0 against round 1's reference, on the four ladders it ran.  That is 27
    fields, the keys both harness versions record.

  Whole logarithms: A1's and A2's trials, relations and oracle word
  operations equal round 2's on all twenty seeds.
- **G2, determinism: not met as registered.**
  - Two of the three rungs reproduce within the gate: chain-holdout-2
    `K_0/2^23 m=3` differs by at most `3.2 × 10⁻⁷`, and R2 `K_0/2^9 m=4` is
    exact in A1 and A2 and `7.4 × 10⁻⁷` in A0.
  - The frozen `K_0/2^13 m=2` does not: `2.5 × 10⁻⁵`, `7.0 × 10⁻⁶` and
    `8.6 × 10⁻⁶` in A0, A1 and A2, against a gate of `10⁻⁶`.
  - The cause, found after the registered runs (`determinism_single_thread/`):
    rayon's scheduling in the parallel Macaulay elimination, which that
    rung's root matrices are large enough to take.  Run single-threaded, the
    rung gives 693,232,587 instructions twice.
  - The noise is at most `2.5 × 10⁻⁵`.  Every ratio reported here is stated to
    two decimals, and the smallest arm difference any conclusion rests on is
    `1.3 × 10⁻⁴`, so no reading below changes.
  - The unit's stated determinism is corrected from `10⁻⁶` to `3 × 10⁻⁵`
    for rungs that reach the parallel path.  A later round that needs better
    can set `RAYON_NUM_THREADS=1`.
- **G3, correctness: met.**  All 60 whole logarithms recover the planted `k`
  and verify `[k]G = Q`, and all 20 rho runs recover and verify theirs.

### 5.2 The stage, in both units

Totals per ladder; `tables.md` has every rung.  Ratios are old / new.

| ladder | rungs | word ops A0 → A1 → A2 | instructions A0 → A1 → A2 | instructions per word op (A2) | round 1: words, **instructions** | round 2: words, **instructions** |
|:--|--:|:--|:--|--:|:--|:--|
| frozen | 6 | 16.70 M → 16.13 M → 15.82 M | 4.19 G → 4.10 G → 4.08 G | 258 | 1.04×, **1.02×** | 1.02×, **1.00×** |
| chain (tuning) | 3 | 13.75 M → 7.55 M → 3.35 M | 1.66 G → 1.37 G → 0.94 G | 280 | 1.82×, **1.21×** | 2.25×, **1.45×** |
| chain-holdout (round 1's holdout) | 4 | 379.3 M → 52.3 M → 32.8 M | 48.3 G → 7.74 G → 6.30 G | 192 | 7.25×, **6.24×** | 1.60×, **1.23×** |
| chain-holdout-2 (round 1's T1′) | 9 | 1,104.5 M → 225.1 M → 136.4 M | 143.5 G → 27.4 G → 23.4 G | 172 | 4.91×, **5.23×** | 1.65×, **1.17×** |
| R2 holdout (round 2's T1) | 19 | 15,410 M → 888.6 M → 591.9 M | 686.5 G → 65.0 G → 58.0 G | 98 | 17.3×, **10.6×** | 1.50×, **1.12×** |

Reading it:

- **Round 1's gains were mostly real.**
  - On its registered holdout T1′ the instruction ratio (`5.23×`) exceeds
    the word-operation ratio.  Its largest cell, `K_1/2^11 m=4`, is `38×` in
    instructions against `33×` in word operations.
  - On the three `K_0/2^31` cells it never ran it is `14–24×`.
  - Where it did little, it did little in both units: the frozen ladder is
    `1.02×`, and `K_0/2^13 m=3` is `1.03×` in instructions against `1.62×`.
- **Round 2's gains were mostly in the phase the unit over-weights.**
  - On the four cells it had never seen, T1 is `1.01×`, `1.01×`, `1.01×` and
    `1.90×` in instructions, against `1.43–1.63×` in word operations.
  - Its **registered `1.3×` threshold is not met in instructions** (`1.12×`).
  - Round 2's status is unchanged, because it was registered in word
    operations (§4).  The page shows both figures.
  - Its genuine instruction gains are on small and middling chains:
    `1.3–2.0×` on the `n ≤ 11` cells and `1.46×` on `K_0/2^13 m=3`.
- **One frozen rung regresses in instructions under round 2:**
  `K_1/2^15 m=2`, 29.5 M → 37.5 M (`0.79×`), with word operations `1.04×`.
  - The rung has 32 targets and no tree: every node is a root.
  - Its generators do not all span every variable, so round 2's builder
    takes the support-local path, which has no layout cache.  The rung's
    layout-cache hits fall from 53 of 55 builds (A1) to 13 (A2), and the
    other 40 builds sort columns (`macaulay_columns`, 3.7 M instructions) and
    hash-pack rows (`pack_rows`, 2.3 M).  The reference used the occurring
    path's cached fused packer instead.
  - Under the convention of §3 this is a per-rung regression.  It is 8 M
    instructions on a 4.1 G ladder.
  - The fix, a layout cache for support-local roots, is listed for the next
    round rather than made here.
- **The unit over-weights elimination.**  Instructions per counted operation
  run from 51 (`K_0/2^31`, where the m4ri kernels dominate) to 586 (the
  `n = 9` cells, where system construction and root builds dominate).

### 5.3 Where the default's stage instructions go

| ladder | system construction | root builds | specialisation | of which `insert` | linear elimination | counted elimination kernels |
|:--|--:|--:|--:|--:|--:|--:|
| frozen | 33% | 23% | 35% | 15% | 0% | 10% |
| chain | 31% | 45% | 12% | 5% | 1% | 15% |
| chain-holdout | 5% | 19% | 61% | 25% | 6% | 6% |
| chain-holdout-2 | 9% | 18% | 46% | 23% | 14% | 6% |
| R2 holdout | 7% | 18% | 54% | 33% | 12% | 9% |
| **all rungs** | **9%** | **19%** | **51%** | **29%** | **12%** | **8%** |

The phases overlap: root builds contain some elimination, and specialisation
contains `insert`.  The shares name where to look.

- **Large cells:** specialisation is the lever.  Its rewrite scatter is
  charged per word but costs per set bit (§1), and `insert` scans for leading
  terms uncharged.
- **Small and quadratic cells:** system construction (31–33%) and root
  builds are the levers.  System construction is target-independent except
  for the last link, and `polynomial_reuse::DecompositionTemplate` already
  implements the split, behind an opt-in serialising cache.

### 5.4 Whole logarithms

`ic run --solver groebner --summands 3`, process instructions, one run per
seed per arm.  The conversion is 923.9 instructions per addition on
`K_0/2^13` and 679.9 on `K_0/2^9` (`add/`).  Rho is `ic_boundary::rho_reference`
on the same instance and seeds (`rho/`).

| cell | runs | rho `S` (instructions; counted) | A0 `S` | A1 `S` | **A2 `S`** (the default) | A2 / rho |
|:--|--:|:--|--:|--:|--:|--:|
| `K_0/2^13`, seeds 201–210, 301–305 | 15 per arm | 19.8; 19.5 | 27,185 | 26,999 | **19,128** (7,508–40,130) | **964×** |
| `K_0/2^9`, seeds 201–205 | 5 per arm | 45.2; 43.9 | 12,395 | 10,464 | **8,674** (7,892–9,721) | **192×** |

- **The conversion checks itself.**  Rho's instruction-derived `S` and its
  counted `S` agree within 3%: its overhead beyond the additions is small, as
  §1 found.
- **The Gröbner route is three orders of magnitude above rho at `n = 13`.**
  The whole-log `S` includes everything: factor base, relation collection
  with the oracle (`89%` of the process at `n = 13`, `64%` at `n = 9`),
  linear algebra and verification.
- **The ratio to rho grows from `n = 9` to `n = 13`.**  Two sizes are not an
  exponent, and none is fitted.
- **The rounds in whole-process instructions** (sums over the seeds, old /
  new):
  - Round 1: `1.01×` at `n = 13` and `1.18×` at `n = 9`.  Its `1.54×` in
    oracle word operations at `n = 13` was almost all in the phase the unit
    over-weights.
  - Round 2: `1.41×` and `1.21×`.  This matches its wall-time `1.41×` from
    its §6.
  - Together: `1.42×` and `1.43×`.

### 5.5 What this does and does not change

- **Round 1 and round 2 keep their registered status.**  Round 2's T1 read
  in instructions is `1.12×`, below its `1.3×`; the scoreboard now shows
  that beside the `1.50×`.
- **The page's verdict is unchanged.**  Nothing measured end to end is below
  rho, and the Gröbner route is now measured at `192–964×` rho on these two
  cells.  The oracle ledger's matrix-F4 figures are word XORs converted to
  GAE.  By §5.2 they undercount the oracle's instructions roughly
  100–500-fold, and the page says so beside them.
- **What it establishes for the next rounds** is the convention of §3:
  - stage instructions beside word operations;
  - a stage gain needs the instruction total to fall on the registered
    holdout;
  - no-regression clauses apply in both units;
  - whole-log claims are made in `S` from process instructions.
- **What it does not establish:**
  - an exponent;
  - any size beyond the ladders;
  - a statement about builds other than this one: instructions are
    build-dependent, which is why every comparison here is on one binary.

### 5.6 The next levers, measured here, none claimed

1. **Specialisation on large cells:** 51% of the default's instructions,
   `insert` 29%.
2. **System construction on small and quadratic cells** (31–33%): keep the
   decoded `DecompositionTemplate` in memory, bypassing the serialising cache
   whose overhead sank it on 2026-09-14.
3. **The `K_1/2^15 m=2` regression:** give support-local roots a layout
   cache.

### 5.7 Reproducing

```bash
cargo build --release --example groebner_stage_bench --example ir_calibration --bin ic
research/groebner_instructions_20260925/run.sh      # 214 callgrind jobs, 4 at a time
python3 research/groebner_instructions_20260925/table.py > research/groebner_instructions_20260925/tables.md
```

`table.py` exits non-zero because G2 failed; `tables.md` records the failure.
