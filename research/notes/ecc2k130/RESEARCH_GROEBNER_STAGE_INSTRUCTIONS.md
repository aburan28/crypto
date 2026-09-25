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

**Status: registered, not yet measured.**
- §0–§4 are committed before any registered run.
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
