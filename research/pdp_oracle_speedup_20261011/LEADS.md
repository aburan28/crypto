# Leads from the pairs-and-solve speedup program

Written 2026-10-12 while PS-0…PS-3 had produced evidence, PS-4/PS-5/PS-7 had stopped early (the authoring agents hit the
organisation's API spend limit; their partial work is on `exp/pdp-ps4-symmetry`, `exp/pdp-ps7-cdcl`, and
`exp/pdp-ps5-gpu-emulation` in cryptanalysis) and PS-3's timing stage was running.  A lead is a thing the evidence points at, with
the evidence, the cheapest test, and what would refute it.  Nothing here is a result.  Status words: **open** (not tested),
**running**, **tested** (with the PR or record), **dropped** (with the reason).

| id | lead | evidence it came from | cheapest test | refuted by | status |
|:--|:--|:--|:--|:--|:--|
| L1 | The repository's own meet-in-the-middle pair table (`ic_boundary.rs` `Oracle::Mitm`, `decompose_mitm`) is the incumbent for `m = 3`; a speedup of pairs-and-solve has to be read against it, and the program has no measurement that does. | PS-1/PS-2: the pairs kernel's candidate `x(P₁ ± R ± P₂)` is the table key computed per target; reading `ic_boundary.rs` (oracle list, `Oracle::Mitm`). | PS-8: a `mitm` arm in the same harness and corpus; probe time, build time, memory, break-even number of targets. | H1…H4 of [`protocols/PS-8.md`](protocols/PS-8.md) failing; in particular a break-even above 100 targets or a probe time above 1/10 of the kernel's. | **open** (arm written, not yet built or measured) |
| L2 | The table's build is the pair loop run once, so PS-3's bit-sliced kernel is a candidate table builder; the memory/time frontier between "no table" (pairs-and-solve) and the full table is a partial table. | L1; PS-3's lane kernel computes exactly `x(Pᵢ ± Pⱼ)`-type quantities. | After L1: build the table with the vector kernel and compare build seconds; measure relation yield against stored fraction on a partial table. | build time not below the scalar build; yield not proportional to the stored fraction. | open, waits on L1 |
| L3 | PS-3's measured gap to its own model has three named causes (product schedule, 30% non-product work, 65% lane utilisation at 512): a lean `stage_c` that computes only planes `l…n`, in-place temporaries, and a diagonal-block lane map would move the bit-sliced arm further. | `ps3_bitsliced/tuning_log.md`, `raw/micro.jsonl`, `raw/utilisation.jsonl`. | Implement the lean `stage_c` (PS-2's program already dead-code-eliminates those planes) and re-run `bs-micro batch8`. | `batch8` cycles do not fall by at least the non-product share. | open |
| L4 | The vector arm's advantage over `group-scalar` is bounded by the host's `vpternlogq` rate, not by the register width: AVX2 (256) and AVX-512 (512) are close. | `tuning_log.md` (1.3 cycles/gate at 512, about 0.9 at 256 per gate), `raw/times.jsonl`. | Repeat the product micro-benchmark on a host with two 512-bit logic ports. | A second host class showing the model's 2 gates/cycle at 512. | open, needs another host class |
| L5 | PS-4's symmetry-reduction experiment (Frobenius/negation orbits on the SAT side) was interrupted after its theorem checks and the panels at `n = 11…39`. | `ps4_symmetry/raw/` snapshots on `exp/pdp-ps4-symmetry`. | Resume the interrupted pipeline (`run_ps4.sh` is resumable), finish README/report. | The pipeline's own stop rules. | **stopped early**, resume |
| L6 | PS-7 (CDCL shell with the algebraic leaf) stopped during code; its frozen hypothesis is a negative. | `examples/pdp_ps7/` snapshot on `exp/pdp-ps7-cdcl`. | Resume and run the 60-input regression corpus and the three rungs. | Success/stop rules of `protocols/PS-7.md`. | **stopped early**, resume |
| L7 | PS-5's device kernel and its CPU emulation are written but unrun: no parity run, no GPU, no allocation. | `suite/cuda/pdp_ps5_*` on `exp/pdp-ps5-gpu-emulation` (cryptanalysis). | Compile the emulation, replay the PS-0 corpus, compare digests. | A digest mismatch. | **stopped early**, device run blocked (no GPU, no owner-approved cap) |
| L8 | PS-6 should price the table (L1) and the kernels in one counted unit, memory included, and let the ratio to the floor say whether the slope moves. | PLAN.md PS-6 finding rule; L1. | After L1 and PS-3 land: an `ecbench-extend` oracle for the kernels beside the existing `Mitm` rows. | The slope of `S/floor` unchanged. | open |

## Not leads (checked and set aside)

* A Gray-code incremental quartic update: the group-law formulation has no quartic to update; the bound on the baseline's own formulation
  is 1.21…1.29× in carry-less multiplies (`ps2_scalar/README.md` §6, arithmetic, not a measurement).  Not built.
* Extending to `ℓ = 13…16` for the table: memory grows as `|B|²` (about 16 M keys at `ℓ = 13`, 64 M at `ℓ = 14`); the 16 GB host fits `ℓ = 14` at
  most.  Whether larger `ℓ` is worth a bigger host is for the owner; this program's open question (`ℓ ≈ 44`) is out of reach for any table.
