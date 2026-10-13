# PS-0: the frozen baseline, host manifest and verdict corpus

Protocol: [`../protocols/PS-0.md`](../protocols/PS-0.md) (frozen by the plan, #1657).
Plan: [`../PLAN.md`](../PLAN.md).  Class: **accounting** (nothing here is faster
than anything; it fixes what every later task is compared with).
Run report with graphs and PDF: [`../reports/run-20261011T2107Z/`](../reports/run-20261011T2107Z/).
Conductor task T-4.

Everything below is an **observation** of one run of
[`run_ps0.sh`](run_ps0.sh) on one host, or arithmetic on those observations
that is shown with its inputs.  No speedup is claimed.

## 1. What was requested and what ran

| PS-0 requirement | status | evidence |
|:--|:--|:--|
| Reference: `semaev_decomp::decompose` (one-word field, `gcd(q, L_V mod q)`), built from `main`, binary kept | **partial** | built from `main` @ `4b6e09b78` plus this PR's diff ([`manifest.json`](manifest.json): `worktree_dirty`); the binary's SHA-256 is recorded, the 2.4 MB binary itself is not committed (it is in the container's `/tmp/ps0_bin` only and dies with it) |
| Counting instrumentation is a separate build; if it changes the timed code, both builds are recorded | **verified complete** | feature `pdp-opcount` (counters in `semaev_decomp::opcount`); two builds, two hashes in `manifest.json`; §7 checks the feature-off build against unedited `main` |
| Rungs `l ∈ {8,…,12}` with the bench's own moduli, `n = 3l` | **verified complete** | [`curves.json`](curves.json), §2 |
| Curve named by its ICV1 slug once PS-0 records it | **verified complete** | §2, registered in `docs/curves/registry.json` in this PR |
| Corpus per rung: 64 uniform + 16 planted targets, SplitMix64 seed `20261011 + l` | **verified complete**, with the clarification in §3 | [`corpus/`](corpus/) |
| Every target decided exhaustively; every decomposition re-added in the group | **verified complete** | §3 and §4 |
| Corpus stores targets, verdicts and full root sets | **verified complete** | [`corpus/corpus_l*.txt`](corpus/) |
| Per pair: instructions retired under `callgrind` at `l = 8, 9` | **verified complete** | §6, [`callgrind/`](callgrind/) |
| Counted field multiplications, squarings, inversions | **verified complete** at `l = 8…12` | §6, [`opcount/counts.jsonl`](opcount/counts.jsonl) |
| ns per pair, single thread, reserved core (`tools/isolated_bench.py`), contended runs reported | **verified complete** | §5, [`raw/`](raw/) |
| A/A (the build against a copy of itself) and ABAB, at least five rounds | **verified complete** | §5 |
| Success: A/A spread recorded; digests of verdict and root-set stream identical in three runs | **verified complete** | §4, §5 |
| ecbench: record which arms already exist | **verified complete** | §8 |

## 2. Host, build, and the five curves

Host ([`manifest.json`](manifest.json), recorded 2026-10-11): Intel Xeon @ 2.80 GHz, KVM guest, 4 vCPU,
one thread per core, no `cpufreq` control, no cgroup CPU limit; `lscpu`: L1d 128 KiB (4 instances), L2 4 MiB
(4 instances), L3 33 MiB; flags `avx2`, `avx512{f,bw,cd,dq,vl,vnni}`, `bmi1/2`, `fma`, `pclmulqdq`, `adx`, `aes`;
**not** reported: `vpclmulqdq`, `gfni`, `avx512vbmi`.  16 GiB, Linux 6.18.44.  `rustc 1.97.0`, default `release`
profile (the repository sets none), `valgrind 3.22.0`.  No GPU.  The field kernel is `pclmulqdq`
(`Gf2::kernel_name()`, recorded in every timing row).

Every timed row ran on CPU 3 through `tools/isolated_bench.py run --cpus 3`.

The curve is `y² + xy = x³ + x² + 1` (a Koblitz-family member) over `F_{2^n}`, `n = 3l`, with the modulus
the bench picks (smallest trinomial, else pentanomial), factor base `V = ⟨1, z, …, z^{l−1}⟩`.  Group orders
come from the Koblitz trace recurrence; the five models are registered from [`curves.json`](curves.json):

| l | n | modulus (hex) | ICV1 slug |
|--:|--:|:--|:--|
| 8 | 24 | `0x100001b` | `icv1-f2m24-tm5983-28fa8b30` |
| 9 | 27 | `0x8000027` | `icv1-f2m27-t7555-cb6a4c9d` |
| 10 | 30 | `0x40000003` | `icv1-f2m30-t10089-8ab21610` |
| 11 | 33 | `0x200000401` | `icv1-f2m33-tm110885-f5e5693b` |
| 12 | 36 | `0x1000000201` | `icv1-f2m36-t473713-f90ff4d7` |

## 3. The corpus, and two clarifications of the protocol

Per rung: 64 *uniform* and 16 *planted* targets from one SplitMix64 stream seeded `20261011 + l`.

* **Uniform** means: a uniformly random `n`-bit word, resampled until it is the abscissa of a rational point,
  with one of its two points chosen at random.  The protocol's "re-added in the group" needs a point `R`;
  the earlier indicative bench used arbitrary words, whose "decompositions" are not points.
* **Planted** means `R = P₁ + P₂ + P₃` for three base points with abscissae in `V`.

| l | n | targets | found | refuted | roots | roots that re-add to ±R in the group |
|--:|--:|--:|--:|--:|--:|--:|
| 8 | 24 | 80 | 25 | 55 | 32 | 32 |
| 9 | 27 | 80 | 28 | 52 | 35 | 35 |
| 10 | 30 | 80 | 25 | 55 | 30 | 30 |
| 11 | 33 | 80 | 25 | 55 | 26 | 26 |
| 12 | 36 | 80 | 23 | 57 | 30 | 30 |

*Roots* are the sorted triples `(x₁ ≤ x₂ ≤ x₃)` in `V³` with `S₄(x₁,x₂,x₃,x_R) = 0`: the full root set of each
target, enumerated without early exit.  *Re-add*: every `xᵢ` lifts to a rational point and some choice of
lifts sums to `±R`, computed with the repository's `binary_ecc::curve::point_add`.  **All 153 roots over the
five rungs are rational decompositions**; no root with a non-liftable abscissa occurs.  (§ PS-1 will explain
why this is expected.)

## 4. Verdict and root-set digests

`verify` replays the baseline `decompose` on every target, requires its verdict to equal the corpus verdict
and its witness to lie in the corpus root set and vanish under `eval_f3`, and recomputes both digests.
Three replays per rung ([`digests.txt`](digests.txt)); each corpus was also regenerated from scratch by a
second run and compared byte for byte (`raw/corpus_gen.log`).

| l | SHA-256 of the verdict stream | SHA-256 of the root-set stream |
|--:|:--|:--|
| 8 | `8e480bd46900c32d2b99ad90cb7aad0f316607dce585aa09e789cea28c8590c0` | `4710adfa9915930971cc01b9987dc38cda8d4bb2bb1be7f1b6801391d1d80374` |
| 9 | `f62c8ab5c1ee5e59e594db9d494536280d6104b4efce704d551e24ed60c0b180` | `ee7bc0a7cafa25ced9ed3101e2cd8fee32899d7e8031226079bc7bf36c207661` |
| 10 | `442e2e80cd0f103b256ac450977b644ea3140ecf1c1a35dc0c7c6b45763c545a` | `6880bbd69bf75aed9f917d62666f1dbf9158c7bf7e4a1d8390fd0904e86c634b` |
| 11 | `eecca3f437196213c68f929859621a554750a8a37ccf3d6186bd08449cecd2c6` | `598d8acb33f5663491db847d443311e176b3790304f04f679d93778c7662daea` |
| 12 | `912b355f5161b92c1b6945befa40d72dfdc07963f9ee2ed9a6a54fd598a784ee` | `6bdcaa39de35e197d6dc104db105497b4199f14ee23d61bbd4bf09fb5a68d3ed` |

The verdict stream is, per target in corpus order, `[planted u8][x_R u64 LE][found u8]`; the root-set stream
appends `[count u32 LE]` and the sorted triples (`u64 LE` each).  All three replays agree on both digests at
every rung, and the replayed verdict digest (built from what the baseline said) equals the corpus digest.

**Independence of the root sets.**  The root sets come from `corpus.rs::all_roots`, which shares only the
quartic's coefficients (`quartic_in_x3`) and the field with the baseline: it normalises each quartic,
forms `t^{2^i} mod f` by generic polynomial squaring, takes a generic Euclidean gcd with `L_V`, and falls back
to evaluating over `V`.  Against a brute-force walk of every sorted triple with the crate's general
`Vec`-backed `S₄` evaluator (`binary_semaev_s4`) the root sets agree on **all 80 targets at `l = 8`** and on
the **first 16 targets at `l = 9`** (`raw/brute_l8.log`, `raw/brute_l9_first16.log`: 0 mismatches).
Brute force at `l ≥ 10` was not run (the cost grows by 8× per rung).

## 5. Timing: A/A with ABAB interleaving

Unit: **wall nanoseconds per pair**, single thread.  The timed region is `decompose` over the first `K`
refuted targets of the corpus, `P` passes (`K, P` = 16×40, 16×10, 16×3, 16×1, 6×1 for `l = 8…12`).  A refutation visits
`span·(span+1)/2` pairs exactly, so `ns/pair = secs / (K·P·span(span+1)/2)`; the one-time table builds
inside `decompose` (`subspace_poly`, Horner coefficients) are inside the timed region and are part of the unit.
One untimed call precedes the region.  Arms: `A` and `B`, **the same build copied under two names**, alternating
`A B A B …`, five rounds per rung (50 timed runs).

| l | rounds | contended runs | A median ns/pair | B median ns/pair | A run-to-run spread (max−min)/median | B/A median | B/A per-round min…max |
|--:|--:|--:|--:|--:|--:|--:|:--|
| 8 | 5 | 0 | 337.7 | 337.0 | 3.24% | 1.0085 | 0.9918…1.0206 |
| 9 | 5 | 0 | 375.4 | 369.7 | **8.81%** | 0.9847 | 0.9364…1.0008 |
| 10 | 5 | 0 | 393.3 | 393.0 | 4.46% | 0.9986 | 0.9783…1.0285 |
| 11 | 5 | 0 | 421.4 | 432.3 | 3.22% | 1.0236 | 0.9717…1.0564 |
| 12 | 5 | 0 | 452.3 | 453.9 | 1.58% | 0.9979 | 0.9917…1.0176 |

([`aa_summary.md`](aa_summary.md) carries the bootstrap interval too; with five rounds it equals the range
and is not independent information.)

* **A/A spread.**  The largest run-to-run spread of one arm is 8.81% (`l = 9`); the largest single-round
  departure of B/A from 1 is 6.4% (`l = 9`, 0.9364).  Both are below the protocol's 10% line, so the
  gates of later tasks are not widened; **a later A/B at rung `l` must be read against that rung's row**.
* **Contention.**  0 of 50 timed runs were marked contended.  One *start* was refused by the tool's pre-flight
  (my own polling used 0.21 CPU-s against a 0.20 limit during the 2 s settle window, before round 4 of `l = 12`);
  the run did not execute, produced no record, and the script was resumed with start retries (`raw/run_ps0_main.log`,
  `raw/run_ps0_resume.log`; the resumed invocation needed no further refusals).
* **Derived, per refuting target** (median A ns/pair × pairs, arithmetic only): `l = 8`: 32,896 pairs, 11.1 ms;
  `l = 9`: 131,328, 49.3 ms; `l = 10`: 524,800, 206 ms; `l = 11`: 2,098,176, 884 ms; `l = 12`: 8,390,656, 3.80 s.
* **Against the plan's indicative numbers** ([`../PLAN.md`](../PLAN.md) §2: 244, 259, 355, 310 ns/pair at `l = 8…11`),
  these medians are 1.4×, 1.4×, 1.1× and 1.4× larger.  The indicative run used arbitrary words as targets, was
  unpinned, and ran once; this run uses curve points, a reserved core and five rounds.  I did not isolate which
  difference accounts for it.  PS-0 replaces the indicative numbers, as the plan says.

## 6. Counted operations per pair

**Instructions retired** (`valgrind --tool=callgrind --toggle-collect='*run_one*'`, one refuted target, collected
only inside the call): `l = 8`: 77,504,893 over 32,896 pairs = **2356 per pair**; `l = 9`: 337,453,547 over 131,328
pairs = **2570 per pair** (arithmetic: Ir / pairs; includes the per-call table builds).

**Field operations** (the `pdp-opcount` build, 2 refuted targets per rung, every call of `Gf2::{mul, sqr, sqr_k,
clmul, sqr_wide, reduce, inv, batch_inv}` counted once at its entry; the folding path only, which the build
asserts).  Per pair:

| l | `mul` | `sqr` | stand-alone products (`clmul`) | squares before reduction (`sqr_wide`) | stand-alone `reduce` | inversions per pair | batch-inverted elements (3 products each) | PCLMULQDQ instructions |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 8 | 28.36 | 4.03 | 74 | 26 | 50 | 0.0078 | 1 | 306.6 |
| 9 | 28.18 | 4.01 | 82 | 30 | 56 | 0.0039 | 1 | 329.8 |
| 10 | 28.09 | 4.01 | 90 | 34 | 62 | 0.0020 | 1 | 353.4 |
| 11 | 28.04 | 4.00 | 98 | 38 | 68 | 0.0010 | 1 | 377.2 |
| 12 | 28.02 | 4.00 | 106 | 42 | 74 | 0.0005 | 1 | 401.1 |

The PCLMULQDQ column is derived: 3 per folded `mul`/`sqr`/`sqr_k` step and per batch-inverted element's 3
products, 1 per stand-alone `clmul`/`sqr_wide`, 2 per stand-alone `reduce`.  The growth with `l` is 8 stand-alone
products, 4 squares and 6 reductions per extra rung (one more step of the Horner loop on `L_V`).  The `inv` counts
are the single-element inversion at the end of each row's batch.  Counts are call counts, not a proof of the
instruction mix the compiler emitted; the `Ir` figures above are the retired-instruction cross-check.

## 7. Does the instrumentation touch the timed code?

The lib edits are: the `opc!` macro (expands to nothing without the feature), the `opcount` module
(`#[cfg(feature = "pdp-opcount")]`), `opc!(…)` statements at the head of eight `Gf2` methods, and `Gf2::folds`
widened from `pub(crate)` to `pub`.  The timing build has the feature off.  PS-0 built the harness a second time
against the unedited `main` sources and compared instructions retired:

| build | `l = 8`, instructions retired in `run_one` | `l = 9` |
|:--|--:|--:|
| harness against the unedited `main` sources (`target/pristine`, `pdp_speedup` SHA-256 `2244cb06…`) | 77,504,893 | 337,453,547 |
| harness against this PR's sources, feature off (the timing build, `9fa12715…`) | 77,504,893 | 337,453,547 |

The two builds retire the **same number of instructions to the instruction** in the timed function at both
rungs.  (Their binaries differ in hash; the hash covers the whole executable, not the timed code.)  This
is evidence that the feature-off lib code generates the same loop; it is not a proof for every call site
the harness does not exercise.

## 8. ecbench arms that already exist (for PS-6)

`ecbench methods` (this checkout) lists: `rho.*` (six variants), `bsgs.*` (three), `kangaroo.vow`,
`claw.pair_table`, `ic.pipeline`, `ic.pipeline_counted`, `ic.large_prime`, `ic.shared_rank`.  The IC arms take an oracle plug-in
(`subtract`, `mitm`, `mitm-frobenius`, `mitm-frobenius-counted`, `descent-algebraic` with `m = 2|3`, the fiber-aware
variants, and `pdp3-koblitz:m=3,engine=inherited-f4|f6-ic`) and a factor-base plug-in that includes
`binary-subspace:dimension=D`.  **No ecbench oracle prices the `S₄` pairs-and-solve loop of
`SubspaceOracle`/`decompose` on a binary subspace**; that loop is priced in the ledger code
(`ic_boundary.rs`, `Oracle::SemaevS4`, `calibrate_s4`) and `ic_oracle_pricing.rs`.  PS-6 therefore adds an
oracle plug-in through the `ecbench-extend` skill.

## 9. Corrections to the plan found while running PS-0

* PLAN §6 says "L2 4 MiB per core".  `lscpu` reports L2 **4 MiB across 4 instances** (1 MiB each) and L1d
  128 KiB across 4 (32 KiB each).
* PLAN §1 and the protocol count `2^{2ℓ−1}` pairs; the loop visits `2^l(2^l+1)/2 = 2^{2l−1} + 2^{l−1}` pairs.
  This PR uses the exact count (§5).

## 10. Not claimed

No kernel is faster than another here.  Nothing transfers to ECC2K-130; the `AGENTS.md` §8a gate is not
touched.  Wall figures are valid for this host class (x86-64 Xeon with AVX-512F/BW, KVM guest) and are
`isolated_bench` runs on a shared container, which the repository's CPU isolation gate does not accept as a
receipt: they are the A/A-qualified reference for this program's A/B comparisons, not a publishable timing.

## Reproduce

```bash
research/pdp_oracle_speedup_20261011/baseline_v1/run_ps0.sh      # on a quiet machine
PS0_STAGES=timing research/pdp_oracle_speedup_20261011/baseline_v1/run_ps0.sh   # resume only the timing stage
```
