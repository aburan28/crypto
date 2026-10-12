# PS-2: the scalar group-law kernel, its gate count and its program

Protocol: [`../protocols/PS-2.md`](../protocols/PS-2.md).  Plan: [`../PLAN.md`](../PLAN.md).  Baseline: [`../baseline_v1/`](../baseline_v1/README.md).
Formulation and exactness: [`../ps1_algebra/`](../ps1_algebra/README.md).  Run report with graphs and PDF:
[`../reports/run-20261011T2255Z/`](../reports/run-20261011T2255Z/).  Class: **engineering**.  Conductor task T-4.
Commands: [`run_ps2.sh`](run_ps2.sh); raw records: [`raw/`](raw/).

## 1. Requested and executed

| PS-2 requirement | status | evidence |
|:--|:--|:--|
| A new module `semaev_gray_pairs` in `src/cryptanalysis/` with a scalar kernel and a gate-count export (program as data, hashed) | **partial**: the kernel and the export exist and are tested, but in `examples/pdp_speedup/` (`kernel.rs`, `bitslice.rs`, `bspairs.rs`, `gatecount.rs`), **not** as the library module the protocol names.  Original statement: a module in `src/cryptanalysis/`.  Reason: development in the example crate rebuilds in about a minute against minutes for the library, and a library edit triggers the repository's frozen-replay workflows; the module is to be promoted when PS-3 settles the interface (PS-3b needs it as a library module in the other repository) | §2, §5 |
| Digests equal PS-0's on every rung | **verified complete**: verdict digest equal on all five rungs for all four arms, and the group kernel's full root set equal on all 400 corpus targets | §3 |
| Scalar ns per pair at most half the baseline's | **verified complete**: 0.12, 0.099, 0.086, 0.071, 0.067 of the baseline's at `l = 8…12` (A/A spread at most 8.15%) | §4 |
| `G` reported with its derivation | **verified complete** | §5 |
| Stop: `G > 1.8·10^5` per pair | **not triggered** (`G` is `1.3·10^4` per task at `n = 33`, batch of 8) | §5 |
| Setup (table build) reported if outside the timed region | **verified complete**: 0.07 ms (`l = 8`) to 1.5 ms (`l = 12`), one table per run | §4 |
| Same compiler profile as the baseline | **verified complete**: one binary, one profile, all arms | §4 |
| The plan's hypothesis as written: *incremental Gray-code coefficient updates* plus the cheaper root test | **partial**: the incremental update is **not built** and is moot in the chosen formulation (no quartic is formed); a bound on what it could give is in §6 | §6 |

## 2. What the scalar kernel is

`kernel.rs`: for each base point `P₁` (a row) the two points `U± = P₁ ± R` are computed once; for every `P₂` of the row
the four candidates `x(U± ± P₂)` come from two batched inversions, four products and four squarings
([`ps1_algebra/README.md`](../ps1_algebra/README.md) §3), and a pair hits when a candidate is below `2^l`.  A zero
denominator or a hit routes the pair to the exact reference (`group::thirds`), which also confirms every hit.  A
target with `x_R ∈ V ∪ V⁻¹ ∪ {0}` is decided by the baseline inside the same call.  Arms of the comparison (one binary, `pdp_speedup time --arm`):

| arm | what it is | candidate set | pairs per refutation |
|:--|:--|:--|:--|
| `baseline` | `semaev_decomp::decompose` (PS-0's reference) | every abscissa of `V` | `2^l(2^l + 1)/2` |
| `so-full` | `SubspaceOracle` over all of `V`: the baseline algorithm through the general-`b` code path | every abscissa of `V` | same |
| `so-liftable` | the same oracle over the **liftable** abscissae of `V` only: the baseline's root test on the group-law kernel's candidate set (**matched control**) | `\|B\|` base points | `\|B\|(\|B\|+1)/2` |
| `group-scalar` | this task's kernel | `\|B\|` base points | same as `so-liftable` |

`\|B\|` is 128, 253, 464, 1007, 1992 at `l = 8…12` (about half of `2^l`: half the abscissae have no rational point).

## 3. Identical outputs

`verify-arm` replays every arm on the whole PS-0 corpus (`raw/verify_arms.jsonl`, 20 records, 0 mismatches).  For each target the arm's
verdict must equal the corpus verdict and the rebuilt verdict digest must equal the corpus digest; the group kernel additionally lists
every root and the sorted set must equal the corpus root set.  80 targets at each of `l = 8…12`.  The digests are PS-0's
([`../baseline_v1/README.md`](../baseline_v1/README.md) §4).

## 4. Timing (interleaved rounds, pinned core)

Unit: wall nanoseconds on one reserved core (`tools/isolated_bench.py run --cpus 3`), **per refuted target** and **per own pair** (an arm's pairs
are its own, so the two columns differ by the pair count).  Five rounds; each round runs `baseline`, `baseline2` (the baseline again: the A/A),
`so-full`, `so-liftable`, `group-scalar` once, in that order, at `l = 8…12` (125 runs).  The corpus targets are PS-0's: the first 16 refuted
targets (6 at `l = 12`); the baseline arms use PS-0's frozen pass counts, the others repeat passes until at least 3 s.  Medians over the five rounds
([`ab_summary.md`](ab_summary.md) has every min…max):

| `l` | baseline, ms/target | `so-full` / baseline | `so-liftable` / baseline | **`group-scalar`, ms/target** | **`group-scalar` / baseline** | baseline / `group-scalar` |
|--:|--:|--:|--:|--:|--:|--:|
| 8 | 10.84 | 1.034 | 0.261 | 0.337 | 0.0308 | 32.5 |
| 9 | 46.99 | 1.026 | 0.250 | 1.137 | 0.0242 | 41.3 |
| 10 | 204.2 | 1.000 | 0.212 | 3.596 | 0.0177 | 56.5 |
| 11 | 876.0 | 1.010 | 0.247 | 15.12 | 0.0173 | 57.8 |
| 12 | 3795 | 1.018 | 0.242 | 60.06 | 0.0158 | 63.3 |

| `l` | baseline, ns/pair | `so-liftable`, ns/pair | **`group-scalar`, ns/pair** | `group-scalar` / baseline per pair | pairs: baseline | pairs: `\|B\|(\|B\|+1)/2` | pair-count ratio |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 8 | 329.5 | 340.1 | 40.8 | 0.1225 | 32,896 | 8,256 | 3.98 |
| 9 | 357.8 | 361.8 | 35.4 | 0.0989 | 131,328 | 32,131 | 4.09 |
| 10 | 389.1 | 400.0 | 33.3 | 0.0860 | 524,800 | 107,880 | 4.86 |
| 11 | 417.5 | 428.2 | 29.8 | 0.0713 | 2,098,176 | 507,528 | 4.13 |
| 12 | 452.3 | 460.3 | 30.3 | 0.0668 | 8,390,656 | 1,985,028 | 4.23 |

* **A/A.**  `baseline` against `baseline2`, per target: medians of B/A 1.044, 1.023, 1.004, 1.009, 0.997 at `l = 8…12`; run-to-run
  spread of the baseline 8.15%, 3.44%, 3.53%, 2.77%, 3.43%.  The `l = 8` row is the least stable (it contains the contended runs below); the
  protocol's 10% line is not crossed anywhere.
* **Decomposition of the per-target ratio** (arithmetic): per-target ratio = pair-count ratio × per-pair ratio.  The matched control
  `so-liftable` isolates the first factor: restricting the baseline's own pair loop to the liftable abscissae is 3.8–4.7× faster per target at a per-pair cost 1.02–1.04×
  the baseline's.  The second factor, the group-law kernel against the baseline's gcd test on the same candidates, is `1/0.1225 … 1/0.0668` = 8.2×…15.0× per pair.
* **Setup.**  The base-point table (rational lifts of `V`) is built once per run outside the timed region: median 0.07, 0.15, 0.32, 0.73, 1.5 ms at `l = 8…12`
  (`raw/times.jsonl`, `setup_secs`).
* **Contention.**  6 of 125 timed runs were marked contended by the tool (`baseline2` at `l = 8`, rounds 2, 4, 5; `so-full` at `l = 8`, rounds 2, 4;
  `group-scalar` at `l = 12`, round 5); in each the dominant "other" process was this session's own agent (`claude[96]`, one a `git` commit).
  9 starts were refused by the pre-flight and retried (`raw/refused_starts.log`).  Contended runs are **kept** in every table above;
  [`ab_summary_excl_contended.md`](ab_summary_excl_contended.md) drops each contended round and the conclusions do not change (for example the
  `group-scalar` / baseline ratio at `l = 12` is 0.0155 instead of 0.0158).
* `so-full` against the baseline measures the cost of the general-`b` code path (1.000–1.034×).

## 5. Gates and the program

`G` is the number of gate evaluations of the straight-line program for one task (one pass over all lanes), counted in fused ternary-logic issues.
PS-1 counted it with the inversion amortised over a batch of 64 (7,410 … 14,608 for `l = 8…12`).  This task records the **program itself**: the batch of
**8** tasks, run on the recording lane, dead operations eliminated, hashed (`raw/programs.jsonl`; the program's numbered inputs are, per task, the planes
`xu, yu, xv, yv, x2, y2`, and its outputs the hit and exceptional masks):

| `l` | `n` | live gates, batch of 8 | per task | recorded ops | SHA-256 of the program bytes |
|--:|--:|--:|--:|--:|:--|
| 8 | 24 | 63,630 | 7,954 | 64,782 | `f34227cade15eb2b3c18852ac5de2a25a4e62a1566c5976ac6e3c6a012bdca98` |
| 9 | 27 | 77,092 | 9,637 | 78,388 | `5faf5360ff9d9fb20ff418166c80d51f8731b5bf56149d14535491cd53461a6e` |
| 10 | 30 | 85,718 | 10,715 | 87,158 | `c2a413f9f136c0a3c585920fa4ee40c5d67d39a260576b153db36e2dda0ac43b` |
| 11 | 33 | 103,030 | 12,879 | 104,614 | `099679633c4d2f0c11f4e8b6024b6aa3d1a8703e59980f712593e90fb7065954` |
| 12 | 36 | 124,886 | 15,611 | 126,614 | `a8bec32c3185d52f68cbd313dbff3840febc9d8c4ef99a826ab521dfc3862b9a` |
| 11 (batch of 64) | 33 | 775,254 | 12,113 | 787,926 | `70faa5b6fb2e740a44faa7918885e902daf31697827607b816df37ef6479517f` |

For every program the interpreter (run on 64 lanes of random plane data) equals the same code executed directly on `u64` lanes.  The counting lane
counts 1.6%, 1.8%, 1.9%, 2.7%, 2.8% and 6.1% more gates than the live program has (`l = 12, 11, 11 (64), 8, 9, 10`: the dead gates the recorder drops; the excess depends on the modulus' squaring maps).  `G` per task is at most `1.6·10^4` for `n ≤ 36`, an
order of magnitude below the `1.8·10^5` of the plan's stop condition.  The program is a function of the field's modulus, `l`, the batch
size and the loop orders in `bitslice.rs`; changing any of them changes the hash.  `export-program --out FILE` writes the program bytes and a JSON sidecar
(`pdp-slp/1`) for other back ends.

Counted field operations per pair (instrumented build, 2 refuted targets; PS-0's definitions): `group-scalar` at `l = 11` does 4.06 products, 4.01
squarings, 0.18 squarings inside inversions, 2 batch-inverted elements (6 products) and **42.75 carry-less multiplies** per pair; the baseline's
loop does 377.2 (PS-0), 8.8× more.  Over `l = 8…12` the ratio of carry-less multiplies per pair is 6.5, 7.4, 8.1, 8.8, 9.5
(`raw/counts.jsonl`).

## 6. The plan's incremental-update hypothesis

The plan expected a Gray-code incremental update of the quartic's coefficients to make the scalar kernel faster.  Two measurements bear on it.

1. **The group-law formulation has no quartic**, so the update has nothing to update; the speedup above comes from a different root test and from the liftable candidate set.
2. **What the update could give in the baseline's own formulation.**  Forming the quartic (`quartic_in_x3`) costs 17 products and 6 squarings, **69 carry-less
   multiplies** per pair (`raw/counts.jsonl`, `what: formation`), against the baseline's 306.6, 329.8, 353.4, 377.2, 401.1 per pair at `l = 8…12`: 22.5%, 20.9%, 19.5%,
   18.3%, 17.2%.  If forming the quartic were free, the carry-less-multiply count (a proxy for time, **not a measured time**) would fall by at most that share, a factor
   of at most 1.29, 1.26, 1.24, 1.22, 1.21 (`1/(1 − share)`).  The arm (incremental update plus the baseline's gcd test) was **not built**; the bound is arithmetic on counts.

## 7. Not claimed

An oracle-stage figure, on this host class only (x86-64 Xeon, AVX-512F/BW, PCLMULQDQ, KVM guest; no Arm64, Metal or FPGA measurement), from `isolated_bench` runs on a shared
container (not a receipt of the repository's CPU isolation gate).  Nothing end to end; nothing for ECC2K-130; the `AGENTS.md` §8a gate is untouched.  The kernel decides the rational question
and requires the baseline for the degenerate targets named in PS-1; **no corpus target is degenerate**, so the fallback was exercised only by the exhaustive small-field tests of PS-1.  The 8.2×…15.0× per-pair
figures compare two different algorithms on different candidate counts: use the per-target column for the total and the control for the split.
