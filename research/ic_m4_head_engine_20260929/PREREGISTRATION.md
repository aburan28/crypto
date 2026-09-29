# The m = 4 exponent audit, rerun on the head engine

This file is registered before any audit cell runs on this engine. It is the
"branch-head engine" successor that
[the first audit](../ic_m4_exponent_audit_20260928/PREREGISTRATION.md) §11 names. That
audit returned **closed for this engine at these sizes** on the frozen `2809b498` engine.
Its fit was `ĉ = 0.985`, band [0.983, 1.018]
([its results](../ic_m4_exponent_audit_20260928/RESULTS.md)).

## 1. Question

After `2809b498`, the Gröbner engine gained the Four-Russians kernel, support-local
multipliers, set-bit row readback and new defaults. At `n = 9` it is 3.1× cheaper than
the frozen one, and at `n = 15` 1.43× cheaper (first audit §5.1).

Do those changes move the **slope**, the exponent in `n`? The decision rule, `c* = 0.25`,
asks about the slope, not a constant.

## 2. Engine

- **Build.** `build.sh` builds commit `4ff512f25813cb66c896860890eb405704f8bf00`, which was
  `origin/main` at registration. It already contains the audit binary
  `examples/m4_exponent_audit.rs` (sha256 `5bfe99f0…ffab1cc8`, the same file the first audit
  ran). The lock file is the first audit's `Cargo.lock.pinned` (`469209e8…`), used with
  `--locked` and unchanged.
- **Binaries.**
  - `m4_exponent_audit-4ff512f2`, sha256
    `307bcd969aec401209b8b7571e9ef3237d0aecc26776fbb8cc0ae6a886826545`;
  - `groebner_stage_bench-4ff512f2`, sha256
    `4af0bb1288e491b0d5ec6689ec3cb575bf166b0d0c3b4fb004d9e7ad72978bc0`.
- **Policy.** The first audit's, unchanged: `KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1
  KIC_F4_DROP=complete`, and the dreg caps for the degree arm.

## 3. Reproduction control: PASSED (exact), before registration

The first audit recorded head-engine totals at `1722bad1`, in `baseline/head-1722bad1`
there. Between `1722bad1` and `4ff512f2`, the engine sources changed only by test code and
one helper, in two commits: `ee9e2ec3` and `32199ca0`.

`baseline.sh` re-ran both ladders with both policies on the new bench, pinned to CPU 3.
Every row matches the recorded head totals exactly
(`baseline/REPRODUCTION-4ff512f2.txt`):

| cell | policy | recorded at 1722bad1 | rerun at 4ff512f2 |
|:--|:--|--:|--:|
| `K_0/2^9` m=4 | reference | 595,511 | **595,511** |
| `K_0/2^9` m=4 | candidate | 294,112 | **294,112** |
| `K_1/2^15` m=4 | reference | 83,790,974 | **83,790,974** |
| `K_1/2^15` m=4 | candidate | 27,620,511 | **27,620,511** |

The match covers word XORs, F4 calls, the verdict digest, the decomposed count, matrix
rows and columns, reductions, splits, infeasible branches and exhaustion.

## 4. Cells, arms, seeds, budgets, metric and decision rule

All of these are **the first audit's, unchanged**:

- the nine cells;
- the arms `semaev`, `enumerate`, `null` and `degree`;
- seed `20260928`;
- 16 targets, 8 for `null`;
- a node budget of 20,000;
- the CPU-seconds limits;
- the lower-median metric;
- the bootstrap and the gating.

The readout is the first audit's `analyze.py` (sha256 `7da0d93d…`), and the decision rule
is its §6:

- **alive** if the band's upper end is below 0.25;
- **closed for this engine at these sizes** if its lower end is above 0.25;
- otherwise **inconclusive**.

Word XORs do not depend on the engine's wall time or the host, so the two audits'
per-cell medians are directly comparable.

## 5. Execution and isolation (AGENTS.md §10)

`run_audit.sh` here is the first audit's with only three changes:

- **Binary.** The binary name.
- **Parallel cells.** They run three at a time, each pinned to its own CPU (1, 2, 3); the
  first audit ran four at a time, unpinned. CPU 0 is left to the system.
- **Degree cells.** They run one at a time, pinned to CPU 3.

The whole run holds the benchmark lock (`tools/isolated_bench.py busy`), so no build,
test or other benchmark overlaps it.

The counted units are deterministic, so the effect of pinning is limited to the CPU-seconds
limits: it keeps them, and so censoring, free of time-slicing.

## 6. Predictions

1. **Verdict: closed for this engine at these sizes.**
2. **The head slope is not lower.** `ĉ_head ≥ ĉ_frozen − 0.10 = 0.885`, and `ĉ_head` lies
   in [0.90, 1.40].
   - The reason, stated before the run: the head engine's advantage shrinks from 3.1× at
     `n = 9` to 1.43× at `n = 15`. On the chain-ladder cells that would put the head slope
     about 0.19 higher.
   - This is a two-cell reading on other cells, so the band is wide.
3. **Controls pass.**
   - enumeration null in [0.60, 1.00];
   - random-system null passes;
   - 0 oracle disagreements;
   - 0 group-check failures.
4. **Cheaper at the small end.** At `n = 9`, each head per-cell lower median is below the
   frozen one on both curves.

A head slope below 0.885 would refute prediction 2. It would mean the post-freeze changes
move the exponent, which is worth following at larger `n` even though it is still far from
0.25. It would not reopen the route unless the band's upper end falls below 0.25.

## 7. What it cannot show

Everything in the first audit's §10 applies unchanged:

- no asymptotic statement;
- nothing beyond `n = 19` or `m = 4`;
- nothing about other solver families.

A *closed* verdict here closes `m = 4` for the head engine at these sizes. It does not
close the route.

## Files

| file | what it is |
|:--|:--|
| `build.sh` | builds the head engine under the benchmark lock |
| `baseline.sh` | the reproduction control against the `1722bad1` head totals |
| `baseline/` | the control's outputs and `REPRODUCTION-4ff512f2.txt` |
| `run_audit.sh` | the rerun (pinned) |
| `RESULTS.md` | written after the run |
