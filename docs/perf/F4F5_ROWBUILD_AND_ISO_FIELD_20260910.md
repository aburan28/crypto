# F4/F5 Macaulay row build and isogeny field: measured negative results

Campaign record, 2026-10-10, Apple M4 Pro (aarch64), rustc 1.93.0,
release profile, single-threaded workloads, interleaved paired A/B with
alternating arm order and A/A noise controls. All wall-clock runs were
**contended** (a shared host running parallel agent sessions at 60–120
load average); each verdict below survived its A/A spread by more than
the stated margin, and every decision-level output was checked identical
between arms before any timing was read.

Host manifest: Apple M4 Pro, 14 logical CPUs (10P+4E), macOS,
`rustc 1.93.0`, release profile, `RAYON_NUM_THREADS=1` for the descent
engines (built into the harness). `tools/isolated_bench.py` cannot run
on darwin (procfs + `os.sched_getaffinity`); interleaving, alternating
arm order and A/A controls stood in, and every number here is labeled
contended. The campaign's follow-up on a quiet host or isolab should
confirm the two surviving margins before they are quoted as clean.

## What was tried, and what each measurement said

| # | change | verdict | evidence |
|:--|:--|:--|:--|
| 1 | shift-budgeted insertion sort for the Macaulay per-multiplier row build (`visit_macaulay_rows_counted`) | **neutral** (kept out) | `ic descent --cells 19:10:2 --targets 16 --repeats 5 --solver matrix-f5`: interleaved A/B medians 5.74 s (base) vs 5.84 s (cand), inside the 3.2 % A/A arm delta; `f4-f2` split over both orders disagreed in sign (6.18→6.78 one order, 6.22→6.02 swapped), stage splits (`build_ms`, `eliminate_ms`) moved together, i.e. load noise. Instrumented row statistics: 813 384 rows, mean 35.0 elements, 23.1 shifts/row, 8.6 % of rows over budget — nearly-sorted as profiled (2.3 adjacent breaks/row) but insertion repair ≈ quicksort in cost. |
| 2 | single decorated `(Reverse(mono_key), mask)` sort in `macaulay_columns` replacing the u64 sort + keyed sort | **regression 25–45 %** (reverted) | same workload, candidate 8.82 s vs baseline 6.03 s median (46 %) under lighter load, 26 % under heavier load, both arm orders; row-only variant (change 1 alone) neutral, isolating the tuple sort. 24-byte tuple swaps lose to two 8-byte sorts. |
| 3 | dedicated SOS squaring + HAC 14.32 reduction for `isogeny_walk::field::Field::sqr` (was `mul(a,a)`) | **regression 19 % walk / 1.5–1.8× the multiply** (reverted) | `isogeny_walk walk --curve p256 --max-ell 13 --max-curves 1200`, 6 rounds alternating order: base median 12.6 s vs candidate 15.6 s, B slower in every round in both positions. Microbench (chained): mul 11.7 ms vs sqr 17.2 ms. Correctness was pinned (6000 randomized + edge inputs equal to `mul(a,a)` on P-256/P-224/P-192; walk records byte-identical) — the loss is codegen, not math: the u128 carry cascades serialize where the CIOS loop pipelines. |

## What landed

- `field_ec/isogeny_field_p256_mul_sqr_pow`: the isogeny-walk field had
  **no perfbench kernel**, so changes to it were invisible to the
  performance index. It now has one: 1 048 576 mul + sqr and 8 pow over
  the P-256-prime four-limb CIOS Montgomery field, fingerprinted
  (`67847ffe72d63c5c` at median 30.4 ms on this host, contended).
  This is the handle any future field change must move.

## What the campaign's profiling found (for the next iteration)

- `ic descent` (matrix-f5, cell 19:10:2) splits roughly: elimination
  (gf2_elim tables) ≈ 36 % of samples, per-multiplier row sort ≈ 15 %,
  `F4ColumnIndex::get` ≈ 5 %, row generation ≈ 5 %; the elimination inner
  loops are already Gray-code table based with x86 SIMD paths; aarch64
  runs the portable unrolled XOR (auto-vectorizing).
- The isogeny walk is ≈ 67 % `Field::mul`; `sqr` is on the point-doubling
  path (5–6 `sqr` per doubling), not just `pow`. A faster squaring must
  beat the CIOS codegen, which this record shows is not free.
- `pow` (5.6 % of walk samples) and `inv` (via `pow(p−2)`) are the next
  candidates: batch inversion where inversions cluster, and windowed
  exponentiation, are algorithmic changes that need their own measured
  iteration before landing.

## Rules applied

- Every candidate was proven decision-identical before timing was read:
  engine tests, `ic descent` four-engine report diff (decision fields),
  and walk-record diff against the baseline binary.
- No speedup is claimed from any of the three tried changes. The e2e
  `ic bench --koblitz-degree 17 ... descent-algebraic + f4-f2` fixture
  (S = 38 564, verified) is the pipeline measurement any future gain
  must move, alongside the frozen WDSat suite for the SAT solver stage
  (untouched here: these changes are in the F4/F5 engines, a different
  solver family; the paired-engine `ic descent` harness is the matched
  protocol for them).
