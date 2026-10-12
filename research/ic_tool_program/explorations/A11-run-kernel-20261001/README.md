# A11 exploration: a run-based least-rotation kernel, rejected before declaration

**What this is.** An exploratory stage timing for the plan's A11, the
scan's canonical key (`research/notes/index-calculus/IC_TOOL_PROGRAM.md`
§8). It ran on 2026-10-01, between R02's timed processes, under the
benchmark lock. It tested whether a candidate is worth declaring as a
round. It is not a round, it claims nothing, and the candidate is not
proposed.

## The candidate

The AVX-512 key takes the least rotation of each key's normal-basis
coordinates. `least_rotation16_avx512` does that by `n − 1`
rotate-and-minimum steps on sixteen words.

The candidate, `least_rotation16_runs_avx512`, vectorises the scalar
method `FrobeniusCanon::least_rotation` already uses:
- grow each word's longest zero runs (`Z_k = Z_{k−1} ∧ rotl^{k−1}(¬c)`);
- compare only the rotations that start at a longest run's top.

It keeps the exhaustive kernel behind `KIC_CANON_RUNS=0`.

The code is in `a11-run-kernel.patch`, on `main` at `465f0088`:
- the kernel;
- its dispatch;
- a test `the_avx512_run_kernel_matches_the_exhaustive_kernel`;
- the example `canon_kernel_bench`.

**The keys are identical.** The test compares the two kernels and the
scalar search on every degree 1–63, on random, sparse, dense and periodic
words. It passes, and so do the 14 other `koblitz_fast` tests.

## The result: slower, so rejected

Each figure is one isolated process (`tools/isolated_bench.py run --cpus 2`,
uncontended). Each process times `canon_in_place` on 2^20 fresh words,
the best of seven passes, in ns a key. Each pair is ABAB. The checksums
of every key agree between the kernels.

| degree | run kernel | exhaustive kernel |
|--:|--:|--:|
| 41 | 8.81, 11.72 | 5.61, 7.07 |
| 53 | 11.15, 12.02 | 8.63, 8.57 |
| 61 | 9.93, 10.18 | 10.19, 10.16 |

**Why it loses.** The run kernel does fewer vector operations in principle.
But its two loops have data-dependent trip counts, and each step runs
through a mask that the next step needs. So it is latency-bound and
mispredicts its exits. The exhaustive kernel is a straight-line,
throughput-bound loop of about 26 vector operations a key at `n = 53`.
It wins at small degrees and ties at 61.

## What it prices instead: the key and the addition, natively

These are stage figures outside the scan, warm. They are not the scan's
own costs.
- **The exhaustive AVX-512 key** is 5.6–10.2 ns a key. Most of that is
  the rotation loop, which grows with `n`.
- **`examples/scan_block_bench.rs`** at `n = 53` (outputs `runs/sb-*`):
  - the 8-lane batched subtraction (`add_many_lazy`) is 8.3–8.4 ns a
    summand, and the scalar path 21–31 ns;
  - the bulk key over 2M fresh abscissae, copy included, is 15.5–15.7 ns.
- **The scan as a whole.** A scanned summand costs about 40 ns at that
  size (R01: 1.49 units). So the key and the addition are about half of
  it. The rest is the filter probe, the admitted candidates and the
  loop's bookkeeping, and is unpriced natively. Pricing it inside the
  scan, with probes compiled into a separate binary, is the next step
  for Track A's scan work. It needs a declared round.

## Files

- `a11-run-kernel.patch`: the code (`git am` on `465f0088`).
- `runs/`: per process, the output, the stderr and the isolation record.
  - `runs-*` and `full-*`: `canon_kernel_bench 41 53 61`, run with the
    run kernel and with `KIC_CANON_RUNS=0`.
  - `sb-simd-*` and `sb-scalar-*`: `scan_block_bench` (`n = 53`), run
    with the defaults and with `KIC_SCAN_SIMD=0`. Both arms ran with
    `KIC_CANON_RUNS=0`, so the exhaustive kernel was used.
- `binaries.sha256`: the two example binaries, built from the patch with
  `cargo build --release --example …`.

The host is R02's (`research/ic_tool_program/rounds/R02-wide-tail-kernel/`,
`runs/host.json` once committed). It is one x86-64 cloud container with
AVX-512 and has no hardware counters.
