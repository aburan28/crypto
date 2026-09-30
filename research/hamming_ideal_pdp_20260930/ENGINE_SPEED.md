# Engine wall time, same Boolean F4

Stage diagnostic of the Hamming-ideal PDP solver. The counted unit of the
protocol (`xor_words`, calls, tame/wild) is unchanged. Wall time is the
practicality note in `PROTOCOL.md`. This is not an `S` measurement, not a
ratio to rho, and not an index-calculus scoreboard row.

**Class: engineering.** The elimination schedule is the same computation
done with less overhead. Nothing moved relative to a counting floor.

## What changed

In `hamming_pdp/src/f4.rs` and `boolpoly.rs`:

- Gebauer–Möller selection keeps the earliest partner for each lcm, then
  drops an lcm only when a lower-degree division-minimal candidate divides
  it. Same-degree lcms are an antichain, so they are not probed as
  divisors. The scratch map is reused across basis elements.
- A propagation restart whose input already contains a linear skips S-pair
  construction. `drain_rows` discards those pairs; field-pair rows,
  redundant-element rows, and retirement still run, because they decide
  which polynomials the restart keeps and in which order.
- Monomial degree, order, and divisibility touch only the limbs the input
  system uses. A guard restores the full width afterwards.
- Pivot substitution tests variables with a bit mask, and the accumulator
  merge reuses one buffer instead of allocating a vector per term.
- Monomial maps use a multiply-add hasher. The keys are monomials this
  process built.

`d`, the XOR budget, the matrix cap, the seeds, and the branching order
are untouched. Frozen `results/main/` is untouched.

## Paired check

Host: Linux 6.12.94+, 4 cores, Intel Xeon (advertised AVX-512), `rustc`
1.83.0. One thread. ABAB, three rounds, encoding `C`, seed 0, the frozen
budget (`2^28` row-word XORs) and matrix cap (`2^26` words). `xor_words`
agreed on every run.

| n | binary | solve median | wall median | xor_words |
|:--|:--|--:|--:|--:|
| 11 | baseline `9f8a1af18` | 17.544 s | 17.560 s | 1233543909 |
| 11 | this engine | 9.275 s | 9.287 s | 1233543909 |

`solve` is the cell's `solve_us` (system build and the exhaustive oracle
sit outside it; together they are a few tens of milliseconds). Median of
three rounds. Ratio of medians: `17.544 / 9.275 = 1.89`.

Phase medians on that cell, milliseconds:

| phase | baseline | this engine |
|:--|--:|--:|
| prep | 1349 | 1436 |
| elim | 1458 | 1576 |
| add | 5708 | 1630 |
| subst | 4915 | 1779 |
| post | 2293 | 1317 |

`n = 7`, encoding `C`, seed 0, one pair on the same host (not interleaved):
solve `0.641 s` to `0.376 s`, `xor_words = 39404234` both.

Correctness on this binary: `hamming_pdp --selftest` seeds 1, 7, and 13
(300 systems each), and `replay_check.py --n 7 --encodings SUB,FC --targets 4`
against `results/main` (zero differences). The `n = 11` seed-0 cell above
matches the frozen raw row on every replay key, including `xor_words`.
