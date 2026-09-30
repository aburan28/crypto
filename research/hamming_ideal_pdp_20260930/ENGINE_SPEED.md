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

The `n = 11` table is the engine at monomial width `W = 8` (512 variables).
The width is now 16 limbs so a larger field fits. A later one-pair check at
`W = 16`, same host, same `n = 11` cell, still has `xor_words = 1233543909`
on both engines (baseline solve 25.184 s, this engine 14.148 s). That pair
is not a replacement for the median above.

## `n = 53`, one root call

`n = 83` does not fit this field. `F_{2^n}` is a `u64` polynomial basis, and
the normal-basis inverse is a `u128` row of width `2n`, so `n ≤ 63`. The
binary refuses `n = 83` before it builds a field.

`n = 53` is the degree that was run. It is prime, so `GF(2^53)` has no
proper intermediate subfield over `GF(2)`, and `ord_53(2) = 52`, the same
cyclotomic shape as `n = 83` and `n = 131`. It is not a substitute for the
`n = 83` gate, and it is not a frozen protocol cell. The frozen set stays
`{7, 11, 13, 17, 19}`. Weight `w = 2` matches the `n = 11` cell; the
protocol's `w = 3` for `n ≥ 17` was not used.

The irreducible is `0x20000000000047` and the smallest normal element is
`3`. The group order comes from the Koblitz recurrence
(`#E = 9007199311360364`); enumerating points stops at `n ≤ 20`, which is
every frozen cell. The weight-`≤ 2` factor base is the normal-basis masks
of that weight (849 abscissae, 1697 points), not a scan of `2^53`.
C-Hamming uses 678 variables (572 of them auxiliaries).

The measured cell is encoding `C`, seed 0, `max-calls 1`, budget `2^28`,
matrix cap `2^26` words. One root Gröbner call. It hits the matrix cap,
returns no witness (`agree = timeout`, `found = false`, label true,
`exhaustive_ordered_pairs = 2`). Both engines, `W = 16`, same field and
factor base. The baseline is `f4.rs` and `boolpoly.rs` from `9f8a1af18`
with only the width and the large-`n` paths changed. ABAB, three rounds,
one thread (`taskset -c 0`). `xor_words = 13017543` on every run, with the
same calls, cap hit, rows and columns.

| n | binary | solve median | wall median | xor_words |
|:--|:--|--:|--:|--:|
| 53 | baseline `9f8a1af18` at `W = 16` | 7.134 s | 8.139 s | 13017543 |
| 53 | this engine | 8.426 s | 9.493 s | 13017543 |

Solve rounds, baseline then this engine: 7.678 / 8.426, 7.005 / 8.167,
7.134 / 9.715 seconds. This engine is slower on the call. Pair insertion
is still cheaper (add median 181 ms to 5 ms) and prep is not (4177 ms to
5925 ms). Prep is where the cap stops the call, so the `n = 11` saving in
substitution never starts. Raising the cap to `2^27` words, one pair, stops
at the same `xor_words`: the next matrix crosses both caps in one step.

**Class on this cell: a wall-time regression of the same computation.**
`xor_words` did not move. Nothing was resolved that the baseline left open.
Not an `S` figure, not a rho ratio, not a scoreboard row.
