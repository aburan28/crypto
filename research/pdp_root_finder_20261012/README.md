# Root-finding decomposition oracles for the point decomposition problem, measured on both tracks

2026-10-12.  Study directory `research/pdp_root_finder_20261012/`.  Protocol:
[`PROTOCOL.md`](PROTOCOL.md) (registered before any timed run).  Requirements,
gates and baseline: [`../pdp_root_finder_requirements_20261012/`](../pdp_root_finder_requirements_20261012/README.md)
(PR #1734).  Driver: [`examples/pdp_rootfind/`](../../examples/pdp_rootfind/).
Run report of the timing stage: [`reports/run-20261012T2135Z/`](reports/run-20261012T2135Z/report.pdf).
Digest with the figures: [`REPORT.pdf`](REPORT.pdf).  **Class: engineering**
(AGENTS.md §3): every arm's counted work sits on a generic floor (the pair
table or enumeration); the constants, the exactness instrumentation and the
measured distance to F4 are the result; no exponent moved.

## 0. Summary

The requirements document asked for a root-finder for the PDP (the rational
zeros of one summation polynomial on the grid `V^m`), specialised to the
tracks, exact, counted, and gated against F4/F5 and the general Gröbner
engines.  This round builds the two mechanisms that meet S1–S8 without a
Macaulay step on either track — the memoryless *pairs-and-solve* and the
*2-sum table* (meet in the middle) — measures them on the frozen PS-0
corpora (binary, `ℓ = 8…14`), on a (23, 8) cell the repository's Gröbner
engines can run, and on the Kummer-tower cells of the PKM rounds (prime,
`m = 2, 3, 4`, the committed F4 rows re-read, not re-run), and reads the
gates.

| what | result |
|:--|:--|
| exactness (G0) | binary: verdict and root-set digests equal to the corpus's on every rung `ℓ = 8…14` for both arms (80 targets a rung at `ℓ ≤ 12`, 20 at `ℓ = 13, 14`); prime: 707/707 brute-force agreements on the committed rows, 0 mismatches on every ladder rung, every planted target found, 479/479 finished F4 verdicts equal to the algebraic count, 0 algebraic solutions without a rational decomposition (486 checked) |
| counted work, binary (`m = 3`, `n = 3ℓ`) | `rf-mitm`: `5·|F| − 3` multiplications, `2|F|` squarings, one inversion and `2|F|` lookups a target (637 mul at `ℓ = 8`, 40,922 at `ℓ = 14`); `rf-pairs`: `≈ 10` multiplications a pair, `|F|²/2` pairs (83,069 at `ℓ = 8`, `3.35·10⁸` at `ℓ = 14`) |
| counted work, prime (`m = 3`, Kummer at `p₁`) | `mitm3` 32 → 20,417 mul at `t = 4 → 13`; `pairs3` 305 → `8.34·10⁷`; `mitm4` 495 → `5.95·10⁵` at `t = 4 → 9`; `subtract2` 22 → 10,392 at `t = 4 → 12` |
| against F4 (G2a) | prime, same cells and targets as the committed rows: F4 multiply-adds ÷ root-finder multiplications `4.6·10¹¹` at `m = 3, N = 18` (`8.99·10¹³` vs 197), `1.0·10¹¹` at `m = 4, N = 16` (`2.71·10¹³` vs 270), `1.0·10⁹` at `m = 2, N = 22`; binary (23, 8): inherited F4 `4.1·10⁷` word XORs and 231 ms a target (median of 8) against 592 multiplications and 3.9 µs; the split-dreg F4 did not finish (memory) |
| growth (G2b, G3) | binary: `rf-mitm` 1.000 ± 0.027 bits per rung of `ℓ`, `rf-pairs` 1.995 ± 0.055 (6 rungs); prime: `mitm3` 1.019 ± 0.020 bits per unit of `t`, `pairs3` 1.990 ± 0.026, `subtract2` 1.052 ± 0.089, `mitm4` 2.028 ± 0.114; F4 on the same kind of cells: 3.214 ± 0.063 (`m = 2`), 6.556 ± 0.864 (`m = 3`) bits per unit of `t`.  The root-finders' slopes are the pair-table and enumeration floors (1, 2 bits per rung at `m = 3`; 2 at `m = 4`) to within the interval; **no slope below a floor, so G3 fails and the class is engineering** |
| wall (isolated core, 5 rounds) | binary, per refuted target (median of the uncontended rounds): `rf-mitm` 3.1 µs at `ℓ = 8`, 57.6 µs at `ℓ = 12`, 482 µs at `ℓ = 14` against PS-0's `baseline` 9.3 ms, 3.07 s, 55.5 s (3,007×, 53,255×, 115,160×) and PS-2's `group-scalar` 276 µs, 48.1 ms, 757 ms (89×, 834×, 1,572×); `rf-pairs` within 0.91–1.18× of `group-scalar` on every rung (P2); prime, per full scan: `mitm3` 0.33 µs at `t = 4` → 306.5 µs at `t = 13`, `pairs3` 2.3 µs → 574 ms (1,873× apart at `t = 13`); A/A spread of the baseline 155–177 ms over 10 runs; 7 of 177 runs flagged contended and excluded, 8 of 185 commands refused by the settle check (fewer rounds on six `ℓ ≤ 12` rows, at least 2 each) |
| memory | the 2-sum table: `|F|²` keys at 12 bytes a key (8 B sorted key + 4 B fingerprint bucket): 198.6 MB at `ℓ = 13`, 804.9 MB at `ℓ = 14`, 150.2 MB at prime `t = 13` |
| not run this round | G2(c) msolve/Singular/PolyBoRi on exported tower systems; G4 (SAT); R8 trace replay; an ecbench end-to-end session (G6 is an extrapolation from the measured floors, §8) |

On the request's wording ("the fastest algorithm you can make"): within the
arms measured here the 2-sum table is the fastest on every cell, in counted
work and in wall time, by the factors in §4–§7; its exponent is the known
pair-table floor, and nothing measured here says anything about the
end-to-end sizes (§8).  Claims stop there.

## 1. What was built

`examples/pdp_rootfind` (Rust, in-tree; the PS harness's rung, corpus,
digest and kernel modules are shared by path).  All counters are
deterministic; every verdict carries its root set when witnesses are kept.

| arm | track | mechanism | per-target work (counted) | memory |
|:--|:--|:--|:--|:--|
| `baseline` | binary | PS-0's `semaev_decomp::decompose` (pairs-and-solve: quartic + `gcd(q, L_V mod q)`), unmodified | not instrumented (ns only) | `O(2^ℓ)` |
| `group-scalar` | binary | PS-2's scalar group-law row kernel, unmodified | not instrumented (ns only) | `O(2^ℓ)` |
| `rf-pairs` | binary | pairs-and-solve by the group law: for each row point `U = R ± P_i` (one batched inversion for all rows), the candidates `x(U ± P_j)` over `j ≥ i` from one batched inversion per row; a candidate decomposes iff it is in the abscissa set of `V` (the field-element test is a branch-free `< 2^ℓ`-span check); a zero denominator routes the row to an exact slow path | `≈ 10` mul, 4 sqr, 4 candidates a pair; `|F|(|F|+1)/2` pairs | none |
| `rf-mitm` | binary | 2-sum table: keys `x(P_i ± P_j)`, `i ≤ j` (one batched inversion a row), sorted; per target one batched inversion over `x_R + x_k`, `x(R ± P_k)` from one λ each, `2|F|` probes | `5|F| − 3` mul, `2|F|` sqr, 1 inv, `2|F|` lookups | `|F|²` keys |
| `subtract2` | prime | `m = 2`: `x(R ∓ P_k)` in the abscissa set of the base | `5|F| − 3` mul | none |
| `pairs3` | prime | `rf-pairs` with the tower's abscissa set (`KeySet`, open addressing) | as `rf-pairs` | none |
| `mitm3` | prime | `rf-mitm` on the Kummer base | as `rf-mitm` | `|F|²` keys |
| `mitm4` | prime | `m = 4`: `R ∓ (P_i ± P_j)` for every pair against the same 2-sum table | `≈ 9` mul a pair, 2 lookups a pair | `|F|²` keys |
| `brute_rational` | prime | general group law over every signed multiset of `m` base points (reference) | — | — |
| `algebraic_count` | prime | the PKM checker's resultant count of the zeros of the summation polynomial on `V^m` (`m ≤ 3`) | — | — |

**The lookup structure** (`fpset.rs`): the sorted key array keeps membership
exact; a probe first reads one 64-byte bucket of 32 sixteen-bit
fingerprints (`entries/16` buckets, overflowing buckets flagged and
answered "maybe"), then confirms a "maybe" by binary search.  False
positives about `16/2^16` a probe; no false negatives.  All `2|F|`
candidates of a target are computed first and probed in one pass with the
bucket of the probe 16 ahead prefetched, so the memory-level parallelism of
the core overlaps the misses.  It replaced a one-bit-per-key filter with a
binary-search fallback after the counts pass and before the timing stage
(§2, deviation 5); the counted units are unchanged and exactness was
re-verified on every corpus.

**Field arithmetic.**  Binary: one-word `GF(2^n)` with `pclmulqdq` folds
(the PS harness's `Gf2`), batched inversion (Montgomery's trick, 3 mul an
element).  Prime: Montgomery `Fp64` with `p₁ = 2013265921`, `FastCurve`
affine law with batched inversion.

## 2. Protocol, and the deviations from it

`PROTOCOL.md` fixed the arms, units, cells, exactness gates, six predictions
and the decision rule before the timed runs.  What differs:

1. **Gröbner reference cell (23, 8), not (24, 8).**  The repository's F4
   drivers refuse `n = 24` ("no K₁ at n = 24"); a (23, 8) corpus (80
   targets) was generated with the PS generator and every binary arm runs
   on it as well.
2. **The split-dreg F4 reference is censored.**  At (23, 8), unsplit `S₄`,
   `g = 0`: the first run was killed by the OS (exit 137), the retry under a
   12 GB address-space cap failed an allocation (exit 134); no solve row.
   The inherited-F4 decision oracle (`f6_ic_probe`, 8 targets) is the
   binary Gröbner reference that finished.
3. **Prime timing mode: full scan on 8 uniform targets at every `t`.**  At
   `t ≥ 12` every uniform target decomposes (no refuted target to time), so
   every rung's timed cost is the full scan — the cost of a refutation —
   on the first 8 uniform targets.  The counts of a refuted target in the
   first-hit mode equal the full-scan counts.
4. **Binary timing targets at `ℓ ≥ 12`: one pass; at `ℓ ≥ 13`: 4 refuted
   uniform targets** (the protocol said 16, as PS-0); the memoryless arms
   cost 2–30 s a target there.  The corpora at `ℓ = 13, 14` have 16 uniform +
   4 planted targets, as registered.
5. **Lookup structure changed before the timing stage** (§1).  The counted
   units (mul, sqr, inv, lookups, candidates) are unchanged; the bytes a key
   went from 9 to 12; every digest and brute-force check was repeated.
6. **Slope intervals are ± 1.96 × the least-squares standard error over
   rungs**, not a bootstrap over targets: the counted work of a table arm is
   the same for every refuted target (`2|F|` lookups), so a bootstrap over
   targets has zero width.
7. **Supersets**: the `m = 2` ladder runs to `t = 12` and the `m = 4` ladder
   to `t = 11` (registered: `t = 4…10`); the committed rows of the pilot
   (`p = 786433`) are read as well as rounds 2, 3 and 7.
8. **Not run** (listed in the protocol as references, not built this
   round): R8 trace replay, the SAT arms (G4), the external Gröbner tools on
   exported tower systems (G2c), the `n = 59, ℓ = 9` `mitm` cell (R2).

## 3. Exactness (G0)

**Binary** (`raw/bin_verify_fp.jsonl`; the pre-change rows in
`raw/bin_verify.jsonl`): for `rf-mitm` and `rf-pairs` on `ℓ = 8…12` (80
targets each: 64 uniform + 16 planted), `ℓ = 13` and `ℓ = 14` (20 each), the
verdict digest and the root-set digest equal the corpus's; 0 per-target
mismatches.  On (23, 8) the verdict digest is equal and 79 of 80 root sets
were compared equal; the 80th target is degenerate (`x_R ∈ V`), is routed to
the baseline decomposer by both arms, and its root set is not listed, so the
corpus's root digest differs on that cell by construction (recorded as
`root_sets_checked: 79, mismatches: 0`).

**Prime** (`raw/prime_rows.jsonl`, `raw/verify_prime_fp_m{2,3,4}.jsonl`,
`raw/prime_ladder_m{2,3,4}.jsonl`):

| check | result |
|:--|:--|
| committed rows read (pilot, rounds 2, 3, 7; Kummer, Dickson, 2-isogeny) | 568 |
| arm root set = `brute_rational` (where it runs) | 707 / 707 |
| finished F4 rows: `inconsistent` = (`algebraic_count` = 0) | 479 / 479; 18 rows unfinished (timed out, oversize or pairs above the degree bound) and excluded |
| algebraic solutions with no rational decomposition | 0 of 486 rows checked (P6's "twist solutions on some uniform targets" did not occur) |
| ladder: planted targets found | 8 / 8 at every rung and arm |
| ladder: brute-force agreement | 0 mismatches at `t ≤ 7` (`m = 2, 3`), `t ≤ 6` (`m = 4`); verdicts of all arms equal at every `t` |

## 4. Counted work: the one table

Medians over the counted targets (`raw/bin_count.jsonl`: 16 refuted uniform
targets a rung, 13 at `ℓ = 14`; `raw/prime_ladder_*.jsonl` and
`raw/prime_time.jsonl`).  `T_enum(m, t) = 2^{(m−1)t}/(m−1)!` candidate
sub-tuples; the pair floor is `2^{⌊m/2⌋ t}` lookups.  Full tables with
squarings, inversions, standard errors and every round: [`ANALYSIS.md`](ANALYSIS.md).

**Binary, `m = 3`, `n = 3ℓ`, curve `y² + xy = x³ + x² + 1`** (and the (23, 8) cell):

| arm | `n` | `ℓ` | `|F|` | mul / target | lookups | candidates | ÷ `T_enum` | ÷ pair floor | table MB |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| rf-mitm | 24 | 8 | 128 | 637 | 256 | 256 | `7.8·10⁻³` | 1.00 | 0.2 |
| rf-pairs | 24 | 8 | 128 | 83,069 | 0 | 33,024 | 1.01 | — | — |
| rf-mitm | 27 | 9 | 253 | 1,262 | 506 | 506 | `3.9·10⁻³` | 0.99 | 0.8 |
| rf-pairs | 27 | 9 | 253 | 322,319 | 0 | 128,524 | 0.98 | — | — |
| rf-mitm | 30 | 10 | 464 | 2,317 | 928 | 928 | `1.8·10⁻³` | 0.91 | 2.8 |
| rf-pairs | 30 | 10 | 464 | 1,080,653 | 0 | 431,520 | 0.82 | — | — |
| rf-mitm | 33 | 11 | 1,007 | 5,032 | 2,014 | 2,014 | `9.6·10⁻⁴` | 0.98 | 12.3 |
| rf-pairs | 33 | 11 | 1,007 | 5,079,305 | 0 | 2,030,112 | 0.97 | — | — |
| rf-mitm | 36 | 12 | 1,992 | 9,957 | 3,984 | 3,984 | `4.7·10⁻⁴` | 0.97 | 48.6 |
| rf-pairs | 36 | 12 | 1,992 | 19,858,245 | 0 | 7,940,112 | 0.95 | — | — |
| rf-mitm | 39 | 13 | 4,052 | 20,257 | 8,104 | 8,104 | `2.4·10⁻⁴` | 0.99 | 198.6 |
| rf-pairs | 39 | 13 | 4,052 | 82,129,985 | 0 | 32,845,512 | 0.98 | — | — |
| rf-mitm | 42 | 14 | 8,185 | 40,922 | 16,370 | 16,370 | `1.2·10⁻⁴` | 1.00 | 804.9 |
| rf-pairs | 42 | 14 | 8,185 | 335,044,787 | 0 | 134,004,820 | 1.00 | — | — |
| rf-mitm | 23 | 8 | 119 | 592 | 238 | 238 | `7.3·10⁻³` | 0.93 | 0.2 |
| rf-pairs | 23 | 8 | 119 | 71,873 | 0 | 28,560 | 0.87 | — | — |

`rf-pairs` tests 4 candidates a pair and `|F|(|F|+1)/2` pairs over the
`|F| ≈ 2^{ℓ−1}` base points (the abscissae of `V` that lift to the curve):
it *is* the enumeration floor (ratio ≈ 1).  PS-0's `baseline` counts its
pairs over the whole span `2^ℓ` of `V` (the quartic method needs no lift):
`2^ℓ(2^ℓ+1)/2` pairs, four times as many, at 337.7–452.3 ns a pair in PS-0's
own measurement; §7 reads it per target.  `rf-mitm` probes `2|F|` keys: it *is* the
pair-table floor (ratio ≈ 1 to the pair floor), `|F|/16`-fold below
enumeration in counted work, exactly P1's prediction (ratio 130× in
multiplications at `ℓ = 8`, 1,994× at `ℓ = 12`, 8,188× at `ℓ = 14`;
`|F|/16` is 8× and 128× — the multiplications a pair of `rf-pairs` are
about 10, not the 5–6 assumed for the base points, so the measured ratio is
`≈ 2 × |F|/16`).

**Prime, Kummer towers at `p₁ = 2013265921`, fresh ladder, uniform targets, full scan:**

| `m` | arm | `t` | `N` | `|F|` | mul / target | lookups | candidates | ÷ `T_enum` | ÷ pair floor | table MB |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 3 | mitm3 | 4 | 12 | 7 | 32 | 14 | 14 | 0.11 | 0.88 | — |
| 3 | pairs3 | 4 | 12 | 7 | 305 | 0 | 112 | 0.88 | — | — |
| 3 | mitm3 | 6 | 18 | 34 | 167 | 68 | 68 | `3.3·10⁻²` | 1.06 | — |
| 3 | pairs3 | 6 | 18 | 34 | 6,083 | 0 | 2,380 | 1.16 | — | — |
| 3 | mitm3 | 8 | 24 | 127 | 632 | 254 | 254 | `7.8·10⁻³` | 0.99 | 0.1 |
| 3 | pairs3 | 8 | 24 | 127 | 81,785 | 0 | 32,512 | 0.99 | — | — |
| 3 | mitm3 | 10 | 30 | 520 | 2,597 | 1,040 | 1,040 | `2.0·10⁻³` | 1.02 | 2.7 |
| 3 | pairs3 | 10 | 30 | 520 | 1,356,677 | 0 | 541,840 | 1.03 | — | — |
| 3 | mitm3 | 12 | 36 | 2,030 | 10,147 | 4,060 | 4,060 | `4.8·10⁻⁴` | 0.99 | 37.2 |
| 3 | pairs3 | 12 | 36 | 2,030 | 20,622,767 | 0 | 8,245,860 | 0.98 | — | — |
| 3 | mitm3 | 13 | 39 | 4,084 | 20,417 | 8,168 | 8,168 | `2.4·10⁻⁴` | 1.00 | 150.2 |
| 3 | pairs3 | 13 | 39 | 4,084 | 83,432,033 | 0 | 33,366,280 | 0.99 | — | — |
| 2 | subtract2 | 4 | 8 | 5 | 22 | 10 | 10 | 0.63 | 0.63 | — |
| 2 | subtract2 | 8 | 16 | 114 | 567 | 228 | 228 | 0.89 | 0.89 | — |
| 2 | subtract2 | 12 | 24 | 2,079 | 10,392 | 4,158 | 4,158 | 1.02 | 1.02 | — |
| 4 | mitm4 | 4 | 16 | 8 | 495 | 128 | 128 | 0.19 | 0.50 | — |
| 4 | mitm4 | 6 | 24 | 34 | 9,608 | 2,312 | 2,312 | `5.3·10⁻²` | 0.56 | — |
| 4 | mitm4 | 8 | 32 | 123 | 127,800 | 30,258 | 30,258 | `1.1·10⁻²` | 0.46 | 0.1 |
| 4 | mitm4 | 9 | 36 | 265 | 595,193 | 140,450 | 140,450 | `6.3·10⁻³` | 0.54 | 0.7 |

(`t = 5, 7, 9, 11` of `m = 3` and the other `m = 2` rungs are in
`ANALYSIS.md`; every rung of every ladder is in the raw records.)  At `m = 2`
the enumeration and pair floors coincide (`2^t`), and `subtract2` sits on
them.  `mitm4` probes about half the pair floor because `R ∓ S` is tested
for each unordered pair-sum `S` once and the table holds both signs.

**Amortised accounting (G5)**, from the recorded build counts
(`raw/bin_count.jsonl`, `prepare`; `raw/prime_ladder_*.jsonl`, `table`):
the 2-sum table costs `≈ 5` multiplications a key, half of what the
memoryless arm spends on one target, so the table pays for itself at `K = 1`.
Binary `ℓ = 13`: `T_prepare` = 41,032,579 mul; `T_target` = 20,257;
`T_prepare/K + T_target` = `4.1·10⁷` (`K = 1`), `2.58·10⁶` (`K = 16`),
`1.81·10⁵` (`K = 256`), against `8.21·10⁷` for `rf-pairs`.  `ℓ = 14`:
`1.67·10⁸` build against `3.35·10⁸` a memoryless target.  Prime `t = 13`:
`4.17·10⁷` build against `8.34·10⁷` for `pairs3`.

## 5. Against F4 (G2)

**Prime — the committed rows, same cells, same targets**
(`raw/prime_rows.jsonl`: each committed row's `p, a, b, g, ζ, t, x_R` rebuilt
into the cell, the arms run on that target; F4's own numbers copied from
the row; finished F4 runs only, medians over the rows of a cell).  Kummer
at `p₁`:

| `m` | `t` | `N` | F4 rows (finished / unfinished) | F4 multiply-adds | F4 ms | F4 degree | arm | mul / target | F4 multiply-adds ÷ mul |
|--:|--:|--:|:--|--:|--:|--:|:--|--:|--:|
| 2 | 6 | 12 | 2 / 0 | `6.50·10⁷` | 53 | 5 | subtract2 | 137 | `4.7·10⁵` |
| 2 | 8 | 16 | 2 / 0 | `6.18·10⁹` | 1,111 | 5 | subtract2 | 632 | `9.8·10⁶` |
| 2 | 10 | 20 | 2 / 2 | `4.84·10¹¹` | 56,800 | 6 | subtract2 | 2,572 | `1.9·10⁸` |
| 2 | 11 | 22 | 2 / 2 | `4.88·10¹²` | 596,937 | 6 | subtract2 | 4,862 | `1.0·10⁹` |
| 3 | 3 | 9 | 2 / 0 | `1.59·10⁸` | 90 | 6 | mitm3 / pairs3 | 17 / 113 | `9.4·10⁶` / `1.4·10⁶` |
| 3 | 4 | 12 | 2 / 0 | `1.06·10¹⁰` | 1,886 | 6 | mitm3 / pairs3 | 37 / 389 | `2.9·10⁸` / `2.7·10⁷` |
| 3 | 5 | 15 | 2 / 0 | `3.19·10¹²` | 370,705 | 7 | mitm3 / pairs3 | 77 / 1,421 | `4.1·10¹⁰` / `2.2·10⁹` |
| 3 | 6 | 18 | 2 / 5 | `8.99·10¹³` | 10,967,967 | 7 | mitm3 / pairs3 | 197 / 8,357 | `4.6·10¹¹` / `1.1·10¹⁰` |
| 4 | 2 | 8 | 2 / 0 | `4.85·10⁷` | 45 | 6 | mitm4 | 113 | `4.3·10⁵` |
| 4 | 3 | 12 | 3 / 0 | `8.55·10¹⁰` | 9,383 | 7 | mitm4 | 113 | `7.6·10⁸` |
| 4 | 4 | 16 | 2 / 3 | `2.71·10¹³` | 3,334,165 | 7 | mitm4 | 270 | `1.0·10¹¹` |

A multiply-add and a field multiplication are compared one to one (the
protocol's unit; it flatters F4).  G2(a) asks for `≤ 0.1×`: the ratios are
`10⁵`–`10¹²`, growing with the cell, on every `m`.  P4 predicted `< 10³`
multiplications for `mitm3` and `< 10⁴` for `pairs3` at `N = 18` (197 and
8,357 measured), `< 10⁴` for `mitm4` at `N = 16` (270), and a ratio above
`10⁹` on every `m ≥ 3` row with `N ≥ 15` (`2.2·10⁹` the smallest): all held.

**Prime — growth (G2b).**  Finished F4 rows on the `p₁` Kummer cells: 3.214
± 0.063 bits per unit of `t` at `m = 2` (6 cells, `t = 6…11`) and 6.556 ±
0.864 at `m = 3` (4 cells, `t = 3…6`); in bits per tower variable that is
1.61 and 2.19, the 1.6 and 2.1 the requirements document quoted.  The
root-finders: 1.052 ± 0.089 (`subtract2`), 1.019 ± 0.020 (`mitm3`), 1.990 ±
0.026 (`pairs3`), 2.028 ± 0.114 (`mitm4`, `m = 4`; no F4 fit at `m = 4`: three
finished cells).  The intervals are disjoint at `m = 2` and `m = 3` over
`≥ 4` sizes: G2(b) passes on the prime track.

**Binary — (23, 8), `m = 3`, 80-target corpus** (`raw/f6_ic_probe_n23_dim8.jsonl`,
`raw/f4_split_3-1-23-8.*`):

| arm | per target | notes |
|:--|:--|:--|
| inherited F4 (`f6_ic_probe`, the default decomposition oracle) | `4.10·10⁷` word XORs, 231.5 ms (medians of 8 targets, 7 refuted at `2.7·10⁷`–`4.7·10⁷` XORs / 225–330 ms, 1 hit at `2.7·10⁶` / 8.6 ms); shared host, not isolated | finishes |
| split-dreg F4 (unsplit `S₄`, `g = 0`) | killed (exit 137); retry under a 12 GB cap: allocation failure (exit 134) | censored: does not finish |
| rf-mitm | 592 mul, 238 sqr, 1 inv, 238 lookups; 3.9 µs (isolated, median of 5 rounds) | finishes, exact |
| rf-pairs | 71,873 mul, 28,798 sqr, 120 inv, 28,560 candidates; 202.2 µs | finishes, exact |
| baseline (PS-0) | 8,698.1 µs | finishes, exact |

Units differ (word XORs against field multiplications), so the binary
comparison is read on wall time: 231.5 ms against 3.9 µs, about `5.9·10⁴` (unlike measurements: a shared host against the reserved core; an order of magnitude, not a ratio to quote to two digits).  The binary F4 slope the
requirements document quoted (`2^{8.6}`–`2^{9.3}` per rung of `ℓ`, two
rungs of the split-dreg grid) was not re-fitted this round (G2b on the
binary track rests on that two-size ratio against the 1.000 ± 0.027 and
1.995 ± 0.055 measured here; it is not a four-size fit on the F4 side).
G2(c), the external tools on exported tower systems, was not run.

## 6. Growth and the generic line (G3, G5)

Least squares of `log₂(work)` against the rung, ± 1.96 s.e. (`ANALYSIS.md`,
"Growth"):

| track | arm | bits per rung | floor (bits per rung) | sizes |
|:--|:--|--:|--:|--:|
| binary `m = 3` (per rung of `ℓ`) | rf-mitm | 1.000 ± 0.027 | pair table 1 | 6 |
| binary `m = 3` | rf-pairs | 1.995 ± 0.055 | enumeration 2 | 6 |
| prime `m = 2` (per unit of `t`) | subtract2 | 1.052 ± 0.089 | both floors 1 | 9 |
| prime `m = 3` | mitm3 | 1.019 ± 0.020 | pair table 1 | 10 |
| prime `m = 3` | pairs3 | 1.990 ± 0.026 | enumeration 2 | 10 |
| prime `m = 4` | mitm4 | 2.028 ± 0.114 | pair table 2 | 6 |

Every root-finder's slope equals its floor within the interval (P3 and P5
exactly).  G3 asks for a slope below the pair table's at polynomial memory,
or a declared table with `log₂(S·T²) ≤ log₂(r/|Aut|)`: the table arm's
`S·T²` gap, with `S = |F|²` keys and `T = (mul + sqr)/6` group-operation
equivalents a target, `r = 2^{n−1}`, `|Aut| = 2`, is +6.9 bits at `ℓ = 8`
and +11.4 at `ℓ = 13`, growing (that line prices a whole DLP; a
decomposition oracle alone cannot pass it, and the row is kept as the G5
accounting).  **G3 fails; nothing here moves an exponent.**

## 7. Wall time (reserved core, five rounds)

Medians of the uncontended rounds, ns per refuted target converted to the unit shown; `raw/bin_time.jsonl`, `raw/isolated.jsonl`; the per-rung rows of every round are in `ANALYSIS.md` and in the run report's graph.  Five rounds were requested for every row; a settle-check refusal (the session's own agent process busy, 8 of 185 commands) leaves fewer rounds on six `ℓ ≤ 12` rows, and 7 runs the isolation tool flagged as contended are excluded (`ANALYSIS.md` lists both counts per row).

**Binary (`n = 3ℓ`; ms or µs per refuted target, median; rounds in parentheses):**

| `ℓ` | baseline (PS-0) | group-scalar (PS-2) | rf-pairs | rf-mitm | baseline ÷ rf-mitm | group-scalar ÷ rf-mitm | group-scalar ÷ rf-pairs | rf-mitm ns per probe |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 8 | 9.3 ms (4) | 276.5 µs (2) | 235.2 µs (4) | 3.1 µs (5) | 3,007 | 89 | 1.18 | 12.1 |
| 9 | 37.7 ms (4) | 809.7 µs (4) | 879.0 µs (5) | 6.4 µs (5) | 5,881 | 126 | 0.92 | 12.7 |
| 10 | 169.7 ms (15) | 2,693.2 µs (5) | 2,662.9 µs (5) | 10.7 µs (5) | 15,802 | 251 | 1.01 | 11.6 |
| 11 | 708.0 ms (5) | 11,290.4 µs (5) | 12,462.9 µs (5) | 23.9 µs (5) | 29,665 | 473 | 0.91 | 11.8 |
| 12 | 3,069.6 ms (5) | 48,094.4 µs (4) | 48,374.1 µs (5) | 57.6 µs (5) | 53,255 | 834 | 0.99 | 14.5 |
| 13 | 13,116.0 ms (5) | 191,078.5 µs (5) | 206,909.0 µs (5) | 156.0 µs (5) | 84,095 | 1,225 | 0.92 | 19.2 |
| 14 | 55,506.1 ms (5) | 757,478.1 µs (5) | 818,311.5 µs (5) | 482.0 µs (5) | 115,160 | 1,572 | 0.93 | 29.4 |
| (23, 8) | 8.7 ms (5) | 213.6 µs (5) | 202.2 µs (5) | 3.9 µs (5) | 2,230 | 55 | 1.06 | 16.4 |

The `ℓ = 10` baseline row pools the A/A pairs (15 runs).  `rf-mitm`'s table build is recorded apart (`setup_secs`: 0.9 ms at `ℓ = 8`, 0.31 s at `ℓ = 12`, 1.7 s at `ℓ = 13`, 8.2 s at `ℓ = 14`) and is not in the per-target figure (§4 charges it in counted operations).  The ns per probe column is the per-target wall divided by the `2|F|` probes, arithmetic included: it rises from 12 ns at `ℓ = 8` (0.2 MB table) to 29 ns at `ℓ = 14` (805 MB, past the 260 MiB L3), which is why `rf-mitm`'s wall grows by 1.196 ± 0.123 bits per rung where its counted work grows by 1.000.

**Prime (Kummer at `p₁`; µs per full scan of a uniform target, median of 5 rounds):**

| `t` | subtract2 (`m = 2`) | mitm3 (`m = 3`) | pairs3 (`m = 3`) | pairs3 ÷ mitm3 | mitm4 (`m = 4`) | mitm3 ns per probe |
|--:|--:|--:|--:|--:|--:|--:|
| 4 | 0.2 | 0.3 | 2.3 | 7 | 5.5 | 23.6 |
| 5 | 0.4 | 0.6 | 10.7 | 17 | 17.5 | 19.3 |
| 6 | 1.1 | 1.2 | 46.4 | 39 | 96.6 | 17.3 |
| 7 | 1.7 | 1.9 | 161.3 | 84 | 248.6 | 15.7 |
| 8 | 2.6 | 3.7 | 549.5 | 149 | 1,257.0 | 14.5 |
| 9 | 5.8 | 7.0 | 2,033.7 | 291 | 5,574.1 | 14.4 |
| 10 | 12.6 | 17.0 | 9,030.2 | 532 | 18,900.8 | 16.3 |
| 11 | 34.3 | 34.7 | 34,147.5 | 984 | 94,566.6 | 17.0 |
| 12 | 61.7 | 90.4 | 135,260.5 | 1,497 | — | 22.3 |
| 13 | — | 306.5 | 574,094.3 | 1,873 | — | 37.5 |

**A/A control** (PS-0's `baseline` against a copy of itself at `ℓ = 10`, one pass over 16 targets, opening every round; ms per target): 171.0 / 159.1, 156.7 / 155.1, 172.9 / 176.6, 170.1 / 161.1, 174.7 / 169.7 (rounds 1–5).  Within a pair the two runs differ by at most 7.5 %; across rounds the spread is 155–177 ms (14 % of the median), which bounds what a wall ratio below about 1.2× can mean on this host.

**Wall growth** (bits per rung, ± 1.96 s.e.): baseline 2.096 ± 0.011, group-scalar 1.935 ± 0.077, rf-pairs 1.973 ± 0.060, rf-mitm 1.196 ± 0.123 (binary, 7 rungs); subtract2 0.998 ± 0.070, mitm3 1.054 ± 0.089, pairs3 1.962 ± 0.031, mitm4 2.009 ± 0.071 (prime).  Predictions: P1's wall ratio to PS-0's baseline (`> 50×` at `ℓ = 8`, `> 1,000×` at `ℓ = 12`) held with 3,007× and 53,255×; P2 (`rf-pairs` within 0.7–1.3× of `group-scalar` in wall) held on every rung (0.91–1.18×).

**Against the Gröbner references in wall.**  Binary (23, 8): the inherited F4 decision oracle's 231.5 ms a target (median of 8, shared host, `raw/f6_ic_probe_n23_dim8.jsonl`) against `rf-mitm`'s 3.9 µs (reserved core): about `5.9·10⁴`, a ratio of unlike measurements read as an order of magnitude; the split-dreg F4 did not finish.  Prime: the committed `m = 3, N = 18` rows took 10,968 s of F4 on their own host against `mitm3`'s 1.2 µs full scan at `t = 6` here (`≈ 9·10⁹`, different hosts; the counted ratio of §5, `4.6·10¹¹`, is the comparable number).

## 8. End-to-end budgets (G6)

Extrapolation only: the measured slopes are the floors, so the per-target
cost at a deployed size is the floor formula with the measured constants.
Binary ECC2K-130 at `m = 4`, `ℓ = 33.9` (`|F| ≈ 2^{33.9}`): the 2-sum-table
oracle probes `≈ |F|²` keys a target, `2^{67.8}`, against the budget
`T_target ≲ 2^{35}` (§6.4 of the requirements), with a table of `2^{67.8}`
keys; at `m = 5`, `ℓ = 27.6`: `≈ |F|³/3` probes, `2^{81}`, against `2^{41}`.
Prime 256-bit at `m = 4`, `t = 64`: `2^{126}` probes against `2^{59}`.  G6 is
not met by any measured arm, as the floors say; an ecbench session was not
started (nothing to price that is not already the floor).

## 9. Gates read

| gate | result | evidence |
|:--|:--|:--|
| G0 exactness | **pass** on every arm and cell | §3 |
| G1 vs enumeration | **pass** for the table arms (ratio `7.8·10⁻³` → `1.2·10⁻⁴` binary, `0.11` → `2.4·10⁻⁴` prime `m = 3`, halving a rung over `≥ 6` rungs) — by memory, not algebra; **fail** (ratio ≈ 1) for `rf-pairs`, `pairs3`, `subtract2` | §4 |
| G2(a) vs F4 counted | **pass** prime: `10⁵`–`10¹²`; binary: the Gröbner arm that finishes is about `6·10⁴` slower in wall (§5, §7), the other does not finish | §5 |
| G2(b) vs F4 growth | **pass** prime (`m = 2, 3`, disjoint intervals over ≥ 4 sizes); binary: not a four-size fit on the F4 side | §5, §6 |
| G2(c) vs msolve | **not run** | — |
| G3 exponent | **fail**: every slope is its floor; `S·T²` gap positive | §6 |
| G4 vs SAT | **not run** | — |
| G5 accounting | reported (`S`, `T_prepare`, `T_target`, `K ∈ {1, 16, 256}`, `S·T²` gap) | §4, §6 |
| G6 end-to-end | **fail** (extrapolation at the floors) | §8 |
| G7 class | **engineering** on every row | — |

Decision rule (`PROTOCOL.md`): G0 held, G1/G2 read from the table, G3 not
passed → engineering result; the constant is recorded; no scoreboard row.

## 10. Figures

- [`figures/binary_cost.svg`](figures/binary_cost.svg) (PNG beside it):
  per-target wall time (isolated rounds, dots) and counted multiplications
  against `ℓ` for the four binary arms, with the fitted bits per rung.
- [`figures/prime_cost.svg`](figures/prime_cost.svg): counted multiplications
  and wall time against `t` for the prime arms, the committed F4 rows
  (multiply-adds) as context.
- [`figures/f4_vs_rootfind.svg`](figures/f4_vs_rootfind.svg): the committed
  F4 rows against the root-finders on the cells they measured.
- Source: `figures/plot.py` (reads `ANALYSIS.json`; lays out recorded
  values only).  Typst tables: `figures/tables.py` → `figures/table_*.typ`.

Canonical graphs checked and left unchanged: `docs/index-calculus-scoreboard.html`,
`docs/ic/progress-timeline.json`, the IC leaderboard, `docs/bounds/FRONTIER.md`
(no end-to-end `S`; an engineering result changes no scoreboard row).

## 11. Reproduction

```sh
cargo build --release --example pdp_rootfind
B=target/release/examples/pdp_rootfind; D=research/pdp_root_finder_20261012
C=research/pdp_oracle_speedup_20261011/baseline_v1/corpus
# exactness
$B bin-verify --l 12 --corpus $C/corpus_l12.txt --arm rf-mitm
$B prime-ladder --p 2013265921 --m 3 --t-min 4 --t-max 8 --seed 20261012 --uniform 16 --planted 8 \
   --arms pairs3,mitm3 --check-max-t 7 --alg-max-t 8 --targets 2 --passes 1 --min-secs 0.05 --label check --out /tmp/m3.jsonl
# counts, the committed rows, the timing stage, the analysis and the figures
$B bin-count --l 12 --corpus $C/corpus_l12.txt --arm rf-mitm --targets 16
$B prime-rows --rows "$(ls research/pkm_tower_{pilot_20260924,round2_20260925,round3_20260925,round7_20261006}/runs/*.jsonl | paste -sd,)" --reps 5 --out $D/raw/prime_rows.jsonl
$D/run_timing.sh 5                    # tools/isolated_bench.py, CPU 3 reserved, A/A then ABAB
$B analyze --dir $D                   # ANALYSIS.json, ANALYSIS.md
/tmp/report-venv/bin/python $D/figures/plot.py && /tmp/report-venv/bin/python $D/figures/tables.py
```

Corpora: PS-0's `corpus_l{8..12}.txt` (frozen, PR #1710); this study's
`corpus/corpus_l13.txt`, `corpus_l14.txt` (seed `20261011 + ℓ`, 16 uniform +
4 planted, exhaustive enumerator: 606 s and 1,731 s on 3 threads) and
`corpus_n23_l8.txt` (80 targets).

## 12. Environment and status

Claude Code cloud container: Intel Xeon @ 2.10 GHz, 4 vCPU (no SMT
siblings), L2 2 MiB a core, L3 260 MiB, 15 GiB RAM, kernel 6.18.44, rustc
1.97.0, base commit `cd0286878` (PR #1734).  Timing through
`tools/isolated_bench.py run --cpus 3` (CPU 3 reserved, other threads
evicted, PSI and other-process CPU recorded a run; runs flagged `contended`
are listed and not pooled), `RAYON_NUM_THREADS = 1`; counts in the same
binary (deterministic).  Conductor's control plane was unreachable in this
container (no Docker daemon or Postgres), so no task card is attached; the
hooks warned and let the edits through.  Python (`figures/*.py`) lays out
recorded values; every number in the tables comes from `ANALYSIS.json`,
written by the Rust `analyze` subcommand from the raw records.

| requested item | status |
|:--|:--|
| build the root-finder on both tracks | **done**: `rf-pairs`, `rf-mitm` (binary); `subtract2`, `pairs3`, `mitm3`, `mitm4` (prime); exact references |
| measure against the baselines (R1 PS-0, PS-2) and the Gröbner references (R3, R7) | **done** on the frozen corpora and the committed cells; R2 (`n = 59`), R4 (SAT), R5 (external tools), R8 (trace replay): not run |
| gates G0–G7 | **read**: G0 pass; G1 pass (table arms) / fail (memoryless); G2(a) pass, G2(b) pass prime, G2(c) not run; G3 fail; G4 not run; G5 reported; G6 fail; G7 engineering |
| "the fastest algorithm you can make" | the fastest arm measured here on every cell, at the pair-table floor; no claim beyond the measured scope |
| full report, graphs, PDF, benchmarks copy | this README, `REPORT.pdf`, `reports/run-20261012T2135Z/`, `benchmarks/20261012-pdp-root-finder/` |

## 13. Files

- `PROTOCOL.md` — the pre-registration.
- `run_timing.sh` — the isolated timing stage.
- `raw/` — `bin_verify.jsonl`, `bin_verify_fp.jsonl` (digests), `bin_count.jsonl`
  (counts), `bin_time.jsonl` (timed rows), `isolated.jsonl` (the isolation
  tool's record of every timed run), `prime_rows.jsonl` (+ `.stderr`; the
  committed rows re-read), `prime_ladder_m{2,3,4}.jsonl` (counts pass),
  `verify_prime_fp_m{2,3,4}.jsonl` (exactness after the lookup change),
  `prime_time.jsonl` (timed rows), `f6_ic_probe_n23_dim8.jsonl` (inherited
  F4), `f4_split_3-1-23-8.*` (the censored split-dreg runs),
  `timing_progress.txt`.
- `corpus/` — the `ℓ = 13`, `ℓ = 14` and (23, 8) corpora with their generator logs.
- `ANALYSIS.json`, `ANALYSIS.md` — written by `pdp_rootfind analyze`.
- `figures/` — `plot.py`, `tables.py`, the SVGs and PNGs, the Typst tables.
- `REPORT.typ`, `REPORT.pdf` — the digest.
- `reports/run-20261012T2135Z/` — the run report of the timing stage.
