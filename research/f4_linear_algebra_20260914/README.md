# Native F4 linear algebra: matched engineering experiment

This round changes only the row-reduction kernel used by the native Boolean F4
and Macaulay profiling paths. The final relation matrix is over the subgroup
modulus and uses a different solver. The frozen contract predates the kernel edit.

For a pivot in column `c`, earlier pivot eliminations and zero skipped columns
ensure that all columns before `c` of the pivot row are zero. The new loop XORs
only words `floor(c/64)..ceil(n_cols/64)`. Disjoint Rust borrows expose that
contiguous suffix to the optimizer. Pivot selection, elimination order, rank and
the full RREF stay identical. No additional allocation, dependency, cache,
configuration flag or CPU instruction-set requirement is introduced.

The operation counter continues to count actual 64-bit XORs performed. It now
excludes the eliminated zero-prefix work; it does not count scans, swaps,
construction, extraction or complete attack cost.

## Boundary and acceptance

With an elimination share `f` of total cost, even free elimination is bounded by
`speedup <= 1/(1-f)`. The generic relation-yield floor is unchanged. This is an
engineering change, with no claim against Pollard rho or the generic-group floor.
Complete calibrated operations, `S`, rho/floor ratios and exponents remain null.

The frozen [contract](contract.json) specifies the numerical kernel target and
correctness gates. The full external WDSat suite is inapplicable because that
executable does not call this native F4 kernel. Instead, the seven native S3
cases from the preceding polynomial reuse experiment are rerun on their original
targets, then on fresh polynomial systems. Reference and candidate must return
identical hashes of reduced rows and complete assignment sets, with exhaustive
evaluation as an independent oracle. Full DLP runs include the previous seeds,
new seeds and planted scalars, and the SAT path as a control.

## Reproduce

Use Rust 1.90.0, the saved `validation/dependencies.lock.txt` as `Cargo.lock` in both
checkouts, no `RUSTFLAGS`, and release builds with `redis-cache` compiled in but
both caches disabled. Check out the contract's base revision separately. Copy
`examples/f4_linear_algebra_bench.rs` into that checkout; leave its production
source unchanged. In each checkout:

```sh
cargo clean --release -p crypto
cargo build --release --locked --features redis-cache --bin ic --example f4_linear_algebra_bench
```

Prefer separate target directories. If sharing a dependency build cache, the
package clean is mandatory: Cargo can otherwise reuse a newer artifact from
the other worktree. Copy the reference executables before cleaning/building the
candidate. The comparator rejects identical reference/candidate binary hashes.

Copy each resulting `ic` and `examples/f4_linear_algebra_bench` executable into
separate reference/candidate directories. In the candidate checkout:

```sh
cargo test --release --locked --features redis-cache --lib koblitz_groebner:: -- --test-threads=1
cargo test --release --locked --features redis-cache --lib --no-run --message-format=json
```

The latter command identifies the executable for the library test harness.
Pass that executable as `--kernel`:

```sh
python3 research/f4_linear_algebra_20260914/run.py \
  --reference-dir /path/to/reference-binaries \
  --candidate-dir /path/to/candidate-binaries \
  --kernel /path/to/crypto_lib-test-executable \
  --output research/f4_linear_algebra_20260914/results/new-run
python3 research/f4_linear_algebra_20260914/compare.py \
  research/f4_linear_algebra_20260914/results/new-run
```

The Linux runner pins every child to one available CPU, clears inherited `IC_`
settings, sets one Rayon thread, alternates pair order, caps each process at
60 seconds and 2 GiB, and preserves stdout, stderr, status, source/binary hashes
and hardware details. Outputs are never overwritten. Kernel timings exclude
input construction and copying, both recorded separately; cold solver timings
include field setup, encoding and solving. Independent exhaustive validation and
an extra diagnostic root reduction are separate from the solver's own work.
Full DLP elapsed times include the entire attempt; reported phase timings are
subsets, not an exhaustive calibrated budget.

The paired percentile bootstrap resamples log candidate/reference ratios 4,000
times with seed 20260914. These intervals are shared-host diagnostics, and the
small repeated DLP corpus does not justify an end-to-end performance claim.

See [RESULTS.md](RESULTS.md) and the canonical
[scoreboard](../../docs/index-calculus-scoreboard.html#f4-linear-algebra-20260914).
