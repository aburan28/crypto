# n83 F6 union-domain width correction

The [registered protocol](PROTOCOL.md) identified a correctness barrier in
the union-of-Frobenius-charts SAT decomposition path. Its finite-coordinate
trie held each allowed abscissa in a `u64`, and the three-summand S4 path
extracted only the first 64-bit word. At field degree 83, coordinates
that differ above bit 63 therefore became indistinguishable to the domain
constraint. A solver's UNSAT result on that encoding could not prove
absence of a decomposition over the actual 83-bit factor base.

The implementation now uses `u128` codes through the trie, packs both
words of each field coordinate, and reports `unsupported` before encoding
fields wider than 128 bits. The existing exhaustive width-1-to-5 test
continues to pass. A new width-83 test accepts two coordinates with
identical low words and different high bits, rejects two neighboring
high-bit values, and checks extraction of field bits 64 and 82. This
fixes the finite-domain predicate; it does not implement a complete
n83 point-decomposition solver.

| Control | Result | Receipt |
| --- | --- | --- |
| Domain trie, widths 1–5 and 83, corrected final source | 2 passed, 0 failed | [domain_tests_clippy_fix.log.gz](domain_tests_clippy_fix.log.gz) |
| Domain trie, before test-only Clippy correction | 2 passed, 0 failed | [domain_tests.log.gz](domain_tests.log.gz) |
| n83 F6 pair geometry | 4 passed, 0 failed | [geometry_tests.log.gz](geometry_tests.log.gz) |
| Small-curve union SAT versus exact enumeration | 1 passed, 0 failed | [union_sat_control.log.gz](union_sat_control.log.gz) |
| rustfmt and `git diff --check` | pass | local command output |

The three original Cargo commands used `--offline --release --lib`,
`RAYON_NUM_THREADS=1`, and the shared release target directory shown
below. Their precise test filters were
`domain_trie_`, `f6_wide_geometry`, and
`frobenius_union_domain_is_invariant_and_sat_matches_search`:

```sh
CARGO_TARGET_DIR=/Volumes/SSD990/crypto/worktrees/f6-n83-geometric-20261004/target RAYON_NUM_THREADS=1 cargo test --offline --release --lib domain_trie_ -- --nocapture
CARGO_TARGET_DIR=/Volumes/SSD990/crypto/worktrees/f6-n83-geometric-20261004/target RAYON_NUM_THREADS=1 cargo test --offline --release --lib f6_wide_geometry -- --nocapture
CARGO_TARGET_DIR=/Volumes/SSD990/crypto/worktrees/f6-n83-geometric-20261004/target RAYON_NUM_THREADS=1 cargo test --offline --release --lib frobenius_union_domain_is_invariant_and_sat_matches_search -- --nocapture
```

GitHub's Rust 1.98 Clippy check then found a redundant `as i32` cast in
the new width-83 test. Removing it changed no production logic. The
corrected final source reran `domain_trie_` with the same command and
passed both tests; its new receipt is in the first table row. The other
two controls were run on the preceding source, before this test-only
correction. The logs are compressed losslessly with `gzip -n`; use
`gunzip -c` to read them. The earlier [initial domain](domain_tests_initial.log.gz) and
[initial geometry](geometry_tests_initial.log.gz) logs are retained from
the first committed version before the direct two-word packing assertion
was added.
The final `koblitz_index_calculus.rs` SHA-256 is
`2caddb1e9cd20df9afa530aa7620a3ba20f7cc2279058a3e441895555a9fd570`;
the preceding test-source SHA-256 for the earlier two controls was
`0f4a2d3e1a7e5cf80524d2bd24aaea6bc6ffa5ecba99961cc1b357b55df879bd`.
The protocol SHA-256 is
`f5a52360b9e72be12d48c535350fa81a4235e7c9294bfee206bbe9b01d69df97`.

This is a correctness prerequisite for n83 F6, not a speed result. The
current S4 SAT encoder still expands a full 83-bit ambient representation
of every summand, and no ordinary n83 relation, complete one-target IC
solve, F4/F5 comparison, or paired rho measurement was performed here.
The next solver gate remains a bounded exact join between small algebraic
blocks and the terminal pair index, with explicit `found`, `refuted`,
`exhausted`, and `unsupported` outcomes.
