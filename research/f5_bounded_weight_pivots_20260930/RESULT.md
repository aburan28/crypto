# Bounded current-weight F5 pivots: structural gate rejected

The [protocol](PROTOCOL.md) was committed as `26534a6a0` before candidate code or measurement, and draft PR #1072 was opened before the candidate was written. The opt-in implementation at `b6c8eae748dc5de0611e6ffbc9984800d7baf97e` selected the lightest current row among only the first 8 or 16 eligible rows for each pivot. It avoided the full remaining-row weight maintenance used by the earlier [unbounded experiment](../gf2_minweight_pivots_20260929/RESULT.md).

Release GF(2) tests passed (7/7, including the new sparse, dense and partial-word row-space test); release F5 tests passed (10/10); and the release benchmark built. The [raw structural receipt](structural_b6c8eae74.json), SHA-256 `8e61d331c45a2926830c50fc8eb573d7875e6ca091aaee80f40dcd982a6f9bfa`, preserves all eight processes: one reference and one candidate for each K and each seed, all seven cases, stdout/stderr, exit status, host, source and binary hashes, exact invariants, active route and observed phase times. Every process succeeded. The release binary had SHA-256 `fba110474d29659c58dd80af06c33ac1c82d165b209b18752042f5abf45c9a64`. The host was Apple ARM64, macOS 26.6, Rust 1.93.1; its load average was about 9–12 during the screen. Timing was explicitly **not** used for this gate and no wall-time speedup is claimed.

| Seed | Reference terms | K=8 terms (candidate/reference) | K=16 terms (candidate/reference) | Canonical row space and rank |
| --- | ---: | ---: | ---: | --- |
| `0` | 13,734,979 | 13,519,606 (0.9843) | 13,342,888 (0.9715) | identical; rank 6,924 |
| `badc0de1` | 13,722,549 | 13,437,669 (0.9792) | 13,229,511 (0.9641) | identical; rank 6,924 |

For all seven cases on both seeds, rank, canonical row-space fingerprint, F5 criterion work, built/pruned row counts and column counts agreed with the reference. The bounded route was active only on the frozen n24 degree-4 case; raw rows, output terms and counted reduction work on inactive cases matched exactly. The baseline process was repeated for each K, and its deterministic outputs agreed across repeats. No failures, timeouts or OOMs occurred.

The frozen structural gate required at most **75%** of reference output terms on **both** seeds. K=8 and K=16 both miss it by a large margin; the best was 96.4% on one holdout. Per the protocol, no local paired timing or x86 speed test follows. This result says that restricting the current-weight choice to the first 16 eligible rows did not capture the earlier unbounded method's roughly 50% term reduction. It does not establish why; useful light rows may occur farther down the strip, or the unbounded method's maintained weights may change choices differently.

The exact measured source change is preserved as [measured_candidate.patch.gz](measured_candidate.patch.gz), gzip SHA-256 `51edf12749506cf034ef2a73fc7db7c63b49fcf5c0a329f81f65405c24afc151`, uncompressed patch SHA-256 `412125a9f38e16d51196a42b90c812bde509243ea2fc699dc5082283327cc108`. On protocol commit `26534a6a0`, replay with `gzip -dc measured_candidate.patch.gz | git apply`. The experimental code was removed after the gate; shared GF(2), F5 and benchmark source are byte-identical to reference commit `4877c7577d7a164e7de0dededa3c1b7486f4ff35`.

The requested further 2× complete-call gain remains unproved. This is a Boolean matrix-F5 solver-stage structural diagnostic, not an IC one-target online DLP or rho speedup.
