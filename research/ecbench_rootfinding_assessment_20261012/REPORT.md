# Root-finding point decomposition: source assessment and ecbench measurement

Status: measurement in progress. Source assessment is bounded to the supplied 13-page PDF and the matching paper branch; this report is not a research-state decision.

## Source binding

The supplied `rootfinding_decomposition.pdf` has SHA-256 `743d3d3e7e3aab6f980813c0944c1d2cbacc906e075e7283de0abc7b923b7c01`, identical to `papers/rootfinding-decomposition/rootfinding_decomposition.pdf` in `crypto` commit `5bb2f20edb` on `feature/semaev-rootfind-nist-paper`. The paper source and end-to-end measurement program and results are present in that commit. Its `README.md` labels the manuscript a research draft and says the novelty audit and a stronger F4 cross-check remain due.

## Algorithm assessed

For two summands, fixing one factor-base point turns the third Semaev polynomial into a quadratic whose roots are the abscissae of `R+F_i` and `R-F_i`. Enumerating the base, subtracting a point and checking a point hash table therefore decides the two-summand membership problem. The paper's prime implementation batches the field inversions across the base. The binary implementation precomputes pair sums and probes the table for two-, three- and four-summand decompositions, charging quadratic setup and memory to the base.

The mathematical identity is exact under the paper's stated generic affine condition. A returned witness can be re-added for soundness. The exhaustive-enumeration completeness argument holds only with the declared sign and factor-base convention; the implementation's binary table stores full point keys and a base containing both signs. The paper's schematic x-only table needs that convention made explicit. The larger-summand enumeration has a growing exponent in factor-base size and is not benchmarked end to end in the paper.

## Paper evidence and its scope

| Scope | Reported observation | Boundary |
| --- | --- | --- |
| Binary per-decomposition, `m=2` | 383× at `n=9`, 1535× at `n=13` versus its matrix-F4 | Pair-table build excluded from the per-query ratio. |
| Five standard NIST P-curves, `m=2`, six-point base | 3021–5802× per query versus its Buchberger implementation | Oracle only; no complete solve or matched rho. |
| Binary complete solve | F4 faster at `n=5,9,13`; pair-table strategy 122× faster at `n=23` | One cold solve includes table construction. |
| Prime complete solve versus `msolve` | 100× at 16 bits, 105× at 20 bits, 50× at 24 bits | 16-bit F4 total directly measured; 20- and 24-bit totals composed from sampled per-query CPU cost and the root-finder's measured query count. |

The NIST benchmark labels its printed baseline "F4", but `examples/semaev_rootfind_nist_bench.rs` calls `groebner_f4::buchberger` at the timed line; those 3021–5802× figures are therefore comparisons to Buchberger. The paper's own classification is engineering. The two-summand oracle has a lower cost than solving the polynomial system for the measured workload, but the relation yield, factor-base size, final matrix, target descent and matched rho determine the cost of an ECDLP method. Its NIST measurements do not price those stages. A direct precomputed pair-sum table also changes the amortization point: one table costs quadratic space and construction, while a probe can be constant time. Consequently the paper's general statement that no implementation can take fewer than linear field operations per ideal is too broad for repeated queries with reusable preprocessing.

## Fresh ecbench protocol

See [PROTOCOL.md](PROTOCOL.md) and the exact two specs under `specs/`. The `ic.pipeline` subtraction oracle implements the same search strategy but uses ecbench's generic group arithmetic, not the paper's batched-inversion kernel. The Koblitz `mitm-frobenius` oracle likewise measures the pair-table strategy in ecbench's pipeline. The primary metric is ecbench's charged operation count and its complete phase accounting. MacOS wall measurements are L0 and descriptive.

## Fresh ecbench results

Pending execution and audit. All statuses will be retained, including errors and budget exits.
