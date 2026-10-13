# Root-finding point decomposition: source assessment and ecbench measurement

Status: local measurement complete and replayed; independent-host validation remains open under the [complete-cost calibration problem](../../problems/rootfinding-complete-cost-calibration.md). Source assessment is bounded to the supplied 13-page PDF and the matching paper branch; this report is not a research-state decision.

Follow-on: the [2026-10-13 prime pair-table size trial](../ecbench_rootfinding_tuning_20261013/REPORT.md) compares 32 versus 64 abscissae on eight fresh public targets. It is a separate counted-cost configuration comparison; the frozen measurements below remain unchanged.

## Source binding

The supplied [`rootfinding_decomposition.pdf`](source-paper.pdf) is frozen here as `source-paper.pdf` with SHA-256 `743d3d3e7e3aab6f980813c0944c1d2cbacc906e075e7283de0abc7b923b7c01`, identical to `papers/rootfinding-decomposition/rootfinding_decomposition.pdf` in `crypto` commit `5bb2f20edb` on `origin/feature/semaev-rootfind-nist-paper`. The paper source and end-to-end measurement program and results are present in that commit. Its `README.md` labels the manuscript a research draft and says the novelty audit and a stronger F4 cross-check remain due.

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

The committed native [PDP root-finder study](../pdp_root_finder_20261012/README.md) at `3745d06911c81227dac551e1dbd4a59f639f82ab` tested the same group-law and two-sum-table mechanisms on separate binary and prime polynomial-grid workloads. It reports 707/707 brute-force agreements on the checked prime rows, binary corpus digest agreements, and slopes at the generic enumeration or pair-table floors. Its large ratios are per-oracle comparisons with different field-operation and wall units; its README explicitly records that an ecbench complete solve was not run in that study. The present sessions supply that different pipeline view on two small public-target curves, so its oracle ratios are not combined with the ecbench ratios.

## Fresh ecbench protocol

See [PROTOCOL.md](PROTOCOL.md) and the exact two specs under `specs/`. The `ic.pipeline` subtraction oracle implements the same search strategy but uses ecbench's generic group arithmetic, not the paper's batched-inversion kernel. The Koblitz `mitm-frobenius` oracle likewise measures the pair-table strategy in ecbench's pipeline. The primary metric is ecbench's charged operation count and its complete phase accounting. MacOS wall measurements are L0 and descriptive.

## Fresh ecbench results

The frozen inputs were committed as `4ffb2c91600a8e82dc7e8790a253e1f0e5b2d523` before execution. The exact release binary has SHA-256 `be956553a5a2cd1227e8f1add8dd24e11dffba5068c799234b55e44ffc3a1b61`. The host capsule reports macOS arm64, 14 logical CPUs, environment class `ECBENV2hd82681268e96`; all measured records earned L0. Both sessions used eight distinct public targets, three measured rounds and one warmup, with matched algorithm seeds and full-rank checked IC linear algebra. Every execution verified the recovered scalar. See `sessions/prime20/` (`ECBS1h40838be1a86e`) and `sessions/koblitz17/` (`ECBS1ha82cdeaf82db`).

| Curve (prime-order subgroup) | Arm | Base points / columns | Verified measured | Mean S, counted GAE / sqrt(r) | Counted S vs matched rho | Bound status |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `icv1-fp20-t727-cd198a38` (r=562427) | rho.negation | - | 24/24 | 1.464 | 1.000 | lower bound |
| same | IC direct subtraction | 64 / 32 | 24/24 | 1524.511 | 1041.344 | lower bound |
| same | IC pair table | 64 / 32 | 24/24 | 28.753 | 19.640 | lower bound |
| `icv1-f2m17-tm101-00378d4e` (r=65587) | rho.signed_frobenius_strong | - | 24/24 | 4.046 | 1.000 | lower bound |
| same | IC direct subtraction | 443 / 14 | 24/24 | 356.894 | 88.219 | fully priced in ecbench.gae |
| same | IC Frobenius pair table | 443 / 14 | 24/24 | 20.360 | 5.033 | fully priced in ecbench.gae |

The paired pair-table/subtraction counted-S quotient is 0.018861 (95% cluster interval 0.018547-0.019242; 24 pairs) on the prime curve and 0.057047 (0.041159-0.092069; 24 pairs) on the Koblitz curve. Both IC arms used the same factor-base digest within each curve and the same public targets. The Koblitz quotient is an ecbench.gae comparison without an unpriced-work flag. The prime quotient is bounded: square roots, Legendre tests, hash probes and row operations remain unpriced, so the counted quotient does not bound a complete-cost ratio. Every comparison against rho is likewise bounded by unpriced native work. These are measured count ratios at the stated curve sizes, not a scaling law or an admitted wall-time gain.

| Curve and IC oracle | Mean factor-base GAE | Mean table-setup GAE | Mean relation GAE | Mean linear-algebra GAE | Mean verification GAE | Mean relation trials |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Prime subtraction | 0 (native work unpriced) | 0 | 1,142,460 | 0 (native work unpriced) | 849 | 17,869 |
| Prime pair table | 0 (native work unpriced) | 1,056 | 19,659 | 0 (native work unpriced) | 849 | 17,869 |
| Koblitz subtraction | 54 | 0 | 91,030 | 14 | 302 | 321 |
| Koblitz pair table | 54 | 3,537 | 1,307 | 14 | 302 | 322 |

The Koblitz table setup is about 68% of that arm's cold counted cost. This session therefore directly charges the reusable table rather than comparing only queries. The small difference in Koblitz mean trial counts comes from when the two runs reach full rank; each independently checks its recovered scalar and base logs. The prime arms attempt the same mean 17,869 relation queries. The count reduction comes from the oracle stage, while verification and matrix work remain.

Peak process RSS from the measured child records was 30,992 KiB for both prime IC arms, 31,328 KiB for Koblitz subtraction, and 31,360 KiB for the Koblitz pair table. The maxima for the matched rho arms were 27,632 KiB prime and 27,968 KiB Koblitz. These are whole-process high-water marks, not isolated pair-table allocations; at these small bases, they do not support a stable table-memory slope. The phase table above charges setup operations separately.

The frozen protocol called for a 12-record `ecbench verify` replay. The release CLI exposes this operation as `ecbench audit`, so that command was used. Its initial 12-replay sample hit only the repeated rho control arm because it selected every eighth record. The subsequent `--replay-all` receipts are the operative local audits: `audit-prime20-all.json` (SHA-256 `d57819e8393132080f01eb653f44606e70f578cf32cbab5876f0df079779b257`) and `audit-koblitz17-all.json` (SHA-256 `bbefb52875cab379ee69d221aab9a68d0b8b26c58e73915311688a3f8351ac07`). Each reports 128 sealed records, 128 verified outcomes, 96 identical measured-run replays across all arms, and zero problems. An independent host class has not replayed them yet; local replay alone does not satisfy an official cross-host claim gate.

The source for the phase summaries is `analysis.jq`, applied to each session's sealed `records.jsonl`; its outputs are `summary-prime20.json` and `summary-koblitz17.json`. The six saved paired comparisons under each session's `comparisons/` directory contain the exact workload-level rows and cluster intervals. The A/A controls were exact in counted operations, but all wall comparisons were descriptive because none reached L2. The observed host load was roughly 80-107 during preparation, and the release build took 24 minutes; neither fact changes the deterministic counts.

## Assessment

The group-law reduction for two summands is sound and useful when a polynomial solver was doing that membership test. The paper's stage measurements give strong evidence that its optimized implementation removes substantial overhead in this narrow oracle. Its complete-solve evidence has two distinct tiers: the 16-bit prime and binary rows directly time both arms, while the larger prime F4 totals are composed estimates. The several-thousand-fold NIST oracle ratios are against Buchberger with a six-point base; the benchmark source's printed `F4` label should not be used to describe that solver. The paper's own draft README leaves a novelty audit and stronger F4 cross-check open.

The fresh ecbench sessions add a separate accounting view on small public-target curves. They show that table preprocessing can reduce counted IC work relative to direct subtraction, yet the matched rho arms still have lower counted S at both sizes. The Koblitz table's setup is material in a cold one-target solve. The prime and rho rows omit native work, the local wall level is L0, and the ecbench subtraction oracle uses generic additions rather than the paper's batched-inversion kernel. These limits prevent transferring the new count ratios to the paper's NIST timings or to larger fields.

## Reproduction and visual record

From the pinned `4ffb2c91600a8e82dc7e8790a253e1f0e5b2d523` source, the measured commands were:

```sh
target/release/ecbench run --wait --allow-busy --cpus none --spec research/ecbench_rootfinding_assessment_20261012/specs/prime20.json --out research/ecbench_rootfinding_assessment_20261012/sessions/prime20
target/release/ecbench run --wait --allow-busy --cpus none --spec research/ecbench_rootfinding_assessment_20261012/specs/koblitz17.json --out research/ecbench_rootfinding_assessment_20261012/sessions/koblitz17
target/release/ecbench audit --session research/ecbench_rootfinding_assessment_20261012/sessions/prime20 --replay-all --json-out research/ecbench_rootfinding_assessment_20261012/audit-prime20-all.json
target/release/ecbench audit --session research/ecbench_rootfinding_assessment_20261012/sessions/koblitz17 --replay-all --json-out research/ecbench_rootfinding_assessment_20261012/audit-koblitz17-all.json
```

The source-derived visual record is [`mean_s.svg`](mean_s.svg), [`pipeline.svg`](pipeline.svg), [`run_report.pdf`](run_report.pdf), and [`report.pdf`](report.pdf), generated by `render_visuals.mjs`, `run_report.typ`, and `report.typ`. The [validation log](VALIDATION.md) records the repository test outcomes and exact reruns. The first chart marks every counted lower bound. The diagram locates the algorithm's setup, oracle, matrix, and verification phases. No canonical frontier graph changes: these local L0 measurements have unpriced work in their rho baselines and do not pass an independent-host claim gate.
