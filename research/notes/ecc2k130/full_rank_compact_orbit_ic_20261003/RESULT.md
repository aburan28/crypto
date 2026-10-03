# Full-rank compact-orbit continuation: n37 gate passes, cold cost remains above rho

**Decision.** The preregistered n37 full-rank gate passed on all 16 measured
fresh public-Q/round pairs. The new mode continued the same relation stream
from target pinning at rank 24–39 to **rank 43/43**, recovered the same target
scalar as the early-pin control and independently passed `[d]G=Q` in every
run. All 48 measured records, including strong rho, reproduced identically
under `ecbench verify --replay-all`; the 24 warmups also verified. The frozen
charged-cost rejection passed: the smallest full-rank IC lower bound was
84,323.502481 GAE, versus the largest strong-rho recorded charge of
5,363.493947 GAE on the same Q/round panel, a **15.722×** separation in
that accounting. The paired aggregate full-rank IC/rho ratio was **22.600×**
[19.505, 26.620]. This is a scoped n37 cold single-target result for this
42-orbit eager folded table. Both arms have unpriced native work, so it is
not a fully priced attack-speed or ECC2K-130 transfer result.

The [protocol](PROTOCOL.md), [spec](SPEC.json) and exact
[`ecbench plan` output](PLAN.txt) were pushed to
[PR #1293](https://github.com/aburan28/crypto/pull/1293) before this
session. Its target seed `202610030741` differs from the prior
`202610030641` cohort; a comparison of the two saved target-coordinate sets
found zero repeated points. The plan fixed eight registered ICV1 n37 public Q,
three arms, one warmup and two measured rounds, 1,000,000 relation trials
per IC attempt and a 600-second per-run timeout. None of these changed
after target outcomes were seen.

## Reproducibility and gates

| Item | Frozen or observed record |
|---|---|
| Curve and boundary | `icv1-f2m37-tm534059-32aad96b`, `r=230603167`, signed-Frobenius order 74, generic floor `S=0.1456948091` |
| Source support | 42 columns, 3,108 signed points, archived `SOURCE42.jsonl` SHA-256 `0a32de24a5680ff46baf9543e8bf8e32447323491e2a2adedce2548c14a25f75`; release source-parity test passed |
| Factor base | `FB1h8dbe08ce02c8` on every measured IC row; source algorithm and parameters identical between the arms |
| `Cargo.lock` | SHA-256 `b28d3c2d81146a40d00df85f45c6460d9e2bc5c26875307a06cfd15efff99365` |
| Binary | SHA-256 `194202898bd7c2705b340c061e0e57bb464e339eae5bbea7f232b3ff2808bf38`, clean source commit `c09894cadae912b5f4131e95942cc74effda6d21` |
| Host | macOS aarch64, 14 logical CPUs, rustc 1.93.1, environment class `ECBENV2hd82681268e96`; all rows L0 |
| Session | [`ECBS1hfaeefe6e296a`](sessions/n37_full_rank_v1/session.json), 72/72 verified including warmups; 16/16 measured per arm |
| Audit | [receipt](AUDIT.json), SHA-256 `95e9eba91b95b3799011fe4d0e81790a03b9ce0be729c7afc5b1a0b7a799f166`; `ok=true`, 48/48 measured replays identical, zero problems |
| Historical compatibility | New binary also replayed all 48 measured records of the preceding early-pin session identically; [receipt](LEGACY_AUDIT.json) SHA-256 `808d93574f07598f93aa9581f9fc9e19b662504cf14ec64a91903d001d77f402` |

The full-rank arm records the first target-pin rank, trial and relation.
All 16 checkpoints exactly match the paired early-pin arm's final rank,
trial and relation count. The algorithm seed, recovered scalar, factor-base
identity, base charge and oracle-setup charge also match in each pair.
This is a strong prefix checkpoint, not a saved row-by-row relation
transcript. The full-rank arm finished after 66–117 trials, with 23–74
dependent rows; every trial returned a valid decomposition and every final
matrix had rank 43. The unit test additionally stops the full-rank arm at
the early-pin trial and confirms that it reports exhaustion rather than a
false full-rank success. Rank 43 makes all 42 base unknowns and the target
unknown mathematically unique in the matrix. This run explicitly extracts
and independently verifies the target scalar; it does not separately emit
and independently point-check each base logarithm.

## One table, one charged unit

`S = charged GAE / sqrt(r)`. The [native ecbench table](TABLE.txt),
[per-pair data](PAIRED.json), [extraction filter](PAIRED.jq), and saved
[comparisons](sessions/n37_full_rank_v1/comparisons/) retain exact costs.
Intervals resample the eight independent workloads by cluster; both IC
and rho rows are lower bounds because native field, hash or table work
remains unpriced. L0 wall times are descriptive only.

| Arm | Measured answers | Mean S [95% interval] | S / generic floor | S / same-Q rho | Class |
|---|---:|---:|---:|---:|---|
| Strong signed-Frobenius rho | 16/16 | 0.250 [0.212, 0.290] | 1.714 | 1.000 | reference, bounded |
| IC, early target pin | 16/16 | 5.433 [5.417, 5.448] | 37.291 | 21.757 [18.738, 25.620] | control, bounded |
| IC, full rank 43/43 | 16/16 | 5.644 [5.612, 5.676] | 38.736 | **22.600 [19.505, 26.620]** | full-rank gate pass; bounded |

The full-rank/early-pin paired charged-cost ratio is **1.0387**
[1.0329, 1.0439]. The added work is online: the common counted factor-base
and oracle setup cost 4,302.554148 and 75,032.630988 GAE respectively,
while mean relation-plus-linear-algebra-plus-verification cost rose from
3,170.114435 to 6,367.044338 GAE. Mean trials rose from 32.3125 to
89.5. Linear algebra itself averaged only 8.222518 GAE in full-rank runs;
the extra relation search, not elimination, drove this increment. The
66,822 eager-table setup additions alone exceed the maximum sampled rho
charge by 12.459×. The minimum full-rank IC charge exceeds that rho
maximum by 15.722×; per-Q/round IC/rho ratios range 15.852–34.242.

The prior [cold-source result](../cold_compact_orbit_ic_20261003/RESULT.md)
listed the same IC unpriced counters: inversions in the base constructor,
raw abscissa scan and lifts, curve-equation checks, several native field
operations, hash probes, point insertion/negation and a storage-byte floor.
The folded table's normal-basis setup and allocation also lack dedicated
prices. Strong rho records unpriced canonicalisations, partition hashes
and distinguished-point table operations. `ecbench` reconstructs the base
to report metadata after the measured pipeline; that metadata work is
outside algorithmic GAE but inside process wall time. No full-speed
ordering follows from comparing two lower bounds, and no wall-time claim
is made from this L0 session.

## Decision for the next experiment

H1 passes: the earlier rank failure was a completion-rule issue for this
source/oracle, not an observed rank ceiling at the frozen cap. H2 passes in
the preregistered charged accounting: carrying the matrix to full rank
raises cold cost by only 3.87%, while the eager table setup is still much
larger than the sampled one-target rho charge. This narrows the immediate
engineering question. Further single-target tuning of elimination cannot
remove the 66,822-addition setup floor; the next informative gate is a
**shared cold table and full-rank batch recovery** against actual
automorphism-aware batched rho on new disjoint public Q. It must freeze
the batch size, charge base/table construction once, separately account
for each online descent, peak memory and unpriced native work, and verify
every scalar. A ratio against independently repeated one-target rho would
answer a different question. n41/n53 capacity, descendant-native/pullback
transport and the n83 confidence gate remain separate before any n131
extrapolation.

Reproduce the derived [paired file](PAIRED.json) with:

```sh
jq -s -f research/notes/ecc2k130/full_rank_compact_orbit_ic_20261003/PAIRED.jq \
  research/notes/ecc2k130/full_rank_compact_orbit_ic_20261003/sessions/n37_full_rank_v1/records.jsonl
```
