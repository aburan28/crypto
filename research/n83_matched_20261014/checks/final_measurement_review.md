# Final measurement and source review

Reviewed 2026-10-15. Outcome: the complete local matched experiment passes its
declared correctness and accounting gates. Rho is cheaper at both measured
complete-process boundaries. The secondary online panel does not establish a
mean-cost winner, and the modeled large-batch crossing remains unmeasured.
This review makes no ledger admission or global-optimality claim.

## Native results and completeness

The authority is `records/postprocess-20261015T0018Z/statistics.json`, generated
by the unchanged frozen Rust analysis. All ratios below are IC/rho, with a
value below 1 favoring IC. Intervals are the source's prospective 20,000-draw
paired percentile bootstrap, not newly fitted intervals.

| Declared boundary | Mean IC seconds | Mean rho seconds | Native IC/rho [95%] | Resampling unit |
|:--|--:|--:|:--|:--|
| Actual cold singleton | 1159.34988210875 | 3.0716487838750006 | 377.43569127919193 [311.8128466078986, 499.01021039372205] | Eight paired fresh-process targets |
| Actual shared-setup 32-target process | 1232.9752171726666 | 78.95740027966667 | 15.615701793694768 [14.743542634654554, 17.614847798270176] | Three paired preparation/process cohorts |
| Secondary online, three-repeat target means | 1.7603639429062503 | 2.4572885118541667 | 0.7163847201556132 [0.4603889573963417, 1.013126948637089] | Thirty-two paired target blocks retaining all repeats |

The crossed target/cohort sensitivity interval is
[0.4521682932867005, 1.0460777825996]. Both online intervals include 1.
The descriptive geometric mean of paired target ratios,
0.3156223546296918, must not replace the predeclared arithmetic cost estimator
to select a favorable conclusion. The actual complete-process comparisons
charge preparation, solving, built-in verification, output and supervisor
observation overhead. They do not use a sum of synthetic cold rho launches as
the batch comparator.

The native complete-panel receipt passes all 32 planned processes, 218 scalar
answers (104 IC, 114 rho), 104 IC four-point witnesses and five same-work rho/rho
A/A pairs, with zero invalid or missing answers. This comprises 192 cohort
answers, 16 actual singleton answers and 10 A/A answers. The extractor retains
32 process rows, 192 cohort target rows and 26 auxiliary rows; no target,
seed, regression or slow case is omitted. All process receipts report exit 0,
completed verification and `timeout:false`. The memory-limit-hit field is
`null` in every receipt: the monitor does not identify allocation refusal from
exit code alone. State the enforced 12 GiB RLIMIT_AS envelope and successful
completion, not a measured zero-cap-hit count.

The main rho panel retains 4,488 internally capped walks and all 1,001,130,567
walk steps. Including A/A, the recorded totals are 5,046 capped walks and
1,123,246,029 steps. These charged internal restarts are distinct from a failed
target or a whole-process timeout. A/A A/B wall ratios range from
0.9587596843493773 to 1.1450331951782278. They are a local descriptive control;
they neither pass an invented 5% gate nor define a universal noise ceiling.

## Costs, tails and extrapolation limits

The exact legacy model is `icv1-f2m83-t6151469093347-9eba6cd5`, with subgroup
order 8569786107849059. Field degree 83 does not mean an 83-bit subgroup. The
conditional target law excludes the 617 declared retained/prior signed orbits;
the eight cold targets are the first eight of this fresh panel. This is not an
unconditional subgroup-wide or cross-curve sample. IC rank seed 29 and both
method policies are fixed before the campaign. There is no new fastest-seed,
policy or target selection on this panel.

Three shared IC preparations do not create 96 independent preparation or host
samples. The batch interval has only three cluster levels and correspondingly
limited uncertainty resolution. Guest isolation records 2,654 samples and zero
observed contention on CPU 3; it does not prove exclusive physical-host use or
independent-host replication. Both methods use one OS thread and the same
2,100-second process cap. Stock IC's online timer omits validation/output
intervals included by rho, so online remains secondary without a silent clock
correction. Cold begins with the retained K600 artifact; its historical
selection/construction cost and OS-cache flushing are outside this boundary.

All native tails, CPU and memory rows remain available. The 96 pooled IC online
draws have median 0.833230443 s, p95 6.67799166125 s and maximum 9.482010347 s;
rho's corresponding values are 2.1688249195000004 s, 5.0169968780000005 s and
6.193311533 s. These are descriptive observations, not population-tail bounds.
The 32 three-repeat IC target means reach 8.590054447 s; the rho means reach
4.016571018666667 s. Slow IC targets remain in the arithmetic estimator.
Mean batch CPU is 1232.8835199999999 s for IC and 78.941576 s for rho.
Maximum batch process peak RSS is 5,811,691,520 bytes for IC and 70,823,936 bytes
for rho; maximum IC cold peak RSS is 5,811,724,288 bytes. Peak process RSS is not
per-target memory, address-space usage or simultaneous memory across serial
runs. CPU, wall and RSS have separate units and are not combined into an S.

The separately versioned native presentation supplement binds the original
producer and execution files and exports all 104 IC phase/residual rows,
114 rho phase/counter rows and 32 process records. IC signed residuals are
retained rather than clipped, assigned to a sixth phase or added twice; their
recorded range is 0.0002999999999993008 to 0.01361100000000004 ms. Rho's phase
identities use exact integer nanoseconds. Nested setup/index/load timers are
not summed twice. This descriptive supplement does not recompute a primary
estimator or change the frozen statistics.

The additive model correctly separates rho's shared curve/basis preparation
from its seed-specific G-only jumps. It charges 0.00122437653125 s of mean jump
preparation per target in the N-dependent slope, rather than retaining all 32
jump tables as a free fixed intercept. A fresh DP/walk state remains per target.
Reusing one common jump table for every Q would be another unmeasured policy.
The model's point threshold is 1,685 targets, beyond the largest measured batch
of 32. Its 814/20,000 nonfinite draws (0.0407) remain in the extended interval:
the lower endpoint is 844 targets and the upper endpoint is no eventual IC
advantage in this model. The finite-only quantiles are explicitly conditional
and cannot be advertised as an unconditional finite confidence interval.
Different online boundaries and excluded residual/audit/base-construction
costs also qualify that model. No actual crossover above 32 was run.
The linear projection does not impose the campaign's 2100-second process cap.
Some projected times exceed that cap; a larger-N validation requires a new
prospective workload and a common budget for both arms.

## Source, correction and verification scope

Reviewed the immutable `PROTOCOL.md` and `provenance/protocol.json`,
`native_runner/src/main.rs` and `src/monitor.rs`, `analysis/src/lib.rs`, the
measured-source manifest and worker sources, and the final input/statistics,
process/target matrices, replay and presentation receipts. The measured-source
hash manifest still matches the current measured source files. The frozen
analysis executable remains
`b1820abad3d17cb4f40f7d4d1a0a6c92ade2dfa056a43a0c4eee2d315c3a8346`.

The original v2 audit failure and frozen `collect` exit 2 are retained. The
v2 receipt's 227 invalid findings concern copied metadata; they are not 227
failed solver answers. The independent diagnosis found 127 differences, all
exactly one decoded binary64 ULP. The separately versioned v3 checker changes
only four authenticated metadata-copy binding contexts, allowing equal or
adjacent finite f64 values at structured allowlisted paths. Integer/u128,
scalar, coordinate, seed, counter, hash, string, key, type, order, matrix,
resource and A/A gates remain exact; the prior derived-arithmetic checks are
unchanged. The native self-test passes four positive cases, 30 negatives and
442 actual copy bindings, including rejection of two ULPs and a one-nanosecond
change. Original v2/frozen sources and executables are not relabeled as passing.

Reviewed `verification/revisions/20261015-v3/verify_matched_campaign.cpp`,
`numeric_binding.hpp`, diagnosis and self-test receipts. The v3 source hash is
`b71469aedf414e48b671675d14721c6a074dcb37200f5c3e57945b0d86efd881`;
its binary hash is
`c26429c20dffa01428709fa206c4e4fb77e4e83d92f352f996fe5c95dc226acc`.

Reviewed the separately versioned `postprocessing/collector-only/src/main.rs`,
config and six passing tests. It pins the full passing v3 receipt, checker,
helper, v2 source, protocol, matrix and original runner; rechecks all 336 unique
snapshot hashes and the exact sealed campaign tree; and hashes/parses the same
cached bytes. It copies already verified observations with arbitrary-precision
number tokens, checks identities/seeds/rank29 and output counts, and never
launches a solver, reconstructs `observe`, changes a false producer verification
flag or recalculates a clock. Its only subprocess is sha256sum. The receipt
retains `original_frozen_collect_succeeded:false`. The extractor source hash is
`c80c62b0da7361bfe2d005a30a4dbae6ae9c014caa43c9ec88d77297671c0177`;
its executable hash is
`32c2d71c5d9215ef154f281d3926b6e31e81c5cdbe6cfa45b856c4a06a3c3c58`.

The final native evidence is committed at
`fb26ccaef5c8f5e97aa81b3503d42e5501b030e4`. Checked content hashes:

| Artifact | SHA256 |
|:--|:--|
| `statistics_input.json` | `8c8b380646e1d1d561492634fc6bdfd10d78e34395483a0f9df7b94c599bf2d9` |
| `statistics.json` | `1ba0dd87653daf737e95add3a37d08ceb905f8d9811b933db61fc156cb83377d` |
| `independent_replay_v3_receipt.json` | `94891f01e158a3a365744a2d349a685d3f1fd844a16e7f5b40e2e38763646c8b` |
| `presentation_supplement.json` | `4804814de5e0fbb7d464b09c26738844314a045cbb60cc33379ab7a9909ef0e8` |

These bindings establish the reviewed public known-answer experiment and its
exported scalar/decomposition witnesses. They do not independently validate a
full unexported IC rank/log transcript. The post-hoc native audit is separately
metered outside solver intervals; no per-arm audit-inclusive production ratio
is inferred. Common calibrated S, sealed legacy-wide ecbench admission,
physical-host admission, independent-host uncertainty, global optimality,
an asymptotic exponent improvement and deployed-curve security conclusions
remain unestablished.

## Canonical publication review

The diagnostic is copied into `docs/ic/progress-timeline.json`, its canonical
embed/table and a source-linked panel in `docs/index-calculus-scoreboard.html`.
Three new context series distinguish actual cold, actual 32-target process and
secondary online; each has `barrier_test:false`. All plotted values and
intervals are copied directly from native fields in IC/rho direction. The date
2026-10-15 is explicitly provisional while PR1828 is unmerged and human review
is pending; the final merge workflow must refresh it if the landing day differs.
The measurement date remains 2026-10-14. Historical observations are retained.

The leaderboard builder adds a separate diagnostic block and requires the full
v3 gate plus receipt/statistics/supplement hash bindings. It preserves full
native binary64 values in that block rather than passing them through the
existing display rounding. All pre-existing non-provenance keys in
`docs/ic/leaderboard.json` equal the saved pre-change JSON, including every
board row, leader, family, floor, oracle, exponent, roster and prior diagnostic.
New diagnostic S, floor-ratio and controlled-speedup fields remain null.
`docs/ic/LEADERBOARD.md` and `docs/ic-leaderboard.html` are regenerated from the
builder. `docs/browser/data.json` changes only its leaderboard source hash.
The required derived panel index in `docs/ic-current-state.html` is regenerated
with the existing, unchanged `scripts/site/panel_index.py`.

The canonical boundary ledger `docs/ic/boundary_targets.json`, boundary
scoreboard `docs/ic/BOUNDARY_TARGETS.md`, ranked recipes and admitted-best
tables are unchanged. No source result is inverted to fit an old axis, no
reciprocal/statistical calculation is performed in the Python presentation
builders, and no existing measurement timestamp is replaced.

Validation: leaderboard regeneration and `--check`, lab-browser regeneration
and `--check`, and the panel-index generator and `--check` pass. The first site
test invocation used the wrong import working directory; the correctly scoped
run exposed the required stale panel index, which was regenerated. The final
full site suite passes all 48 tests after that correction (4.806 seconds).
The 12-file measured-source hash check also passes. Conductor was
attempted before these scoped edits and before the derived index expansion;
the executable is unavailable (`command not found`). Root announced and
approved the single generated index scope expansion. No measurement or
experiment was rerun for this review.
