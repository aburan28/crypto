# Rank-243 S3 preparation across 33 matched N41/N53 targets

The source-identical rank-243 preparation reduces online S3 root calls by a
geometric mean of **402.84× on N41** and **378.65× on N53** across 33 matched
single-target indices per curve. Every one of the 264 paired runs verified;
independent Sage replay recovered all 264 scalars and checked 756 distinct
relation witnesses. The matched N41/N53 root-call gain ratio is 1.064
(95% target-index bootstrap interval 0.540–2.054), so these 33 pairs do not
resolve which curve gains more from preparation.

[`FREEZE.json`](FREEZE.json) selects the previously unused T025–T033 indices
on both curves: **18 public points and 72 fresh-process solves** in the same
ABBA/BAAB order and with the same 14 workers and 180-second per-run cap as
the [first 24-pair pilot](../koblitz-s3-crosscurve-20261009/README.md).
The original Rust snapshot `b562cdaed14254d4f6c5f2a4640b302592e8882b`
was rebuilt offline. Both release executable SHA-256 values exactly match the
first pilot: `c4ff99ace090551747b2a1083b9546754e0ce1814f309a47e46d2ab0bd132ca7`
for the reference and
`50d3ad29fbd28c146ed60cbd6a255638fc6d48733ce08bd8097b08e03acbd195`
for rank-243 preparation. The source, toolchain, binary, original archive,
runner, and analyzer identities are all bound in the freeze.

The [runner](run.py) checks those hashes before launching the original frozen
runner. The [analyzer](analyze.py) binds the old 192-run ledger to the new
72-run ledger and retains every run. After all 264 solves verify, it
bootstraps 20,000 matched target-index resamples for each curve's exact
online S3 root-call ratio and exploratory online wall ratio, plus the ratio
of N41 to N53 gains. The operation boundary and seed were fixed before
execution. A paired 95% interval excluding 1 supports a difference between
curves for that metric. Independent Sage replay of the pooled ledger checks
the scalar and relation equations.

Target-independent rank preparation remains outside the single-target online
interval; its cost is retained separately. Local wall timing is exploratory
without the repository's host-level CPU isolation and noise receipt. The
root-call count is a deterministic target-dependent stage diagnostic, not an
end-to-end wall-time substitute.

## Measured result and fidelity

The extension completed **72/72 verified fresh-process solves** for T025–T033.
Together with the hash-bound first panel, it contains 33 targets per curve,
two variants and two repetitions per target. The operation boundary is the
target-dependent S3 root calls in each complete online solve, including
fallback attempts. Ratios below are reference divided by prepared; intervals
resample the 33 target indices 20,000 times using the frozen seed 20261010.

| Curve | Root-call gain, geometric mean [95% interval] | Online wall gain, exploratory [95% interval] | Minimum per-target root-call gain | Online fallback attempts |
| --- | ---: | ---: | ---: | ---: |
| N41 | 402.84× [248.05, 656.93] | 165.75× [111.64, 241.71] | 6.68× | 84 |
| N53 | 378.65× [266.65, 544.11] | 118.69× [83.93, 170.10] | 67.95× | 0 |

N41 target T028 accounts for all 84 fallback attempts: 42 in each
repetition. Those attempts are included in both the operation counts and
online wall intervals. The exploratory matched N41/N53 **wall-gain ratio** is
1.396 [0.789, 2.385]; it also does not resolve a curve winner.

The prepared A/A 95th percentile absolute log wall ratios are 0.961 on N41
and 0.738 on N53. Both exceed the 10% noise gate, and this host lacks an
exclusive physical CPU partition receipt. The
[`analysis.json`](local/analysis.json) therefore records
`controlled_speedup: null`. The online wall ratios above describe this local
run; a ≥2× controlled single-target wall-time result needs the frozen panel
run on an auditable isolated host.

## Artifacts and reproduction

The [extension raw archive](local/raw_extension.tar.gz) contains the 72-run
ledger and per-process stdout, stderr, and report files (218 files total;
SHA-256 `6aa0c63b68a4e8673c4d24e92e5a79e88b39b5d36d096fec1be9859341e07973`).
The [completion receipt](local/completion.json) binds its ledger hash
`7fa08df44a7d044810c6876b953b06c3dea843cc899f60ebbfaa51964b8e5926`.
The [pooled stage analysis](local/prior_analysis.json),
[33-pair bootstrap](local/analysis.json), and
[independent replay](local/independent_sage_replay.json) bind the 264-run
ledger hash `da63c4f19710fc36b9612d4588322c90bd24aa4619151151044333c400ce6c92`.
The checked Sage [runtime receipt](local/sage_runtime_info.json) was saved
before replay.

The extension was launched with the rebuilt original binaries using:

```sh
python3 experiments/koblitz-s3-crosscurve-extension33-20261010/prepare.py \
  --binaries /private/tmp/s3-crosscurve-b562-build/release/examples --verify
python3 experiments/koblitz-s3-crosscurve-extension33-20261010/run.py \
  --binaries /private/tmp/s3-crosscurve-b562-build/release/examples \
  --output /private/tmp/s3-crosscurve-extension33-20261010
```

After extracting the raw archive, `analyze.py --extension <extracted-dir>
--output <fresh-output-dir>` reconstructs the pooled ledger and analysis from
the first pilot's frozen raw archive and this extension. The original
[`replay_sage.py`](../koblitz-s3-crosscurve-20261009/replay_sage.py), launched
through `/Volumes/SSD990/cryptanalysis/sage -python`, rechecks the pooled
ledger against the saved Sage runtime receipt.
