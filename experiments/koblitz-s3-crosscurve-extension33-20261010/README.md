# Nine-pair extension of the frozen N41/N53 rank-243 panel

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

At freeze time the extension has **0/72 launched runs**. The intended command
uses the rebuilt binaries and writes raw data outside the checkout:

```sh
python3 experiments/koblitz-s3-crosscurve-extension33-20261010/prepare.py \
  --binaries /private/tmp/s3-crosscurve-b562-build/release/examples --verify
python3 experiments/koblitz-s3-crosscurve-extension33-20261010/run.py \
  --binaries /private/tmp/s3-crosscurve-b562-build/release/examples \
  --output /private/tmp/s3-crosscurve-extension33-20261010
```
