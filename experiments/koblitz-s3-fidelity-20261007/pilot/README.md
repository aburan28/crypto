# Frozen pilot inputs

These are **inputs**, not timing results. There are twelve distinct public
subgroup points for each of the registered Koblitz instances at field degrees
41 and 53. `fixtures.json` holds the known scalars for independent checking;
each solver invocation receives only one `Tnnn/public_target.json` file.
`workload_id` is null until the complete Linux curve/candidate and workload
records are frozen. No measured `IC1` identity is claimed by this directory.

The exact generation commands were:

```sh
./target/release/s3_pilot_targets 41 12 s3-pair-pilot-20261007-v1 \
  experiments/koblitz-s3-fidelity-20261007/pilot/n41
./target/release/s3_pilot_targets 53 12 s3-pair-pilot-20261007-v1 \
  experiments/koblitz-s3-fidelity-20261007/pilot/n53
```

The seed law is implemented in [`s3_pilot_targets.rs`](../../../src/bin/s3_pilot_targets.rs):
`1 + int_be(SHA256("s3-pair-target-v1:" || seed_label || ":" || u32be(n)
|| u32be(index))) mod (r−1)`. Inputs are immutable; a repeat invocation
refuses to replace its output directory. Verify the file list from this
directory with `shasum -a 256 -c SHA256SUMS`. The two fixture manifests have
SHA-256 values `16fa23c01e6e6b0ae133079c84777ddfb49f6683415c01dc8c67ce8906cafbc3`
and `5f18351f40dd13e97dcf8921b7784f110abe963c418601df78a1cceac18c5112`,
respectively.

The pilot must be measured on the qualifying Linux host before the
confirmatory target count is chosen. The confirmatory seed label and count
remain unset so no result-dependent target selection can occur now.
