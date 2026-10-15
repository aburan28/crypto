# N83 root Bloom profile records

This directory seals the six-arm comparison in
[`PROTOCOL.md`](../n83_bloom_profile_20261015/PROTOCOL.md). Three candidate
manifests distinguish the exact Bloom size and hash count on the same K3300
factor base and root11 index. Six measurement rows retain two fresh processes
per profile and the common 4,000-point public panel.

The executed source commit is `42fc71945cd7087f78464abc65c4bae7c52a93e8`.
The public-point SHA-256 is
`35e8b43966d15065cf0186ace762f823a602aa1be7d7c0a9039e020550d01491`;
the known-scalar sidecar SHA-256 is
`6e1181ba958708803316903a645463425bc79130074041d39c8094778e2f4379`.
The sidecar is used for independent replay outside the measured intervals.

The sealer checks the exact raw ZIP and independent analysis digests, source,
fixture, base, Bloom profile, resource envelope, and analysis identities.
Repeating the command verifies byte equality of each canonical record. The
[result report](RESULT.md) records the selected profile. Artifact `11969951389`
from run `38718569391` has ZIP SHA-256
`d4f7c9f0b773eaa05e0f18043c8d11503a456d5aa7a5a426e4e891bd9959c43d`;
the independently replayed `analysis.json` SHA-256 is
`f15dc07f8b816038d824ef56962caf46e062a568e1ca998134cfd08f4f2119b6`.

```sh
python3 research/n83_k3250_sparse_fresh_20261015/generate.py
python3 research/n83_k3350_sparse_fresh_20261015/generate.py
python3 research/n83_k_abba_20261015/generate.py
python3 research/n83_bloom_profile_20261015/verify.py /path/to/results
python3 research/n83_bloom_profile_results_20261015/materialize.py \
  --artifact-zip /path/to/raw.zip \
  --artifact-sha256 d4f7c9f0b773eaa05e0f18043c8d11503a456d5aa7a5a426e4e891bd9959c43d \
  --analysis-sha256 f15dc07f8b816038d824ef56962caf46e062a568e1ca998134cfd08f4f2119b6 \
  /path/to/results
```
