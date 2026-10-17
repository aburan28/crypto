# Hosted N83 fixed-witness SAT-access result

[GitHub Actions run 39051579293](https://github.com/aburan28/crypto/actions/runs/39051579293)
completed at exact source/workflow commit
`0891c2e062ff95e923302e570a8101f83a59b6a5`. Its artifact ID is
`12107379138`, named `n83-hosted-witness-39051579293-1`. The input is the
earlier independently checked, **known-answer** K1182/K2048 ordered
five-point witness: five summand x-coordinates and three affine partial-sum
x-coordinates were fixed before SAT solving. The exact retained domain input
has SHA-256 `63accc5b9dc9de6327c5b07e93b3e62717f8376a96c3f91f7bb70b77095f23bb`.

The worker reported **`SAT_verified_witness_x`** after **787 of 5,000**
allowed conflicts. Its decoded model passed original chained-S3 verification
and matched every fixed summand x-coordinate. The separate [retained audit](audit.json)
also recomputes the ordered point sum and all four original S3 residuals from
the pinned fixture. This establishes a reachable valid model for the exact
fixed domain, first identity pattern and first-pair order. It does not measure
how often natural query targets admit or yield a partial relation.

| Measurement | Hosted ARM64 result |
| --- | ---: |
| Circuit built | 23.632 ms after worker start |
| Exact orbit domains and first-pair order installed | 82.035 ms |
| Eight fixed x-blocks installed | 82.098 ms |
| SAT conflicts | 787 |
| Worker wall | 214.870 s |
| External solve wall | 215 s |
| Variables / clauses / XOR rows after fixed units | 101,636 / 714,736 / 4,980 |
| Orbit-domain clauses | 465,738 |
| Resource envelope | 1 CPU, 4 GiB memory, no swap, no network |

The container exited 0 without OOM and was removed. The GitHub job ran from
17:03:03 to 17:11:12 UTC, about 489 seconds including checkout, image/input
freeze, dependency fetch, build, solve and artifact upload. The conservative
hosted charge is the full **2,400-second job envelope**; observed phase walls
are reported separately. This control used a different host and execution
shape from the earlier capped natural-search trials, so its 214.870-second
worker wall is not a matched factor-base runtime comparison.

The [protocol audit](protocol-audit.json) replays the GitHub artifact's 22
files, source/input hashes, binary, image, resource configuration, phase
markers, raw solve result, exit code and cleanup. Its SHA-256 is
`ec2f5f2c8bc7429d343be521bc97feeea0ca39bd2b38aabb9e1518c1c17140e9`.
The [independent replay](replay.py) adds source-pinned reference group/S3
arithmetic, exact result counters, and a complete archive-member comparison.
The [S3 roundtrip receipt](s3-roundtrip.json) records a checksum-protected
upload and fresh GET of the content-addressed 474,794-byte artifact at
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/a0/evidence/hosted-fixed-witness-run-39051579293/75362a2d12c48a6e7854eaf070d9430a5fbb9d3f484229f167fcff44a621135d.tar.gz`.
The downloaded bytes have SHA-256
`75362a2d12c48a6e7854eaf070d9430a5fbb9d3f484229f167fcff44a621135d`;
their extracted protocol audit is byte-identical to the GitHub-artifact audit.
Two corruption controls reject a changed solver status and a changed S3
archive byte.

To replay with the downloaded GitHub artifact directory and S3 tarball:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/hosted-fixed-witness-run-39051579293/replay.py \
  --artifact ARTIFACT_DIR --s3-tar DOWNLOADED_S3_TARBALL --check
N83_WITNESS_ARTIFACT=ARTIFACT_DIR N83_WITNESS_S3_TAR=DOWNLOADED_S3_TARBALL \
  python3 -m unittest discover \
  -s research/koblitz_n83_factor_base_sweep_20261008/verification/hosted-fixed-witness-run-39051579293 \
  -p test_replay.py -v
```

This control resolves the narrow question of fixed-witness encoding access.
The uncapped domain search, natural partial-relation yield, large-prime graph
closure, useful rank, linear algebra, individual-log extraction and complete
single-target IC runtime remain unmeasured. The factor-base winner is still
unselected; `selected_best_total_runtime` remains null in the saved receipts.
