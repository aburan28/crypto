# Input preflight for the pending N83 known-witness control

The prepared [v7 control](../mixed-orbit-control-reuse-19800-20261013/RESULT.md)
can use the pinned K=1,182 and K=2,048 signed-Frobenius archives without a
new Linux worker build. This [preflight receipt](preflight.json) verifies the
two local archive files against both SHA-256 and BLAKE3 from their frozen
manifests and verifies the decompressed BLAKE3 through each gzip stream. It
also rechecks the prior 19,800-second evidence chain. **No control SAT call
was launched.**

| Input | Compressed bytes | SHA-256 | Compressed BLAKE3 |
| --- | ---: | --- | --- |
| K=1,182 | 5,956,114 | `eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380` | `76e1cbe4ea5a9d26e5b24a2d27f9401271dc471846529190a0416d8e09837a92` |
| K=2,048 | 10,313,647 | `aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a` | `82990cda5d57b1c4c4b635673a881ad233c24f707e7ac644a11478aeddddf851` |

Read-only S3 HEAD calls on 2026-10-13 returned those two lengths. The
K=1,182 ETag was `ee87502a43c5775af9041666473e6699`; the K=2,048
multipart ETag was `96006317b77bf3a271eb8330ab209ce5-2`. HEAD metadata
is an availability snapshot, not archive-content proof. The K=1,182 object
was then downloaded from the manifest's S3 URI to a local cache and matched
both pinned compressed hashes. The existing K=2,048 local object matched
the same checks. The saved receipt stores the full S3 URIs, source hashes,
plain BLAKE3 values and sizes without a machine-specific cache path.
The raw [K=1,182 HEAD](s3-head-k1182.json),
[K=2,048 HEAD](s3-head-k2048.json) and
[Docker image-ID output](docker-image-id.txt) are retained beside the
receipt, and the generator pins their bytes by SHA-256.

`reuse_control_supervisor.prior_evidence_error` returned no error for the
prior audits and the reused Linux ELF. A read-only Docker image inspection
found the pinned local image
`sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31`.
The [preflight generator](preflight.py) pins both manifests, the v7
supervisor and its independent auditor by SHA-256; its two focused tests
reject modified compressed and decompressed archive data. Reproduce the
local-byte check with a Python environment containing `blake3`:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/control-v7-input-preflight-20261013/preflight.py --check --small-object "$N83_K1182_OBJECT" --envelope-object "$N83_K2048_OBJECT"
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/control-v7-input-preflight-20261013 -p 'test_*.py' -v
```

The frozen v7 supervisor still requires a separate fresh 1,200-second
authorization record before creating an output directory. This preflight
does not supply that record, does not validate solver access, and leaves the
conservative local retained-run ledger at 19,800 seconds. Natural useful
yield and the minimum complete cold index-calculus runtime remain unknown.
