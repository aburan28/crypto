# Retained-N83 mixed-orbit build and import gate

**Result: `PASS_import_only`.** From clean source commit `d68cee586ae4b59b18ac5885ef064fe3fe709772`, the pinned Linux arm64 offline build produced the mixed-orbit worker. A separate no-network, one-CPU, 4-GiB, zero-swap container imported the exact retained K=1,182 small base and K=2,048 envelope. The worker checked 196,212 small points in 1,182 signed-Frobenius orbits and 143,756 residual points in 866 signed-Frobenius orbits. Its receipt states `solver_search_executed=false`.

| Stage | Outer wall | Disposition | Conservative charge |
| --- | ---: | --- | ---: |
| Offline Linux build | 177.906 s | Exit 0; Linux ELF; cleanup verified | 600 s |
| Exact nested-base import | 7.172 s | Exit 0; worker receipt PASS; cleanup verified | 240 s |
| Receipt and cleanup reserve | Included in 186.157 s supervisor wall | Frozen input and source audit PASS | 60 s |
| Total | 186.157 s | `PASS_import_only` | **900 s** |

The full 900-second additional allowance is conservatively charged, bringing the cumulative local charge to **12,600 of 12,600 seconds**. The 186.157-second elapsed supervisor wall is a stage diagnostic, not a complete cold index-calculus runtime. No relation search, partial yield, graph completion, rank, sparse linear algebra, individual logarithm, or matched-rho comparison ran. A best base is still unset.

The [outer receipt](outer.json), [budget](budget.json), [authorization](authorization.json), [worker import receipt](verify.stdout), and raw [build stderr](build.stderr) and other stage logs retain the run. The separate [audit](audit.py) passed against the exact committed source blobs, recorded input and log SHA-256 values, stage commands, image, cgroup settings, container exits/cleanup, and import receipt. It independently checked the two local compressed object SHA-256 values: `eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380` and `aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a`. Both retained S3 keys were also present when checked with `aws s3api head-object`, at 5,956,114 and 10,313,647 bytes. The head responses establish current object presence and lengths; the earlier upload/download receipts and local SHA-256 check establish byte identity for the copies used here.

This is an implementation and input-admission gate for the [separate search protocol](../mixed-orbit-supervisor-source-20261012/PROTOCOL.md). That protocol was frozen when the cumulative charge was 11,700 seconds. Before any later search launch, its authorization accounting must be updated to the now-charged 12,600-second baseline, with a new run-specific grant and clean source receipt. The current 900-second grant cannot authorize relation search.
