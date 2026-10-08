# Matched n37 native six-sum batch versus signed-Frobenius rho

Status: protocol committed before this paired cost run. The public b01 points
and the successful native six-sum recovery are already published, so this is
a **retrospective engineering diagnostic**, not held-out selection evidence.
It answers whether the complete descendant-native 1,024-target batch is
competitive with the strongest available same-point signed-Frobenius batched
rho on this host. The primary one-target online IC/rho question is separate.

Native-driver amendment, committed before the repeat panel: the first
20-arm pilot used a Python process driver and audit script. Repository
`AGENTS.md` disallows that execution path for research comparisons. Preserve
that pilot's raw archive and source commit as historical provenance, but do
not use its timings for this decision. Repeat the same frozen schedule and
unchanged producer binaries with the Rust
`examples/n37_native_batch_rho_panel.rs` driver. The Rust analyzer and
published independent replay must validate the repeat panel. The new raw
directory and receipt remain separate from the pilot; no pilot rows are
selected or substituted into the repeat.

Freeze repository parent `80a05cc22ade46e5ddd5918bb5a01da9d66c0d3c`
and the three unchanged example sources:

| Executed component | SHA-256 |
| --- | --- |
| `examples/n37_native_m6_residual.rs` | `086a6fbdd12e2b104ec18f809d4048453d747ad76c5a74f695befd08d4efffa8` |
| `examples/n37_native_m6_residual_replay.rs` | `68a4ce03029b82821fdabb9854a14559dde243bd5ac0f409de884df032f7394a` |
| `examples/koblitz_rho_batch_ks_v3.rs` | `98b6e8a27d821ebfc7ae716410184dffe89aadf27dfe36bcac8ce834c39cf04c` |
| `Cargo.toml` | `f88ac8c9527bc2403cc550b8b00aaa778630911a4c9d1de47f67dbbff055d6a2` |
| pinned `Cargo.lock` from `n37_native_m6_mitm_20261002` | `b28d3c2d81146a40d00df85f45c6460d9e2bc5c26875307a06cfd15efff99365` |

Use Rust/Cargo 1.93.1 and build both methods once in release mode with
`--locked --offline`. The pinned lockfile is copied to the checkout root
before the build. Archive build output and binary hashes. The native solver
reads its own SHA-pinned base, degree-73 archive, freeze, and b01 public
point file; the rho arm reads only the same b01 point JSONL. That file has
SHA-256 `78553fdff5ae66521d3a1052978962258e78d48c92285df6e1df962034b43717`.
The separate fixture scalar file has SHA-256
`3453995c0423e4911ad4a6afa7cfe50ec1d6ed680df` and is available only
to post-run verification. Neither timed producer may read it.

Rho command: `koblitz_rho_batch_ks_v3 37 0 signed_frobenius 1024 2026100110101`,
with `KIC_RHO_POINT_INPUT` set to that exact public file,
`KIC_RHO_BATCH_CORPUS=compact-disjoint-cold-v2-n37-L1024-b01-20261001`,
`KIC_RHO_CANON_BACKEND=normal_basis`, and `KIC_RHO_DP_BITS=8`.
This uses 32 parallel walks and a shared distinguished-point table across
the 1,024 targets. Native command: `n37_native_m6_residual <fresh-output>`
with its published 42-column, bounded 16-shift, full-rank policy unchanged.

Run five fresh-process blocks on the same host. Blocks 0, 2, 4 execute
`IC, rho, rho, IC`; blocks 1, 3 execute `rho, IC, IC, rho`. Each of the four
arms gets a fresh process, empty in-process caches, a 180-second wall cap,
and a 4 GiB peak-RSS acceptance cap. Record every arm's command, environment,
exit/timeout, wall, user/system CPU, peak RSS, stdout/stderr and complete
output, including failures. Compute each block's IC and rho geometric mean
of the two complete child CPU costs and its IC/rho ratio. Also report each
method's within-block duplicate ratio and all five paired ratios. A duplicate
drift over 10% flags that block's timing as noisy; never discard it or select
replacement runs. This macOS host lacks an auditable exclusive CPU partition,
so **all CPU and wall ratios remain exploratory even if drift passes**.

After all timed arms, run the published independent full-rank native replay
on every successful IC output. For each rho output, require 1,024 indexed
public-point rows in input order, every `verified=true`, recovered scalars
equal to the independently published b01 fixture, and `published_q` equal
to both the public input and fixture. Cross-check the IC and rho recovered
scalars on each shared Q. Preserve raw failures, replay failures and timeouts;
a missing or unverified answer is not a completed batch. Report phase counts
and native addition requests and rho charges separately. They are not yet a
common calibrated operation unit, so `S=operations/sqrt(r)` and any
algorithmic speedup remain unset. The batch result must never replace the
single-target online comparison.

Report complete-process cost and uncertainty for this exact fixed policy.
If IC loses, call it a diagnostic fixed-cell result, not an ECC2K-130 or
degree-131 no-go. If it wins, require disjoint-Q confirmation and isolated
benchmark evidence before promotion. Archive exact source/input/binary
hashes, every timed arm, independent replay, analysis and decision in the
follow-on commit and PR.
