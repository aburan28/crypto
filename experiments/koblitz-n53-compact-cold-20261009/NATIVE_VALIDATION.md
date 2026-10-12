# Native n53 comparison path: correctness gates

The native protocol reconstructs the frozen n53 target point
`[6322155974735900, 5110849281364254]` from its SHA-256 seed and subgroup
generator, yielding scalar `7948768810114`. The native relation verifier
replays the archived n13 certificate to the same JSON bytes as the original
archive, including the base, four signed-Frobenius orbits, four rank-increasing
relations, S3 pair roots, recovered representative logs, and target scalar.

## Local validation on 2026-10-09

| Check | Result |
| --- | --- |
| `koblitz_n53_compact_protocol verify-workload` | `PASS`, frozen Q and scalar above |
| `koblitz_n53_compact_replay --run-dir validation/n13-smoke` | `PASS`; output byte-equal to archived `replay.json` |
| `cargo test --release --example koblitz_n53_compact_replay` | 4 passed: archive and tampered modular row, S3 root, target scalar |
| `cargo test --release --example koblitz_rho_fixture` | 8 passed, including the strong backend over its fixture rungs |
| `koblitz_rho_fixture 13 0 signed_frobenius 1 strong 20261009 7` | exit 0, verified scalar 7 on `[384,1476]` |
| `cargo clippy --release --example koblitz_n53_compact_replay --example koblitz_n53_compact_protocol -- -D warnings` | exit 0; repository toolchain emits the pre-existing unknown `chunks_exact_to_as_chunks` lint warning |
| `rustfmt --check` on the three native files | exit 0 |
| `koblitz_n53_compact_protocol preflight` in approved process-inspection context | `PASS`, sampled own RSS `1802240` bytes |
| Committed-freeze `run` with deliberate ambient `KIC_SENTINEL=1` | exit 2, preflight failure receipt; no IC or rho arm files created |

The ordinary sandbox denied `/bin/ps` with `Operation not permitted`; the
native preflight rejected that context. This is the intended fail-closed
behavior. The approved process-inspection context passed the same preflight.
The deliberate ambient-variable attempt verified the committed freeze and
binary hashes before writing
[`validation/native-freeze-preflight-20261009/receipt.json`](validation/native-freeze-preflight-20261009/receipt.json).
These checks establish input, replay, and control-path correctness. The n53
timed comparison still requires the committed source/binary freeze, clean
prerequisite PRs, full-rank producer result, independent replay, and paired
same-point strong-rho result under the frozen resource limits.

The archived n13 receipt and original Python scripts remain historical
evidence. New comparisons use `native_run.rs`, `native_replay.rs`, and
`reference.rs`; the latter implements independent polynomial-basis group
arithmetic rather than calling the producer's curve operations.
