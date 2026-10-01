# Prepared F5 native source and runtime

This follow-up implements the missing source/build/registration and frozen-audit
transport for the target-only F5 adapter accepted in
[PR #1093](https://github.com/aburan28/crypto/pull/1093), merge
`5f0cec874b5b627b6e9c4bb598c7d04fe13e4eb3`. It does not consume a new production
invocation or generate a fresh point. Use the [control protocol](PROTOCOL.md).
Implementation acceptance, separately frozen F5/SAT development executions and
the fresh paired incumbent/rho comparison remain distinct gates.

The native solver and online observer are unchanged. Its Rust unit fixture now
contains only mathematical state, ordered geometry and column logs; the old
entire F5 certificate included historical preparation evidence in the
conservative build manifest. Tests independently compare the new mathematical
fixture to both accepted certificates. The new manifest includes that fixture
and excludes both whole certificates. All historical inputs, source hashes,
validators, archives and consumed registrations remain unchanged.

`prepared_f5_inputs_v1.py` verifies a new controlled worker build and retains all
root and registry dependency source bytes before sealing assets.
`prepared_f5_runtime_v1.py` freezes canonical math/code identity and a separate
one-use invocation containing the certificate, supplied point, seeds and limits.
Its entrypoint uses the existing isolated Python/native wrappers, registered
stdin and controller watchdog. Only disclosed n17 controls are admitted.
The independent auditor checks preparation, source gates, exact native argv,
all failed attempts, observed F5 dispatch, group witnesses, scalar replay and
exclusive online phases. Timeout/error rows retain unknown verified costs.

The audit CLI uses the same auditor source and interpreter frozen **before**
the native execution. It starts a fresh `-I -S -B` process in the read-only
extraction and checks loaded modules before/after replay. It executes no native
solver, changes no execution evidence and never restarts a failed invocation.
The exact public-input and five-phase online interval remains the primary
contract; whole-process resource diagnostics are separate and uncalibrated.

The new retained [macOS ARM64 assets](native-inputs-macos-arm64/manifest.json)
are 12,228,690 compressed bytes and contain the binary, build receipts, complete
root sources and complete registry dependencies. These are local source/build
receipts, not remote attestation or a hermetic toolchain claim. They were built
and staged without running an IC query. Native input identity is:

| Item | SHA-256 |
|---|---|
| Canonical Rust source manifest | `09321707259fea2336bca8b81525473e9fa3f1632974e20996c41bbcbb9f972b` |
| Canonical controlled build policy | `b1b4dd927a1376d773dc9edd74efb6ba0bdcac6c00e98be994e68debee317348` |
| Worker bytes | `2728c99756835b270362b7237747d354276b2ad2ec120fbdd934bcd0e8917b83` |
| Canonical native asset manifest | `4a5b0215aedcc7ca26d2b0e537a913fc2653bf0231b53b1dfb04e14bdfb79a3d` |
| Asset archive bytes | `54565f75696b5811c87d5d24a7bf1bbcb764398ffe0d2326493ef4b7f4f7ce29` |

Candidate identity binds the mathematical state and executed code/binary;
certificate history and native build/transport receipt hashes are in the sealed
run. Native assets contain no target, scalar answer, query stream or certificate.
Compiler policy and source archives stay reproducible beside the binary. Each
execution must use its own build platform. Linux staging is a separate local
build; these macOS bytes cannot establish Linux correctness or performance.

Nine focused Python controls pass. They check mathematical-only build inputs,
complete native/dependency retention, identity separation, one-use registration,
unexecuted transport rejection, retained negative/witness replay and tampering,
and unknown-cost timeout rows without retries. Native/source calls in the
mathematical execution tests are explicitly mocked; their returned admission
objects are test controls, not new execution receipts. A real fresh isolated
import control loads both new F5 and SAT adapters from the frozen source surface
and passes before/after module gates without executing either solver.

Run controls from the repository:

```sh
PYTHONDONTWRITEBYTECODE=1 python3.12 -m unittest discover \
  -s research/ic_candidate_tournament_20260915 -p 'test_prepared*.py' -v
RAYON_NUM_THREADS=1 cargo test --locked --offline --release --no-default-features \
  --example ic_tournament_worker prepared_target_tests -- --ignored --test-threads=1
```

Build/stage a **new** platform input rather than reusing the old worker:

```sh
PYTHONDONTWRITEBYTECODE=1 python3.12 research/ic_candidate_tournament_20260915/generic_build.py \
  --out NEW_BUILD_DIRECTORY
PYTHONDONTWRITEBYTECODE=1 python3.12 research/ic_candidate_tournament_20260915/prepared_f5_inputs_v1.py \
  --build NEW_BUILD_DIRECTORY --root /absolute/path/to/accepted/checkout --out NEW_ASSET_DIRECTORY
```

A new measured development invocation still requires its separately versioned
panel and externally recorded source/candidate/workload/run/whole-invocation
seals before `register`, `execute` and independent `audit`. None has executed
in this implementation follow-up. Natural ordinary-query evidence remains in
the accepted preparation histories. Fresh exposure/reference/hardware/observer/
calibration/resource/order reconciliation, fresh comparison, uncertainty,
operation totals, S/floor ratios, speedup and promotion remain pending/unknown.
