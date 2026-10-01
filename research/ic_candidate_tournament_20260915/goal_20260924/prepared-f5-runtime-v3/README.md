# Corrected prepared F5 runtime gate

PR #1127 repaired the actual prepared serializer and accepted native release
interface controls. The corrected worker needs its own source-admission version:
the immutable input v1 pin and consumed runtime v2 remain historical.
[PROTOCOL.md](PROTOCOL.md) was committed at `43b76b342` before this build and
staging. This change adds input v2, runtime v3 and frozen audit transport v2.
No native scientific registration or solver execution is made by this gate.

The new controlled macOS ARM64 release build and
[retained assets](native-inputs-macos-arm64/manifest.json) pass admission and
relocation replay. The archive retains the worker, complete root/dependency
source archives, source manifest and controlled build receipt/policy/log:

| Retained input | Recorded identity or count |
| --- | --- |
| Corrected worker source | `8683bfd2ffebbada30397e5ced9825ee625f85c6838ca56ada6665f6a43c90ca` |
| Root source files | 533 |
| Registry packages / source files | 67 / 3,011 |
| Complete native asset archive | 12,278,251 bytes |
| Asset archive SHA-256 | `f28711697074edd65c6b3d37010d7e1986525348836d7b9b291fb5ec999e4400` |
| Rust source manifest SHA-256 | `08c759415b2346d0d41c49ce2e958bf0c147a32796d0ac44f907ea7118095b4c` |
| Worker executable SHA-256 | `5d113a46f67c982bfdb6b81dc49aec7203ca990bd5b7574b4c2feec4c5011624` |

The [relocation/admission receipt](controls/native-asset-replay.json) checks every
asset and both nested source archives against their manifest, proves byte equality
after copying to a temporary location, and records rejection by the unchanged
old input gate. Staging launches no worker. The controlled builder's native
`--build-identity` invocation reads metadata and performs no mathematical query.
The new executable is not substituted into any historical result or runtime.

`prepared_f5_inputs_v2.py` changes only the accepted worker-source pin and its
version description relative to input v1. Runtime v3 uses that input and its own
registered entrypoint. Registration rejects missing frozen-helper coverage before
deriving candidate identity or consuming an execution claim. Its audit CLI
delegates to transport v2, which requires
its helper to have been frozen before execution, runs the registered auditor
under `-I -S -B`, verifies the interpreter and loaded sources before/after, and
checks the complete/incomplete claim boundary. Transport v1 remains unchanged
and cannot admit the new entrypoint. There is no mutation or retrofit of a
historical registration. The complete bounded Python source surface is retained
at registration; interpreter executable and stdlib contents are bound by hashes,
not packaged as a portable Python installation.

[local-tests.log](controls/local-tests.log) retains all 36 focused old/new runtime
and transport test results. New controls cover real retained-source/build/asset
admission, source/binary/dependency/recipe substitution, cross-version rejection,
candidate versus workload identity, target and resource policy, one-use
registration and rejection of an unexecuted audit, real isolated import-only
source gates, and mocked native calls around independent mathematical replay.
The transport controls use a clearly synthetic auditor to exercise real isolated
children, retained-source selection despite changed live files, external import
rejection, artifact/claim tampering, timeout and mandatory preexecution helper
coverage. They execute no native solver or fresh mathematical target. Require
the linked PR's exact-head applicable checks before treating this implementation
as accepted; local controls are not full native scientific execution evidence.

The original implementation is accepted in [PR #1138](https://github.com/aburan28/crypto/pull/1138)
at merge `6dc9a7c661f28ae4bc38de8be6e4920a56c7a56b`; its four applicable
workflows passed on head `2ffbf6cbf518be91feb29fa0a62db9b6ce1b0500`, including
459 harness and 16 boundary tests. The subsequent early registration-helper
guard has [37 focused local controls](controls/registration-guard-tests.log),
including a specific rejection before native-job or candidate construction.
Its acceptance still requires its own exact-head CI; the original native assets,
source pins, archived tests and scientific-execution status are unchanged.

For a future disclosed control, first publish a new frozen protocol, derive its
candidate/workload/run from the accepted full runtime and these assets, and retain
the external execution seal. Execute that registration once; any failure consumes
it. Replay with its own retained `prepared_runtime_transport_v2.py` helper. The
old F5/SAT prepared registrations and all three confirmation sets stay closed.
Neither this implementation nor its tests authorize rerunning them.

The next scientific gate is a new complete source-bound native prepared control.
SAT budget/encoding diagnostics on its retained source-valid query remain
separate. Natural ordinary-query yield and failed-attempt evidence are the
original independently audited preparation histories. Current exposure union,
reference bindings, physical comparison-host build, calibration/resources/order
and fresh paired one-target qualification remain pending. Scientific online cost,
matched rho ratio, promotion and full goal completion remain unestablished.
