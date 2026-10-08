# GPU tau-adic arithmetic: compiler candidate, device execution pending

2026-09-29. Follow-up to [PR #913](https://github.com/aburan28/crypto/pull/913).
The CUDA reference backend compiles and agrees with independently verified
known-scalar fixtures. **No device execution or GPU speedup was measured.**
The readiness follow-up to [PR #960](https://github.com/aburan28/crypto/pull/960)
found and removed compiler-reported spilling in an opt-in candidate. This is
compiler evidence only; see the paired resource table below.

## Optimization candidates

| Variant | Recoding | Squaring | GPU result |
|---|---|---|---|
| Control | Binary NAF | Generic field multiplication | Not measured |
| Tau candidate | Reduced tau-NAF | Generic field multiplication | Not measured |
| Squaring candidate | Binary NAF | Precomputed linear basis map | Not measured |
| Combined candidate | Reduced tau-NAF | Precomputed linear basis map | Not measured |

Reduced tau-NAF previously reduced scalar arithmetic CPU cost, but its many
squarings, digit-dependent additions, and register use need a GPU comparison.
The new squaring candidate uses linearity over GF(2): XOR selected precomputed
images of polynomial-basis vectors. It removes the multiplier's shift/reduce
work but adds table reads. Whether that trade pays is unmeasured.

All four arms use the same portable field multiplier and Itoh-Tsujii inverse.
This is an isolated reference backend, **not** the existing optimized CLMAD,
normal-basis, batched walk backend. It does not modify a rho walk, initialization,
collector, production worker, or challenge computation. A known-scalar result
cannot be reported as a walk-speed improvement.

## What actually ran

| Check | Result | Evidence |
|---|---|---|
| Independent Python integer-polynomial oracle versus archived Sage output hashes | All 96 m=83/m=131 training and holdout fixtures agree | `check_host.py`, parent `results/run-02.json` |
| C++ compilation of the same arithmetic header; field squaring and inversion | 536 checks passed | [host-01.json](results/host-01.json) |
| Both scalar recodings and both squaring variants, including infinity/torsion/zero/negative-scalar edges | 504 point checks passed | [host-01.json](results/host-01.json) |
| NVRTC 12.8 to PTX, m=83/m=131 x both squaring variants x compute_89/compute_120 | All 8 compilations passed, empty compiler logs | [nvrtc-01.json](results/nvrtc-01.json) |
| Modal 1.6.0 app-definition import | Passed locally; no server call | `cloud_modal.py` |
| Ephemeral Modal RTX-PRO-6000 attempt | Failed: could not connect to Modal server | [cloud-attempt-01.json](results/cloud-attempt-01.json) |
| Authenticated RunPod read-only pod inventory | HTTP 403 Forbidden | [runpod-access-01.json](results/runpod-access-01.json) |

The initial Modal authentication probe timed out after 25 seconds without
configuring a profile. The subsequent bounded CLI attempt exited 1 before a
device result was obtained. RunPod access was tested using the documented
Bearer authorization header; no pod mutation was sent. These errors do not
establish whether the credentials are valid or whether the failure is a
network/access restriction. No provider workaround or permission bypass was
attempted. No cloud execution is confirmed.

The parent CPU receipts remain byte-for-byte unchanged. The archived successful
receipt is pinned to SHA-256
`f108dfa628f141b7da5c2885472aa960a25403a8f157e84a74cac2e58be8b797`.
Host and NVRTC receipts bind the exact CUDA arithmetic source hash. PTX hashes
are compilation evidence, not SASS instruction counts or timing measurements.

## What the runner measures

The fixed workload repeats each 24-case panel 1,024 times in a deterministic
shuffle: 24,576 evaluations per launch, **24 unique scalar/point fixtures**.
It records exact EC1 aliases and curve UIDs, device and host details, compiler
and register metadata, every paired sample, and verified output hashes.

Five A/A pairs precede seven alternating rounds. The control selects a common
launch-repetition count aiming for 50 ms samples, capped at 16; shorter samples
remain flagged. Every timed batch's final output is checked. Kernel-only
figures explicitly amortize recoding, allocation, and transfers.

Preparation, table construction, JIT/module reuse, allocation/upload, kernel,
download and verification are also reported. `accounted_stage_seconds` is a
sum of recorded arithmetic stages, with module/fixture reuse flags; it is not
a fresh-process cold-start measurement or an ECDLP pipeline measurement.
This first port does not measure maximum GPU occupancy or production throughput.
`rho_speedup` and `walk_iterations_per_second` remain null.

The runner now checks all 96 unique fixtures in all four arms first: **384
device outputs must verify before timing starts**. `--mode smoke` stops after
that correctness gate and produces no throughput summary. Device smoke also
warms the modules; their compile costs are recorded under `smoke`, and the
later timing panels explicitly record module reuse. Atomic checkpoint writes
preserve the last complete JSON on a process interruption. A submitted launch,
completed execution, and verified output are recorded separately. Missing
remote receipts mean unknown execution, not proof that nothing ran.

## Compiler-resource follow-up

The frozen extension is [READINESS_PROTOCOL.md](READINESS_PROTOCOL.md). One
candidate adds `#pragma unroll 1` to the linear squaring loop, guarded by
`LINEAR_SQUARE_NOUNROLL`. It changes no arithmetic. The paired device runner
can now select it with `--study square-unroll`; `original` remains the default.
Both recodings share a kernel, so these are kernel-level
resource counts, not measurements of either recoding's runtime.

| Target | m | Squaring | Registers: default / candidate | Candidate/default registers | Default spill store/load bytes | Candidate spill store/load bytes |
|---|---:|---|---:|---:|---:|---:|
| sm_89 | 83 | Generic multiply | 56 / 56 | 1.000 | 0 / 0 | 0 / 0 |
| sm_89 | 83 | Linear map | 255 / 48 | 0.188 | 1840 / 1852 | 0 / 0 |
| sm_89 | 131 | Generic multiply | 66 / 66 | 1.000 | 0 / 0 | 0 / 0 |
| sm_89 | 131 | Linear map | 96 / 56 | 0.583 | 0 / 0 | 0 / 0 |
| sm_120 | 83 | Generic multiply | 56 / 56 | 1.000 | 0 / 0 | 0 / 0 |
| sm_120 | 83 | Linear map | 255 / 48 | 0.188 | 5680 / 5908 | 0 / 0 |
| sm_120 | 131 | Generic multiply | 80 / 80 | 1.000 | 0 / 0 | 0 / 0 |
| sm_120 | 131 | Linear map | 96 / 70 | 0.729 | 0 / 0 | 0 / 0 |

All 16 matched NVRTC/PTXAS compilations passed using version 12.8.93 and the
same source hash. The candidate meets the preregistered **static** screen:
zero reported spill bytes at m=83 on both targets, none introduced at m=131.
Classification: engineering candidate. It has not established higher occupancy,
lower latency, or a runtime speedup. These spill figures are compiler reports,
not measured dynamic traffic; loop overhead and memory access still require
device profiling. Host regression checks again passed 536 field checks and
504 point checks ([host-02.json](results/host-02.json)). Host execution cannot
validate a CUDA-only pragma.

Evidence: [paired comparison](results/resource-comparison-01.json),
[control](results/ptxas-control-02.json),
[candidate stdout](results/ptxas-nounroll-01.stdout.json), and the initial
[compiler audit](results/ptxas-01.json). The candidate process exited 0 and
printed a complete eight-configuration receipt, but the surviving file was a
seven-configuration checkpoint with status `started`. Its cause is undiagnosed.
Both are retained: the complete stdout was copied byte-for-byte, and its first
seven entries exactly match the [partial file](results/ptxas-nounroll-01.json).
The comparison uses complete stdout; the partial file is not called a pass.

Five regression tests exercise atomic replacement failure, overwrite refusal,
missing/malformed receipts, smoke-only behavior, and abort-before-timing on a
mismatched result. They use a **CPU stub**, not a CUDA emulator or GPU.
Their output, Python compilation, Modal import, and parent receipt replay are
retained in [readiness-validation-01.json](results/readiness-validation-01.json).
The [local smoke attempt](results/local-smoke-attempt-01.json) failed before
launch because CuPy was absent; this runtime also has no NVIDIA device.
Provider HEAD probes returned Modal 503 and RunPod 403
([access recheck](results/access-recheck-01.json)). No new paid job was launched.

## Reproduce local validation

From this directory in an isolated Python 3.11 environment with g++:

```bash
python -m pip install numpy==2.2.6 nvidia-cuda-nvrtc-cu12==12.8.93
python check_host.py --output results/host-reproduction.json
python check_nvrtc.py --output results/nvrtc-reproduction.json
```

To reproduce the compiler-resource comparison (no GPU required):

```bash
python -m pip install nvidia-cuda-nvcc-cu12==12.8.93
python check_nvrtc.py --assemble --output results/control-reproduction.json
python check_nvrtc.py --assemble --square-no-unroll --output results/candidate-reproduction.json
python compare_resources.py --control results/control-reproduction.json --candidate results/candidate-reproduction.json --output results/comparison-reproduction.json
python -m unittest -v test_readiness
```

Output paths must be new. The NVRTC check needs no GPU; the host arithmetic
tests need no Sage installation because the oracle uses independent integer
polynomial arithmetic and checks archived Sage hashes.

## Run on Modal when server access works

Authenticate privately through Modal's normal setup or credential environment;
never place credentials in this repository. From this directory:

```bash
python -m pip install modal==1.6.0
modal run cloud_modal.py --mode smoke --output results/modal-smoke.json
```

After a verified smoke receipt, request the full original four-arm comparison
with `--mode benchmark` and a new output path. It repeats the smoke gate.

The no-unroll comparison now has its own frozen
[device protocol](UNROLL_DEVICE_PROTOCOL.md). To exercise the candidate and
its matched controls in the same invocation:

```bash
modal run cloud_modal.py --study square-unroll --mode smoke --output results/unroll-smoke.json
modal run cloud_modal.py --study square-unroll --mode benchmark --output results/unroll-benchmark.json
```

Run the second command only after verified smoke. It retains default and
no-unroll kernels for each recoding, separate compiled modules, five A/A pairs
for each control, and seven alternating rounds. Each candidate is compared
with its own recoding's control. Samples shorter than 50 ms are ineligible
for the screening result. The overall screen requires both degrees and both
training/holdout panels to pass. No new device measurements exist yet.

Modal and Docker use the shared `launch_gpu.py` entrypoint. Timed mode invokes
the current repository `tools/isolated_bench.py run`, reserves a complete
visible SMT sibling group while leaving another CPU free, and retains its
conditions record. A busy machine, unavailable isolation, contention, or user
threads left on reserved CPUs prevents an accepted timing result. The outer
receipt's `timing_eligible` field governs acceptance; inner kernel summaries
are diagnostics and must not override it. Smoke mode is correctness-only.
CPU isolation does not establish exclusive GPU use or isolate other VM tenants.
Dedicated GPU allocation and the recorded hardware/noise checks remain necessary.

Timeout rejection is unconditional: a worker that exits cleanly during the
termination grace period may retain a complete inner receipt, but the outer
receipt always sets `timing_eligible: false` when its deadline was exceeded.
The regression test first reproduced the contradictory timeout/eligible
receipt, then passed after the gate correction. All 16 local tests passed;
see [timeout-gate-validation-01.json](results/timeout-gate-validation-01.json).
This is an accounting correction with CPU-stub validation, not a new device
measurement or arithmetic optimization. Existing archived receipts are unchanged.

RunPod/local-container users can invoke the same study by appending
`--study square-unroll --mode smoke --output /results/unroll-smoke.json` to
the Docker command below, then use benchmark mode with a new output path.
Within an existing GPU environment, use
`python launch_gpu.py --study square-unroll --mode smoke --output results/unroll-smoke.json`.
The raw `run_gpu.py` is the worker; use the launcher for timed execution.

Validation of this wiring is retained in
[paired-runner-validation-02.json](results/paired-runner-validation-02.json).
Fifteen tests cover the cache, matched denominators/noise, complete alternating
loop, short samples, isolation policy, timeout/interruption cleanup, and prior receipt gates.
CuPy calls and event durations in these tests are explicit CPU stubs, with no
saved device-performance receipt. A real subprocess test confirms termination
of a parent/child group and cleanup on timeout. The
[local candidate smoke attempt](results/paired-smoke-local-01.json) stopped
before device execution because CuPy is absent; its failure is retained.

This uses one RTX-PRO-6000, at most one container, zero retries, a 600-second
function timeout and a 540-second benchmark-process timeout. It is an
ephemeral invocation, not a persistent deployment. Failure/partial receipts
are retained when the entrypoint or remote process starts far enough to write
them. A provider-level failure before the entrypoint starts needs the CLI log,
as in the captured attempt above.

The pinned CUDA 12.8.1 image installs CuPy 13.6.0 and NumPy 2.2.6. Image-build
and actual CuPy device integration remain unvalidated because server access
failed. Successful PTX compilation alone does not discharge that gate.

## RunPod portability

The same `run_gpu.py` can run inside an already authorized GPU pod with
CUDA 12.8 and `cupy-cuda12x==13.6.0`, `numpy==2.2.6`. A one-shot Dockerfile is
provided; build from the repository root, not this subdirectory:

```bash
docker build -f research/notes/ecc2k130/tau_adic_arithmetic_20260929/gpu/Dockerfile -t tau-arithmetic .
mkdir -p gpu-results
docker run --rm --gpus all -v "$PWD/gpu-results:/results" tau-arithmetic
```

The process is bounded by 540 seconds. **Container exit does not necessarily
stop RunPod pod billing**; the pod lifecycle is separate. The image was not
built or published in this environment, and no RunPod worker was launched.
The container defaults to smoke mode and `/results/gpu-smoke.json`. For the
full original timing protocol, append
`--mode benchmark --output /results/gpu-benchmark.json` to the `docker run`
command. The timeout sends TERM at 540 seconds, then KILL after a 10-second
grace period if necessary.

## Decision

- [x] Freeze the protocol before execution (local commit `afe0630`).
- [x] Implement the four-arm arithmetic comparison and independent oracle.
- [x] Validate host arithmetic and compile Ada/Blackwell CUDA variants.
- [x] Attempt cloud access and preserve the unsuccessful outcomes.
- [x] Compare default/candidate compiler resources on Ada and Blackwell.
- [x] Add the device smoke gate and test receipt/failure handling locally.
- [ ] Build the remote image and complete device correctness checks.
- [ ] Collect paired GPU timing/metadata receipts and apply the screening rule.
- [x] Freeze and wire the matched no-unroll device comparison and isolated launcher.
- [ ] Execute that comparison on a GPU and retain uncontended device results.
- [ ] Establish an eligible runtime fraction or any end-to-end rho benefit.

The next concrete requirement is a working authorized Modal connection or
accessible RunPod GPU environment. Runtime tuning and selecting the no-unroll
candidate require the first verified device comparison.
The initial compiler protocol and its single-candidate extension were frozen
locally in commits `50ef214` and `4ad9010`, respectively, before those runs.

Provider references checked 2026-09-29:
[Modal GPU support](https://modal.com/docs/guide/gpu),
[Modal CUDA images](https://modal.com/docs/guide/cuda),
[Modal pricing](https://modal.com/pricing),
[RunPod GraphQL specification](https://graphql-spec.runpod.io/).
Compiler statistics: [NVIDIA CUDA compiler documentation](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/).
See [PROTOCOL.md](PROTOCOL.md) for the preregistered accounting and stop rules.
