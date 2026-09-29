# GPU tau-adic arithmetic: validated package, cloud execution blocked

2026-09-29. Follow-up to [PR #913](https://github.com/aburan28/crypto/pull/913).
The CUDA reference backend compiles and agrees with independently verified
known-scalar fixtures. **No device execution or GPU speedup was measured.**

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

## Reproduce local validation

From this directory in an isolated Python 3.11 environment with g++:

```bash
python -m pip install numpy==2.2.6 nvidia-cuda-nvrtc-cu12==12.8.93
python check_host.py --output results/host-reproduction.json
python check_nvrtc.py --output results/nvrtc-reproduction.json
```

Output paths must be new. The NVRTC check needs no GPU; the host arithmetic
tests need no Sage installation because the oracle uses independent integer
polynomial arithmetic and checks archived Sage hashes.

## Run on Modal when server access works

Authenticate privately through Modal's normal setup or credential environment;
never place credentials in this repository. From this directory:

```bash
python -m pip install modal==1.6.0
modal run cloud_modal.py --output results/modal-reproduction.json
```

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

## Decision

- [x] Freeze the protocol before execution (local commit `afe0630`).
- [x] Implement the four-arm arithmetic comparison and independent oracle.
- [x] Validate host arithmetic and compile Ada/Blackwell CUDA variants.
- [x] Attempt cloud access and preserve the unsuccessful outcomes.
- [ ] Build the remote image and complete device correctness checks.
- [ ] Collect paired GPU timing/metadata receipts and apply the screening rule.
- [ ] Establish an eligible runtime fraction or any end-to-end rho benefit.

The next concrete requirement is a working authorized Modal connection or
accessible RunPod GPU environment. Additional kernel optimization should wait
for the first verified device result.

Provider references checked 2026-09-29:
[Modal GPU support](https://modal.com/docs/guide/gpu),
[Modal CUDA images](https://modal.com/docs/guide/cuda),
[Modal pricing](https://modal.com/pricing),
[RunPod GraphQL specification](https://graphql-spec.runpod.io/).
See [PROTOCOL.md](PROTOCOL.md) for the preregistered accounting and stop rules.
