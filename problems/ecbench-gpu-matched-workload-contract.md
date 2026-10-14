# EC Bench GPU matched-workload contract and device validation

- **Status:** `open`
- **Kind:** implementation gap; access/tool limit
- **Affected work:** [EC Bench GPU stage runner](../docs/ecbench/README.md#13-gpu-stage-runner), `ecbench compare` CPU/GPU wall or operation-cost claims
- **Next action owner:** EC Bench maintainer on a CUDA host
- **Last updated:** 2026-10-13

## Problem

EC Bench cannot yet make an admitted CPU/GPU comparison of an ECDLP method: the existing CUDA driver accepts its own fixed curve, target and options, while `ecbench run` requires a frozen workload and seed, an exact cost ledger, an independently replayed answer, and paired isolation evidence. The `gpu-stage` command only runs the device self-test and a bounded stage measurement. This blocks a CPU/GPU confidence interval for the same target and method. The CPU isolation score does not grade GPU contention.

## Evidence and uncertainty

- **Observed:** `src/cryptanalysis/ecbench/gpu.rs` snapshots `gpu/ecc/bench`, runs `selftest`, and saves hashed stdout, stderr, arguments, status and timeout. `gpu/ecc/bench.cu` has fixed curve/target stage modes. `docs/ecbench/README.md` describes the stage receipt and its comparison boundary. The local Apple M4 Pro host has no NVIDIA CUDA device or `nvcc`, so device execution has not been validated here. The GPU wrapper's fake-driver unit test checks control flow and hashes, not CUDA arithmetic or throughput.
- **Derived or proved:** none.
- **Suspected explanation:** no claim about device correctness or performance follows from a local missing device.
- **Missing capability:** a CUDA input/output adapter for EC Bench's exact curve, subgroup, generator, challenge point and seed; an operation and phase ledger; independently checkable recovered answer; GPU device identity and contention telemetry; an accessible CUDA host for validation.
- **What this evidence does not establish:** any CPU/GPU speed ratio, GPU benchmark fidelity level, device failure, or impossibility of implementing the adapter.

## What would unblock it

- [ ] On one public known-answer EC Bench workload, the CUDA method accepts the exact frozen workload and seed, reports the recovered scalar and exact cold-start and online phase counts, and the existing independent `ecbench verify` path checks the answer and ledger. A failed, timed out or wrong-answer run remains a sealed record.
- [ ] On a CUDA host, run the device self-test and interleaved CPU/GPU repeats of that workload with recorded device identity, CPU placement, GPU activity/clock/memory pressure, background GPU activity and A/A controls. Define and enforce GPU admission thresholds before interpreting a wall-time interval.
- [ ] Produce an `ecbench compare` result with the same workload identity, matched seeds, complete cost accounting and a paired interval; keep any unmeasured transfers or setup outside an admitted end-to-end result.

These conditions would establish a comparable result for the tested curve, target and host. They would not establish other curves, scaling, or performance on other GPUs.

## Next discriminating step

- **Action or proof obligation:** implement one workload adapter for the existing CUDA BSGS path, then run its device self-test and one known-answer fixture on a CUDA host.
- **Inputs and controls:** one frozen EC Bench prime-field workload supported by the CUDA driver, identical CPU reference workload and seed, device self-test, independent scalar check, and A/A repeats. If the driver cannot represent that workload, record the exact representation mismatch before changing the fixture.
- **Execution limits and stop condition:** one fixture and bounded repeats under a frozen protocol; stop on wrong answer, self-test failure, missing cost fields, missing contention telemetry, or watchdog expiry. No broad GPU campaign is implied.
- **Outcomes and next actions:** a verified fixture permits a predeclared paired measurement; a bounded mismatch directs adapter correction; incomplete telemetry keeps wall results descriptive; infrastructure failure calls for a CUDA host retry without a performance conclusion.
- **Recheck trigger:** a CUDA workload adapter and host are available.

A planned step is not a completed experiment. Apply the existing execution, review, budget, and archival rules before acting on it.

## Resolution and history

- 2026-10-13: Opened after the first bounded GPU stage runner was added. Unresolved.
