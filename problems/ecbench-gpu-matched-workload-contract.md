# EC Bench GPU matched-workload contract and device validation

- **Status:** `open`
- **Kind:** implementation gap; access/tool limit
- **Affected work:** [EC Bench GPU driver](../gpu/ecc/README.md#ec-bench-matched-workload), `ecbench compare` CPU/GPU wall or operation-cost claims
- **Next action owner:** EC Bench maintainer on a CUDA host
- **Last updated:** 2026-10-13

## Problem

EC Bench cannot yet make an admitted CPU/GPU comparison of an ECDLP method. The new `gpu-match` path passes a frozen prime-curve target to CUDA BSGS and independently checks the answer, but it has not run on a CUDA device. It writes a separate receipt: CPU and GPU executions are not yet interleaved in one session or priced in a common complete-cost unit. The CPU isolation score does not grade GPU contention; the GPU receipt has its own conservative timing gate. This blocks a calibrated CPU/GPU confidence interval for the same target and method.

## Evidence and uncertainty

- **Observed:** `src/cryptanalysis/ecbench/gpu.rs` snapshots `gpu/ecc/bench`, saves raw streams and hashes, rebuilds the workload, checks the compiled curve and replays `[k]G = Q`. `gpu/ecc/bench.cu` now accepts an explicit target for full-subgroup BSGS and emits phase, host and transfer counters. `ecbench gpu-match-verify` rechecks a receipt without a device. `gpu/ecc/ecbench_spec.py` emits deterministic 36- and 40-bit prime-curve specs. Local tests use a simulated executable; they establish wrapper behavior only. The Apple M4 Pro host has neither an NVIDIA device nor `nvcc`.
- **Access limit:** the requested RunPod launcher was inspected read-only. `python3 cloud/runpod_pod.py list` stopped with `RUNPOD_API_KEY is not set`; no pod was created and no device test ran.
- **Derived or proved:** none.
- **Suspected explanation:** no claim about device correctness or performance follows from a local missing device.
- **Missing capability:** CUDA compilation and device validation of the adapter; an interleaved CPU/GPU session and comparable complete-cost model; empirical A/A interval coverage and contention-gate calibration on an isolated Linux host; a wider-than-one-word EC Bench group representation for broader curves. RunPod credentials or an authorized existing runner connection are needed for the device step.
- **What this evidence does not establish:** any CPU/GPU speed ratio, GPU benchmark fidelity level, device failure, or impossibility of implementing the adapter.

## What would unblock it

- [ ] On one public known-answer EC Bench workload, compile and run the CUDA method with the exact frozen target, report the recovered scalar and exact cold-start and online phase counts, and have the independent GPU receipt audit check the answer and ledger. Include failed, timed out and wrong-answer receipts. The current BSGS method is deterministic and has no algorithm seed; the target seed is frozen in the workload.
- [ ] On a CUDA host, run the device self-test and interleaved CPU/GPU repeats of that workload with recorded device identity, CPU placement, GPU activity/clock/memory pressure, background GPU activity and A/A controls. Define and enforce GPU admission thresholds before interpreting a wall-time interval.
- [ ] Produce an interleaved same-host CPU/GPU comparison with the same workload identity, complete cost accounting, A/A controls and a paired interval. Keep any unmeasured transfers or setup outside an admitted end-to-end result.
- [ ] Demonstrate wider curve coverage with an independently replayed answer and versioned raw receipts; do not transfer the small-curve result to a 256-bit curve by assumption.

These conditions would establish a comparable result for the tested curve, target and host. They would not establish other curves, scaling, or performance on other GPUs.

## Next discriminating step

- **Action or proof obligation:** supply RunPod access to the execution environment, compile the adapter, run its device self-test and one frozen known-answer fixture, then audit the fetched receipt independently.
- **Inputs and controls:** one frozen EC Bench prime-field workload supported by the CUDA driver, identical CPU reference workload and seed, device self-test, independent scalar check, and A/A repeats. If the driver cannot represent that workload, record the exact representation mismatch before changing the fixture.
- **Execution limits and stop condition:** one fixture and bounded repeats under a frozen protocol; stop on wrong answer, self-test failure, missing cost fields, missing contention telemetry, or watchdog expiry. No broad GPU campaign is implied.
- **Outcomes and next actions:** a verified fixture permits a predeclared paired measurement; a bounded mismatch directs adapter correction; incomplete telemetry keeps wall results descriptive; infrastructure failure calls for a CUDA host retry without a performance conclusion.
- **Recheck trigger:** RunPod access becomes available to this execution environment, or an authorized CUDA runner is identified.

A planned step is not a completed experiment. Apply the existing execution, review, budget, and archival rules before acting on it.

## Resolution and history

- 2026-10-13: Opened after the first bounded GPU stage runner was added. Unresolved.
- 2026-10-13: Added the matched-target adapter and separate GPU receipt/audit implementation; device validation and complete paired comparison remain open. RunPod API access was unavailable to this session.
