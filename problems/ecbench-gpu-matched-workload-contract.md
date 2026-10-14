# EC Bench GPU matched-workload contract and device validation

- **Status:** `open`
- **Kind:** implementation gap; device telemetry limit
- **Affected work:** [EC Bench GPU driver](../gpu/ecc/README.md#ec-bench-matched-workload), `ecbench compare` CPU/GPU wall or operation-cost claims
- **Next action owner:** EC Bench maintainer for GPU admission and paired runs
- **Last updated:** 2026-10-13

## Problem

EC Bench cannot yet make an admitted CPU/GPU comparison of an ECDLP method. The `gpu-match` path has solved one frozen prime-curve target on an RTX 4090 and independently checked the answer, but that device's `nvidia-smi` does not expose the PCIe traffic counters required by the GPU wall-time admission rule. The verified receipt therefore refuses wall time. CPU and GPU executions are not yet interleaved in one session or priced in a common complete-cost unit. The CPU isolation score does not grade GPU contention. These gaps block a calibrated CPU/GPU confidence interval for the same target and method.

## Evidence and uncertainty

- **Observed:** `src/cryptanalysis/ecbench/gpu.rs` snapshots `gpu/ecc/bench`, saves raw streams and hashes, rebuilds the workload, checks the compiled curve and replays `[k]G = Q`. `gpu/ecc/bench.cu` accepts an explicit target for full-subgroup BSGS and emits phase, host and explicit `cudaMemcpy` transfer counters. `ecbench gpu-match-verify` rechecks a receipt without a device. `gpu/ecc/ecbench_spec.py` emits deterministic 36- and 40-bit prime-curve specs. On 2026-10-13 PDT, [the RunPod RTX 4090 receipt](../gpu/ecc/validation/20261013-runpod-4090/README.md) for `W673c598bf1d4` recorded a passing device self-test, recovered scalar `34790300172`, independent scalar replay, and an offline audit with `issues: []`. The wall-time admission was refused because `pcie.rx_util` is unsupported by that driver. The full binary snapshot and raw streams are in the linked archive; failed setup and earlier wrapper attempts remain in the local fetched directory named there.
- **Access/tool limit:** RunPod access was available for this device run, and its pod was deleted after receipt fetch. The tested driver's PCIe traffic query is unavailable; no alternative counter has yet been qualified against a contention control.
- **Derived or proved:** none.
- **Suspected explanation:** the unavailable PCIe query reflects the tested driver interface, not an observed competing workload. No general device-timing conclusion follows from it.
- **Missing capability:** qualified GPU traffic telemetry or a separately validated admission rule for this driver; an interleaved CPU/GPU session and comparable complete-cost model; empirical A/A interval coverage and contention-gate calibration on an isolated Linux host; a wider-than-one-word EC Bench group representation for broader curves. The transfer byte counters cover explicit `cudaMemcpy` calls, not all implicit runtime traffic.
- **What this evidence does not establish:** any CPU/GPU speed ratio, admitted GPU wall-time interval, calibrated GPU fidelity level, performance on other devices, or wider-curve result.

## What would unblock it

- [x] On one public known-answer EC Bench workload, compile and run the CUDA method with the exact frozen target, report the recovered scalar and phase counts, and have the independent GPU receipt audit check the answer and ledger. The completed scope is `W673c598bf1d4` on one RTX 4090; see the linked receipt archive. Earlier failed setup and wrapper attempts are retained separately. The current BSGS method is deterministic and has no algorithm seed; the target seed is frozen in the workload.
- [ ] On a CUDA host, run the device self-test and interleaved CPU/GPU repeats of that workload with recorded device identity, CPU placement, GPU activity/clock/memory pressure, background GPU activity and A/A controls. Define and enforce GPU admission thresholds before interpreting a wall-time interval.
- [ ] Produce an interleaved same-host CPU/GPU comparison with the same workload identity, complete cost accounting, A/A controls and a paired interval. Keep any unmeasured transfers or setup outside an admitted end-to-end result.
- [ ] Demonstrate wider curve coverage with an independently replayed answer and versioned raw receipts; do not transfer the small-curve result to a 256-bit curve by assumption.

These conditions would establish a comparable result for the tested curve, target and host. They would not establish other curves, scaling, or performance on other GPUs.

## Next discriminating step

- **Action or proof obligation:** qualify a GPU traffic measurement on this RTX 4090 driver, or select a device exposing the required counter; then run a predeclared, interleaved CPU/GPU comparison of the exact same workload with A/A controls and full phase accounting.
- **Inputs and controls:** the frozen `W673c598bf1d4` workload and verified GPU receipt, an independently verified CPU record for that workload, idle and deliberately contended GPU controls, CPU placement and contention samples, and repeated A/A pairs. Check any substitute traffic measurement against injected GPU transfers before using it for admission.
- **Execution limits and stop condition:** bounded interleaved repeats under a frozen protocol; refuse wall admission on wrong answer, self-test failure, missing cost fields, missing or unqualified contention telemetry, or watchdog expiry. No broad GPU campaign is implied.
- **Outcomes and next actions:** a validated counter and successful controls permit paired wall-time analysis; continued telemetry absence keeps wall results refused; a correctness mismatch directs adapter correction; infrastructure failure calls for a host retry without a performance conclusion.
- **Recheck trigger:** supported GPU traffic telemetry or a qualified substitute and a paired runner are available.

A planned step is not a completed experiment. Apply the existing execution, review, budget, and archival rules before acting on it.

## Resolution and history

- 2026-10-13: Opened after the first bounded GPU stage runner was added. Unresolved.
- 2026-10-13: Added the matched-target adapter and separate GPU receipt/audit implementation; device validation and complete paired comparison remain open. RunPod API access was unavailable to this session.
- 2026-10-13: RunPod RTX 4090 device self-test and frozen `W673c598bf1d4` solve passed; independent receipt audit passed. The device wall-time gate refused because the driver lacks the queried PCIe traffic counter. Source-coordinate parsing and telemetry failure classification were corrected after the first device attempts; paired CPU/GPU comparison remains open.
