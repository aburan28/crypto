# Generalized GPU rho workers for the nine Certicom binary challenges

- **Status:** `open`
- **Kind:** implementation gap; measurement and hardware access gap
- **Affected work:** the nine Certicom binary challenge records in `src/binary_ecc/curve.rs`, the `ecdlp-nist` rho CLI, and the curve-specific `ecc2k130/` GPU campaign; see also [the matched GPU workload problem](ecbench-gpu-matched-workload-contract.md)
- **Next action owner:** Conductor task `T-3352` for the curve and CPU path; GPU worker task unassigned
- **Last updated:** 2026-10-16

## Problem

The repository cannot yet run a persistent, independently replayable GPU Pollard rho campaign on all nine Certicom binary challenges because the existing high-throughput client is specialized to ECC2K-130 and the generic binary rho solver retains its distinguished-point table in process memory.

The affected groups are ECC2K-130 and ECC2-131 over GF(2^131); ECC2K-163 and ECC2-163 over GF(2^163); ECC2-191 over GF(2^191); ECC2K-238 and ECC2-238 over GF(2^239); and ECC2K-358 and ECC2-353 over GF(2^359). Each has its own published field polynomial, coefficients, subgroup and target point. Random binary curves cannot inherit a Koblitz Frobenius quotient without a separately valid automorphism. The requested acceptance criterion is an end-to-end worker, durable distinguished-point ingest and collision verification for each exact challenge, with the fastest validated arithmetic selected by matched measurements on its intended hardware.

## Evidence and uncertainty

- **Observed:** `ecc2k130/src/main.cu` dispatches its configured degrees 131, 97 and smaller validation fields; it has no dispatch for all nine challenge records. `src/cryptanalysis/ecdlp_nist/rho.rs` uses an in-memory `HashMap` for distinguished points. The current development host is arm64 and has no `nvcc` or NVIDIA device utility. Certicom's parameter source is https://www.certicom.com/en/curves-list.
- **Derived or proved:** A walk rate alone omits durable DP transfer, merge, collision replay and verification. End-to-end throughput and a fastest-kernel choice therefore require measurements covering those phases. No mathematical obstruction is established.
- **Suspected explanation:** Unknown whether one field representation and kernel geometry will be best at all four degrees; this is a measurement question.
- **Missing capability:** GPU field arithmetic and exact walk/coefficients for every curve, durable per-curve identity and DP storage, independent replay, device safety validation and matched performance receipts.
- **What this evidence does not establish:** It gives no speed ranking for the nine curves and no conclusion about the feasibility of their discrete logarithms.

## What would unblock it

- [ ] Bind all nine exact curve and target records to independently checked generator/target subgroup tests and stable curve IDs.
- [ ] For each curve, pass CPU-reference equivalence of field operations, point operations, coefficient bookkeeping, distinguished-point records, restart and collision verification across ordinary and exceptional cases.
- [ ] Persist and resume DP collection with exact curve/workload identity, zero unexplained drops and independent collision replay.
- [ ] On each intended GPU class, compare candidate kernels with the current valid baseline using the same frozen workload, A/A controls, device safety checks, full phase accounting and uncertainty; retain the measured winner per curve and hardware profile.

These conditions would establish runnable, validated implementations and qualified performance choices at the measured devices and workloads. They would not establish an end-to-end attack-cost advantage over another algorithm without a separate complete cost and success-probability analysis.

## Next discriminating step

- **Action or proof obligation:** After the nine CPU parameter/subgroup checks pass, implement a bounded ECC2-131 random-curve GPU worker using the published target and a seeded planted control, then compare every emitted DP and checkpoint with an independent CPU replay.
- **Inputs and controls:** the exact ECC2-131 field/curve/subgroup record; an ordinary binary affine reference; a planted scalar with known answer; the published target as an answer-free collection input; matched kernel baseline and A/A timings.
- **Execution limits and stop condition:** freeze the seed, step count, DP predicate, hardware and maximum process memory before running; stop on any arithmetic mismatch, record drop, sanitizer fault or resource watchdog. These limits bound the validation run, not the broader campaign budget.
- **Outcomes and next actions:** passing replay permits measured tuning and then the remaining eight curve ports; a mismatch localizes a correctness defect; a bounded no-collision run says nothing about the target logarithm; missing hardware or a failed launch is infrastructure evidence and is rerun only after the dependency is restored.
- **Recheck trigger:** the first exact ECC2-131 GPU replay receipt or a generalized worker design with a verified backend.

## Resolution and history

- 2026-10-16: Opened while integrating the nine published parameter records. GPU and durable campaign support remain unimplemented for the full set.

Unresolved. Completing this record does not complete the requested integration.
