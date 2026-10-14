# EC Bench CUDA matched-workload validation, 2026-10-13 PDT

This is a device-correctness receipt for one frozen public known-answer workload. It makes no CPU/GPU speed or wall-time fidelity claim.

| Item | Observation |
| --- | --- |
| Source | PR branch based on `49fb057e8293936189a5e56276acbfe4f71ef5ba`; final `src/cryptanalysis/ecbench/gpu.rs` SHA-256 `91e5c23dfc4648d160ad54aa237ca3c395676d80eddec43059a0db3a9ef4fc12` |
| Device | NVIDIA GeForce RTX 4090, CUDA 12.8, `sm_89`; RunPod pod `0hp47p30opj3tt`, deleted after receipt fetch |
| Workload | `W673c598bf1d4`, target `(0x151248995, 0x5aa04ecdf)`, subgroup order `39615274243` |
| CUDA binary | SHA-256 `2f1a6c24afbc03079770adfd6d2a1886a006517466eb519e4b4bd47aa7050c08` |
| Device self-test | Exit 0; host/device field, point, rho, and BSGS checks passed |
| Matched solve | `status: verified`; recovered scalar `34790300172`; independent `[k]G = Q` check passed |
| Offline receipt audit | `valid: true`, `issues: []` |
| GPU wall-time admission | `refused`: this driver's `nvidia-smi` does not expose `pcie.rx_util`/`pcie.tx_util` |

The receipt archive includes the compiled CUDA binary snapshot, exact workload, raw command streams, device samples, receipt JSON, audit JSON, build logs, and environment hashes. Its SHA-256 is `a80d6abb6c1a059966b7f15827198e7cbb2c35beae27ad2588d9419e44a206f5`.

To replay the audit from this branch on a host with a built `ecbench` binary:

```sh
mkdir -p /tmp/ecbench-runpod-replay
tar -xzf gpu/ecc/validation/20261013-runpod-4090/receipt.tar.gz -C /tmp/ecbench-runpod-replay
ecbench gpu-match-verify --dir /tmp/ecbench-runpod-replay/gpu-match-v4
```

The full fetched directory, including earlier failed setup and wrapper attempts, is retained locally at `/private/tmp/ecbench-runpod-1809/ecbench-1809`. The first `gpu-match` attempt rejected the workload's normal `0x` coordinates; the parser fix is part of this PR. An earlier telemetry receipt incorrectly cascaded an unsupported PCIe counter into identity and clock blockers; the final receipt reports the missing counter directly. A RunPod RTX 5090 request returned no capacity before the RTX 4090 allocation. The first Rust build lacked system OpenSSL headers and `pkg-config`; its log was retained before the successful retry.

One verified device solve does not supply an admitted timing interval or a same-host CPU/GPU comparison. Interleaved CPU/GPU trials, A/A calibration, complete cost accounting, and wider-curve coverage remain open in [the problem record](../../../../problems/ecbench-gpu-matched-workload-contract.md).
