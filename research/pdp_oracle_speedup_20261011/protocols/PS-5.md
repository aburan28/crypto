# PS-5: GPU kernel

Part of [`../PLAN.md`](../PLAN.md); frozen when the plan merges.  Repo:
cryptanalysis (C/CUDA core, cloud runners, isolated-benchmark service).
Needs PS-3 and an owner-approved allocation.  The companion protocol is
`experiments/pdp-oracle-speedup-20261011/README.md` in that repository.

**Hypothesis.**  Executing PS-3's straight-line program on a CUDA device,
with logic operations mapped to `lop3`, one pair per lane and the Gray
position per warp, reaches at least 50× PS-3's single-core pairs per second.
The expected ratio is the device's 32-bit logic throughput over one core's
512-bit logic throughput, to be filled from the device datasheet and a
microbenchmark at the start of the task; the plan does not guess it.

**Frozen inputs.**  PS-3's program (hash), PS-0's corpus and digests, `ℓ`
from 11 to 15 (and 16 if time allows), and the device's `gpu` block.

**Gate before any timing.**  The whole PS-0 corpus is replayed on the device
with identical digests.

**Unit.**  Pairs per second per device and per core-equivalent, counted bit
operations per pair, A/A on the device, ABAB against PS-3 on the same host.

**Success.**  ≥ 50× PS-3 on one core, digests identical, run records carrying
the hardware block.  **Stop.**  Below 20×.

**Spend.**  Requires the owner's explicit approval and a cap (planning cap:
8 GPU-hours).  Every pod or instance started is deleted; a run whose
`gpu.devices` is empty is labelled a CPU-only measurement.  Nothing is
launched by merging this plan.

**Layout.**  A device kernel with a CPU emulation of it (the pattern of
`suite/cuda/f4_gf2_kernel.cu` and `f4_gf2_emulate.c` in `cryptanalysis`), so
digest parity is checked on a CPU-only host before any device is rented.

**Paired change.**  `fes_gpu.rs` is byte-identical in both repositories; any
change lands in `cryptanalysis` first and is mirrored in `crypto` in a paired
PR.

**Class.**  Engineering.  **Effort.**  7–10 days plus GPU hours.
