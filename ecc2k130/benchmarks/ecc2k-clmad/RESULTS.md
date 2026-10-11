# gpu/ecc2k native carry-less product on the RTX PRO 6000

One NVIDIA RTX PRO 6000 Blackwell Server Edition (driver 580.95.05, 188 SMs,
2.43 GHz max), CUDA 13.3 (`compiler.38244171_0`), `-arch=sm_120`, through
`ecc2k130/modal_job.py` with [`gpujob.sh`](gpujob.sh) (sources staged by
[`prepare.sh`](prepare.sh) from `gpu/ecc2k` at the revision in
`results/host.txt`). Two binaries from the same source: the default device
build, which multiplies with PTX `clmad.lo/hi.u64` (`F2M_CLMAD=1`), and the
`-DF2M_CLMAD=0` control, which uses the software product the host harness
runs. `cuobjdump -sass` counts 894 `CLMAD` instructions in the native binary
and 0 in the control.

**Correctness.** `./bench selftest` passed for both binaries
(`results/selftest-*.txt`): every device kernel against the host
implementation, no `FAIL` line.

**Rates**, three alternating repetitions each (`./bench field`,
`./bench rho --iters 2048`, 385,024 class walks, lowmem variant):

| build | mul throughput (Gop/s) | mul dependent (Gop/s) | inv (Mop/s) | rho (Gstep/s) |
|---|---:|---:|---:|---:|
| software, rep 1/2/3 | 97.29 / 97.31 / 97.33 | 97.22 / 96.87 / 97.22 | 2764.9 / 2774.6 / 2766.1 | 1.415 / 1.415 / 1.436 |
| native, rep 1/2/3 | 127.26 / 127.30 / 127.29 | 127.09 / 127.19 / 127.20 | 3430.9 / 3424.7 / 3434.0 | 1.684 / 1.685 / 1.685 |
| ratio (medians) | 1.308 | 1.308 | 1.238 | 1.190 |

Observed: the native product raises the F(2^97) multiply rate by 1.31x, the
Itoh–Tsujii inversion by 1.24x and the ECC2K-95 class walk by 1.19x on this
card. Squaring (spread-based, unchanged) and the class-weight read are flat.
This is the first GPU execution of any kernel in `gpu/ecc2k`; the
multiplication is no longer the walk's whole cost, so the remaining items in
`gpu/ecc2k/README.md` (per-thread batch width, in-walk reseeding, Frobenius
divergence) are now the larger share. `ECC2K-95` itself is a retired
challenge; no claim about ECC2K-130 follows from these rates.

Raw outputs: `results/` (host line, per-run logs, `ptxas -v` for both builds,
source hashes, Modal launch receipt).
