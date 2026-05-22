# SHA-1 chosen-prefix collision: FPGA vs GPU cost model

Companion to `hdl/sha1/`. Numbers are estimates; the dominant uncertainty is the **effective** (not raw) compression rate, which depends on how well the search kernel is implemented. See "Sensitivity" at the bottom.

## Workload

Leurent & Peyrin's 2020 chosen-prefix collision on SHA-1 ("SHA-1 is a Shambles") requires roughly:

| Phase | Work | Cost driver |
|---|---|---|
| Birthday / near-collision search | 2^61.2 compressions | dominant, parallel, GPU/FPGA-friendly |
| Path conditioning (4-block chain) | ~2^33 compressions | runs on CPU once per chain |
| **Total** | **~2^63.4 ≈ 1.3e19 compressions** | |

This document assumes the attacker is paying for the birthday phase.

## Hardware reference points

### AWS f2.6xlarge (FPGA)
- 1× AMD Virtex UltraScale+ **VU47P** (2.85M LUTs, 5.5M FFs, 9024 DSP, 84 Mb BRAM)
- 24 vCPU, 256 GiB RAM, PCIe Gen4 x16 to FPGA
- **$2.65/hr** on-demand (us-east-1, mid-2026), **$1.40/hr** 1-yr reserved

### AWS p5.48xlarge (GPU)
- 8× NVIDIA H100 80 GB
- 192 vCPU, 2 TiB RAM
- **$98.32/hr** on-demand, **~$55/hr** 1-yr reserved

## Per-instance effective throughput

"Effective" = compressions actually advancing the search per second (raw hashcat numbers overstate by 2–3× because of conditions, masking, and message-modification overhead).

| Instance | Raw SHA-1 ceiling | Search-effective | Notes |
|---|---|---|---|
| f2.6xlarge (1× VU47P) | ~100 GH/s (150 × 0.5 GH/s pipelines @ 500 MHz) | **~50 GH/s** | from `hdl/sha1/sha1_pipe.vhd` × 150 with the structures in `sha1_search.vhd` |
| H100 (single GPU) | ~25 GH/s (hashcat) | **~10 GH/s** | Shambles authors hit ~35% of peak; H100 better at staying high-occupancy → ~40% |
| p5.48xlarge (8× H100) | ~200 GH/s | **~60 GH/s** | minus host-side coordination overhead |

## Cost per effective GH/s-hour

| Instance | $/hr (on-demand) | Effective GH/s | **$/GH/s-hour** |
|---|---|---|---|
| f2.6xlarge | $2.65 | 50 | **$0.053** |
| p5.48xlarge | $98.32 | 60 | **$1.64** |

**FPGA is ~31× cheaper per unit of search work**, *assuming* the HDL stack from `hdl/sha1/` is implemented to the level sketched there.

With 1-yr reserved pricing the ratio narrows slightly: f2 $0.028 vs p5 $0.92, still ~33×.

## Wall-clock vs cluster cost

Total work: **1.3e19** compressions.

### Single instance
| Instance | Wall-clock |
|---|---|
| 1× f2.6xlarge | 8.25 years |
| 1× p5.48xlarge | 6.88 years |

### Same wall-clock, increasing cluster size
| Wall-clock | f2 count | f2 cost (on-demand) | p5 count | p5 cost (on-demand) | FPGA cheaper by |
|---|---|---|---|---|---|
| 90 days | 34 | $187K | 28 | $5.94M | 32× |
| 47 days | 64 | $192K | 53 | $5.93M | 31× |
| 12 days | 256 | $195K | 209 | $5.91M | 30× |
| 6 days | 512 | $195K | 419 | $5.93M | 30× |
| 1 day | 3037 | $193K | 2517 | $5.94M | 31× |

(The per-job total cost is roughly invariant in cluster size — the ratio is set by $/effective-GH/s and not by parallelism.)

### Equal-budget comparison
At a $200K budget:
- **f2**: 75K instance-hours → 12 days × 260 instances → **finishes in ~12 days**
- **p5**: 2034 instance-hours → 12 days × 7 instances → **finishes in ~4 years**

### Reproducing Shambles' actual 2019 spend (~$45K)
- At $45K on-demand f2: 17K instance-hours, ~10 instances over 71 days
- At $45K on-demand p5: 458 instance-hours, ~1 instance over 19 days — only 1/30th of the work done

## Sensitivity

The ratio is dominated by **effective rate**, not pricing. Plausible-bound table:

| FPGA effective | GPU effective | f2 $/GH/s-h | p5 $/GH/s-h | Ratio |
|---|---|---|---|---|
| 30 GH/s (pessimistic — routing/IO limit) | 12 GH/s (optimistic GPU) | $0.088 | $1.37 | 16× |
| **50 GH/s** (base) | **10 GH/s** (base) | **$0.053** | **$1.64** | **31×** |
| 80 GH/s (optimistic — VU47P fully utilized) | 8 GH/s (Shambles-like) | $0.033 | $2.05 | 62× |

Floor on the FPGA advantage: ~10×, even under pessimistic FPGA / optimistic GPU assumptions.

## What this model does NOT account for

1. **Engineering cost.** Shambles took ~6 person-years of crypto + software work for the GPU implementation. An equivalent FPGA stack does not exist publicly. For a one-shot attack, engineering dominates; for sustained or repeated workloads (research labs, internal red teams), amortization tips back to the per-GH/s ratio above.
2. **Memory-bandwidth bottleneck** at the distinguished-point storage stage. Shambles used ~500 GB of DP storage; f2.6xlarge's 256 GiB host RAM + PCIe Gen4 x16 (~32 GB/s) can support the DP throughput from `sha1_search.vhd`, but a longer chain or smaller DP bit-count would hit the PCIe ceiling.
3. **CPU-side path conditioning** (~2^33 work per chain) is small in compressions but tricky to parallelize; budget a single c7i.4xlarge running for ~hours per chain. Negligible vs. the FPGA/GPU costs.
4. **Spot pricing** would cut both columns by ~70%; the ratio is unchanged.
5. **Future hardware**: an MI300 or B200 would compress the GPU column but not change the structural argument — FPGAs win on $/Joule for fixed-function hashing.

## Recommendation

For a research lab or red team doing repeated CPC-class work on SHA-1 (or other Merkle-Damgard hashes), the FPGA path is **>10× cheaper at the floor of plausible assumptions**, and ~30× at the central estimate. The break-even engineering effort vs. one-shot cost is roughly **one full collision attack** — if you only run the attack once, GPUs are still cheaper because you don't pay for the HDL.
