# Pursuing 22 B/s on one RTX PRO 6000

Target declared before this round: exceed **22.0 billion complete scalar
updates/s** with a correct table walk on one RTX PRO 6000. The reference is
the verified 18.110 B/s composed result in
[`benchmarks/throughput-18b-composed`](benchmarks/throughput-18b-composed/README.md).

Result: **20.134 B/s**, a verified **1.1084×** gain over the equal-work
18.165 B/s paired reference, but only **0.9152 of the 22 B/s target**. The
target was not met.

## Boundary and falsification condition

The unit is billions of complete scalar updates per second from `finished:`.
The derived one-affine-addition practical boundary remains approximately
22--25 B/s on this device. This is walk-kernel engineering, not a change to
Pollard rho's generic-group operation count.

Success required median `>22.0`, 300/300 device reports reproduced by the
scalar reference, zero dropped reports, and unchanged walk accounting.
Changing GPUs, omitting selection/replay work, counting field operations, or
reporting a scout as a result was inadmissible.

## Single table

| variant | median B/s | / 22 target | / paired reference | products/update | correct | class |
|---|---:|---:|---:|---:|---|---|
| 18 B composed table walk, B16, 256×2 | 18.164748 | 0.8257 | 1.0000 | 5.3125 | 300/300, 0 dropped | reference |
| inline inversion wrapper | 18.041 | 0.8200 | 0.993 | 5.3125 | arithmetic path unchanged | engineering, rejected |
| one 512-thread block, B16 | 19.208 | 0.8731 | 1.058 | 5.3125 | scout | engineering |
| selection/product pipeline, 120 registers | 19.210 | 0.8732 | 1.058 | 5.3125 | 300/300 | engineering, neutral |
| half-width top CLMAD correction | 18.626 | 0.8466 | 1.026 | 5.3125 | arithmetic differential pass | engineering, rejected |
| denominator recomputation, B16 | 19.816 | 0.9007 | 1.091 | 5.3125 | 300/300 | engineering |
| B17, all persistent fields in L2 | 20.131 | 0.9150 | 1.108 | 5.2941 | scout | engineering |
| **B17 + shared-bank row padding** | **20.134326** | **0.9152** | **1.1084** | **5.2941** | **300/300, 0 dropped** | **engineering, target not met** |
| warp-cooperative LUT gather | 17.944 | 0.8156 | 0.988 | 5.2941 | 300/300 | engineering, rejected |
| B18, 5% beyond persist window | 19.956 | 0.9071 | 1.099 | 5.2778 | scout | engineering, rejected |
| one-addition practical boundary | 22--25 | 1.00--1.14 | 1.21--1.38 | derived | — | boundary |

The reference and selected candidate each processed exactly
53,619,982,336 updates per confirmation sample: 34 launches at B16 and 32
launches at B17. Six alternating pairs all favored the candidate. Its range
was 20.090--20.185 B/s; the paired geometric mean ratio was 1.10848 with an
exploratory 95% log-ratio t interval `[1.10537, 1.11161]`.

## What moved

1. **One 512-thread block per SM.** It retains 16 active warps but gives ptxas
   one launch-bounds block instead of two 256-thread blocks. The 576/640
   experiments forced allocation down to 96 registers and regressed.
2. **Recompute table denominators.** The latest table tag is already in
   history. The reverse pass rebuilds `x+x_T` with four shared reads instead
   of writing and rereading a fourth 17-byte field. The persistent blob falls
   from 104.7 MB to 78.5 MB at B16.
3. **Batch 17.** Three persistent fields occupy 83,453,952 bytes, just below
   the 83,886,080-byte persisting-L2 cap, while inversion amortization improves
   from 5.3125 to 5.2941 products/update.
4. **Shared-bank padding.** The addend phase-row stride changes from 66 words
   (16 reachable start banks) to 67 (all 32). Polynomial-to-normal sign rows
   similarly change from four to five words. This is a small measured gain.

The denominator lever changes storage, not arithmetic: reconstruction uses the
same selected table point, and full device reports replay identically. The
batch change slightly lowers the applicable arithmetic floor; the table keeps
the product count visible so that movement is not mislabeled as a same-floor
instruction gain.

## Why 22 remains open

The selected Nsight profile reports 56.38% issue slots busy, 70.4% FP64/CLMAD,
42.81% cycles with no eligible warp, and 1.31 eligible warps per scheduler.
L2 hits 97.57%. Random shared addresses still produce 58% excessive
wavefronts. A warp-cooperative gather removed those conflicts but its shuffle
cost was larger, dropping throughput to 17.944 B/s.

The remaining gap is 9.27%. No measured row in this round prices that gain.
Reaching 22 now requires a cheaper selector/table gather or a reduction/product
schedule that raises eligible warps without moving enough work onto the
70%-utilized CLMAD pipe. The target is not extrapolated from the 20.134 result.

## Evidence

[`benchmarks/throughput-22b-attempt/result.json`](benchmarks/throughput-22b-attempt/result.json)
freezes configuration, samples, ratios, profile counters and binary hashes.
`balanced-comparison.log`, verification logs, build logs, profile output and
all retained scout/regression logs are in the same directory.
