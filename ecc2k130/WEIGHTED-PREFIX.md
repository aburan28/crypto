# Weighted batch prefixes and paired Frobenius

`PACKED_WEIGHTED_PREFIX` is an optional packed-kernel schedule. Its default
is zero. Nonzero modes require polynomial state, polynomial chains, the
denominator cache and paired products. Mode 2 also requires the walk
permutation network.

| Mode | Forward pass |
|---:|---|
| 0 | Original denominator prefixes and separate coordinate permutations |
| 1 | Numerator-weighted prefixes with separate coordinate permutations |
| 2 | Numerator-weighted prefixes with a shared x/y permutation-mask stream |

To validate the optional mode through the public RTX preset:

```sh
ECC_PACKED_WEIGHTED_PREFIX=2 make audit-rtx-pro6000
```

The environment value is retained in the image, baked-build identity,
rebuild command and benchmark metadata. Arithmetic and timed records must
contain one matching schedule marker. The audit also requires the direct
paired-Frobenius GPU check to complete, with 6,240 pairs compared against
independent routing. The general default remains mode 0; the RTX preset
selects mode 2 after the controlled comparison below.

Let `D_i` and `E_i` be the existing coordinate differences, and let `P_i` be
the product of denominators before slot `i`. The new schedule stores
`W_i = E_i * P_i` in the existing prefix buffer. Slot zero directly copies
`D_0` and `E_0`. Later forward slots use a paired product to obtain `W_i`
and the next denominator product. The reverse pass pairs `inverse * W_i`
with the next inverse update, yielding the slope directly.

Every weighted prefix slot is initialized, including the formerly unused
last slot. Coordinates, reports, seed/restart/guard behavior and checkpoint
representation are retained. Prefix contents are transient scratch and are
not part of a checkpoint. The existing whole-batch zero-denominator behavior
is preserved.

At batch 16, the source schedule changes from 47 ordinary plus 15 paired
product calls to 17 ordinary plus 30 paired calls. Both retain the eight
normal-basis products in the inverse chain: **85 field products per batch**.
The change reduces field-function boundaries from 70 to 55; it does not
remove mathematical field products.

Mode 2 applies each loaded walk-permutation mask to both coordinates. The
generator owns the paired helper and retains the existing scalar and inverse
networks. The CPU compilation screen observed 56 global read-only load
instructions in both the scalar helper and the paired helper. The paired
helper handles two coordinates with one mask stream.

Host validation passed 18,720 paired-permutation calls across ordinary and
emulated CUDA preprocessing branches. Actual kernel host simulation passed
45 runs across modes 0/1/2, batches 1/4/16/32 and partial/tiled worker
boundaries. Persistent coordinates, seed/start/dead arrays, reports and
guards matched; canaries and undefined-behavior checks passed. The tests
covered all weighted slots and whole-batch zero poisoning. These are host
checks, not device execution or checkpoint serialization tests. The separate
GPU comparison below supplies device and checkpoint evidence.

The CPU-only CUDA screen retained zero stack/local/shared bytes and used
122/101/100 walk registers for modes 0/1/2. Its examined live, non-reporting,
unguarded path counted 2,375.125 / 2,345 / 2,242 static instruction visits
per scalar update. The native carryless count was unchanged. These are
compiler observations; measured throughput is reported separately below.

The schedule adds an original-Y read and uses the last prefix slot. The
common-step logical state-traffic increase is `20 + 40/B` bytes per scalar,
or 22.5 bytes at batch 16. This is not measured DRAM traffic.

## Complete-walk GPU comparison

The [retained comparison](benchmarks/weighted-prefix/comparison.json) ran
all three modes on one RTX PRO 6000 Blackwell Server Edition with CUDA
13.3.73, native sm_120 code and disabled PTX JIT. All timed rows completed
**201,863,462,912 scalar updates** with B16/T256/minBlocks2, 385,024 workers
and 32 launches of 1,024 steps.

Mode 2 passed the predefined screening threshold and three alternating
paired confirmations. Mode 1 did not qualify. Warmups and screening rows
are excluded from the following confirmation statistics:

| Workload | Control median B/s | Mode-2 median B/s | Ratio-of-medians gain |
|---|---:|---:|---:|
| Complete scalar benchmark | 13.110335 | **13.323276** | **1.6242%** |
| DP34 collection | 12.851526 | **13.054029** | **1.5757%** |

Benchmark ranges were 13.109104–13.113692 B/s for the control and
13.318348–13.338980 B/s for mode 2. Collection ranges were
12.847556–12.856025 and 13.050171–13.055579 B/s. Every pair favored mode 2.
All six collections produced the same sorted multiset: 5,149 records,
164,768 bytes and zero drops, with SHA-256
`ab237b6352380547fd37fcdc9e2aa83f6ae1a718b842590b9a19e4224e5d336b`.

Every mode passed the five established arithmetic suites and the direct
6,240-pair Frobenius test, full client replay/restart checks and common-state
checks at 64 and 16,384 walks. All 36 checkpoint children completed with
their expected outcomes, including three worker-mismatch rejections and
preserved incompatible input files. Bidirectional resumes across schedules
retained the same logical checkpoints.

The [independent retained-result audit](benchmarks/weighted-prefix/comparison-review.json)
passed all 19 timed rows, source/binary/native-code bindings, mode/resource
identities, exact counters and paired statistics. It verifies recorded
corpus hashes; transient corpus payloads were not retained. The raw result
SHA-256 is `1ad29dc750ce1b1a634e3f850d34e5a9bb1a67e40db8a0aea619e6fb6c08b8c2`.

This is a modest measured improvement; the 15 B/s target remains
unachieved. A separate audit of the updated public Make command is pending.
