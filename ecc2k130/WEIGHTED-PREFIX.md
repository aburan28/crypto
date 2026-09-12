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
independent routing. Mode 0 remains the preset default while the matched
GPU comparison is pending.

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
checks, not device execution or checkpoint serialization tests.

The CPU-only CUDA screen retained zero stack/local/shared bytes and used
122/101/100 walk registers for modes 0/1/2. Its examined live, non-reporting,
unguarded path counted 2,375.125 / 2,345 / 2,242 static instruction visits
per scalar update. The native carryless count was unchanged. These are
compiler observations; no throughput gain is established yet.

The schedule adds an original-Y read and uses the last prefix slot. The
common-step logical state-traffic increase is `20 + 40/B` bytes per scalar,
or 22.5 bytes at batch 16. This is not measured DRAM traffic. GPU arithmetic,
state/replay/checkpoint checks and matched timing are still pending.
