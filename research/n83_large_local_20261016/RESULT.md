# N83 K4096 packed root-state limit

The frozen current-source K4096 dynamic one-root diagnostic stopped at the
producer's state-ID assertion before base or index allocation. Its exact
state count is 696,256,512, while the four- and five-byte packed slots
encode fewer than 2²⁹ = 536,870,912 states. The maximum folded-column count
for that encoding is K3596: K3596 has 536,647,262 states and K3597 has
536,945,772. The Lean theorem
`N83OnlineRange.packed_state_id_excluded_above_3596` proves the exclusion
for every natural K ≥ 3597 from the source's state-count formula.

The local producer exited 101 after 1.547 seconds, with sampled peak RSS
52,445,184 bytes. Its exact stderr says `one-root complete key bound exceeds
fixed indexed slots`. The independent classifier recorded
`UNSUPPORTED_STATE_ID`; it is a format result, not a timing comparison or
a candidate selection. The raw receipts are retained at
`/Volumes/SSD990-2/n83-archive-20261016/local-large/k4096-run1/`.

| Receipt | SHA-256 |
| --- | --- |
| `analysis.json` | `0eb4ee7953802302757a5292c3ce05ff551fb5fb9a4817b7ddef8d1ff9fbd9fd` |
| `runs/k4096root11dynamic/stderr.log` | `ed34b34725dd15e77025586ab1decbb354bdc7fd5be580e31090a4f8b45073d2` |

The current cloud range screen tests K3584 near the encoding ceiling. The
recorded K3596 four-byte layout from the earlier cross-layout panel is the
next exact-format endpoint for a same-point online comparison on the current
source. A larger K requires a different state-ID encoding.
