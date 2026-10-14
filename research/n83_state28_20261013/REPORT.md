# N83 compact states: verified K512 run and resource boundary

**The 28-byte N83 state layout completed a 512-column, 84,992-point
prefix-base run with full rank, a four-point target relation and independent
scalar replay.** The archived base's distinct points, subgroup
representatives, signed-Frobenius orbits and coefficient labels also passed
an independent generic GF(2^83) checker. The one public target needed 650,715
probes and 1,533.656333 ms online after preparation. This is an algorithmic
correctness and memory result; the host was shared.

| Control arm | Usable points | Folded columns | State bytes each | Allocated state bytes | Peak child RSS bytes | Wall s | Status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| K512 wide | unknown for incomplete run | 512 planned | 48 | unknown | 2,143,174,656 | 300.110 | TIMEOUT |
| K512 packed | 84,992 | 512 | 28 | 609,222,656 | 1,696,169,984 | 291.115 | PASS_INDEPENDENT_REPLAY |

The paired peak-RSS difference is 447,004,672 bytes. The source-level
element-width saving at the measured 21,757,952 states is 435,159,040
bytes; the remaining RSS difference reflects other allocations and allocator
behavior. The wide control reached the five-minute process cap before
producing a target row, so this pair supplies no completed wall-time ratio.
The packed arm spent 36.444 s in setup, 252.863 s in rank and 0.000885 s in
final matrix algebra. Its target online interval is separate from these
reusable stages.

## Larger-base resource guard

The K2202 compact-state attempt was stopped after 653.425 s when macOS
reported 6.5 GiB of its pages swapped out and the system data volume fell
to 2.2 GiB free, then 1.2 GiB after the stop. The direct-child wrapper
recorded signal `-15`; the shell recorded its modulo-256 status `241`.
The next K1798 arm began automatically and was stopped after 36.019 s while
swap use remained above 8 GiB. Neither arm emitted a summary or target row.
The [independent status analysis](analysis_primary.json) classifies both as
`STOPPED_RESOURCE_GUARD`. The raw [post-guard snapshot](host_after_guard.txt)
records disk, swap and memory state. The K2202 `vmmap -summary` observation
before stopping showed a 6.6 GiB physical footprint, including 6.5 GiB
swapped-out writable pages; its terminal output is retained in the task
record. The wrapper's 1.60 GiB peak resident-set reading consequently does
not represent the process's swapped working set.

The larger factor-base comparison has no completed primary row. The K1798
dense layout from the earlier [64-point study](../n83_dense_root_20261013/REPORT.md)
remains the lowest-mean-probe completed prefix base in that frozen panel.
The new five-byte indexed-slot design will spend two additional GiB of
table storage at K2202 to lower the bound on table occupancy while retaining
exact root checks against the compact state array. It remains an unmeasured
implementation lead in this report.

## Method and verification

The [frozen protocol](PROTOCOL.md) fixes the public point, source and binary
hashes, controls, primary arm order, rank seed, threads, caps and accounting.
The accepted target point is
`(6156598173116761880680889,2569387936113339852934176)` on the
Koblitz `a=b=1` curve over `F_2^83`, polynomial
`z^83+z^7+z^4+z^2+1`, subgroup order `8569786107849059`.
The packed control recovered scalar `2927171373054142`; the
[checker](analyze.py) replays that scalar from the generator independently,
checks all 84,992 archived base points against their orbit construction,
and sums the four relation points to the public target.

The [Lean proof map](PROOF_MAP.md) identifies the three built conditional
numeric declarations for the 11-byte root encoding and 28-byte state bound.
The source-loop and encoding bridges remain open proof obligations. Thirteen
focused release example tests passed, including root round trips at the
83-bit boundary. [TESTS.md](TESTS.md) lists the commands and evidence.
