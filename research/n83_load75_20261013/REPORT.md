# N83 75% root-table request: K2000 completes in 8 GiB

**The 75% request rule extends the completed N83 prefix base to K2000 with
332,000 distinct subgroup-usable points and 2,000 signed-Frobenius columns
in the same 8 GiB packed root-table bucket as K1798.** On the frozen single
public target, both candidates reached full rank and independently recovered
scalar `2927171373054142`. K1798 used 111,004 target probes and 59.985167
ms online; K2000 used 334,563 probes and 150.499459 ms online. The same-point
strong-rho receipt recorded 1,173.29325 ms online. These shared-host CPU
ratios are exploratory; the [cloud allocation receipt](cloud_tier_attempt.json)
records why a Tier 0/1 rerun was unavailable.

| Primary arm | Usable points | Folded columns | Regular states | Exact root keys | Packed slots | Rank | Target probes | Verified online ms | Setup / rank / LA s | Status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| K1798 dense control | 298,468 | 1,798 | 268,322,732 | 268,316,606 | 8 GiB | 1,798 | 111,004 | 59.985167 | 381.113 / 21.253 / 0.010 | PASS |
| K2000 load75 | 332,000 | 2,000 | 332,000,000 | 331,992,947 | 8 GiB | 2,000 | 334,563 | 150.499459 | 579.381 / 16.640 / 0.011 | PASS |
| K2202 load75 | unknown | 2,202 planned | unknown | unknown | 8 GiB bound | unknown | unknown | unknown | unknown | INTERRUPTED_RESOURCE |

The public point is
`(6156598173116761880680889,2569387936113339852934176)` on the
Koblitz `a=b=1` curve over `F_2^83`, polynomial
`z^83+z^7+z^4+z^2+1`, subgroup order `8569786107849059`.
The [independent primary checker](analysis_primary_v2.json) validated each
completed scalar with generic curve multiplication, the exact online phase
sum, same-point rho identity, complete rank, and K1798 first-witness parity
against the prior frozen binary. The K1798 and K2000 local arithmetic
rho/IC ratios are 19.5597 and 7.7960 respectively. Their online clocks
exclude reusable base/index/log preparation. Peak RSS for these primary
rows is missing because the corrected macOS `time -p` wrapper reports wall
and CPU time without a peak-memory counter.

## Secondary 64-point panel and resource receipts

The frozen reverse-order panel attempted K2000 then K1798 on the identical
64 public points. Both reached the 900-second process cap before emitting
target rows. The [independent panel checker](analysis_panel.json) retains
both as TIMEOUT, marks K2202 as not run under its resource guard, and leaves
the panel selection null. No mean probe comparison is computed from these
censored rows.

| Panel arm | Exit | Elapsed s | Child CPU user+system s | Peak RSS B | Involuntary switches | Target rows |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| K2000 load75 | 124 | 901.027 | 336.555 | 7,796,932,608 | 3,636,481 | 0 |
| K1798 dense | 124 | 900.936 | 277.456 | 8,409,251,840 | 2,544,226 | 0 |

The [direct-child wrapper](run_with_usage.py) records the kernel peak RSS,
elapsed and CPU time, exit code, faults and context switches in each
`runs/panel/<arm>/usage.json`. The large wall-versus-child-CPU gap and switch
counts are observed host conditions. The raw outputs and zero-row files
remain under `runs/panel/`.

## Capacity boundary and next implementation

The producer asks its existing twofold, power-of-two allocator for
`ceil(2*regular_states/3)` entries under `KIC_PACKED_ROOT_TABLE_LOAD75=1`.
The merged [Lean 4 capacity proof](../n83_dense_root_20261013/Load75Capacity.lean)
checks that this request fits `2^29` slots at K2202 under the source count
bound `states <= 83*K^2`, while the K2203 maximum crosses that bucket. The
K2000 run measured 331,992,947 unique keys in 536,870,912 allocated slots.

K2202 entered macOS uninterruptible I/O after about 9.5 minutes, with all
13,056 MiB of configured swap used and no producer target or summary row.
The owned process was stopped to release host pressure; its partial primary
directory, process snapshot and swap snapshot are retained. The first
outer-timer failure and all run-plan amendments are in [AMENDMENT.md](AMENDMENT.md).

The next implementation stores each N83 state in 28 bytes, encoding its two
83-bit roots as 11 bytes each. The current struct is 48 bytes on this ARM64
build. At the K2202 maximum, the element-width difference is over 7 GiB;
the separate state28 proof map records the numeric Lean results and the
remaining source bridge. A new frozen run will
recheck first-witness identity and the complete one-target solve before any
online comparison.

## Reproducibility

[PROTOCOL.md](PROTOCOL.md) fixes the public points, scalar sidecar, source
and binary SHA-256, arm order, rank seed 29, 12 rank threads, one target
thread and process caps. The first producer source commit is `63264cc2de`;
the corrected and resource-guarded panel plans were committed before their
corresponding attempts. The source and release binary hashes remain
`d4bf0f34d9d57eb26ed795a610be4baa469b3835455d69ff902ef46a20036fdd`
and `fc4b62da8e5dc7980c3e2a27e774166d594b5661641833252e8ec8ae6894fa40`.
All raw records, including the failed timer attempt, interrupted K2202,
two panel timeouts, exact resource receipts, checker outputs and command logs
are retained. [TESTS.md](TESTS.md) lists the actual verification commands.
