# ecbench rootfind actual-speed loop (2026-10-14)

**Decision: engineering.** Extends PR #1774. After cycle-1, inherited-F4
prices pinned word XORs into admitted S, so **complete** rootfind-negfold
beats **complete** F4 on claim-scale costing (smoke: 22.2 ≤ 540.7). Wall
column added to `ecbench table`. Rho still wins cold S. No security claim.

| | |
|---|---|
| protocol | [`PROTOCOL.md`](PROTOCOL.md) |
| cycles | [`cycles/`](cycles/) |
| binary | `ecbench` release SHA-256 `27fd2653b706b6e82e8da1c39a8c26c8f686e3f6e4831953b85eb9798610b76c` |
| host | Apple M4 Pro; `--cpus none --allow-busy` → **L0** |
| classification | **engineering** |

## Cycle log

| cycle | change | outcome |
|---|---|---|
| 1 | word-XOR op_unit prefix match | **win**: F4 complete S~541; negfold 22.2 ≤ F4 |
| 2 | `ecbench table` med wall ms column | **win**: wall visible beside S |
| 3 | negfold d=4 | (pending) |

## Before (parity claim-n13-v2)

| arm | mean S | lower bound |
|---|---:|---|
| rho-strong | 14.359 | yes |
| ic-c2-frob-counted-m2 | 15.552 | yes |
| ic-f4-m3-d5 | 16.749 | yes |
| ic-negfold-m2-d5 | 23.220 | no |

## After cycle-1 smoke (XOR priced)

| arm | mean S | lower bound | med wall ms |
|---|---:|---|---:|
| rho-strong | 13.961 | yes | 0.073 |
| ic-c2-frob-counted-m2 | 15.513 | yes | 0.284 |
| ic-negfold-m2-d5 | 22.228 | no | 0.244 |
| ic-f4-m3-d5 | 540.708 | **no** | 224.530 |

## Reproduce

```bash
export CARGO_TARGET_DIR=/Volumes/SSD990/wt-ecbench-rootfind-improve/target-faster
cargo build --release --bin ecbench
ECBENCH=$CARGO_TARGET_DIR/release/ecbench bash research/ecbench_rootfind_faster_20261014/run.sh cycle1
```
