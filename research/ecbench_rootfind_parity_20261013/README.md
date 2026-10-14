# ecbench rootfind IC row-parity loop (2026-10-13)

**Decision: engineering.** Extends PR #1774 / improve loop. New engines
`rootfind-negfold` and `rootfind-window-negfold` cut PairSumTable set-up by
building `PairTable::build_negation_folded` (~|F|²/4 adds). Best cold rootfind
S on the claim panel: **24.462** (~1.70× rho) vs F4 **16.480** (~1.14× rho).
**Cold IC row parity: not achieved.** No security claim; no ECC2K-130 inference.

| | |
|---|---|
| protocol | [`PROTOCOL.md`](PROTOCOL.md) |
| cycles | [`cycles/cycle-1.md`](cycles/cycle-1.md), [`cycles/cycle-2.md`](cycles/cycle-2.md), [`cycles/cycle-3.md`](cycles/cycle-3.md) |
| engines | `rootfind-negfold` (new); `rootfind-window-negfold` + required `window_div` (new); prior `rootfind` / `rootfind-window` unchanged |
| binary | `ecbench` release SHA-256 `1622e3fa8e638031c8406b651a8e2e596573189bd2c78aa8b83d5e979fa372dc` |
| host | Apple M4 Pro; env `ECBENV2h59034d4cd61b`; `--cpus none --allow-busy` → **L0** |
| claim audit | SHA-256 `6bf26b3be576ab7490419d6ffcbdcd499198db3ea43362902359ecbb66771399` |
| classification | **engineering** |

## Cycle log

| cycle | change | cold S outcome |
|---|---|---|
| 1 | `rootfind-negfold` | **win**: m2/d5 mean S 23.5 vs full 29.2 (smoke) |
| 2 | `window_div` ∈ {2,3,4} | **dead end** vs m2 negfold |
| 3 | claim panel 8×3 | negfold-m2-d5 **24.462**; parity **N** vs F4 16.480 |

## Claim-n13 table (this loop)

| arm | mean S | S / strong-rho | lower bound |
|---|---:|---:|---|
| rho-strong | 14.414 | 1.000 | yes |
| ic-f4-m3-d5 | 16.480 | 1.143 | yes |
| **ic-negfold-m2-d5** | **24.462** | 1.697 | no |
| ic-negfold-m2-d4 | 25.992 | 1.803 | no |
| ic-rootfind-m2-d4 (#1774 best full) | 27.802 | 1.929 | no |
| ic-rootfind-m2-d5 | 30.182 | 2.094 | no |

vs improve-loop best (full rootfind m2/d4 mean S 28.461): negfold m2/d5 is
lower on this panel's seeds; vs same-panel full m2/d5 ops ratio **0.810**.

## Sessions

| session | id | records | audit |
|---|---|---:|---|
| [`sessions/cycle1`](sessions/cycle1) | `ECBS1hcd05f1dc18a3` | 48/48 | OK |
| [`sessions/cycle2`](sessions/cycle2) | `ECBS1hd0002abcd05a` | 36/36 | OK |
| [`sessions/claim-n13`](sessions/claim-n13) | `ECBS1h5834884bc326` | 256/256 | OK |

## Secondary amortized (not a cold claim)

| arm | cold S | amS@8 | amS@64 |
|---|---:|---:|---:|
| ic-negfold-m2-d5 | 24.3 | 18.3 | 17.5 |
| ic-f4-m3-d5 | 16.4 | 16.3 | 16.3 |

## Caveats

- F4 `S` is a lower bound (unpriced Boolean XORs); rootfind/negfold `S` is complete.
- Do not quote `negfold/F4` ops as a speedup.
- Wall figures are L0 / descriptive only.
- Multi-target amortization is secondary analysis only.
- Conductor daemon was down at start (`127.0.0.1:8443` refused).

## Reproduce

```bash
export CARGO_TARGET_DIR=/Volumes/SSD990/wt-ecbench-rootfind-improve/target-parity
cargo build --release --bin ecbench
ECBENCH=$CARGO_TARGET_DIR/release/ecbench bash research/ecbench_rootfind_parity_20261013/run.sh cycle1
ECBENCH=$CARGO_TARGET_DIR/release/ecbench bash research/ecbench_rootfind_parity_20261013/run.sh cycle2
ECBENCH=$CARGO_TARGET_DIR/release/ecbench bash research/ecbench_rootfind_parity_20261013/run.sh claim
```
