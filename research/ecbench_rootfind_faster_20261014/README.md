# ecbench rootfind actual-speed loop (2026-10-14)

**Decision: engineering.** Extends PR #1774. Inherited-F4 now prices pinned
word XORs into admitted S, so **complete** `rootfind-negfold` beats
**complete** F4 on claim-n13-v3 (22.97 ≤ 421.37). Wall column on
`ecbench table`. Rho still wins cold S (hard stop ~1.55× for complete IC).
No security claim; no ECC2K-130 inference.

| | |
|---|---|
| protocol | [`PROTOCOL.md`](PROTOCOL.md) |
| cycles | [`cycles/`](cycles/) (1–3 + claim-v3) |
| binary | `ecbench` release SHA-256 `27fd2653b706b6e82e8da1c39a8c26c8f686e3f6e4831953b85eb9798610b76c` |
| host | Apple M4 Pro; `--cpus none --allow-busy` → **L0** |
| claim audit | `audits/audit-claim-n13-v3.json` |
| classification | **engineering** |
| Conductor | T-2575 |

## Cycle log

| cycle | change | outcome |
|---|---|---|
| 1 | word-XOR op_unit prefix match | **win**: F4 complete; negfold ≪ F4 |
| 2 | `ecbench table` med wall ms | **win**: wall beside S |
| 3 | negfold d=4 smoke | modest; claim keeps d=5 best complete |
| claim-v3 | 8×3 panel | complete-S parity **Y**; vs rho **N** |

## Before (parity claim-n13-v2 @ `3f4a87aae0`)

| arm | mean S | lower bound |
|---|---:|---|
| rho-strong | 14.359 | yes |
| ic-c2-frob-counted-m2 | 15.552 | yes |
| ic-f4-m3-d5 | 16.749 | yes |
| ic-negfold-m2-d5 | 23.220 | no |

## After (faster claim-n13-v3)

| arm | mean S | lower bound | med wall ms |
|---|---:|---|---:|
| rho-strong | 14.822 | yes | 0.085 |
| ic-c2-frob-counted-m2 | 15.648 | yes | 0.169 |
| **ic-negfold-m2-d5** | **22.972** | **no** | **0.262** |
| ic-negfold-m2-d4 | 25.729 | no | 0.327 |
| **ic-f4-m3-d5** | **421.369** | **no** | 154.013 |

Paired: `negfold-d5 / F4` = 0.0545 [0.047, 0.064] (fully priced).
Wall median ratio ≈ 0.0017 (L0 descriptive).

## Hard stop vs rho

Cold complete IC does not beat rho-strong (~1.55×). Best LB compact-frob
is ~1.06× rho but remains a lower bound (FB-construction natives unpriced).
Amortized secondary amS@64 for counted-c2 ≈ 11.2 (beats rho) — not a cold
claim per hygiene.

## Reproduce

```bash
export CARGO_TARGET_DIR=/Volumes/SSD990/wt-ecbench-rootfind-improve/target-faster
cargo build --release --bin ecbench
ECBENCH=$CARGO_TARGET_DIR/release/ecbench bash research/ecbench_rootfind_faster_20261014/run.sh cycle1
ECBENCH=$CARGO_TARGET_DIR/release/ecbench bash research/ecbench_rootfind_faster_20261014/run.sh cycle3
ECBENCH=$CARGO_TARGET_DIR/release/ecbench bash research/ecbench_rootfind_faster_20261014/run.sh claim
```
