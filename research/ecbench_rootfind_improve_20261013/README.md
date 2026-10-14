# ecbench rootfind improve loop (2026-10-13)

**Decision: engineering.** Extends PR #1765
(`research/ecbench_rootfind_engine_20261012/`). New engine
`pdp3-koblitz:m=3,engine=rootfind-window` and a smaller factor-base
dimension cut cold `S` without changing the asymptotic group-operation
count of the decomposition. Against inherited-F4, rootfind `S` remains a
complete count beside a lower bound and is **not** a speedup; descriptive
L0 wall still favours rootfind. No security claim; no ECC2K-130 inference.

| | |
|---|---|
| protocol | [`PROTOCOL.md`](PROTOCOL.md) |
| cycles | [`cycles/cycle-1.md`](cycles/cycle-1.md), [`cycles/cycle-2.md`](cycles/cycle-2.md) |
| engine ids | `pdp3-koblitz:m=2\|3,engine=rootfind` (unchanged); **`pdp3-koblitz:m=3,engine=rootfind-window`** (new); `rootfind-subtract:m=2\|3` (unchanged) |
| binary | `ecbench` release SHA-256 `394e0c1bfe705f2bb0795425099beb02c50f04abe5027fe8b4020b91d40ec973` |
| host | Apple M4 Pro; env class `ECBENV2h59034d4cd61b`; `--cpus none --allow-busy` → **L0** |
| classification | **engineering** |

## What improved / what did not

| change | cold `S` | notes |
|---|---|---|
| `dimension=4` + `m=2` vs #1765/`d=5` `m=2` | modest ↓ (~5%, CI includes 1 on claim panel) | setup GAE 190 vs 561 |
| `rootfind-window` vs full `m=3` at `d=5` | **↓ ~17%** (ops ratio 0.828, CI excludes 1) | same table build; shorter third-summand scan |
| `rootfind-window` at `d=4` | worse | cycle-2 **dead end** |
| `dimension=3` + `m=2` | worse than d=4/d=5 | cycle-1 dead end for cold S |
| vs inherited-F4 in `S` | still higher | complete count vs lower bound |
| vs F4 in solve wall (L0) | ~250–300× lower (descriptive) | not admitted |
| multi-target amortization | secondary only | d=5 benefits most when setup is shared |

## Sessions

| session | id | records | audit | replay cert (SHA-256 of audit JSON) |
|---|---|---:|---|---|
| [`sessions/select`](sessions/select) | `ECBS1h5cb3bc5c9cd5` | 42/42 | OK | `8421ef064ead9419c7ab8a520ac8f968fdd0f9fe37c75e1ef4474f249460401a` |
| [`sessions/window-d4`](sessions/window-d4) | `ECBS1hbac35523d128` | 30/30 | OK | `2c1bf6fd6517b4cd79263eb67721bef32d370e500db3cac163b1941364a71506` |
| [`sessions/claim-n13`](sessions/claim-n13) | `ECBS1h5dc382c51ae9` | 256/256 | OK, 192/192 replays identical | `fd0572768b02a45defe7d5a70afc638c1c686eca184c95c477d31309f49778aa` |

## Before / after (Koblitz `E_0`, `n=13`)

### #1765 claim-n13 (before; different seeds)

| arm | mean S | S / strong-rho | lower bound |
|---|---:|---:|---|
| `rho-strong` | 14.118 | 1.000 | yes |
| `ic-f4-m3` (d=5) | 15.928 | 1.128 | yes |
| `ic-rootfind-m2` (d=5) | 28.383 | 2.010 | no |
| `ic-rootfind-m3` (d=5) | 38.822 | 2.750 | no |

### This improve claim-n13 (after; 8 public × 3 rounds)

| arm | mean S | S / strong-rho | lower bound |
|---|---:|---:|---|
| `rho-strong` | 14.342 | 1.000 | yes |
| `ic-f4-m3-d5` | 16.322 | 1.138 | yes |
| **`ic-rootfind-m2-d4`** | **28.461** | 1.985 | no |
| `ic-rootfind-m2-d5` (#1765 baseline) | 29.927 | 2.087 | no |
| **`ic-rootfind-window-m3-d5`** | **32.671** | 2.278 | no |
| `ic-rootfind-m3-d5` (full scan) | 39.462 | 2.752 | no |

Paired ops (claim-n13):

- `ic-rootfind-m2-d4 / ic-rootfind-m2-d5` = 0.951 [0.855, 1.061]
- `ic-rootfind-window-m3-d5 / ic-rootfind-m3-d5` = **0.828** [0.783, 0.863]
- `ic-rootfind-m2-d4 / ic-f4-m3-d5` = 1.744 [1.541, 1.939] (**bounded**; not a speedup)
- Descriptive wall median `rootfind-m2-d4 / f4-m3-d5` ≈ 0.0032 (L0)

### Secondary amortized reading (setup shared across N targets)

From [`tables/amortized-claim-n13.txt`](tables/amortized-claim-n13.txt)
(`analyze_resources.py`). Not a cold one-target claim.

| arm | cold S | amS@8 | amS@64 | solve wall ms (L0) |
|---|---:|---:|---:|---:|
| ic-rootfind-m2-d4 | 27.7 | 24.0 | 23.5 | 0.90 |
| ic-rootfind-m2-d5 | 29.6 | 18.6 | 17.2 | 1.11 |
| ic-rootfind-window-m3-d5 | 32.8 | 21.8 | 20.4 | 1.09 |
| ic-f4-m3-d5 | 16.4 | 16.3 | 16.3 | 228 |

## Reproduce

```bash
cargo build --release --bin ecbench
ECBENCH=target/release/ecbench bash research/ecbench_rootfind_improve_20261013/run.sh select
ECBENCH=target/release/ecbench bash research/ecbench_rootfind_improve_20261013/run.sh claim
python3 research/ecbench_rootfind_improve_20261013/analyze_resources.py \
  --records research/ecbench_rootfind_improve_20261013/sessions/claim-n13/records.jsonl
```
