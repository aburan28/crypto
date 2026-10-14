# ecbench rootfind vs-rho chase (2026-10-14)

**Decision: engineering.** Extends PR #1774 after complete-S ≤ F4.
Two measure cycles probed denser/sparser bases, window-negfold, and
plain rootfind. **Cold complete IC still loses to rho-strong**
(~1.54× on cycle-5 smoke; claim-v3 was ~1.55×). Best complete arm
remains `rootfind-negfold` m=2 d=5. No security claim; no ECC2K-130
inference.

| | |
|---|---|
| protocol | [`PROTOCOL.md`](PROTOCOL.md) |
| cycles | [`cycles/`](cycles/) (4–5) |
| binary | `ecbench` release SHA-256 `27fd2653b706b6e82e8da1c39a8c26c8f686e3f6e4831953b85eb9798610b76c` |
| host | Apple M4 Pro; `--cpus none --allow-busy` → **L0** |
| classification | **engineering** |
| Conductor | T-2575 (daemon unreachable at start; prior lease) |

## Cycle log

| cycle | change | outcome |
|---|---|---|
| 4 | d=6 negfold; window-negfold div=2\|3 | **dead end** — all worse than d=5 negfold |
| 5 | d=4/d=7 negfold; plain rootfind m=2 d=5 | **dead end** — d=5 negfold still best complete (~1.54× rho) |

## Hard stop vs rho (unchanged)

| arm | cycle-5 mean S | LB? | notes |
|---|---:|---|---|
| rho-strong | 14.145 | yes | unpriced canons/table ops |
| ic-c2-frob-counted-m2 | 15.456 | yes | FB natives unpriced |
| ic-negfold-m2-d5 | 21.803 | **no** | best complete; ~1.54× rho |
| ic-f4-m3-d5 (claim-v3) | 421.369 | **no** | XOR-priced; still ≫ negfold |

Negfold GAE split (claim-v3): oracle_setup ≈ 305 (pair table),
relations ≈ 700 (walk + trials). Dimension sweeps move setup vs hit
rate the wrong way for cold S. Secondary amS@64 can dip near/below
rho (d=6 ≈ 13.5) — not a cold claim per hygiene.

## Reproduce

```bash
export CARGO_TARGET_DIR=/Volumes/SSD990/wt-ecbench-rootfind-improve/target-faster
cargo build --release --bin ecbench
ECBENCH=$CARGO_TARGET_DIR/release/ecbench bash research/ecbench_rootfind_vs_rho_20261014/run.sh cycle4
ECBENCH=$CARGO_TARGET_DIR/release/ecbench bash research/ecbench_rootfind_vs_rho_20261014/run.sh cycle5
```
