# ecbench rootfind IC row-parity loop (2026-10-13)

**Decision: engineering.** Extends PR #1774. After cycles 4–6 the admitted
cold candidate is `compact-orbit-scan:columns=2` +
`mitm-frobenius-counted:m=2` with mean S **15.552** vs F4 m=3/d=5
**16.749** on claim-n13-v2. **Cold IC row parity: Y under harness admitted
S** (both arms remain lower bounds). Best *complete*-S rootfind on
standard-subspace is still negfold **23.220**. No security claim; no
ECC2K-130 inference.

| | |
|---|---|
| protocol | [`PROTOCOL.md`](PROTOCOL.md) |
| cycles | [`cycles/`](cycles/) (1–6) |
| engines / oracles | prior `rootfind-negfold`; this loop uses existing `mitm-frobenius` / `mitm-frobenius-counted` on `compact-orbit-scan` (E₀ has no `koblitz-orbit:divisor=0;1`) |
| binary | `ecbench` release SHA-256 `1622e3fa8e638031c8406b651a8e2e596573189bd2c78aa8b83d5e979fa372dc` |
| host | Apple M4 Pro; env `ECBENV2h59034d4cd61b`; `--cpus none --allow-busy` → **L0** |
| claim audit | `audits/audit-claim-n13-v2.json` |
| classification | **engineering** |

## Cycle log

| cycle | change | cold S outcome |
|---|---|---|
| 1 | `rootfind-negfold` | **win**: m2/d5 mean S 23.5 vs full 29.2 (smoke) |
| 2 | `window_div` ∈ {2,3,4} | **dead end** vs m2 negfold |
| 3 | claim panel 8×3 | negfold-m2-d5 **24.462**; parity **N** vs F4 16.480 |
| 4 | orbit divisor + subtract | **dead end** (no invariant subspace; subtract S~120) |
| 4b | compact-orbit columns sweep | **win** vs negfold; c2-frob ~14.1 (LB) |
| 5 | counted c2/c3 | counted-c2 **15.271** ≤ F4 15.950 (smoke) |
| 6 | claim-n13-v2 | counted-c2 **15.552** ≤ F4 **16.749**; parity **Y** (LB/LB) |

## Claim-n13-v2 table

| arm | mean S | S / strong-rho | lower bound |
|---|---:|---:|---|
| rho-strong | 14.359 | 1.000 | yes |
| **ic-c2-frob-counted-m2** | **15.552** | 1.083 | yes |
| ic-f4-m3-d5 | 16.749 | 1.166 | yes |
| ic-c3-frob-counted-m2 | 18.643 | 1.298 | yes |
| ic-negfold-m2-d5 | 23.220 | 1.617 | no |

## Sessions

| session | id | records | audit |
|---|---|---:|---|
| [`sessions/cycle1`](sessions/cycle1) | `ECBS1hcd05f1dc18a3` | 48/48 | OK |
| [`sessions/cycle2`](sessions/cycle2) | `ECBS1hd0002abcd05a` | 36/36 | OK |
| [`sessions/claim-n13`](sessions/claim-n13) | `ECBS1h5834884bc326` | 256/256 | OK |
| [`sessions/cycle4`](sessions/cycle4) | `ECBS1h6718277fab42` | 30/48 | OK (18 orbit errors) |
| [`sessions/cycle4b`](sessions/cycle4b) | `ECBS1h5ae06f07ab10` | 48/48 | OK |
| [`sessions/cycle5`](sessions/cycle5) | `ECBS1h3600bd37dbee` | 42/42 | OK |
| [`sessions/claim-n13-v2`](sessions/claim-n13-v2) | `ECBS1h2ef2bb46ad6c` | 256/256 | OK |

## Secondary amortized (not a cold claim)

| arm | cold S | amS@8 | amS@64 |
|---|---:|---:|---:|
| ic-c2-frob-counted-m2 | 15.57 | 11.65 | 11.16 |
| ic-f4-m3-d5 | 16.64 | 16.57 | 16.56 |
| ic-negfold-m2-d5 | 23.28 | 17.25 | 16.49 |

## Caveats

- F4 `S` is a lower bound (unpriced Boolean XORs); compact-frob `S` is
  also a lower bound (unpriced inversions + FB-construction natives).
  Negfold `S` is complete.
- Paired ops ratios involving F4 or compact-frob are **diagnostic only**.
- Do not quote `negfold/F4` or `frob/F4` as a speedup without matching
  priced work.
- Secondary pinned-XOR sensitivity: F4 mean S → ~421
  (`tables/f4-sensitivity-claim-n13-v2.txt`).
- Wall figures are L0 / descriptive only.
- Multi-target amortization is secondary analysis only.
- Conductor: `conductor check` allow_with_warning (overlap T-861);
  daemon claim timed out; task T-2514 created.

## Reproduce

```bash
export CARGO_TARGET_DIR=/Volumes/SSD990/wt-ecbench-rootfind-improve/target-parity
cargo build --release --bin ecbench
ECBENCH=$CARGO_TARGET_DIR/release/ecbench bash research/ecbench_rootfind_parity_20261013/run.sh cycle4
ECBENCH=$CARGO_TARGET_DIR/release/ecbench bash research/ecbench_rootfind_parity_20261013/run.sh cycle5
ECBENCH=$CARGO_TARGET_DIR/release/ecbench bash research/ecbench_rootfind_parity_20261013/run.sh claim-v2
```
