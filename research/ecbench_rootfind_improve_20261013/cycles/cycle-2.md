# Cycle 2 — window at d=4 (dead end) + claim panel kickoff

Session: `sessions/window-d4` (`ECBS1hbac35523d128`), 30/30 verified, audit OK.

## Mean cold S

| arm | mean S | solve wall (ms, L0) |
|---|---:|---:|
| rho-strong | 14.6 | 0.13 |
| ic-f4-m3-d5 | 16.0 | 192 |
| **ic-rootfind-m2-d4** | **27.8** | 0.72 |
| ic-rootfind-window-m3-d5 | 33.0 | 0.91 |
| ic-rootfind-window-m3-d4 | 35.5 | 1.15 |

## Decision

`rootfind-window` at `dimension=4` is a **dead end** for cold S: the smaller
table does not offset the thinner m=3 yield. Keep cycle-1 selection:

- cold candidate: `pdp3-koblitz:m=2,engine=rootfind` + `dimension=4`
- m=3 candidate: `engine=rootfind-window` at `dimension=5`
- honest F4 baseline: `engine=inherited-f4` m=3 d=5 (wall favourite for
  showing the constant-factor engine win; S still a lower bound)

## Next

Run `claim-koblitz-rootfind-improve-n13.json` (8 public targets × 3 rounds),
write before/after vs #1765, amortized secondary panel, open PR.
