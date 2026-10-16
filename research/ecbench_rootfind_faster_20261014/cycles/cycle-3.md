# Cycle 3 — negfold m=2 at dimension=4

**Change:** smaller standard-subspace `|F|` (d=4 vs d=5) on complete
`rootfind-negfold`, aiming to cut PairTable set-up. Not a repeat of
window@d4 / window_div.

Session: `sessions/cycle3` (`ECBS1h627cfb7e3fb1`), 30/30 verified, audit OK.

## Mean cold S

| arm | mean S | lower bound | med wall ms |
|---|---:|---|---:|
| rho-strong | 14.536 | yes | 0.085 |
| ic-c2-frob-counted-m2 | 15.837 | yes | 0.161 |
| **ic-negfold-m2-d4** | **25.680** | no | 1.033 |
| ic-negfold-m2-d5 | 28.652 | no | 0.326 |
| ic-f4-m3-d5 | 446.754 | no | 205.950 |

## Outcome

**Modest / keep on claim:** d=4 mean S below d=5 on this smoke (CI
overlaps); wall favours d=5. Both complete arms remain ≪ XOR-priced F4.
Carry both to claim-n13-v3. Not a dead end; not a clear cold win.
