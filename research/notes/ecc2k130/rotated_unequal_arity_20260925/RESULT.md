# Exact n131 unequal-slot census: m10 alone passes the frozen two-axis gate

The [protocol](PROTOCOL.md), four allocations, source and inputs were frozen at
`c798909` in draft [PR #784](https://github.com/aburan28/crypto/pull/784), with
hash-only CI green before the first outcome. `FROZEN.json` SHA-256 is
`ba93ce394df2c57f498f62bd59b74ed1faac8b748119d3e15594659d197f09cf`.
The first and only four-arm run completed successfully. Its independent verifier
rebuilt every one of 1,056,768 low/high masks with separate bit-serial field
arithmetic, checked all row and chunk hashes and exact signed-column sets, and
passed 64 fixed sample ordinals in each of eight spaces. Of those 512 samples,
269 had rational lifts and passed a full curve-equation and two-doubling check.
Every arm has a normal-basis rank and combined slot-coordinate rank of 131;
the low normalized column set is a subset of the high one. This is **field
coordinate rank**, not relation-matrix rank.

| Arm (high slots first) | F low | F high | Nonzero signed C low | Global nonzero signed C high | Exact N/q | Raw affine-chain variables | Frozen `N≥q, C≤16,384` |
|:--|--:|--:|--:|--:|--:|--:|:--|
| m7: 19×5, 18×2 | 261,709 | 523,721 | 130,854 | 261,860 | 3.965224835223 | 786 | fail C |
| m8: 17×3, 16×5 | 65,747 | 131,331 | 32,873 | 65,665 | 4.088952050082 | 917 | fail C |
| m9: 15×5, 14×4 | 16,743 | 33,247 | 8,371 | 16,623 | 4.690573928758 | 1,048 | fail C by 239 |
| m10: 14×1, 13×9 | 7,977 | 16,125 | 3,988 | 8,062 | 3.098750509013 | 1,179 | **pass** |

Here `F=1+2L` is the physical F0 point count in one normalized space,
including `(0,1)`, and `N=F_high^r F_low^(m−r)` is the exact ordered physical
tuple count. The nonzero signed projected columns are unique within each
space in this census; the identity adds one further column. Every low/high
space had exactly one `x=0` mask and **zero `x=1` masks**. The low F and C
values agree exactly with the corresponding balanced #778 arms. The high
slots occupy i=0..r−1, so the slot conjugate indices partition 0..130; after
inverse Frobenius the low normalized space is nested in the high space.

All four tuple counts exceed `q`, but `N/q` is only a necessary uniform-target
count ceiling. It clips at 100% and is neither an observed PDP hit rate nor a
solver yield. The m9 threshold miss is **239 columns, or 1.4587% of the frozen
cap**; it does not prove m9 mathematically impossible. Under this preregistered
gate the smallest admitted unequal arm is m10. A later cost or Pareto study
could reconsider m9 only under a new frozen protocol and criterion.

Against [#778](https://github.com/aburan28/crypto/pull/778)'s balanced
m10,d13, unequal m10 has 3.098750509013 versus 1.532944670412 ordered
tuples per `q`, a 2.02144× count increase. Its normalized signed-column count
is 8,062 versus 3,988, a 2.02156× increase, and its raw affine-chain floor
is 1,179 versus 1,178 variables. Thus the count gain brings an almost
proportional column burden. Neither arm dominates on these two axes, and the
gate pass selects a **feasibility arm**, not a winning index-calculus attack.
There is no branch-complete implicit m10 PDP exporter, measured target
support, relation yield or independent relation rank for these n131 spaces;
SAT/Gröbner cost, full ECDLP cost `S`, and matched rho ratio remain unset.

The complete first run began 2026-09-25 14:36:44 UTC and ended 14:54:05 UTC.
Producer/independent-verifier wall times were 76.907/706.144, 15.986/168.856,
3.920/44.765 and 2.136/21.232 seconds for m7–m10. All eight children exited
zero under the frozen 600/1,200-second caps. The largest recorded child RSS
upper bound was 79,659,008 bytes, below the 256-MiB cap. These are census and
replay times, **not** PDP solver or attack times. The [receipt](evidence/receipt.json)
(SHA-256 `f22041008f4cb3fef622f7381150dc198fe085758aacfb17e6e0d011c5869e47`)
preserves all commands, UTC intervals, raw results, hashes and stream files;
[EVIDENCE.md](EVIDENCE.md) explains the independent replay.
