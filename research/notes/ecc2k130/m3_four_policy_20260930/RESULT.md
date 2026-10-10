# Degree-7 four-policy m3 result: native selection misses the cold-cost gate

The preregistered degree-7, eight-useful-point, three-summand panel is a
verified **negative for native-base cost advantage** on this toy curve. All
16 same-Q cells reached rank nine and independently recovered `[360]G=Q`,
but descendant-native bases used fewer cold field multiplications than the
transported controls in only one of four paired seed/holdout cells. The
original source policy was cheapest in all four. This says nothing about
ECC2K-130 natural PDP yield or a full attack crossover.

The [original protocol](PROTOCOL.md) was merged in
[PR #1065](https://github.com/aburan28/crypto/pull/1065), before the new
target streams were evaluated. The producer/source lock was merged in
[PR #1066](https://github.com/aburan28/crypto/pull/1066). Its one-shot
[hosted run 36722040881](https://github.com/aburan28/crypto/actions/runs/36722040881)
on `5548b284c1c5a95538f05e23a2a7a19322300ebc` produced all 16 cells
in 8.856 wall seconds, 8.851 process CPU seconds and 26,705,920 bytes
peak RSS. The initial verifier failed **before reading any cell** because an
imported pilot hook replaced its bare curve constructor with a metered one.
That `AttributeError` remains in [the original replay receipt](evidence_run_36722040881/replay.json).
The [repair protocol](REPLAY_REPAIR_PROTOCOL.md) and
[repair source lock](FROZEN_REPLAY.json) were committed before replaying
the same archived result. The repaired independent verifier passed all
8,192 target cases, 16 cell hashes, 228 separate bit-polynomial point checks,
all first witnesses and rank trajectories, recovered scalars, map covariance
and the charged ledger. No producer rerun or target regeneration occurred.

The table uses **cold field multiplications** through first verified rank
nine as its operation unit. Its ratio is to the original source policy in
the same seed/holdout. Each cell also preserved all 512 post-rank outcomes;
the first-rank cost excludes that audit tail. `Eligible` is exact distinct
three-sum support inside the assigned 210-point orbit holdout, computed
without looking at sampled hits. All four policies use eight useful points.

| Seed | Holdout | Policy | Eligible / 210 | Hits / 512 | First rank | Cold field mul | / original | Correctness |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| 2026093001 | A | Original | 42 | 109 | 79 | 139,226 | 1.000 | rank 9; `[360]G=Q` |
| 2026093001 | A | Transported | 42 | 109 | 79 | 416,852 | 2.994 | rank 9; `[360]G=Q` |
| 2026093001 | A | Native | 46 | 100 | 40 | 361,049 | 2.593 | rank 9; `[360]G=Q` |
| 2026093001 | A | Pullback | 46 | 100 | 40 | 339,835 | 2.441 | rank 9; `[360]G=Q` |
| 2026093001 | B | Original | 52 | 126 | 45 | 86,844 | 1.000 | rank 9; `[360]G=Q` |
| 2026093001 | B | Transported | 52 | 126 | 45 | 346,100 | 3.985 | rank 9; `[360]G=Q` |
| 2026093001 | B | Native | 50 | 131 | 36 | 352,931 | 4.064 | rank 9; `[360]G=Q` |
| 2026093001 | B | Pullback | 50 | 131 | 36 | 334,643 | 3.853 | rank 9; `[360]G=Q` |
| 2026093002 | A | Original | 46 | 114 | 33 | 88,330 | 1.000 | rank 9; `[360]G=Q` |
| 2026093002 | A | Transported | 46 | 114 | 33 | 342,658 | 3.879 | rank 9; `[360]G=Q` |
| 2026093002 | A | Native | 44 | 101 | 46 | 356,377 | 4.035 | rank 9; `[360]G=Q` |
| 2026093002 | A | Pullback | 44 | 101 | 46 | 332,215 | 3.761 | rank 9; `[360]G=Q` |
| 2026093002 | B | Original | 46 | 122 | 28 | 84,700 | 1.000 | rank 9; `[360]G=Q` |
| 2026093002 | B | Transported | 46 | 122 | 28 | 335,464 | 3.961 | rank 9; `[360]G=Q` |
| 2026093002 | B | Native | 48 | 118 | 56 | 374,153 | 4.417 | rank 9; `[360]G=Q` |
| 2026093002 | B | Pullback | 48 | 118 | 56 | 346,361 | 4.089 | rank 9; `[360]G=Q` |

The exact full-group support is 92–96 distinct points per base, below the
unordered-multiset capacity `C(10,3)=120` and the uniform-subgroup
one-shot ceiling `120/421`. Within the preregistered disjoint holdouts,
eligible support is 42–52 of 210 points. Each eligible support set's
first-witness base coefficient rows have rank eight, so this panel's bases
have no structural base-row-rank obstruction. More exact native support did
not reliably yield more sampled hits: seed 1/A has 46 versus 42 eligible
points but 100 versus 109 hits. That is why a target-independent exact
support score and disjoint held-out streams are both needed.

The native-versus-transported gate fails three ways. Native cold field
multiplications were lower only at seed 1/A (361,049 versus 416,852);
native had fewer hits in three of four pairs, later first rank in two,
and higher modular-row operation burden in seed 2/A. The charged kernel
search alone cost 236,057 field multiplications; it accounts for about
63–67% of a native cell's cold field-multiplication total. The smaller
native scan-to-rank work in seed 1/A (6,226 versus 12,826 field
multiplications) is a stage observation, not a whole-method gain.
Even subtracting the shared kernel setup as a warm diagnostic leaves no
consistent native advantage over the original source policy.

The next factor-base test should preregister a **target-blind** score built
from exact eligible three-sum support, first-witness row rank/diversity,
and charged base/orbit construction cost. It should compare selected bases
with random quota-matched controls on new seeds and disjoint target streams,
using the same four-policy transport checks. This is a new hypothesis from
the toy panel, not an established improvement. A Koblitz-family m=31
exploratory test and the required m=83 confidence gate would still precede
any transfer claim; the actual degree-263 leaf's changed endomorphism order
does not supply a cheap native orbit action by itself. The review-gated
m10 capacity PR #937 remains a separate prerequisite for its own n=131
dispatch. Solver comparisons (FES, SAT, crossbred, Gröbner/F4/F5) require
the same explicit complete positive/negative PDP corpus and charged
encoder/oracle costs before timing can rank them.

The complete [raw archive and SHA-256 manifest](evidence_run_36722040881/MANIFEST.json),
[repaired replay](evidence_run_36722040881/replay_repaired.json),
[analysis source](analyze.py), and [derived analysis](ANALYSIS.json) are
committed together. From the repository root, independently regenerate
the receipts with:

```sh
python3 research/notes/ecc2k130/m3_four_policy_20260930/verify.py \
  --evidence research/notes/ecc2k130/m3_four_policy_20260930/evidence_run_36722040881 \
  --out /tmp/m3-replay-check.json
cmp /tmp/m3-replay-check.json research/notes/ecc2k130/m3_four_policy_20260930/evidence_run_36722040881/replay_repaired.json
python3 research/notes/ecc2k130/m3_four_policy_20260930/analyze.py \
  --evidence research/notes/ecc2k130/m3_four_policy_20260930/evidence_run_36722040881 \
  --out /tmp/m3-analysis-check.json
cmp /tmp/m3-analysis-check.json research/notes/ecc2k130/m3_four_policy_20260930/ANALYSIS.json
```

This is an **accounting and policy-selection negative** at `GF(2^21)`,
not a measured n=131 PDP yield, full-ECDLP `S`, rho ratio or method
crossover. Those quantities remain null in the machine-readable analysis.
