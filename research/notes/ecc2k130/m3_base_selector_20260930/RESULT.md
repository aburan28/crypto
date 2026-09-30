# Charged degree-7 m3 base selector: verified negative

**Decision: `NO_CHARGED_SELECTOR_ADVANTAGE` under the frozen two-candidate
protocol.** The independent verifier passed all 64 cells, all 32,768 case
records and all recovered `[260]G=Q` scalars. Every candidate's first-witness
base rows had rank eight, and source/transported and native/pullback witness
and rank trajectories agreed. Yet the selected arm used fewer cold field
multiplications in only **2/8 source** and **0/8 native-leaf** seed/holdout
pairs. Source selection had two hit regressions and three later rank-nine
solves; native selection had two hit regressions and two later solves. The
predeclared all-pairs gate therefore fails. This is a negative for this exact
selector at `GF(2^21)`, not for factor-base selection in general or for
ECC2K-130 index calculus.

The [protocol](PROTOCOL.md) was merged in
[PR #1086](https://github.com/aburan28/crypto/pull/1086) before new inputs
were generated. The [27-file source lock](FROZEN.json) and its read-only CI
were merged in [PR #1089](https://github.com/aburan28/crypto/pull/1089), at
`991b904720d0b08b6cbc171f195a433930cd4f5e`. A first local command
stopped in `checked_config()` because the sparse checkout omitted a pinned
workflow file. Its [preflight record](evidence_run_20260930_2125Z/preflight_failure.json)
shows that no candidate, Q or target was generated. After materializing that
tracked file, the sole completed producer run generated the new challenge
`Q=(1073582,1159790)` with audit-only scalar 260 and all frozen cells in
11.324 wall seconds, 11.089 process CPU seconds and 35,176,448 bytes peak
RSS. The [raw manifest](evidence_run_20260930_2127Z/EVIDENCE.json) pins
`result.json`; the result pins every one of the 64 cell JSON files. The
[independent replay](evidence_run_20260930_2127Z/replay.json) regenerated
both candidates per seed and role, exact support and scores, new target
streams, witnesses, ranks, scalar logs and charged ledgers. It checked
selected full points and target edges with a separate bit-polynomial law.

The table is the **cold field-multiplication cost to first verified rank
nine**. Each control uses candidate 0 and pays for that one scan. Each
selected arm pays for both scans and both 324-addition exact-support scores,
even when it chooses candidate 0. All policies have eight useful points,
two from each of the same four source signed-Frobenius orbits. Hit counts
include all 512 accepted targets; later audit targets are excluded from the
cold rank cost. Ratios compare selected/control on the same Q and role.

| Seed | Holdout | Role | Chosen candidate | Control mul | Selected mul | Ratio | Hits control/selected | First rank control/selected |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2026100101 | A | source | 0 | 86,025 | 114,363 | 1.329 | 125/125 | 46/46 |
| 2026100101 | A | leaf | 1 | 346,410 | 439,881 | 1.270 | 98/109 | 45/71 |
| 2026100101 | B | source | 0 | 186,873 | 215,211 | 1.152 | 96/96 | 127/127 |
| 2026100101 | B | leaf | 1 | 342,340 | 368,227 | 1.076 | 125/95 | 43/38 |
| 2026100102 | A | source | 1 | 176,776 | 147,695 | 0.835 | 79/98 | 81/36 |
| 2026100102 | A | leaf | 1 | 327,897 | 364,054 | 1.110 | 118/111 | 41/44 |
| 2026100102 | B | source | 1 | 156,272 | 149,697 | 0.958 | 89/102 | 70/40 |
| 2026100102 | B | leaf | 1 | 345,035 | 354,484 | 1.027 | 80/129 | 51/40 |
| 2026100103 | A | source | 1 | 82,504 | 142,021 | 1.721 | 100/92 | 45/48 |
| 2026100103 | A | leaf | 0 | 427,999 | 475,783 | 1.112 | 122/122 | 85/85 |
| 2026100103 | B | source | 1 | 67,764 | 152,845 | 2.256 | 140/103 | 36/62 |
| 2026100103 | B | leaf | 0 | 351,175 | 398,959 | 1.136 | 104/104 | 46/46 |
| 2026100104 | A | source | 1 | 57,933 | 92,526 | 1.597 | 123/140 | 28/30 |
| 2026100104 | A | leaf | 0 | 404,920 | 447,865 | 1.106 | 106/106 | 56/56 |
| 2026100104 | B | source | 1 | 104,485 | 126,098 | 1.207 | 80/107 | 67/59 |
| 2026100104 | B | leaf | 0 | 374,142 | 417,087 | 1.115 | 111/111 | 43/43 |

Candidate scores on the whole 420-point nonzero subgroup ranged from 64 to
96 distinct three-sums, each with base-row rank eight. The selected arm's
extra second scan plus two exact scores cost **28,338–51,245 source** and
**31,647–47,784 native-leaf** field multiplications per seed before any
online improvement. When the selector chose candidate 0, this overhead
raised cold cost without changing hits or rank. Seed 2026100102 shows a
real, narrow counterexample to a blanket dismissal: choosing source
candidate 1 saved 29,081 and 6,575 cold multiplications in A and B,
respectively, after paying that overhead. The other source seeds and every
native-leaf seed failed the gate. Native-leaf cold cost also exceeded the
original source policy in all 16 control/selected holdout pairs; this
continues the earlier toy accounting result, not a degree-263 result.

The [analysis](ANALYSIS.json) includes all 16 primary pair comparisons,
selector choices, operation counts, rank/hit regressions and a **post hoc**
support breakdown across the two 210-point holdouts. This breakdown was not
the frozen selector's score and cannot be promoted as a tested alternative.
For example, source seed 2026100103 increased all-group support from 90 to
92, but moved A/B support from 38/52 to 46/46 and yielded fewer observed
hits in both held-out streams. Exact global support and base-row rank alone
did not predict the charged first full-rank cost on this stream.

The next testable lead is to obtain multiple quota-matched eight-point bases
from **one modest superset scan** and reuse a single table of its triple
sums when scoring subsets. This might avoid paying two independent cofactor
scans and two repeated exact-support enumerations; it still must charge the
larger scan, the superset table, subset ranking, row operations and memory.
Preregister a cost-aware, target-blind selection rule, fresh seeds and
held-out Q/target labels before opening that experiment. A positive toy
result would justify a Koblitz-family m31 pilot, then the required m83
confidence gate. Neither this panel nor the proposed direction measures
natural n=131 PDP yield, a full ECDLP cost, calibrated `S`, or a matched-rho
crossover; those fields remain null.

Reproduce the deterministic replay and analysis from the repository root:

```sh
python3 research/notes/ecc2k130/m3_base_selector_20260930/verify.py \
  --evidence research/notes/ecc2k130/m3_base_selector_20260930/evidence_run_20260930_2127Z \
  --manifest research/notes/ecc2k130/m3_base_selector_20260930/evidence_run_20260930_2127Z/EVIDENCE.json \
  --out /tmp/m3-selector-replay-check.json
cmp /tmp/m3-selector-replay-check.json \
  research/notes/ecc2k130/m3_base_selector_20260930/evidence_run_20260930_2127Z/replay.json
python3 research/notes/ecc2k130/m3_base_selector_20260930/analyze.py \
  --evidence research/notes/ecc2k130/m3_base_selector_20260930/evidence_run_20260930_2127Z \
  --out /tmp/m3-selector-analysis-check.json
cmp /tmp/m3-selector-analysis-check.json \
  research/notes/ecc2k130/m3_base_selector_20260930/ANALYSIS.json
```
