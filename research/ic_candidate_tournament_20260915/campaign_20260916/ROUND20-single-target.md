# Round 0020 — the replication round

Pre-registered before any round-0020 stage ran. Nothing below reads round 0020.

## 1. Why this round exists

Round 0017 held `beats_rho_strict` on the eight-cell panel. Round 0018b ran the
same executable under a fresh seed and did not. Round 0019 established why —
twelve fixtures cannot resolve `n23a1`, whose per-case log spread is 0.329 —
fixed the instrument with a per-cell allocation, and promoted `both` with
`beats_rho_strict` on all eight cells in both metrics on both final stages.

**That is one seed.** By this campaign's own demonstrated failure mode, a
single-seed strict win is not a decisive one, and round 0019 said so in its own
scope paragraph. This round is the test that round 0017 never got.

Round 0019 also left a specific question. `round19_seed_variance.py` combines
three rounds' per-cell log ratios by inverse variance and finds three of eight
cells scattering beyond within-round sampling — `n13a0` (p = 0.043), `n23a1`
(0.044) and `n31a0` (0.0001). The test over-rejects there, because rho's cost
is right-skewed and two of the three rounds carry only twelve fixtures a cell,
so those flags were recorded as questions rather than findings. **A fourth seed
at round 0019's allocation is the cheapest thing that answers them.**

## 2. The round

Identical to round 0019 in every respect but the seed. That is the point: a
replication that changes anything else is not a replication.

| | |
|:--|:--|
| seed | `2026092120` (rounds 0017, 0018b, 0019 used 2026091717, 2026091818, 2026092119) |
| panel | `--cells 13a0,17a1,19a0,23a0,23a1,31a0 --holdout-cells 19a1,29a1` |
| allocation | `--confirmation-cases n23a1=40` — 124 cases, unchanged |
| profile | pilot, one target, objective `rho`, `--require-native-progress` |
| arms | `incumbent`, `both` — unchanged |
| baseline | `runs/round-0019/source` (the round-0019 incumbent, worker `e9f263b8…`) |
| trials | 2,646 |

The baseline stays the **incumbent**, not round 0019's promoted winner. A
replication must re-measure the same comparison, and promoting `both` into the
baseline would make `both`/incumbent identically 1 and destroy the thing being
replicated. Round 0019's promotion stands on its own record; this round asks
whether that record reproduces.

## 3. Predictions

1. **`beats_rho_strict` holds again.** All eight cells below one in both
   metrics on both final stages. This is the whole point of the round; a
   failure here is the most informative outcome it can produce and would mean
   round 0019's win was seed-luck after all.
2. **`both`/rho lands near round 0019's value.** Panel instruction ratio in
   [0.62, 0.74] (round 0019: 0.6751), native in [0.85, 0.92] (0.8848).
3. **`n23a1` clears with margin.** `both`/rho at `n23a1` below 0.95 in
   instructions. Round 0019 read 0.7498 on forty fixtures; the cross-seed
   combined value for that cell is about 0.78. Predicting 0.95 rather than 0.78
   is deliberate: prediction 2 of round 0019 was falsified by predicting a
   point value from a pooled two-seed mean, and this round does not repeat that
   mistake. The claim under test is that the cell clears, not where it lands.
4. **`both`/incumbent reproduces.** Instruction ratio in [0.95, 0.98] with both
   upper limits below one on both final stages. Round 0018b read 0.9652 and
   round 0019 read 0.9657, so this is the most stable quantity the campaign
   measures and a miss would indicate something wrong with the build or the
   harness rather than with the fixture draw.
5. **The three flagged cells are tested, not assumed.** After this round,
   `round19_seed_variance.py` over four seeds either keeps `n13a0`, `n23a1` and
   `n31a0` above p = 0.05 or does not. **No outcome here is pre-registered as
   the expected one** — that is what it means for round 0019 to have recorded
   them as open questions. The round records the four-seed chi-square for each
   cell and states which way it moved.
6. **The incumbent still fails somewhere.** Round 0019's incumbent missed the
   strict gate at `n29a1` natively (1.0007). Predicting it fails *somewhere*
   again is weak; predicting the cell is not. `n29a1` is where the cofactor is
   12,646 and where `column` does the most work, so: the incumbent's worst
   native cell is `n29a1`, whether or not it exceeds one.

## 4. Falsification and scope

Any bad, missing or unmatched certificate; any arm returning a different
logarithm or factor-base hash from the incumbent on the same fixture; a `both`
report failing to declare `column_convention: representative`; the A/A control
outside [0.95, 1.05] at any cell; a prepare whose contract differs from §2 in
seed, panel, allocation, profile, objective, arms or limits.

What two seeds would and would not establish: two independent strict wins at
this allocation make the result robust to the fixture draw **on this panel, at
these subgroup orders, for this executable**. They say nothing about larger
`r`, where round 0019's boundary analysis puts this collector at `Θ(r^{2/3})`
against rho's `Θ(r^{1/2})` and its same-degree contrast measures the ratio
rising 16.6% per doubling at degree 23. The classification stays
**engineering** under AGENTS.md §3 whichever way this round goes.
