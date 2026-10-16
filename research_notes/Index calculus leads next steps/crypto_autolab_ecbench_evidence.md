# Empirical index-calculus evidence in aburan28/crypto: autolab, ecbench, benchmarks, experiments

Scope: local repository `/home/user/crypto` only, read on 2026-10-14. Sources are given as repository-relative paths. Every number below is quoted from a committed record. Arithmetic I did myself is labelled as an inference. Abbreviations: S = counted cost / sqrt(r), where the counted cost is usually ecbench GAE (group-addition equivalents) and sometimes instructions (Ir); r = prime subgroup order; "cold" = one target with all setup charged; "online" = one target after reusable setup is excluded.

## 1. What the IC autolab does, what it tried, and whether any variant beats rho on a complete solve

### Takeaway
The autolab (`research/ic_candidate_tournament_20260915/`) is a heavily gated tournament. It treats the factor base, the decomposition solver (pair_table, enumerate, F4, F5, inherited F4, SAT-XOR, SAT-CNF) and dense or sparse linear algebra as interchangeable parts, and compares them against a qualified incumbent and matched rho. Its three registered improvement rounds (10,203/10,203 pairs verified) found no challenger that met the 20% gate, so the `pairinv` three-summand pair-table incumbent was kept. The generic F4/F5/SAT backends never qualified: two campaigns were operationally censored and v3 ended in NEGATIVE_FAMILY_QUALIFICATION. On its own toy Koblitz panel (n17–n29, r about 2^18), the incumbent IC is cheaper than rho in cold instructions and about 7.4x faster online, but statistically tied with rho in cold native time. No result anywhere in the repository shows IC beating rho on a cold, setup-inclusive complete solve. IC wins exist only on the online boundary, where retained setup is excluded.

### Cited Findings
**Architecture and components**
- Pipeline: public curve and targets → factor-base recipe and orbit closure → PDP encoding and solver → relation verification → rank and linear algebra modulo the subgroup order → certified descent → independent checks, with matched optimized rho controls. — [AUTOLAB.md](research/ic_candidate_tournament_20260915/AUTOLAB.md)
- Interchangeable parts:
  - factor bases: factor, divisor, Frobenius union, subgroup orbits, pruned, torsion-saturated;
  - solvers: `pair_table`, `enumerate`, `f4`, `f5`, `inherited_f4`, `sat_xor`, `sat_cnf`;
  - linear algebra: `dense`, or `sparse` with filtering and Wiedemann;
  - other knobs: batch trials and collection window; rho at 1, 8 or 32 walks.
  External SAT solvers, Sage/Magma F4/F5 and GPU kernels are "not silently admitted"; the protocol is CPU-only. — [AUTOLAB.md](research/ic_candidate_tournament_20260915/AUTOLAB.md)
- The skill says native Rust is now mandatory. The Python tournament scripts are historical, and new cross-method work must use native ecbench. — [.claude/skills/ic-autolab/SKILL.md](.claude/skills/ic-autolab/SKILL.md)

**The three bounded improvement rounds**
- Results by round (gate: at least 20% lower cold Ir and cold native time, i.e. ratio ≤ 0.8):

  | round | selected challenger | cold-Ir ratio | cold-time ratio |
  |---|---|---|---|
  | 1 | `stop6` | 0.957 | 0.965 |
  | 2 | `stop5_word` | 0.923 | 0.987 |
  | 3 | `stop7_word` | 0.990 | 0.986 |

  In round 2 the best cold-Ir improvement came with an online regression. The incumbent was retained, and the note calls this "an audited negative result on the registered synthetic toy panel". — [BOUNDED-GOAL-RESULT.md](research/ic_candidate_tournament_20260915/goal_20260924/BOUNDED-GOAL-RESULT.md)
- Scope of the rounds: they exercised "optimized three-summand pair-table point decomposition, Frobenius-orbit factor bases, relation collection, scalar-field row kernels and direct target descent". The F4/F5/SAT and sparse relation-LA adapters exist, but these rounds did not qualify them competitively. — [BOUNDED-GOAL-RESULT.md](research/ic_candidate_tournament_20260915/goal_20260924/BOUNDED-GOAL-RESULT.md)

**Round 3 confirmation panel (216 rows per arm)**
- Cells: Koblitz n17a1, n19a0, n23a0, n23a1, n29a1.
- Incumbent base: 272–464 usable points folded to only 8 columns (rank 8). `stop7_word` uses 7 columns. — [round3/README.md](research/ic_candidate_tournament_20260915/goal_20260924/improvement/round3/README.md)
- Online: incumbent 0.0339 ms against rho_online 0.2526 ms, so rho/IC = 7.45x.
- Cold Ir: incumbent 1.80e6 (S = 2,822) against rho 2.50e6 (S = 3,912), so IC/cold rho = 0.721.
- Cold native time: incumbent 0.888 ms against rho 0.860 ms, so rho/IC = 0.968 [0.915, 1.027]. — [round3/README.md](research/ic_candidate_tournament_20260915/goal_20260924/improvement/round3/README.md)

**Generic-backend campaigns**
- Campaign 1 (seed 2026092901): hit the six-hour job cap and the artifact upload failed, so it is censored.
- Campaign 2 (seed 2026092902): measure step timed out at 300 minutes and zero bundle was uploaded, so it is censored. A static audit found that all four F4/F5 arms exceed the 64-variable encoder cap on every cell.
- Dimension-6 F4/F5 pilot: zero accepted relations. — [ic-autolab SKILL.md](.claude/skills/ic-autolab/SKILL.md)
- Campaign 3 (Actions run 37399021219) returned `NEGATIVE_FAMILY_QUALIFICATION`: 45 trial slots, 34 verified, 11 failures, `qualified_f4_f5_arms` empty, SAT out of scope. — [generic-backend-qualification-v3/RESULT.md](research/ic_candidate_tournament_20260915/goal_20260924/generic-backend-qualification-v3/RESULT.md)

**Native one-target controls (n17)**
- F5 v2: rank 29/29, 61 verified relations in 216 ordinary queries, scalar 24886 after three target attempts, 10.5 s online on an uncalibrated Mac. It had no rho arm. — [ic-autolab SKILL.md](.claude/skills/ic-autolab/SKILL.md)
- Prepared one-target F5 and SAT controls: both closed at eight attempts with no scalar. The SAT attempts were conflict-budget inconclusive. — [ic-autolab SKILL.md](.claude/skills/ic-autolab/SKILL.md)
- STATUS (Oct 3): "zero new ordinary queries, no fresh paired performance qualification and no measured speedup". — [STATUS.md](research/ic_candidate_tournament_20260915/goal_20260924/STATUS.md)

**Separate IC1 "autolab" series on a 256-bit prime curve**
- About 24 round directories `research/icv1-fp256-…_round1…round24_20261004–06` exist (Dickson, S3/S5/S9/S17 image, factor-base symmetry, root-orbit FB, negation cross-colour…). I did not read them in depth. — [research/ directory listing](research/)

### Inferences
- The autolab's IC "beats rho in cold Ir" result rests on the tiny 8-column regime (r about 2^18.6, my inference from S = Ir/sqrt(r)). The tie in cold native time shows that instruction counts and native time disagree at this scale. This is not evidence of a scalable win.
- The autolab's main output so far is infrastructure and negative qualification results, not a better decomposition method.

### Gaps
- I did not inventory each of the 16 candidates per round or the contents of the fp256 IC1 rounds 1–24.
- No autolab run is at a cryptographic size.

## 2. What the ecbench bounds and frontier say about IC's constant and exponent, and which arms, sizes or curves are missing

### Takeaway
There is no admissible sealed IC bound. The three `ic.pipeline` bounds are inadmissible (48/96 runs unverified, only 2 sizes). The canonical cold frontier is rho with signed Frobenius (Koblitz medium tier: 1.089x floor, α = 0.451 over 144 runs at log₂ r 32.1–47.2) and BSGS-negation or rho-negation for prime toy curves. Every complete IC arm measured in ecbench costs more than matched rho:
- binary ladder: 2.9x to 25,116x;
- prime PKM tower: 29–141x;
- Koblitz one-target: 5–21x;
- fiber IC: 8–63x, growing like #E^{1/6}.

IC beats rho only online, and then by a lot: n61 13.9x admitted; N83 18.5x exploratory; N83 54.1x on 5 fresh points.

### Cited Findings
**Sealed bounds and frontier**
- The frontier is built from 46 bounds. No IC method is on any frontier.
- The inadmissible list holds three `ic.pipeline` bounds "on 2 koblitz curves, toy tier: 48 of 96 measured runs did not verify… 2 size(s)". — [docs/bounds/FRONTIER.md](docs/bounds/FRONTIER.md)
- Koblitz medium-tier ops leader: `rho.signed_frobenius_strong` at 1.089x floor [0.927, 1.261], α 0.451 [0.399, 0.512], over sizes 32.1–47.2. `claw.pair_table` sits at 2.390x. — [docs/bounds/FRONTIER.md](docs/bounds/FRONTIER.md)
- Prime toy leaders: `bsgs.negation` at 1.100x floor and `rho.negation` at 1.453x (15.7–25.3 bits). The prime challenge is described as "the question every rho variant and every index-calculus candidate on prime curves is asked first". — [FRONTIER.md](docs/bounds/FRONTIER.md); [prime-toy-reference.json](docs/bounds/challenges/prime-toy-reference.json)
- The only two challenge verdicts (epoch 1) are not IC ("Nothing here is index calculus"):
  - bsgs.negation against rho.negation: trade;
  - rho.signed_frobenius against strong: advances at the constant level, 0.299×. — [ecbench_bounds_challenges_20261006](research/ecbench_bounds_challenges_20261006/README.md)

**Frontier synthesis (2026-10-13)**
- Decision: "The cold, one-target, setup-inclusive frontier is still generic rho. None of the measured index-calculus candidates beats its matched rho arm." — [ecbench_frontier_20261013/README.md](research/ecbench_frontier_20261013/README.md)
- Factor-base arms (mean S, ratio to matched rho):

  | curve | factor base | S | S / rho |
  |---|---|---|---|
  | Koblitz n17 | compact orbit (272 points / 8 columns) | 9.114 | 9.985x |
  | Koblitz n17 | full orbit | 17.081 | 18.7x |
  | Koblitz n17 | binary or standard subspace, dim 4 | – | about 47x |
  | Koblitz n23 a=0 | compact orbit | 2.996 | 8.36x |
  | prime 20-bit | abscissa (64 points / 32 columns) | – | 12.3x |
  | prime 20-bit | GLV orbit | – | 24.3x |

  The native catalog holds ten sealed bases; only seven were exercised. — [same](research/ecbench_frontier_20261013/README.md)
- Online boundary results:
  - n=61 compact-orbit, K=500: rho/IC 13.9x [9.8, 19.1] over 64 unseen targets, about 28 s setup excluded, cold IC/rho median 83.4.
  - N83 dense K1798: 18.53x on one point (63.3 ms against 1,173 ms); complete IC process 339 s; derived what-if crossover q* = 306 targets. — [same](research/ecbench_frontier_20261013/README.md); [docs/ic/LEADERBOARD.md](docs/ic/LEADERBOARD.md)
- Replay audit: the current-binary root-finder replay-all matched 120/168. All 48 inherited-F4 rows differ in solver word-XOR counts across hosts, a failed audit that was preserved. — [same](research/ecbench_frontier_20261013/README.md)

**Newest N83 result (2026-10-14, not yet in the frontier synthesis)**
- K3240 root21 (537,840 points, 3,240 columns), Tier-1 GitHub Actions run on 5 fresh points:
  - median rho/IC = 54.10 [17.18, 58.90];
  - means: IC 58.9 ms, rho 2,597 ms;
  - IC setup 580–640 s per point, peak RSS about 13.9 GiB.

  Setup is excluded from the ratio. — [n83_tier1_final_results_20261014/RESULT.md](research/n83_tier1_final_results_20261014/RESULT.md)
- The 100-point selection panel chose K2940 root22 (22.548 ms mean online). The runner-up difference [-5.68, 6.00] ms is not resolved. A reserve five-point fixture "is the next same-point IC/rho validation for K2940". — [n83_factor_base_100_results_20261014/RESULT.md](research/n83_factor_base_100_results_20261014/RESULT.md)

**Koblitz one-target (§23, leaderboard Table A)**
- Six curves, log₂ r 36.6–47.2, 64 points each.
- Online speedup: 8.84–17.5x. Cold: rho faster at all six (cold ratio 1.98–31.6x).
- At baseline v3, IC costs 5.23x to 21x the reference, and 73.5–220x the floor.
- Reusable setup is 98.5–99.7% of IC S. — [ic_single_target_20260930/README.md](research/ic_single_target_20260930/README.md); [docs/ic/LEADERBOARD.md](docs/ic/LEADERBOARD.md)

**Exponent and constant fits**
- Ledger §20, Koblitz m=3 against batch rho at k=32, nine sizes from 2^18 to 2^47.2:
  - "No crossing"; best 4.86x [4.21, 5.55] at 2^39;
  - fitted β = 0.137 [−0.10, 0.37], so the exponent is undecided;
  - the r^{1/6} law was falsified at the small end (12.7–24.8x measured against 0.8–1.2x predicted). — [ic_exponent_20260926/README.md](research/ic_exponent_20260926/README.md)
- Fiber-closed MITM IC against rho.negation: envelope ratio grows from 8.4 to 37.0 (cofactor 4) and from 13.4 to 63.1 (cofactor 72) between 18 and 30 bits. Slopes are 0.177 [0.145, 0.209] and 0.168, i.e. #E^{1/6}, a negative result. — [fiber_ic_scaling_20261011/README.md](research/fiber_ic_scaling_20261011/README.md)
- PKM Kummer-tower prime IC: 29–141x rho-negation on five curves of 20–28 bits. Slope 0.26 ± 0.28 bits per bit, against 0.17 predicted. — [pkm_tower_ecbench_20261013/README.md](research/pkm_tower_ecbench_20261013/README.md)
- Binary E_0 ladder: S_IC/S_rho = 2.9, 15, 61, 3,053 and 25,116 at n = 13, 19, 23, 31, 37. The jump at n31 is mostly the 1,492 cofactor. — [pdp_root_finder_ecbench_20261013/README.md](research/pdp_root_finder_ecbench_20261013/README.md)

**Leaderboard Tables B and C**
- Table B (prime ladder 9.7–23.4 bits): 5.5–13.5x reference. Random-binary ladder: 10.6–38.2x.
- Table C (Koblitz, 32-target batch): 3.62x to 14.9x. — [docs/ic/LEADERBOARD.md](docs/ic/LEADERBOARD.md)

### Inferences
- Measured exponent evidence for linear-scan IC on prime and binary curves points the wrong way, with the ratio to rho growing at about #E^{1/6}. The only regime with theoretical sub-rho scaling, E(F_{p^k}) with k ≥ 3, has not been put into ecbench (see §5).
- The online-boundary wins (13.9x to 54x) are precomputation trades. Using the reported means (N83: about 600 s setup, about 2.54 s saved per target), break-even would be roughly 236 targets. This is my arithmetic, not a measured batch frontier.

### Gaps / missing arms
- IC with four or more sizes and all runs verified (needed for an admissible exponent) does not exist.
- Missing arms:
  - torsion-symmetrised bases have no ecbench plug-in;
  - the N83 dense-root engine is not registered in `ecbench methods`;
  - no ecbench IC on E(F_{p^4});
  - no kangaroo or BSGS comparison against IC.
  — [ecbench_frontier_20261013/README.md](research/ecbench_frontier_20261013/README.md) "What would move the frontier"
- The m=83 IC gate in ecbench is only a plan: [ic_wide_shared_rank_plan_20261005/PLAN.md](research/ic_wide_shared_rank_plan_20261005/PLAN.md) says "no code, no measurement"; [ecbench_m83_ic_20261009](research/ecbench_m83_ic_20261009/PROTOCOL.md) is PROTOCOL only.
- The F6 ladder at m=31 hit 3,600 s on the first target and stopped without a decision. Matrix F4/F5 at m=19 timed out (1,800 s). — [ecbench_frontier_20261013/README.md](research/ecbench_frontier_20261013/README.md)

## 3. benchmarks/README.md: rows on IC, Semaev, root finding, Gröbner, GPU and sparse LA, and the recorded speedups

### Takeaway
The index has only 16 rows, but `benchmarks/` holds about 50 directories. Many IC-relevant directories have no index row; they are listed below with their own headlines. Every recorded speedup is a stage or engineering gain: root finding against F4, F5 against F4, D₄ symmetry on C₄, solve-speed and eager rewrite. None is a cold IC-over-rho gain. The one sub-rho number recorded is linear algebra alone at k=4.

### Cited Findings
Index rows, all from [benchmarks/README.md](benchmarks/README.md):

**Root finding in ecbench, end to end**
- `20261013-ic-rootfind-ecbench`: whole IC solve with root finder against F4 is 350x faster in wall at n=13 (admitted) and about 10,000x at n=23 (descriptive). IC is still 2.9–25,116x the strong rho (binary) and 29–141x rho (prime).
- `20261014-rootfinding-lanes-coverage`: 1,088/1,088 verified. The B48 train selection did not resolve a holdout winner (B56/B48 0.979 [0.919, 1.038]). Field operations are missing from sealed records, so costs are lower bounds.
- `20261013-prime-rootfinding-batchtable`: batch against scalar inversions gives the same cold GAE (ratio 1.000). Batch/rho = 16.229.
- `20261013-prime-rootfinding-arity`: m3/m2 = 1.431 [1.225, 1.716]. m2/rho = 20.03.
- `20261013-rootfinding-pair-size`: B64/B32 = 0.830 [0.697, 0.962]. B64/rho = 18.87.
- `20261012-rootfinding-decomposition`: pair-table/subtraction S = 0.0189 (prime) and 0.057 (Koblitz). Matched rho is smaller on both curves.

**Root-finding oracles as a stage**
- `20261012-pdp-root-finder`:
  - `rf-mitm` against PS-0 is 3,007x to 115,160x faster at ℓ = 8–14;
  - prime m=3: 197 multiplications against F4's 8.99·10¹³ multiply-adds;
  - "every slope at its generic floor — engineering, no exponent moved".

**Gaudry k=4 (C₄)**
- `20261013-gaudry-quartic-d4`: C₄ falls from 1.21·10¹² to 2.54·10¹⁰ F_p multiplications (47.5–48.4x); a solve takes 5 s instead of 3 min. The model crossover with rho moves from 2^150 to 2^135 (h = 2).
- `20261013-gaudry-quartic-solve-speed`: 1.9x / 4.1x / 6.4x faster at d = 4/5/6; C₄ is unchanged.

**Frontier, N83 and conductor gap**
- `20261013-ecbench-frontier`: synthesis (see §2).
- `20261013-n83-dense-root`: K1798 dense against K1500 gives a paired probe difference of 130,952 [52,058, 209,991]. Online 18.53x exploratory.
- `20261012-conductor-gap-counts`: IC floor/crater ratios 0.975 and 1.127. "no full-cost or wall advantage".

**Gröbner**
- `20261011-groebner-vs-msolve-singular-polybori`:
  - F_p: msolve (1 thread) is 12–24x faster than the repository's `f4_fp` on square systems;
  - Boolean n=17: the repository's oracle decides a target in 7–24 ms, against 40–52 s for msolve and PolyBoRi to compute a full basis.
- `20261011-inherited-f4-eager-rewrite`: lazy/eager 1.40–1.62x.

**GPU**
- `20261013-nist-gpu-rho-ladder`: rho on P-256/P-384 only, 1.632 and 0.548 billion iterations/s. Not IC.

Directories without an index row (headline from each directory's README):
- `20261011-fiber-ic-scaling`: envelope IC/rho grows from 8.4–8.8 to 35.1–37.0 (h=4) and from 13.4–20.1 to 60.4–63.1 (h=72); negative. — [README](benchmarks/20261011-fiber-ic-scaling/README.md)
- `20261005-ic-semiregular-law`: registered verdict "neither". Heavy degree-10 solves of 239k–247k rows × 199k columns took 42–51 min. — [README](benchmarks/20261005-ic-semiregular-law/README.md)
- `20261011-lit2026-w4-early-abort`: no change or negative (word-XOR ratios 0.997–1.004). — [README](benchmarks/20261011-lit2026-w4-early-abort/README.md)
- `20261011-lit2026-w5-gate-first-lane-roots`: gate-first word-XOR ratios 1.5675–1.7788. Lane roots cost 2.6–4.0x the sequential control in the three batch cells, recorded as "(negative)". — [README](benchmarks/20261011-lit2026-w5-gate-first-lane-roots/README.md)
- `20261012-lit2026-w5-gate-first-rref`: ratios 1.64–2.41. — [README](benchmarks/20261012-lit2026-w5-gate-first-rref/README.md)
- `20261011-lit2026-w8-quotient-s3`: F4 word-op ratios 1.34–1.64. — [README](benchmarks/20261011-lit2026-w8-quotient-s3/README.md)
- `20261011-lit2026-w8-summation-models`: convolution splitting uses 1.21–2.07x more operations; Frobenius auxiliaries 57–701x more. — [README](benchmarks/20261011-lit2026-w8-summation-models/README.md)
- `20261012-lit2026-w12-zero-minor`: zero-minor/rho = 260–6,715x; negative. — [README](benchmarks/20261012-lit2026-w12-zero-minor/README.md)
- `20261011-pdp-ps2-scalar-group-law`: group-scalar/baseline per pair 0.067–0.123. — [README](benchmarks/20261011-pdp-ps2-scalar-group-law/README.md)
- `20261011-pdp-split-hybrid`: closest relation-phase projection is 22.8x strong rho. — [README](benchmarks/20261011-pdp-split-hybrid/README.md)
- `20261011-conductor-gap-round3`: E5 combined/blind F4 is 234–529x; negative. — [README](benchmarks/20261011-conductor-gap-round3/README.md)
- Also present without index rows: `20261012-f4-f5-*` (dense packing, fused readback, native hashing, history), `20261012-f4-fp-vs-msolve`, `20261011-f4-f5-row-allocation`. — [benchmarks/ listing](benchmarks/)

GPU and Semaev outside `benchmarks/`:
- GPU root-find decomposition: only host-emulation validation exists ("built without RF_HAVE_CUDA"). It passed 8,720 targets at 61 bits with fb4096, but `collect_relations` found 0 relations in 200,000 targets scanned. — [validate_61bit_fb4096_seed20261012.json](research/gpu_rootfind_decomposition_20261012/validate_61bit_fb4096_seed20261012.json)
- Semaev root finding against Gröbner F4 at |F| = 6: 34,602–46,911x at 16–24 bits. Prime root finding at |F| = 120–160 where F4 is infeasible: 28-bit run in 2.83 s. — [semaev_root_vs_groebner_20261012/console.txt](research/semaev_root_vs_groebner_20261012/console.txt)

Sparse linear algebra:
- The only status is an applicability review: "No ECDLP speedup or blocked-workload result is reported here." — [SPARSE_ALGEBRA_FOLLOWUP.md](SPARSE_ALGEBRA_FOLLOWUP.md)
- lit2026 W7 (Wiedemann pruning, Galois k² saving) is "not started". — [workstreams.json](research/lit2026_speedup_program/workstreams.json)

### Inferences
- The benchmark index understates IC work: several negative IC results (fiber, W4, W5, W8, W12) appear only as un-indexed directories.
- Recorded speedups cluster in the decomposition oracle (10²–10⁵x over F4), which moves constants but leaves every slope "at its generic floor".

### Gaps
- No GPU IC run on real hardware. No sparse-LA benchmark against dense at sizes where linear algebra matters.

## 4. Experiments that were frozen or preregistered but not run, and reports that say "next" or "TODO"

### Takeaway
There is a large backlog of preregistered work that never ran:
- the six PDP protocols D-1 to E-3;
- the m=83 IC gate in ecbench;
- lit2026 W1, W2, W6, W7 and W9;
- the Magma F4 comparison;
- dedicated two-socket hardware qualification;
- PDP τ-slots;
- k=4 steps 3–5 (E(F_{p^4}) in ecbench, the compiled S₅ solver);
- C2 sparse-spectrum and solving-degree census.

Several runs were censored or inconclusive.

### Cited Findings
Preregistered, never run:
- `pdp_preregistrations_20261007`: six protocols, all "PENDING": D-1 Artin–Schreier coordinates, D-2 linearisation reach, D-3 relation-matrix filtering against the LA cap, E-1 packed root index, E-2 batch verification, E-3 rank from index self-collisions. Instruments "not built". No other research note cites them, which I checked by grep. — [README.md](research/pdp_preregistrations_20261007/README.md)
- `ic_wide_shared_rank_plan_20261005`: m=83 IC arm in ecbench; "plan… no code, no measurement". — [PLAN.md](research/ic_wide_shared_rank_plan_20261005/PLAN.md)
- `ecbench_m83_ic_20261009`: PROTOCOL only (full-width coefficient sampler repair). — [PROTOCOL.md](research/ecbench_m83_ic_20261009/PROTOCOL.md)
- `ecbench_hardware_20261009`: "preregistered; no host launched and no hardware measurements taken". — [PROTOCOL.md](research/ecbench_hardware_20261009/PROTOCOL.md)
- `magma_f4_comparison_20261013`: "Registered… before any Magma output exists". It contains bundles and our-side jsonl but no Magma results. — [PROTOCOL.md](research/magma_f4_comparison_20261013/PROTOCOL.md)
- `pdp_split_slots_20261012`: registered; contains specs and `run_all.sh` but no session or result files. — [PROTOCOL.md](research/pdp_split_slots_20261012/PROTOCOL.md)
- `ic_exponent_top_20260930`: §23 at k=32 was withdrawn: "Nothing of it ran". — [ic_single_target_20260930/PROTOCOL.md](research/ic_single_target_20260930/PROTOCOL.md)
- `f5_projected_sparse_basis_20260930`: PROTOCOL and BASELINE only. — [directory](research/f5_projected_sparse_basis_20260930/PROTOCOL.md)
- lit2026 workstreams "not started":
  - W1 solving-degree diagnostics;
  - W2 multi-graded regularity;
  - W6 crossbred;
  - W7 sparse LA (Wiedemann);
  - W9 F_p F4 Las-Vegas echelon.

  W4, W8, W10, W11 and W12 are marked "PR open", although benchmark directories exist for W4, W5, W8 and W12, so the JSON may be stale. — [workstreams.json](research/lit2026_speedup_program/workstreams.json)
- PDP follow-ups blocked by the environment:
  - Xorricane XOR-CDCL;
  - WDSat/SATIC;
  - Bosphorus;
  - GPU pairs-and-solve (no GPU; CPU reference only to ℓ = 11).

  — [pdp_followups_20261010/README.md](research/pdp_followups_20261010/README.md)
- Exponent falsifiers "Still open": C2 sparse-spectrum zero finding on the tower "Not run yet"; solving-degree census with symmetries "Not run yet"; the k=4 symmetric rewrite "not built" (later done as D₄, see §5). — [ic_exponent_falsifiers_20261013/README.md](research/ic_exponent_falsifiers_20261013/README.md)
- IC scope programme:
  - planned levers L1–L7 (choose d by cost model; small-cofactor curves; lazy ordinates; Frobenius-folded pair table; block size; interleaved trials; threads);
  - five steps; step 3 (E(F_{p^4}) in ecbench) and step 4 (compiled straight-line S₅ solver, success gate C₄ at least 10x lower) are planned;
  - I found no ecbench F_{p^4} result. — [ic_speed_and_rho_scope_20261013/README.md](research/ic_speed_and_rho_scope_20261013/README.md)
- N83 next step: the K2940 reserve five-point IC/rho validation is named as "next". — [n83_factor_base_100_results_20261014/RESULT.md](research/n83_factor_base_100_results_20261014/RESULT.md)
- Generic-backend autolab: "next useful campaign is a new preregistered… complete-DLP comparison" qualifying F4/F5/SAT and larger sparse relation-LA. — [BOUNDED-GOAL-RESULT.md](research/ic_candidate_tournament_20260915/goal_20260924/BOUNDED-GOAL-RESULT.md)

Inconclusive, censored or negative:
- Generic-backend qualification v1 and v2: censored. v3: negative. — [ic-autolab SKILL.md](.claude/skills/ic-autolab/SKILL.md); [v3 RESULT.md](research/ic_candidate_tournament_20260915/goal_20260924/generic-backend-qualification-v3/RESULT.md)
- Prepared F5 and SAT one-target controls: eight attempts, no scalar. — [ic-autolab SKILL.md](.claude/skills/ic-autolab/SKILL.md)
- Semi-regular degree law: registered verdict "neither". — [ic_semiregular_law_20261005/RESULTS.md](research/ic_semiregular_law_20261005/RESULTS.md)
- Ledger §20 exponent: β undecided, [−0.10, 0.37]. — [ic_exponent_20260926/README.md](research/ic_exponent_20260926/README.md)
- PKM round 2: m=4 N=16 ran out of memory after 29 minutes. — [pkm_tower_round2_20260925/README.md](research/pkm_tower_round2_20260925/README.md)
- F6 m=31 stopped; matrix F4/F5 at m=19 timed out. — [ecbench_frontier_20261013/README.md](research/ecbench_frontier_20261013/README.md)
- Rootfind against rho chase (2026-10-14): cycles 4–5 were "dead end". Cold complete IC is still about 1.54x rho-strong at n=13 (negfold m2 d5, S 21.8 against 14.1). Recorded dead ends: orbit divisor subtract, window_div, d=3 subspace, d=4/6/7. — [ecbench_rootfind_vs_rho_20261014/README.md](research/ecbench_rootfind_vs_rho_20261014/README.md); [ecbench_rootfind_faster_20261014/README.md](research/ecbench_rootfind_faster_20261014/README.md)
- SMT ladder to n=23: "nothing decided at n = 23 by any arm". Hamming ONB SAT at n=23, m=3: neither solver decides any of 8 targets in 180 s. — [pdp_followups_20261010/README.md](research/pdp_followups_20261010/README.md)
- Root-finding lanes: the holdout winner was not resolved. — [benchmarks/README.md](benchmarks/README.md)

### Inferences
- The highest-leverage unrun items are the k=4 E(F_{p^4}) ecbench port with the compiled S₅ solver, and D-3 (relation filtering against the LA cap). These are the only places where theory predicts an exponent change or where linear algebra becomes binding.

### Gaps
- I did not verify, PR by PR, whether some "PR open" workstreams were later merged.

## 5. Bottlenecks in the measured pipeline: relation collection, decomposition or linear algebra?

### Takeaway
On every measured binary, Koblitz and prime ladder, relation collection (trials × per-trial scan or decomposition cost) dominates: 70–95% of S. Linear algebra is below 3% everywhere and often 0.0–0.2%. Gröbner-based decomposition was the historical bottleneck, but root finding removed it at small sizes. In one-target Koblitz runs, reusable setup is 98.5–99.7% of S. At k=4 over F_{p^4}, linear algebra alone already beats rho (0.52x), and C₄ (the cost of one decomposition) is the only blocker.

### Cited Findings
- Binary E_0 rootfind ladder (d = 5…10):
  - relation phase is 82%, 71%, 81% and 94% of S at n = 13, 19, 23, 31; oracle setup 6–29%;
  - "Factor base and linear algebra are under 2 % of S at every size";
  - hit rate falls from 0.51 to 0.085 at n31, set by the whole group: about (base points)³/(6·#E).
  - "Index calculus with any linear-scan oracle has total cost about 2^{2n/3}… against rho's 2^{n/2}." — [ic_speed_and_rho_scope_20261013/README.md](research/ic_speed_and_rho_scope_20261013/README.md)
- Leaderboard phase shares:
  - Table B prime/binary: Relations (PDP) 47–96%, factor base 4–52%, linear algebra 0.0–0.4%.
  - Table C Koblitz batch: linear algebra 0.1–2.7%; descent 7.8–38%; relations rising to 78.8% at 2^47.2.
  - Table A one-target: reusable setup 98.5–99.7%. — [docs/ic/LEADERBOARD.md](docs/ic/LEADERBOARD.md)
- Rootfind negfold GAE split at n13: oracle setup about 305 (pair table), relations about 700. — [ecbench_rootfind_vs_rho_20261014/README.md](research/ecbench_rootfind_vs_rho_20261014/README.md)
- Ledger §20: at small r the descent per target exceeds batch rho's cost below about 2^37. Big-integer workflow constructions add 7–70% of S. — [ic_exponent_20260926/README.md](research/ic_exponent_20260926/README.md)
- Decomposition cost per oracle:
  - F6-IC charged S is 64.6–68.8x inherited F4 at m13 and 1,303–1,353x at m19;
  - matrix F4 does 198.6–222.1x inherited F4's work. — [ecbench_frontier_20261013/README.md](research/ecbench_frontier_20261013/README.md)
  - Gaudry cubic (k=3) over F_{p^3}, p = 2083 (#E about 2^33): about 9.2·10⁵ F_p multiplications per residual, decomposition rate about 0.16–0.19. — [experiments/21_gaudry_cubic_la.log](experiments/21_gaudry_cubic_la.log)
- Extension field k=4:
  - C₄ was 17x its floor (1.2·10¹²); linear algebra alone is 0.52x rho on eight curves; the crossover was about 2^150;
  - "every 10× off C₄ moves the crossover down by about 2^13";
  - k=3 crosses only past 2^230. — [ic_speed_and_rho_scope_20261013/README.md](research/ic_speed_and_rho_scope_20261013/README.md)
  - D₄ symmetry then cut C₄ 48x, but the base decomposes about 1/8 as often (0.133 pooled), with half as many relations needed; net model crossover is 2^135 (h = 2). — [gaudry_quartic_d4_20261013/README.md](research/gaudry_quartic_d4_20261013/README.md); [benchmarks/README.md](benchmarks/README.md)
- Online regime: the N83 IC needs 580–640 s of setup and about 14 GiB RSS (an 8 GiB table in the K1798 variant), against about 60 ms per target online. — [n83_tier1_final_results_20261014/RESULT.md](research/n83_tier1_final_results_20261014/RESULT.md); [ecbench_frontier_20261013/README.md](research/ecbench_frontier_20261013/README.md)
- The `store/index.json` remote store holds only 10 factor-base dumps from the 2026-10-13 ecbench FB protocol (Koblitz n11/n17 cyclotomic, standard subspace, compact orbit). — [store/index.json](store/index.json)
- `taskq`, `isolab` and `spotlab` are infrastructure ("makes no performance claim"). They contain no IC measurements. — [taskq/README.md](taskq/README.md); [isolab/README.md](isolab/README.md); [spotlab/README.md](spotlab/README.md)

### Inferences
- For prime and binary curves at measured sizes, the bottleneck is the relation hit rate (a property of the whole group and the cofactor) times the per-trial scan cost. Faster linear algebra or sparse methods would change essentially nothing until base sizes grow by orders of magnitude. That explains why sparse-LA work (W7, D-3) was deprioritised, and why it would only matter at m=83-scale bases (3,000+ columns) or at k=4.
- The single most promising empirical lever recorded is driving C₄ down by a further 10–100x on E(F_{p^4}): a compiled record-and-replay S₅ solver combined with D₄ invariants. Then port it to ecbench to obtain the first harness-measured IC-below-rho ratio. The repository's own plan says this; the measurement is missing.

### Gaps
- No measured linear-algebra cost at N83 or other large bases appears in the files I read. The N83 RESULT reports online phases and setup, not a separate LA share.
- No GPU decomposition timing on a real device. No evidence on memory bandwidth as a bottleneck for the N83 root index beyond peak RSS.
