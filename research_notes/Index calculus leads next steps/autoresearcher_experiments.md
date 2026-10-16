# Index calculus experiments in aburan28/crypto-autoresearcher (executed runs and outputs)

Scope: local clone `/home/user/aburan28/crypto-autoresearcher` (repo root abbreviated `AR/`). Read only. The sources are local file paths, not URLs. `experiments/` holds 59,238 files. A grep for IC terms hits about 400 `EXP-*` directories, but most of those are draft or "closed-form calibrator" specifications that have no runs. This note covers the directories that have runs, analyses or reports. `outputs/` (10 files) holds only ledger digests, and `examples/` holds a single crypto-ml stub. Neither contains IC results. I did not read `harness/` or `evals/` in depth (see Gaps).

## Q1. Which IC experiment families exist (and were actually run)?

### Takeaway
About 15 executed families exist. All are toy-scale, at most 36-bit prime fields or GF(2^53) for binary fields. They cover prime-field Semaev S_3/S_4 decomposition with sympy Groebner, multi-target amortization, factor-base geometry batteries, binary Weil-descent d_reg/Macaulay rank, EcGFp5 n=5 symmetrized PDP, function-field xedni, BKK sparse checks, and complete cold-start IC vs rho on GF(2^53) Koblitz. Every executed comparison against Pollard rho came out negative or inconclusive. The formal/ Lean work proves only structural lemmas, not complexity bounds.

### Cited Findings
**Prime-field Semaev / multi-target amortization**
- EXP-SEMAEV-001 is an S_3 decomposition by sympy Groebner at 8/10/12 bits with 2 seeds each. Groebner took 0.157–0.173 s, flat across bit sizes, and the max-degree proxy was ≤4. Decompositions were found in 3 of 6 runs. Matched rho needed 6–30 group ops. The analysis itself says "This does not establish a scaling law — the bit range is too narrow and timing is sympy-bound." — [AR/experiments/EXP-SEMAEV-001/analysis.md](AR/experiments/EXP-SEMAEV-001/analysis.md)
- EXP-IC-001 measured the multi-target crossover K*. It covered 8–36 bits, seeds 1–3, factor base B=14, m=2 (S_3) and K=10 targets, with 528 runs. T_desc came to about 32K–55K group-op equivalents: Groebner takes about 0.04 s and the calibration is about 1M rho gops/s. The mean T_desc/sqrt(N) ratio fell from 10,055 (8 bits) to 32 (36 bits), and only 1 of 3 seeds at 36 bits had a finite K* (ratio 0.17). The extended 36-bit set had 13 seeds, and 4 of them had a finite K*, all with N > 3×10^9. The result was rated inconclusive. A key caveat: "all Groebner solves return trivial ideals (no decompositions)". — [AR/experiments/EXP-IC-001/analysis.md](AR/experiments/EXP-IC-001/analysis.md)
- EXP-MTIC-001 measured the multi-target K* against the S·T² = N frontier at 16, 20 and 24 bits. At 16 bits (N=17,623): S_rel 0.022 s, S_LA 0.011 s, and the T_desc median hit the 60 s cap on all 25 targets. Rho needed a median of 1,720 group ops (0.001 s). K* was infinite and the frontier ratio was 3.0×10^16. At 20 bits, T_desc was capped on 13 of 13 attempted targets (37 cancelled), K*=∞ and the frontier ratio was 1.2×10^16. At 24 bits (N=11,000,719): S_rel 9.86 s, S_LA 8.64 s, T_desc capped on 12 of 12, K*=∞, frontier ratio 1.3×10^16. The 20-bit ablations at B=52 and B=115 were also all capped. — [AR/experiments/EXP-MTIC-001/analysis.md](AR/experiments/EXP-MTIC-001/analysis.md)
- EXP-FCP-001/002/003 tested fixed-curve factor-base preprocessing at 8–36 bits with B=14, m=2, K=3–5, over 72+144+72 runs. The amortized fixed/naive cost ratio was 1.125 (16 bits), 1.169 (20), 1.067 (24), 1.061 (28), 1.023 (32) and 1.040 (36), and never below the 0.95 target. The precompute/Groebner share plateaued at about 1.6%. Groebner stayed at about 0.04–0.05 s regardless of p. "No decompositions found at any bit size" for bits ≥16. EXP-FCP-001 had no frozen spec and all its runs were marked `dirty: true`. — [AR/experiments/EXP-FCP-001/analysis.md](AR/experiments/EXP-FCP-001/analysis.md), [AR/experiments/EXP-FCP-002/analysis.md](AR/experiments/EXP-FCP-002/analysis.md), [AR/experiments/EXP-FCP-003/analysis.md](AR/experiments/EXP-FCP-003/analysis.md)
- EXP-MTBK-306bdb tested the BKK sparse-check speedup on the p=1009 family over 36 cells (m=3,4,5). The mean descent ratio was 1.170, 1.037 and 0.957 against the theorem's 4, 8 and 16 (2^(m-1)). The relation ratio was 1.822, 0.968 and 1.427. Only the m=3, b=0.4 cells with B=13–16 reached 3.35–4.00. The conclusion was "WEAKENED", and the rescue-window claim was called "UNFALSIFIABLE at toy scale (empty corridor)". The spec status is still `draft`. — [AR/experiments/EXP-MTBK-306bdb/analysis.md](AR/experiments/EXP-MTBK-306bdb/analysis.md)

**Factor-base structure (symmetric / invariant / structured bases)**
- EXP-FB3-001 compared six structured bases against matched random ones at group sizes 2^14–2^18. 0 of 6 geometries met the success criterion. Three (high_bit_interval, small_height, coset_union) sit exactly at ratio 1 with slope 0. Three (mixed_two_base, asymmetric_sizing, greedy_optimized) show deficits below 1 at every size. The analysis notes that the primary metric was "mathematically incapable of showing a matched-size geometry effect" and that this was known before execution (conservation.md). — [AR/experiments/EXP-FB3-001/analysis.md](AR/experiments/EXP-FB3-001/analysis.md)
- EXP-PFDR-1b78f7 is a coincidence census of subgroup/small_x/dickson bases against random ones at 12–32 bits, m=3,4,5, five curves per rung, over 17 runs. It returned four "excursions" among 166 cells (z 3.39–4.75). None replicated at Stage R (z 0.298, −1.989, −2.959, −0.431). In relation units the count drops to 0 of 133. The known-null random arms "excursed" more often (5.8%) than the structured arms (3.6%). The m=4 small_x harvest-on exponent is 0.6076 [0.601, 0.614]. All 14 A2 slopes are negative (−0.22 to −0.11). The bootstrap intervals under-cover: 0.864–0.916 against a nominal 0.95. — [AR/experiments/EXP-PFDR-1b78f7/analysis.md](AR/experiments/EXP-PFDR-1b78f7/analysis.md)
- EXP-PFDR-011cd0 is the relation-level census follow-up. R16 returned outcome_ids [O-NULL-BOUNDED, O-HEUR], so R17 and R18 were not triggered. EXP-PFDR-0b3699 replicated the small_x TT4 excess at 30/32 bits and returned outcome O-A3. — [AR/experiments/EXP-PFDR-011cd0/execution-report-236691.yaml](AR/experiments/EXP-PFDR-011cd0/execution-report-236691.yaml), [AR/experiments/EXP-PFDR-0b3699/execution-report-9456e8.yaml](AR/experiments/EXP-PFDR-0b3699/execution-report-9456e8.yaml)

**Weil descent / binary fields / degree of regularity**
- EXP-DREG-001 used block-m4ri exact GF(2) rank on the boolean chained Semaev m=3 system (t=3, D=5). The deficit vs semi-regular was 1,322, 1,862 and 1,999 at n=12, 15 and 18, with ranks 28,096, 69,073 and 143,882 reproduced by three engines. The null equals sr_pred at n≤15. d_ff at n=12 was 3 (ti=0) and 2 (ti=1), with only 2 of 8 targets resolved. Every past-wall cell (n≥17) was censored. "Exponent fit over three points is meaningless." — [AR/experiments/EXP-DREG-001/analysis.md](AR/experiments/EXP-DREG-001/analysis.md)
- EXP-DREG-004 is the decisive n=21 cell: a 279,048 × 778,394 Macaulay matrix with sr_pred 268,674. It stopped at 192,000 of 778,394 columns processed (24.67%), with a partial rank of 186,574, peak RSS 12.0 GB and done=false. It lost data in two exFAT payload deletions (INC-20260724/25-EXFAT-01/02). Status is "running"/"IN PROGRESS". — [AR/experiments/EXP-DREG-004/analysis.md](AR/experiments/EXP-DREG-004/analysis.md)
- research/dreg-linear-law, FINDING_v2 (independent session, not a ledger record): semi-regular d_reg for the t=3 family is Θ(n) with density c* ≈ 0.2374790709. For the crypto-scale t=7, n=161 instance (nb=966), d_reg=150, which puts Groebner at 2^1194–2^1415 against rho at 2^80.5. The red team retracted claim 5 (d_reg(sem) > d_reg(null)), and it was re-established with a graded-HF instrument at n=6 and n=9 only. First fall is low (D4 at n=6, D3 at n=9), but the relation supply there is O(1)–O(n). — [AR/research/dreg-linear-law/FINDING_v2.md](AR/research/dreg-linear-law/FINDING_v2.md), [AR/research/dreg-linear-law/RED-TEAM.md](AR/research/dreg-linear-law/RED-TEAM.md)
- The EXP-SEMBIN-* family covers binary Semaev experiments. 7e1371 locates the boundary of Semaev's Assumption 1 (d_F4 = 4 vs ≥5) using degree-4 Macaulay rank certificates. f4a17b and 2c40bb run a memory-charged crossover against van Oorschot–Wiener rho (status: analyzed). 92724f localizes the Semaev-2015 exponent constant to the m! factor. 992e73 (Certicom/FIPS binary curves) is only "specified". 79a02d is marked "invalid". I read the titles and statuses only, not the numbers. — spec titles under [AR/experiments/EXP-SEMBIN-*/specification.yaml](AR/experiments/)
- EXP-FROB-30006a compares the WDSat conflict ratio against 2^{ml}/m! on Frobenius-stable vs random factor bases (n=41 primary, n=43 capacity-gated), with verdicts "basis_blind" or "pruning". I did not locate an outcome. — [AR/experiments/EXP-FROB-30006a/specification.yaml](AR/experiments/EXP-FROB-30006a/specification.yaml)
- EXP-ICPERF-66fd51 (status: analyzed) ran "every installed engine on Trimoska's 60 shipped S_4 Weil-descent instances (n,l) ∈ {(15,5),(17,6),(19,6)}". I did not extract its numbers. — spec title, [AR/experiments/EXP-ICPERF-66fd51/specification.yaml](AR/experiments/EXP-ICPERF-66fd51/specification.yaml)

**Extension field F_{p^5} (EcGFp5) symmetrized descent**
- EXP-GFPN-05ff43 is an n=5 Gaudry PDP ladder on EcGFp5 with raw, S_5 and 2-torsion⋊S_5 arms. Of 48 run packages, 34 completed validly, 5 failed on infrastructure, 3 on implementation and 6 on resource exhaustion (11 GiB RLIMIT_AS). No m=5 cell completed, so the success criterion was met for no arm. The m=4 system had 4 variables, 5 equations of degree 8 and 495 monomials, and msolve reproduced the F4 trace. The Joux–Vitse timing anchor ratio was 0.470, which fails the [0.5, 2.0] band. The "torsion" arm as built is a norm with the same support and F4 traces as the S_5 arm (20,349 = 20,349 at m=5), so it gives no extra quotient. The frozen degree prediction was "mis-anchored by 120" (5!). Re-anchored model bands are 2^134.4–2^152.7, "model only". — [AR/experiments/EXP-GFPN-05ff43/analysis.md](AR/experiments/EXP-GFPN-05ff43/analysis.md)

**Xedni (function-field lift)**
- EXP-XEDN-002 gives the exact lift probability 2p^3·P_lift = 1 − 2/p + p^{-3}. The finite-size exponents were 2.9858, 2.9932 and 2.9965 (OLS 2.9917), so the limiting α = 3, not < 1/2. P[≥9 sections] is 0 at p=5, 7, 11 and 1.35e-6 at p=13. The success criterion is "not met, on both clauses", and the falsification conditions are met. That is a scoped rejection of H-XEDN-001. XEDN-003 to 007 (Mordell–Weil relation-coefficient censuses) and many review_required audits follow. — [AR/experiments/EXP-XEDN-002/analysis.md](AR/experiments/EXP-XEDN-002/analysis.md)

**Binary Koblitz complete-solve IC vs rho (GF(2^53))**
- research/cold_single_ic_20260921 used the fixed Koblitz curve over GF(2^53) with subgroup order 21,044,858,204,113, one worker and ARM PMULL arithmetic in both arms. Signed-Frobenius rho took a median of 1.324 s. The selected 75-column IC took 4.788 s, a paired wall ratio of 3.511, with 439.9 MB RSS against 174.6 MB. "The cold-start rho target has not been met." — [AR/research/cold_single_ic_20260921/REPORT.md](AR/research/cold_single_ic_20260921/REPORT.md)
- cold_capacity_ic (doubling the hash slots to 33,554,432) gave an IC/rho paired median of 3.561 [2.622, 4.396], against 4.068 at the original capacity, and failed. cold_target_stop_ic (RUN-KIC-56c897, a target-certificate early stop) saved 6.4% against the IC baseline (ratio 0.9358), but its IC/rho ratio was still 4.153 [2.900, 5.707] and it failed. — [AR/research/cold_capacity_ic_20260921/REPORT.md](AR/research/cold_capacity_ic_20260921/REPORT.md), [AR/research/cold_target_stop_ic_20260921/REPORT.md](AR/research/cold_target_stop_ic_20260921/REPORT.md)
- The N19/N7/N23 affine-flat SAT/XOR factor-base work (EXP-KIC-*, research/golden_factor_base_sat, n19_affine_sat, admissible_affine_factor_base) consists mostly of design and derivation drafts. One example: "The cost tradeoff is unresolved" for the orbit-canonical pair transport. — [AR/research/n19_affine_sat_20260921/NEXT_STEP_PAIR_TRANSPORT.md](AR/research/n19_affine_sat_20260921/NEXT_STEP_PAIR_TRANSPORT.md), [AR/research/admissible_affine_factor_base_20260921/COFACTOR_PROJECTION_DERIVATION_DRAFT.md](AR/research/admissible_affine_factor_base_20260921/COFACTOR_PROJECTION_DERIVATION_DRAFT.md)

**Infrastructure and formal**
- groebner_compare/ is "solver-stage infrastructure" for Boolean systems with 1–12 variables. It "does not implement an ECC attack ... or establish a speedup", and its smoke data are "correctness controls, not performance evidence". The ordinary-relation campaign lives in an external repo (cryptanalysis PR #161). — [AR/groebner_compare/README.md](AR/groebner_compare/README.md), [AR/groebner_compare/ORDINARY_RELATIONS.md](AR/groebner_compare/ORDINARY_RELATIONS.md)
- formal/ Lean contains Semaev/Symmetrization.lean (esymm homogeneity), Semaev/FrobeniusCollapse.lean (X^q − X in the field equations, Frobenius twist mod the field equations), Semaev/DegreeTwoFall.lean (a trace witness annihilates the degree-two part, degeneracy root uniqueness), and FallDegree/FirstFallBounded.lean and LastFall.lean (monotonicity of V_i filtrations). The README says nothing containing `sorry` lands there. None of these is an IC complexity bound. — [AR/formal/README.md](AR/formal/README.md), [AR/formal/CryptoResearch/Semaev/](AR/formal/CryptoResearch/Semaev/)
- ML-guided decomposition exists only as specs that are review_required or approved: EXP-RELN-82f487 (a GNN/GBT decomposability predictor), EXP-RELN-141a86 (PySR symbolic regression) and EXP-RELN-715aa1 (a learned search-policy gain, Stage-0 zero-compute). I found no executed ML results. — spec titles in [AR/experiments/EXP-RELN-*/specification.yaml](AR/experiments/)

### Inferences
- The executed evidence is concentrated in two settings. One is toy prime-field S_3 (m=2) with sympy Buchberger and B=14, where the systems are mostly trivial ideals. The other is binary-field Macaulay-rank and d_reg work. Very little completed data exists for m≥4 in the prime-field setting with a real F4 or msolve solver.
- A large share of the "experiments" are protocol and audit artifacts: calibrators, forced-value audits and review rounds. The executed run count is far smaller than the directory count suggests.

### Gaps
- I did not extract numbers from the SEMBIN family, ICPERF-66fd51 (Trimoska instances), FROB-30006a, ICINV-4d33aa/55c2d8 (within-isogeny-class over-dispersion), YIELD-001/004, DTREE-001, DS-001, ICI-001, LPF/SMTH (Dickman smoothness), ICEX-aaccfc/153c34 (m=5 membership), or STR-003/004 (phi-invariant factor base displacement rank). These need a follow-up read.
- I did not read harness/ (135 files) or evals/ (5 files) for embedded results. The tool-call budget ran out.

## Q2. What scaling trends were measured; are fits recorded and reliable?

### Takeaway
A few exponents are recorded. They are the m=4 small_x harvest exponent 0.6076, slopes of −0.22 to −0.11, the xedni α ≈ 2.99 (exact limit 3), and the analytic d_reg density c* ≈ 0.2375. Most "trends" are 3–8 point toy series with 2–5 seeds, and the system itself flags its own bootstrap intervals as under-covering.

### Cited Findings
- EXP-PFDR-1b78f7 rests on five curves per rung. stats.bootstrap_slope covers 0.887–0.909 and the A2 curve-index bootstrap covers 0.864–0.875, against a nominal 0.95, which the analysis explains by the sqrt((n−1)/n) deflation at n=5. — [AR/experiments/EXP-PFDR-1b78f7/analysis.md](AR/experiments/EXP-PFDR-1b78f7/analysis.md)
- EXP-DREG-001 has a three-point deficit series. Its own analysis says "Exponent fit over three points is meaningless". — [AR/experiments/EXP-DREG-001/analysis.md](AR/experiments/EXP-DREG-001/analysis.md)
- EXP-IC-001's T_desc/sqrt(N) ratio falls about 300× over 8–36 bits. The crossover "depends on N (the subgroup order), not the field bit size", and the Groebner cost is constant because the systems are trivial ideals. — [AR/experiments/EXP-IC-001/analysis.md](AR/experiments/EXP-IC-001/analysis.md)
- EXP-FCP-003: the Groebner time is p-independent because the polynomial system "(degree 14, 2 variables) is structurally fixed". — [AR/experiments/EXP-FCP-003/analysis.md](AR/experiments/EXP-FCP-003/analysis.md)
- EXP-XEDN-002's exponents are exact rational computations, not statistical fits, which makes them reliable. — [AR/experiments/EXP-XEDN-002/analysis.md](AR/experiments/EXP-XEDN-002/analysis.md)
- dreg-linear-law's c* is analytic for the semi-regular model only. The structured system's d_reg was checked at n=6 and n=9 only. — [AR/research/dreg-linear-law/FINDING_v2.md](AR/research/dreg-linear-law/FINDING_v2.md)

### Inferences
- Wall-time trends from sympy Buchberger, as in IC-001, FCP and MTIC, measure implementation overhead and not algebraic complexity. They should not be extrapolated.
- The only robust scaling statements are the exact or analytic ones: the xedni α=3, and the d_reg Θ(n) under the semi-regular assumption.

### Gaps
- No prime-field degree-of-regularity series for m≥3 with a real F4 or msolve solver over a meaningful bit range was found among the analyses I read.

## Q3. Aborted, timed-out or half-finished experiments

### Takeaway
Several key cells were never completed. The decisive cells (n≥17 and n=21 d_reg, the EcGFp5 m=5 cells, every MTIC descent) were censored by memory caps, timeouts or data-loss incidents.

### Cited Findings
- EXP-DREG-001: all past-wall cells (n=17–24) were censored. EXP-DREG-002 (the n=17 cells) is still "running". — [AR/experiments/EXP-DREG-001/analysis.md](AR/experiments/EXP-DREG-001/analysis.md)
- EXP-DREG-004 (n=21) stopped at 24.67% of columns. Two exFAT incidents deleted checkpoints, and the countermeasure was mirroring carries into git objects. — [AR/experiments/EXP-DREG-004/analysis.md](AR/experiments/EXP-DREG-004/analysis.md)
- EXP-GFPN-05ff43 completed no m=5 cell. It had 6 resource-exhaustion failures under the 11 GiB cap, though the cap recorded in the record (12 GB) was not the cap actually in force. msolve -g 1 "crashes inside F4 at the degree-9 round at p' = 65551". The implementation was never committed before the runs. — [AR/experiments/EXP-GFPN-05ff43/analysis.md](AR/experiments/EXP-GFPN-05ff43/analysis.md)
- EXP-MTIC-001: every descent hit the 60 s cap, and 25+37+38 targets were cancelled across the 16/20/24-bit cells. — [AR/experiments/EXP-MTIC-001/analysis.md](AR/experiments/EXP-MTIC-001/analysis.md)
- EXP-IC-001: 40-bit instance generation took about 262 s per instance and was not pursued. — [AR/experiments/EXP-IC-001/analysis.md](AR/experiments/EXP-IC-001/analysis.md)
- Statuses: SEMBIN-509d41 "running", SEMBIN-79a02d "invalid", SEMAEV-a097b5 "cancelled", SMTH-9d04ba "cancelled", SMTH-71b1b0 "frozen_awaiting_execution", ICEX-146ff5 "approved_execution_withheld", MTBK-306bdb still "draft" despite having runs. — spec status fields, [AR/experiments/](AR/experiments/)

### Inferences
- The scaling questions the program cared most about are the ones that hit infrastructure walls: memory at n≥17 Macaulay and m=5 F4, and disk loss. The negative conclusions at those scales therefore rest on models, not measurements.

### Gaps
- I could not tell whether EXP-DREG-004 was ever resumed after the July 2026 incidents. Its analysis is the latest record I read.

## Q4. Surprising or anomalous results never followed up

### Cited Findings
- EXP-IC-001 at 36 bits: 4 of 13 seeds had a finite K*, all those with N > 3×10^9. The analysis calls the conditional "definitive", but the "crossover" arises because the Groebner calls return trivial ideals cheaply. It is not a real descent. EXP-IC-002 (a yield-charged control) has only a spec. — [AR/experiments/EXP-IC-001/analysis.md](AR/experiments/EXP-IC-001/analysis.md), [AR/experiments/EXP-IC-002/](AR/experiments/EXP-IC-002/)
- EXP-DREG-001 found d_ff(ti=0) = 3 at n=12, against the prior series' flat 2. It was flagged as surprising and stayed unresolved, since only 2 of 8 targets resolved. — [AR/experiments/EXP-DREG-001/analysis.md](AR/experiments/EXP-DREG-001/analysis.md)
- EXP-PFDR-1b78f7: (subgroup, TT, m5, 16 bits) showed z 4.752 with kappa 4.412, the largest excursion. Stage R had power 0.89–1.00 for this cell and still did not replicate it. In the subgroup m=5 cell, delta was −0.0375 [−0.058, −0.017]. These remain detection-limited non-detections. TT m=4 and all TB tests were "unresolved". — [AR/experiments/EXP-PFDR-1b78f7/analysis.md](AR/experiments/EXP-PFDR-1b78f7/analysis.md)
- EXP-MTBK-306bdb: the m=3, b=0.4 relation-collection cells reached 3.35–4.00×, the full sweep bound. They were attributed to near-exhaustive enumeration at tiny B and were not followed up at larger B. — [AR/experiments/EXP-MTBK-306bdb/analysis.md](AR/experiments/EXP-MTBK-306bdb/analysis.md)
- EXP-GFPN-05ff43: the cross-arm ratio of 39.63 becomes 13.0–148.7 depending on the accounting convention. "No committed F3 verdict exists." — [AR/experiments/EXP-GFPN-05ff43/analysis.md](AR/experiments/EXP-GFPN-05ff43/analysis.md)
- cold_target_stop_ic: a real 6.4% IC-internal saving (paired ratio 0.9358) was measured. It is too small to change the IC/rho gap of about 4×. — [AR/research/cold_target_stop_ic_20260921/REPORT.md](AR/research/cold_target_stop_ic_20260921/REPORT.md)

### Inferences
- The most interesting unclosed anomaly is the PFDR subgroup TT m=5 cell. The others look like artefacts: trivial ideals, a pair-unit counting convention, tiny B.

## Q5. What conclusions did the system draw, and do they hold up?

### Takeaway
The system's conclusions are mostly appropriately cautious and scoped to toy settings: "inconclusive", "weakened", "scoped rejection" or "detection-limited non-detection". Its review rounds caught several real defects: the mis-anchored 5! factor, the miscalibrated nulls, the pair-unit artefact, and the retracted d_reg claim. A few headline readings overreach the data. The strongest well-supported negatives are xedni (α=3 exactly) and GF(2^53) cold IC at about 3.5–4.2× slower than rho.

### Cited Findings
- PFDR: the coordinator explicitly refused the frozen O-NULL closure: "The measured band is miscalibrated, so it is not a valid obstruction record." — [AR/experiments/EXP-PFDR-1b78f7/analysis.md](AR/experiments/EXP-PFDR-1b78f7/analysis.md)
- dreg-linear-law concludes that binary Weil-descent Semaev is a "clean negative control — no sub-rho signal at any tested scale". This was hardened after the red team, but the analysis notes it is "not an official ledger record". — [AR/research/dreg-linear-law/FINDING_v2.md](AR/research/dreg-linear-law/FINDING_v2.md)
- GFPN: the success criterion was met for no arm. The ratio predictions survive re-anchoring (D_raw/D_torsion 2^10.9), but the band prediction "does not survive in its frozen form". — [AR/experiments/EXP-GFPN-05ff43/analysis.md](AR/experiments/EXP-GFPN-05ff43/analysis.md)
- FB3: "no geometry produced any advantage at any size on the primary metric", with a pre-known blind spot in that metric. — [AR/experiments/EXP-FB3-001/analysis.md](AR/experiments/EXP-FB3-001/analysis.md)
- IC-001 claims its trend "strongly supports the hypothesis direction". The same file admits the systems were all trivial ideals and the results are implementation-bound. — [AR/experiments/EXP-IC-001/analysis.md](AR/experiments/EXP-IC-001/analysis.md)

### Inferences
- The IC-001 "supports the hypothesis direction" reading does not hold up. A falling T_desc/sqrt(N) with constant-time trivial-ideal Groebner calls reflects sqrt(N) growing, not descent getting cheaper. MTIC-001, where real descents all time out at 60 s, contradicts any favourable crossover.
- The crypto-scale binary-field verdict (2^1194+ against 2^80.5) rests on the semi-regular model applied to a structured system. The structured system was validated only at n ≤ 18 and was measured to have a *higher* d_reg than the null. The direction of the claim is robust, but the exact exponent is model-dependent.
- No executed experiment shows IC beating rho at any size. The closest approach is the GF(2^53) Koblitz complete solve, at about 3.5× slower.

### Gaps
- Coordinator decision records (DEC-*, EV-*) live in ledger/. I did not read them, so the final official status transitions for many of these experiments are unverified here.
