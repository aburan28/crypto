# Index calculus in aburan28/crypto-autoresearcher: ideas, directions, claims, ledger and knowledge base

All links are paths relative to the local clone `/home/user/aburan28/crypto-autoresearcher` (read-only survey, 2026-10-14). Scope: direction documents, ideas/, ledger/ (hypotheses, proposals, decisions, evidence records, goals), knowledge/, claims/, focus/, notes/, papers/, research_digest/. The repo's experiments/ and outputs/ were left to the other researcher. Size check: `ideas/` 1,174 files (873 match IC terms), `ledger/` 16,118 files (4,523 match), `knowledge/` 8,485 files (2,671 match, mostly literature stubs), `kb/` is a code service (79 files, no research content), `claims/` holds only a README, `ledger.json` has no IC content.

## Q1. What IC directions appear in research_directions_20260717/18, research_goals_*, ROADMAP.md, plan.md and PUBLICATION-CANDIDATES.md?

### Takeaway
The two July direction rounds proposed 24 IC-adjacent candidates. The 2026-07-17 round dispatched all 12 at once and every attack candidate was scoped-rejected within a day, while its barrier candidates were confirmed. The 2026-07-18 round ranked monodromy census (A1, 27), imprimitive monodromy (C2, 26), S-unit budgeting (A3, 25) and function-field xedni (B2, 24) highest, and ran only harness smoke tests. The 30-goal program of 2026-07-23 then reorganised IC into cost-ruler, degree, monodromy, factor-base, linear-algebra and m≥4 goals, and explicitly deferred xedni, Coppersmith, Buium and Heegner. PUBLICATION-CANDIDATES.md and plan.md contain no IC directions beyond the 07-17 dispatch table.

### Cited Findings
**Round 1, 2026-07-17. Twelve candidates, all dispatched (plan.md, EXP-*-001). Outcomes are from ledger decisions:**
- A1 tangent-split / dual-number jet summation → **rejected_scoped**. The survival rate was σ=1.0000 at p=101..1009 (9,984 pairs), so the ε-system is implied by the zeroth-order solution set — [research_directions_20260717.md](research_directions_20260717.md) §A1; [ledger/DEC-20260718-001.yaml](ledger/DEC-20260718-001.yaml)
- A2 output-sensitive incidence reporting for chord harvesting → **rejected_scoped**. Candidate-test exponent was 2.078 = Θ(B²) and no pairs were avoided — [ledger/DEC-20260718-002.yaml](ledger/DEC-20260718-002.yaml)
- A3 displacement-rank relation matrices from arithmetic-progression supports → **rejected_scoped**. α/B was 0.56–0.77 (generic), and the AP relation-stage penalty was 17.5x–4,128.6x against a 1.5x gate — [ledger/DEC-20260718-003.yaml](ledger/DEC-20260718-003.yaml)
- B1 elliptic-net (EDS/Somos) representation → **rejected_scoped**. Collision enrichment of 1.05–1.20 matched birthday plus the tautology class — [ledger/DEC-20260718-004.yaml](ledger/DEC-20260718-004.yaml)
- B2 BKK / Newton-polytope sparse elimination of Semaev → **rejected_scoped**. Supports are box-saturated, so MV = Bézout_box exactly at m=3,4,5 — [ledger/DEC-20260718-005.yaml](ledger/DEC-20260718-005.yaml). Later **proved** for all m, with MV_6=125,829,120 certified — [ledger/DEC-20260722-003.yaml](ledger/DEC-20260722-003.yaml)
- B3 equivariant / isotypic index calculus on the Semaev fibre product → **rejected_scoped**. Six of 10 isotypic blocks are identically zero and rank never splits — [ledger/DEC-20260718-006.yaml](ledger/DEC-20260718-006.yaml)
- C1 Koopman / transfer-operator spectral channel → **rejected_scoped**. Localisation L=O(1), with fitted δ=0.0195 — [ledger/DEC-20260718-007.yaml](ledger/DEC-20260718-007.yaml)
- C2 tensor-network contraction of the Semaev tensor → **rejected_scoped** on the truncation gate (α=1.161/1.322) — [ledger/DEC-20260718-008.yaml](ledger/DEC-20260718-008.yaml). An exact bond-rank law χ = 3, 6, 15, 45 at d = 2, 4, 8, 16 was later certified — [ledger/DEC-20260722-002.yaml](ledger/DEC-20260722-002.yaml)
- C3 noncommutative path-algebra syzygy search → **rejected_scoped**, then **proved** as a theorem: every word relation descends to commutative subset-sum — [ledger/DEC-20260718-009.yaml](ledger/DEC-20260718-009.yaml), [ledger/DEC-20260722-005.yaml](ledger/DEC-20260722-005.yaml)
- D1 jet barrier, D2 MV growth law and D3 chord-incidence ceiling → all **supported**, and later **proved** — [ledger/DEC-20260718-010..012](ledger/DEC-20260718-010.yaml); [ledger/DEC-20260722-004.yaml](ledger/DEC-20260722-004.yaml); [ledger/DEC-20260722-006.yaml](ledger/DEC-20260722-006.yaml)

**Round 2, 2026-07-18. Twelve candidates with a ranking table** — [research_directions_20260718.md](research_directions_20260718.md) §5:
- Scores: A1 Chebotarev / monodromy census 27; C2 imprimitive-monodromy resolvent 26; A3 Mason–Stothers / S-unit relation budget 25; B2 function-field xedni (Mordell–Weil lattices over F_p(t)) 24; A2 isogeny-class advice 23; B1 AG-code / MacWilliams relation enumerator 23; B3 Coppersmith small-root windowed relations 22; C1 Buium δ-characters 20 (rejected as winner); C3 Heegner / Drinfeld heights 19 (rejected). D1 log-space GAP barrier, D2 full-monodromy barrier and D3 cd(F_p)=1 barrier are barrier candidates exempt from the score.
- The document says no candidate scored 3 or more on "exponent chance", because "761 negative rows have closed every constant-level channel tested" — same file §5
- Winners A1, B2 and C2 ran smoke tests only. A2's cross-curve advice was null at p≤431. MONO phase 1 passed controls, with "the substantive Chebotarev gate is phase 2 (m=3 Semaev census) — not yet run". The xedni scarcity test found 0 sections in 5,760 slots, which the document calls "consistent with JKSST 2000" — same file §9.1–9.4

**30-goal program, 2026-07-23** — [research_goals_20260723.md](research_goals_20260723.md):
- IC goals in the program:
  - G30-001: group-op-equivalent cost ruler
  - G30-004/005/006: true d_reg backend, D6@n=12 cell, prime-field first-fall census m∈{2,3,4} on an 8–20-bit ladder
  - G30-007/008/009/012: Chebotarev census, imprimitive resolvent, wreath barrier, exceptional-curve sieve
  - G30-010: structured factor-base geometry battery
  - G30-011: Edwards / Jacobi / Hessian representation census
  - G30-013: cyclic-cover label-conditioned smoothness
  - G30-014: AG-code enumerator
  - G30-015: split-Jacobian projected smoothness
  - G30-016/017: sparse and structured LA benchmark; displacement-rank supports net of penalty
  - G30-018: 14/16 rank-saturation origin
  - G30-022/023: D1 and D3 barriers
  - G30-024: Ritt indecomposability
  - G30-026: m=3/m=4 decomposition with a certified 40–64-bit ladder
- Explicit deferrals (§4): CLM-ECDLP-BLIND-DESCENT, the end-to-end integration owner, and round-2 B2 xedni, B3 Coppersmith, C1 Buium and C3 Heegner, "unstaffed and deferred pending theory formalization" — same file §4
- [research_goals_20260831_coverage_gap_batch.md](research_goals_20260831_coverage_gap_batch.md) records an external "do not fund" list: "further F_p Semaev tuning for a rho win, xedni, canonical-height lifting, pairing inversion, first-fall-degree over F_p". The batch adopts it as "independent concurrence" (lines ~85, 259). The other research_goals_* files (PQC, AES, hash and NGCC batches) have 0 IC matches.

**ROADMAP, plan and publication candidates:**
- [ROADMAP.md](ROADMAP.md) Phase 4 names four programs: an isogeny-neighbourhood audit with matched Semaev/Gröbner; representation search (coordinate models, factor bases, symmetrisations); pair-selection learning (heuristic, imitation and RL Buchberger policies); and predicate/decomposition channel tests. Its ranking criterion is Wesolowski-style "target the asymptotic exponent … against explicit numbered heuristics".
- [plan.md](plan.md) is the 07-17 dispatch table (EXP-JET/INC/STR/NET/BKK/EQJ/TRA/TTN/NCP/JETB/BKKMV/INCB-001).
- [PUBLICATION-CANDIDATES.md](PUBLICATION-CANDIDATES.md) covers only the AES campaign. It has no IC content.
- IC papers drafted elsewhere: [papers/semaev-conservation-specialization/paper.tex](papers/semaev-conservation-specialization/paper.tex) (the Galois group of S_m is (Z/2)^{m-2}, plus a factorisation dichotomy and yield conservation), [papers/nagao-lemma4-note/](papers/nagao-lemma4-note/note.md) (counterexample to Nagao 2015/984 Lemma 4), and [papers/binary-ic-blog/](papers/binary-ic-blog/DRAFT.md).

### Inferences
- Round 1 was "dispatch everything, kill fast". Every attack candidate died on the named fatal obstruction, and the most durable output was barrier theorems (jet simulability, MV = Bézout, chord ceiling, NC to commutative collapse).
- Round 2's three winners all centred on one object, the monodromy of the summation cover. The program's own red team admitted that A1, C2 and D2 are "one object with three stances" (§6.1). That object was later resolved mathematically (see Q5).

### Gaps
- No document re-ranks the round-2 non-winners (A3 S-unit budget at 25, B1 AG-code at 23, B3 Coppersmith at 22) after G30. I found no evidence record for A3 or B3. G30-014 (B1) has proposals but I did not confirm an EV record.

## Q2. What do the claims and the ledger record about IC, both verified and refuted?

### Takeaway
No IC speedup claim exists. `claims/` is empty apart from its README, and the three campaign-critical claims (relation collection, blind descent, sub-rho end-to-end) stay `not_attempted`. The ledger's main IC result is a negative of record: the membership-constrained Semaev route costs ℓ^{2.05} (90% CI [1.86, 2.29]) against rho's ℓ^{0.5}. The only charged end-to-end exponent run was inconclusive at about 9% power and cost 1,739–9,789x rho. Several "positive" knowledge findings conflict with proved negatives and look unreliable.

### Cited Findings
**Claims register:**
- `claims/` holds only README.md: "absent evidence means no qualified claim" — [claims/README.md](claims/README.md)
- The claim matrix lists CLM-ECDLP-RELATION-COLLECTION, CLM-ECDLP-BLIND-DESCENT and CLM-ECDLP-SUBRHO-END-TO-END as `not_attempted`, not independently verified. Of about 30 P15xx mechanism claims, most are `open`, P1511/12/13 are `not_reproduced`, and only P1509/P1510 (Hasse-section, global compiler) are `reproduced` — [focus/current_plan.md](focus/current_plan.md). The same file lists "Critical Experiments: none".

**Refuted, negative of record:**
- [ledger/FINDING-PF-IC-001.md](ledger/FINDING-PF-IC-001.md) gives per-PDP d_reg = 2, solve cost ~ d^γ with γ ∈ [4.57, 5.86], and total IC cost ~ ℓ^{(1+γ)/m}, an exponent of 2.05 with CI [1.86, 2.29] at m=3. None of the levers moved it: curve model (Edwards vs Weierstrass, EXP-REP-001/002), isogeny class (EXP-ISO-001), factor-base structure (EXP-FB-001) or arbitrary explicit membership (R6). It is labelled toy scale (p ≤ 2^16, m ≤ 3), "not a proof".
- R6 → reject_scoped: α = 2.074 with CI90 (1.565, 2.621), 2.0e4–1.5e5x above rho — [ledger/DEC-20260718-016.yaml](ledger/DEC-20260718-016.yaml)
- ICI (crossbred / MITM) → reject_scoped: totals 0.891 and 0.583, both above 0.5 — [ledger/DEC-20260718-013.yaml](ledger/DEC-20260718-013.yaml)

**Charged exponent goal and other IC goals:**
- GOAL-ICEX-001: "RUN-ICEX-0ad4d8 … frozen verdict inconclusive (exponent 0.62, CI [0.11, 2.32]) … about 9% power at run noise. Complete cost is 1,739-9,789x rho on all six fixtures; one successful enumerative B0 membership costs 50-306 rho solves … enumerative B0 runs stopped." The successor EXP-ICEX-153c34 was designed with "execution not admitted" — [ledger/goals/GOAL-ICEX-001/](ledger/goals/GOAL-ICEX-001/goal.yaml)
- GOAL-SDEG-001 (Semaev solving-degree law, feeds ICEX) has been on HOLD since 2026-08-18 on unmet verifier-hash residuals — [ledger/goals/GOAL-SDEG-001.yaml](ledger/goals/GOAL-SDEG-001.yaml)
- GOAL-RELN-001 (P(decomp) audit) is in design only, with "no runs" — [ledger/goals/GOAL-RELN-001/goal.yaml](ledger/goals/GOAL-RELN-001/goal.yaml)
- GOAL-ICLIFT-001 is completed: "both routes closed within narrowed scope". A Bose–Chowla B_3 base reached coverage ratio 1.1071, and "it is NOT shown that no factor-base geometry beats a matched random base". GOAL-XEDN-001/002/003 are completed, all finding max|coeff| = 1 Mordell–Weil relations at p ≤ 31 — [ledger/goals/GOAL-ICLIFT-001.yaml](ledger/goals/GOAL-ICLIFT-001.yaml), [ledger/goals/GOAL-XEDN-00*.yaml](ledger/goals/)

**Confirmed or proved structural results:**
- MV_m = (m−1)!·2^{(m−1)(m−2)} = box Bézout for all m — [ledger/DEC-20260722-003.yaml](ledger/DEC-20260722-003.yaml)
- Jet simulation theorem: Shoup's bound lifts — [ledger/DEC-20260722-004.yaml](ledger/DEC-20260722-004.yaml)
- Chord ceiling: B = p^α with α<1/3 implies T3 = 0 w.h.p. — [ledger/DEC-20260722-006.yaml](ledger/DEC-20260722-006.yaml)
- Syzygy residual law, residual = 2n/3 + 1 at n = 12..24 — [ledger/DEC-20260720-001.yaml](ledger/DEC-20260720-001.yaml)
- DREG: n=17 sem rank 125,099 vs prediction 126,922, deficit 1,823, still `inconclusive` / expand — [ledger/DEC-20260718-021.yaml](ledger/DEC-20260718-021.yaml), [ledger/DEC-20260720-002.yaml](ledger/DEC-20260720-002.yaml)
- Semaev cover is E^{m-1}/(diagonal −1) with monodromy (Z/2)^{m-2}, and the factor-base locus is totally split at every m — [knowledge/findings/KN-FIND-a8990a.md](knowledge/findings/KN-FIND-a8990a.md), [KN-FIND-a1f3c2](knowledge/findings/KN-FIND-a1f3c2.md), [KN-FIND-6edecd](knowledge/findings/KN-FIND-6edecd.md)
- Yield conservation: factor-base geometry cannot change mean yield — [KN-FIND-007](knowledge/findings/KN-FIND-007.md)
- Weil bound cannot close the yield gap — [KN-FIND-c93d45](knowledge/findings/KN-FIND-c93d45.md)
- Semaev chain length t = m is optimal — [KN-FIND-9643e7](knowledge/findings/KN-FIND-9643e7.md)
- Nagao Lemma 4 false over F_2 — [KN-FIND-936151](knowledge/findings/KN-FIND-936151.md)
- GGM closure of incidence and endomorphism oracles — [KN-FIND-b7e091](knowledge/findings/KN-FIND-b7e091.md)
- "Betti-Yield duality: H-PSEUDO is the exact condition for sub-rho combinatorial ECDLP" — [KN-FIND-9d2f56](knowledge/findings/KN-FIND-9d2f56.md)

**Suspect or conflicting knowledge entries:**
- [KN-FIND-2a8b7e](knowledge/findings/KN-FIND-2a8b7e.md) ("BKK sparse Semaev gives growing speedup 2^{m-1}·γ_m … ~3-5x at crypto scale", preliminary_empirical) and [KN-FIND-c7d31e](knowledge/findings/KN-FIND-c7d31e.md) ("BKK Speedup Theorem (m+1)/2", proved). Both are a half-factor-base enumeration trick that shares the BKK name. Neither concerns the proved MV = Bézout saturation in DEC-20260722-003.
- [KN-FIND-f8c290](knowledge/findings/KN-FIND-f8c290.md) says the "Semaev m=2 complexity exp(c·sqrt(log N·log log N)) is ACHIEVABLE at toy scale" on the basis of a single p=1009 solve (45/54 factor-base logs). This contradicts FINDING-PF-IC-001.
- C(p) ~ p^{0.055} (KN-FIND-d4f820) conflicts with p^{0.079} (KN-FIND-e7a3b1). [REVIEW-UNADDRESSED-20260807.md](REVIEW-UNADDRESSED-20260807.md) §6 says this is covered by IDEA-20260805-c4f675, "which … argues the null was never subtracted".

**Data integrity:**
- [knowledge/open-problems/KN-OPEN-5b3a08.md](knowledge/open-problems/KN-OPEN-5b3a08.md) records "harness/semaev.py's confirmed s4_expr bug, or build_factor_base's whole-curve (not target-subgroup) sampling". It says no audit of whether EXP-MTIC-001 / EXP-FIB-001 results were affected has been done. `s4_expr` is still present at harness/semaev.py:52.

### Inferences
- The program's verified IC content is almost entirely negative or structural. The positive-sounding findings (2a8b7e, c7d31e, f8c290) are constant-factor or single-instance and should not be cited as leads without re-derivation.
- The highest-value unknown is still the charged exponent. ICEX ran but was underpowered (CI [0.11, 2.32]), so neither "IC beats rho" nor "IC loses" is established at the program's own standard. FINDING-PF-IC-001's 2.05 is for one specific membership-constrained route.

### Gaps
- I did not open the evidence for EXP-MTIC-001 / EXP-FIB-001 to judge the s4_expr impact. That falls under the experiments/outputs remit.

## Q3. What does REVIEW-UNADDRESSED-20260807.md say about IC?

### Takeaway
The review never uses the literal words "index calculus". Its §4.1 "Prime-field ECDLP" register is nonetheless largely an IC gap list, with 13 rows. Its central charge: IC's dominant cost term (the relation matrix / linear algebra) and its most promising axis (arity m ≥ 4) have never been measured, 18 ECDLP contracts were frozen and never run, and several closures rest on a mis-specified object or a one-sided gate. As of 2026-10-14 every one of its IDEA-20260807 proposals is still `status: proposed`, and most of their follow-on hypotheses have no evidence.

### Cited Findings
**Overall shape of the corpus:**
- Claim tiers are toy 188, theory 44, medium 10, crypto 2. Proof status is "certificate 1" — [REVIEW-UNADDRESSED-20260807.md](REVIEW-UNADDRESSED-20260807.md) §2
- "ECDLP core: 18 frozen contracts with zero run directories, including EXP-IC-002 (the yield-charged descent control that DEC-20260727-008 says has 'NO MEASURED REPLACEMENT')", plus EXP-TLD-0a6b4d and EXP-AOM-e057f7 against KN-OPEN-019 — same file §3

**§4.1 gap rows directly about IC:**
- #6: "The relation matrix has never been measured as an object, although FINDING-PF-IC-001 says the linear-algebra stage dominates the cost. The only end-to-end run records one datum ('45/54 FB discrete logs solved'), never followed up." → IDEA-20260807-291968
- #7: "Arity m ≥ 4 has zero measured cost cells anywhere. The cost model is decreasing in m; EV-R6-001 censored every m=4 arm … EXP-ICEX-146ff5 sits in a three-decision routing deadlock." → IDEA-20260807-32f96d
- #3: "Bond tensors were measured DRAMATICALLY low-rank and only ever gated on lossy truncation … χ(m) = C(2^{m−3}+2, 2) exactly … The exact χ(m)-term factorization was never costed as a construction algorithm. KN-OPEN-007's 'Current state' says 'bond ranks near full' — the opposite of the committed evidence." → IDEA-20260807-049ed4
- #11: "KN-OPEN-009's exploitable branch was disposed of in one clause. The monodromy is C_2^{m−2} — smaller than symmetric … No cost model, no measurement." → IDEA-20260807-5126f4
- #1: path-algebra negative "forced before it ran": a one-vertex quiver with no isogeny or correspondence arrows. #2: transfer-operator barrier is "a property of the partition, not the walk", because x-interval cells identify P with −P. #9: a 152-σ EQJ distinct-row collapse was "filed as 'not the kill' and dropped". #10: INCB found EC 15–35% *below* prediction, outside the one-sided gate. #8: every banded statistic is calibrated on one curve instance per cell ("THE SINGLE LARGEST HOLE", a ten-minute fix, orphaned by a batch-label collision). #12: seven barrier claims exist "only as unreviewed decision prose".

**Follow-through as of 2026-10-14** (my trace through ledger/hypotheses and evidence):
- All 13 §4.1 IDEA records still read `status: proposed` — [ledger/proposals/IDEA-20260807-*.yaml](ledger/proposals/)
- 291968 (relation matrix) became H-ICEX-267ebc ("THE DOMINATING TERM OF THE PROGRAM CENTRAL NEGATIVE HAS NEVER BEEN MEASURED") and H-ICEX-760fd5. Both are `proposed` with 0 evidence records — [ledger/hypotheses/H-ICEX-267ebc.yaml](ledger/hypotheses/H-ICEX-267ebc.yaml)
- 32f96d (m ≥ 4) became H-ICEX-5acf1c (`proposed`, 0 EV). It also seeded five H-BINSTD hypotheses and EV-ICPERF-784b25 on the binary side, but there is no prime-field m ≥ 4 cost cell.
- 049ed4 (exact TTN factorisation) became H-TTN-806b84, `proposed`, 0 EV.
- 5126f4 (C_2^{m-2} monodromy exploitation) became H-MONO-6adf3c, which is "refuted_reject_branch_mathematically_unreachable…", and H-MONO-39b1f0, which is "…main_partial_locus_claim_remains_untested…successor_construction_named_not_yet_designed". H-ECDLP-a0d3ea is `proposed`.
- 347a83 (one-curve calibration) became H-ECDLP-fe99ef (`proposed`) and H-ECDLP-75b566 (`specified`). Both have 0 EV.

**Register staleness:**
- KN-OPEN-004 to 010 all still read `status: open` — [knowledge/open-problems/](knowledge/open-problems/). This holds even though the ledger *proved* the BKK, jet and NC results (DEC-20260722-003/004/005) and resolved the monodromy group (KN-FIND-a8990a).
- KN-OPEN-007 still says the tensor-train result was a "SCOPED NEGATIVE" — [KN-OPEN-007.md](knowledge/open-problems/KN-OPEN-007.md)

### Inferences
- Rows #3, #6, #7 and #11 are reviewer objections that were acknowledged (turned into hypotheses) but never executed. They are the clearest "unaddressed" IC objections. Row #8 (n=1 curve per cell) silently scopes several existing negatives, so some "closures" are weaker than the register suggests.

### Gaps
- I did not check whether EXP-IC-002 or EXP-ICEX-146ff5 later gained run directories (experiments/ is out of scope). The other researcher should confirm.

## Q4. Which IC ideas in ideas/ and kb/ (and ledger/proposals) are untested, by cluster?

### Takeaway
`ideas/` holds about 436 numbered ECDLP-IDEA records: 20 active `proposed_unapproved` (none run), 28 deferred pending a theorem, and about 370 rejected. These are dominated by "source-lift / source-router" variants of direct-source relation finding: 235 of 418 filenames contain "source". The much larger live pool is `ledger/proposals/`, with 1,072 IC-related proposals, 892 of them `proposed`. Only about 20–30% in any cluster ever reached an evidence record. Linear algebra, lifting/xedni, multi-target and ML clusters are the least executed. `kb/` is infrastructure code with no idea content.

### Cited Findings
**ideas/ corpus:**
- The registry has 20 active ideas, all "proposed, unapproved … novelty_unverified" — [ideas/idea_registry.tsv](ideas/idea_registry.tsv), [ideas/README.md](ideas/README.md). IC-relevant actives:
  - IDEA-002 split-Jacobian projected smoothness (contract 002)
  - IDEA-006 elliptic-net short annihilator
  - IDEA-007 Miller S-unit descent
  - IDEA-009 nonequivariant trace-zero transfer
  - IDEA-010 torsor deck-orbit descent
  - IDEA-012 aggregate complement-divisor compression
  - IDEA-056 block-Krylov transition-intersection extractor
  - IDEA-059 Cremona-shrunk toric decomposition
  - IDEA-434 isogeny-class variation of Semaev yield
  - IDEA-435 function-field lift target-coefficient audit
  - IDEA-436 torsion-lift valuation profile
- README states: "The corpus is 15 active proposed/unapproved, 12 deferred, and 69 rejected records, not 96 validated algorithms" (later rounds grew it) — [ideas/README.md](ideas/README.md)
- 28 deferred records each name a missing theorem. Examples: 049 bounded-root transducer, 069 Massey secondary relation, 098 squarefree Gröbner shelling, 102 elliptic dynamical-R transfer, 103 tropical source atlas, 119 elliptic Bloch-polylog regulator descent, 121 shared bivariate common norm, 133 nonlinear apolar flat extension ("structured/adaptive routes remain open") — [ideas/deferred/README.md](ideas/deferred/README.md)
- The rejected records keep exact no-go boundaries. Example: IDEA-001 incidence flattening gives max(T, M) ≥ N^{0.6} — [ideas/rejected/README.md](ideas/rejected/README.md)

**ledger/proposals, 1,072 IC-matching records.** Status: 892 `proposed`, 164 unparsed, 2 `proved`, 1 withdrawn. The top question IDs are RQ-ECDLP-002 (131), RQ-BINSTD (82), RQ-JINV (39), RQ-FROB (35), RQ-SDEG-001 (25), RQ-MONO-001 (25), RQ-RELN-001 (21), RQ-ICINV (20), RQ-ICEX-001 (17) and RQ-FBG-001 (13). Cluster counts from regexes over titles, with the number reaching evidence directly or through an evidenced hypothesis:

| Cluster | Proposals | Reached evidence |
|---|---|---|
| Factor-base geometry / P(decomp) | 207 | 54 |
| Cost accounting / charging / rho baseline | 198 | 53 |
| Symmetry / monodromy / Galois | 196 | 50 |
| Summation polynomial / Gröbner / degree | 185 | 64 |
| Transfer / cover / descent | 121 | 30 |
| Linear algebra / relation matrix | 116 | 23 |
| Binary / extension fields | 89 | 47 |
| Lifting / xedni / function field | 79 | 12 |
| Multi-target / preprocessing | 34 | 5 |
| Large prime / sieve | 29 | 8 |
| Hardware / GPU / SAT | 26 | 15 |
| ML / learning | 13 | 4 |

(My regex classification, run over [ledger/proposals/](ledger/proposals/) and [ledger/evidence/](ledger/evidence/).)

**Hypothesis families.** Each line gives count, how many have evidence, and the status mix:
- H-ICEX: 62, 2 with evidence. 52 `proposed`.
- H-DREG: 30, 3. 18 `approved`.
- H-XEDN: 23, 5.
- H-RELN: 22, 5.
- H-ICINV: 19, 3. 14 `specified`.
- H-SATIC: 13, 0.
- H-SDEG: 10, 1.
- H-PFDR: 52, 9. 28 `proposed`.
- H-ALPF: 32, 1. All `analyzed`.
- Better-executed families: H-MONO 33/19, H-SEMBIN 26/16, H-BINSTD 148/48 (100 `approved`), H-CERTBIN 31/10.

Source: [ledger/hypotheses/](ledger/hypotheses/).

**Untested IC ideas by requested cluster** (all `proposed` with 0 evidence unless stated):
- *Summation polynomials / degree:*
  - H-ICEX-3bccf3: solving degree vs d_reg vs first fall on S_3
  - H-ICEX-357d9d: exponent stability across monomial orders
  - H-ICEX-604b5f: singular locus of the m=3 Semaev hypersurface = 2-torsion, (0, 6)
  - H-ICEX-14a46a: n=15, D=6 Macaulay memory wall is a representation cost
  - KN-OPEN-02200b: H-TOP top-form conjecture at m ≥ 5
  - KN-OPEN-d218ec: where Semaev's Assumption 1 fails
  - G30-006: prime-field first-fall census; SDEG is on hold
- *Factor bases:*
  - H-ICEX-12d125: P(decomp | target) concentration
  - H-ICEX-c273c2: decomposition probability as an exact count
  - H-ICEX-ab1da2: "THIS REPOSITORY HOLDS ZERO MEASURED RANK VALUES FOR ANY FACTOR-BASE GEOMETRY"
  - H-ICEX-87ad66: LP-1 large-prime rate S_m/S_{m−1} = (B+m−1)/m
  - KN-OPEN-020: universal algebraic-factor-base no-go, open
  - Bose–Chowla B_3 base at 1.1071 coverage, not followed up (GOAL-ICLIFT-001)
- *Symmetry:*
  - H-ICEX-70ddd4: the only reproducible endomorphism effect is the automorphism-quotient r = 3; "linear-algebra claim is open in both directions"
  - KN-OPEN-025: automorphism-order speedups of TNFS linear algebra
  - The 2-torsion halvability (h+, h−) structure at m=4 (KN-FIND-19f5ea, KN-FIND-edd62c) is established but not costed
- *Transfer:*
  - G30-015 / IDEA-002 / CLM-P1546 split-Jacobian projected smoothness: `open`, never run
  - IDEA-009 trace-zero transfer
  - H-ICEX-50e403 Kummer/Picard-bridge charged transfer IC: `approved`, no EV
  - KN-OPEN-ed922d: an advantage hidden behind a subgroup-preserving correspondence
  - KN-OPEN-cc1988: non-isogeny-invariant weakness
- *Lifting / xedni:*
  - B2 function-field xedni: only coefficient controls were done (XEDN-001..003). The P_lift trend gate over three sizes was never run.
  - IDEA-005 / CLM-P1543 height-compressing global lift: `open`
  - B3 Coppersmith windowed relations: deferred, never staffed
- *ML-guided:*
  - ROADMAP Phase 4.3 Buchberger pair-selection learning (heuristic, imitation, RL)
  - Only 13 IC proposals, 4 reached evidence
  - KN-OPEN-7f0d85 and KN-OPEN-875d43 concern learned feature search, not IC
- *Linear algebra:*
  - H-ICEX-267ebc: relation matrix as an object
  - H-ICEX-5666b2: "THE SECOND HALF OF EVERY INDEX-CALCULUS EXPONENT IS CITED, NOT MEASURED", block Wiedemann/Lanczos on real relation matrices
  - H-ICEX-fe2944: LA share of cost across a fourfold |V| range
  - H-ICEX-c35cf3: charge LA memory
  - G30-016/017: an LA benchmark module and displacement-rank supports net of penalty
  - IDEA-056: block-Krylov extractor
- *Hardware / SAT:*
  - GOAL-SATIC-c49b77 (SAT on summation-polynomial systems) is `partially_cleared`, with 0 of 13 H-SATIC hypotheses having evidence and two tasks blocked — [ledger/goals/GOAL-SATIC-c49b77.yaml](ledger/goals/GOAL-SATIC-c49b77.yaml)
  - No GPU/FPGA IC idea reached evidence
- *Cost and baseline:*
  - H-ICEX-835c4f: "The program has no measured Pollard-rho baseline"
  - H-ICEX-9d75dd / a2e95b: charging of right-censored attempts
  - H-ICEX-5f425f: a four-rung ladder cannot identify the exponent statistically
  - H-ICEX-b8308b: zero-compute go/no-go with the most favourable inputs
  - H-ICEX-71bce3 / 9e54c2: closed-form arity threshold

### Inferences
- The ideas/ corpus shows heavy mode collapse. The "source"-centric direct-source family was generated hundreds of times and rejected by exact no-go boundaries. It is a poor place to look for leads.
- The live untested backlog is concentrated in ICEX: measurement and accounting hypotheses written mostly as critiques of the program's own cost model, cheap and decisive but never dispatched.

### Gaps
- The cluster counts come from title regexes, so they overlap and are approximate. 164 proposals had unparsable status fields.
- I did not read the 20 active ideas' contracts in full.

## Q5. Which ideas recur many times without execution?

### Takeaway
Three IC themes recur across four or more independent documents yet still have no evidence: (1) measuring the relation matrix / linear-algebra stage, (2) arity m ≥ 4 cost cells, and (3) split-Jacobian / cover-based projected smoothness. A fourth, monodromy-based relation rates, recurred most of all and was executed. It was resolved mathematically as (Z/2)^{m-2} with totally split factor-base fibres, but GOAL-MONO-001 is still `active` under unresolved reviewer dissents.

### Cited Findings
**Relation-matrix / LA measurement.** It appears in:
- 0717 A3 (displacement rank; executed but rejected on supply)
- G30-016 and G30-017 — [research_goals_20260723.md](research_goals_20260723.md)
- review gap #6 — [REVIEW-UNADDRESSED-20260807.md](REVIEW-UNADDRESSED-20260807.md)
- IDEA-20260807-291968
- H-ICEX-267ebc, H-ICEX-5666b2, H-ICEX-fe2944, H-ICEX-ab1da2, H-ICEX-c35cf3 — [ledger/hypotheses/](ledger/hypotheses/)
- the red team's "PO68/PO71 closures show linear algebra and descent dominate complete cost" — [research_directions_20260718.md](research_directions_20260718.md) §6.4

None of these has an evidence record. The only datum is "45/54 FB discrete logs solved" (KN-FIND-f8c290).

**Arity m ≥ 4.** It appears in:
- G30-006 (m ∈ {2,3,4}) and G30-026 (m=3/m=4 plus a 40–64-bit ladder)
- review gap #7
- IDEA-20260807-32f96d
- H-ICEX-5acf1c, H-ICEX-372bec ("arities 3–5", `specified`), H-ICEX-9e54c2 (arity threshold, `specified`), H-ICEX-d8e858 (m ∈ {5,6,8} batched sieve, `approved`, no EV)

EV-R6-001 censored every m=4 arm.

**Split-Jacobian / genus-2 projected smoothness.** It appears in:
- ideas/ECDLP-IDEA-002 (registry rank 1, with contract)
- G30-015 ("the program needs at least one relation engine whose relations are not x-line Semaev summations")
- claim CLM-P1546 (`open`) — [focus/current_plan.md](focus/current_plan.md)

No evidence record exists.

**Function-field xedni / lifting.** It appears in:
- 0718 B2 (round-2 winner)
- the G30 deferral (c)
- GOAL-ICLIFT-001 and GOAL-XEDN-001..003 (coefficient controls only)
- IDEA-005, IDEA-435, CLM-P1543
- the 0831 "do not fund" list

The P_lift trend gate was never run.

**Monodromy / Chebotarev.** It appears in:
- 0718 A1, C2 and D2; G30-007/008/009/012
- GOAL-MONO-001 (33 hypotheses, 19 with EV, 28 EV-MONO records)
- KN-OPEN-009, review gap #11
- KN-FIND-a8990a / a1f3c2 / c41ea9 / 6edecd, and the paper [papers/semaev-conservation-specialization/paper.tex](papers/semaev-conservation-specialization/paper.tex)

GOAL-MONO-001's latest next action (2026-10-06) still says "NO completion criterion of this goal has been adopted as met. The J2, J3 and J4 DISSENTs … persist" — [ledger/goals/GOAL-MONO-001/goal.yaml](ledger/goals/GOAL-MONO-001/goal.yaml)

**Matched rho baseline and cost charging.** H-ICEX-835c4f ("no measured Pollard-rho baseline"), H-ICEX-9d75dd, H-ICEX-a2e95b, H-ICEX-4ccbaa and H-ICEX-5f425f all argue that the ICEX exponent is not identifiable. GOAL-ICEX-001's 9%-power verdict confirms this — [ledger/goals/GOAL-ICEX-001/goal.yaml](ledger/goals/GOAL-ICEX-001/goal.yaml)

**Exact tensor factorisation of S_m.** It appears in:
- 0717 C2
- EV-TTN-001/002 (χ law certified)
- review gap #3, IDEA-20260807-049ed4
- H-TTN-806b84 (`proposed`, 0 EV)

It was never costed as a construction algorithm against the dense-resultant exponent of about 1.979.

**High-ranked round-2 candidates never executed:**
- A3 S-unit / Mason–Stothers relation budget (score 25, third-highest overall)
- B1 AG-code enumerator (23; G30-014 drafted)
- B3 Coppersmith windowed relations (22; deferred)

### Inferences
- Recurrence here tracks *self-criticism of the measurement apparatus* (LA term, m ≥ 4, rho baseline, censoring) more than new mechanisms. These are cheap, decisive measurements that the governance machinery (holds, quorum dissents, routing deadlocks) has kept from running. They are the most defensible "next steps".
- The monodromy line is effectively answered: (Z/2)^{m-2}, and the factor-base locus is totally split. Further monodromy census work has low expected value. The remaining open sub-question is the "partial locus" claim of H-MONO-39b1f0.
- The tensor χ(m) ~ d²/2 law against the (d+1)^m dense array is the one structural *positive* (an exponential size gap) that nobody has costed. It is the strongest unexecuted mechanism-level lead in this corpus, with the caveat that construction and specialisation costs may erase it.

### Gaps
- Recurrence counts are from my targeted greps across direction docs, goals, hypotheses and the review. I did not do a full-corpus dedup.
- I could not verify whether experiments/ contains unrecorded runs for split-Jacobian (G30-015) or m ≥ 4. The other researcher should cross-check.
