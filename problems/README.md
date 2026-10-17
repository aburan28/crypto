# Research problems and obstructions

Track the major named blockers that prevent a research claim, proof step,
implementation, experiment, or audit conclusion from being established.
The [agent rule](../AGENTS.md#name-the-obstruction) requires each record to
explain the gap and what would let the work move past it.

## Index

Keep resolved and superseded records in the index.

| Problem | Status | What it blocks | What would unblock it |
| --- | --- | --- | --- |
| [Prime61 exhaustive horizontal class memory](prime61-horizontal-class-memory.md) | open; disk-backed pilot verified | certified full-class coverage and all-class workload execution; retained-vector requirement removed on the new path | certified degree-set coverage, provisioned disk/compute envelope and replayed remaining workloads |
| [Screenshot candidate driver is not located](screenshot-factor-base-driver-unlocated.md) | open | faithful 1,000/2,000-candidate isogeny benchmark at a requested large field | original file/branch/PR or source bundle, native fixture replay and validated larger-field map |
| [Generalized GPU rho workers for nine Certicom binary challenges](certicom-binary-gpu-rho-integration.md) | open | persistent, independently replayable GPU rho and qualified fastest-arithmetic choices for all nine exact challenges | exact CPU bindings, per-curve GPU replay and durable DP/collision path, then matched device performance gates |
| [Genus-3 IC-2LP study 04 missing](ic-improve-04-genus3-2lp-study-missing.md) | investigating (04 has since landed on `main`) | pairing 04b's arms against the original IC-2LP arm (same instances, seeds, unit); 04b pairs against a rebuilt baseline instead | the original 04 tool and records at a pinned revision, then a matched baseline comparison |
| [p61 selected-class census and adjacency completion](p61-full-isogeny-class-census.md) | investigating | complete selected-class graph and its use as an admitted large-field corpus | exact 358,094,304-vertex ledger, verified selected-degree adjacency, and adequate durable storage |
| [SEA Rust-to-group refinement proof](sea-rust-group-refinement-proof.md) | open | complete formal verification of the Rust SEA/counting implementation | connect concrete field/point operations and congruences to the checked abstract search and uniqueness proofs |
| [cairn certificate checkers are Python-only](cairn-certificate-checker-python-only.md) | open | objective mode for backlog cards with native checkers (the upgrade path of `XP-20261015-p192-isogeny-aux-search-t48`, which runs as lab tasks meanwhile) | cairn accepting a pinned native checker, or the `workspace` verifier running a zero-dependency Rust checker under its jail |
| [ECC arithmetic cross-platform performance](ecc-arithmetic-platform-measurement.md) | open | controlled multi-platform and worldwide arithmetic comparisons | identified isolated ARM/x86 hosts with A/A controls, validated backend implementations and matched comparators; new X25519 Modal bundle awaits explicit source-export approval |
| [ECC arithmetic PR merge gates](ecc-arithmetic-pr-merge-gates.md) | investigating | delivery of ECC arithmetic PR #1798 with a passing current-head check suite | resolve the recorded registry, evidence, lint and workflow failures; obtain final-head review and checks before a merge commit |
| [Root-finding pipeline complete-cost calibration](rootfinding-complete-cost-calibration.md) | open | qualified full-cost or wall-time comparison for the assessed root-finding IC pipeline; 4/8/16-lane, B train/holdout and coverage panels are verified and replayed; a fresh-target equal-prefix component control exists, while native field counts are absent from sealed records and wall data are L0 | serialize and price native work; repeat matched public targets under L2 isolation and on an independent host |
| [N83 cloud online A/A noise](n83-online-aa-noise.md) | investigating | selecting the exact-root factor-base finalist for the five sealed Tier 1 points after two K2543/root22 controls exceeded 10% | frozen mirrored K3072/K2543 adjacent A/A screen with exact base, relation and scalar replay, then an exact IC1 selection record |
| [EC Bench GPU matched workload and validation](ecbench-gpu-matched-workload-contract.md) | open | admitted CPU/GPU comparison of the same ECDLP workload | one RunPod device solve audited; qualify GPU traffic telemetry, interleave full-cost CPU/GPU trials with A/A controls, and extend curve coverage |
| [N83 recorded-query sampler contract](n83-recorded-query-sampler-contract-mismatch.md) | open | green whole-library validation on the pinned baseline | traced and versioned query law, point checks and passing regression |
| [GHS → Smith characteristic mismatch](ghs-to-smith-pipeline-characteristic-mismatch.md) | open | the GHS → Smith → Diem pipeline (char 2 as stated: blocked; odd-char n = 3 descent onto genus-3 hyperelliptic curves exists, 2026-10-14) | a char-2 hyperelliptic→quartic isogeny, or Smith Props. 2/4 passing on the odd-char family (cycle-type screen done: #S(H) = 3 or 1) |
| [Smith 18.57% walk extension](smith-kernel-isogeny-class-walk.md) | open | costing walks that lift Smith's coverage toward the measured class ceiling (43.0% at `p = 7`, 47.5% at `p = 11`; Prop. 2/4 tests landed 2026-10-14) | explicit isogenies between hyperelliptic class members and their walk length and cost; the empty-class deficit is not iso-class sparsity (2026-10-15); the Prop. 4 sign follows #H(F_p), explaining ~3/4 of it (2026-10-16); that bias shrinks like 1/√p (2026-10-16) nor a `J[2]`/2-Sylow/`J[4]` invariant of `J_H` (2026-10-16); next: quotient-curve classes |
| [No class arithmetic on the Jacobian of a plane quartic](plane-quartic-jacobian-class-arithmetic-missing.md) | resolved (toy smokes) | executing the transfer homomorphism, relation collection and any cost comparison on plane-quartic Jacobians (blocked until the linear-equivalence routine landed, #1768) | done for the toy discriminating checks: `quartic_ic.rs` linear-equivalence checks and a pure IC target path; deg-3 canonical forms and any cost statement remain open |
| [Binary decomposition systems capped at 64 Boolean unknowns](binary-boolean-system-64-variable-cap.md) | open | Gröbner/SAT decomposition arms on chain systems with `md + (m−2)n > 64` unknowns (the `m = 4`, n ≥ 23 cells of `research/ic_m3_algebraic_vs_mitm_20261014/`); the study itself was completed with the `u128` `wide_groebner` solver (amendment 1) | a ≥128-variable monomial type in `build_decomposition_system`, inherited F4 and the SAT encoder, agreeing with the floor on 12 targets at `m = 4`, n = 23 |

## Add or update a problem

1. Check this index and the affected work for an existing record. Reuse it when
   the obstruction is the same; link the canonical record across repositories
   rather than duplicating its evidence.
2. Copy [TEMPLATE.md](TEMPLATE.md) to a stable, descriptive filename such as
   `subgroup-preserving-transfer-not-established.md`. Name the mechanism or
   missing obligation, not a vague "hard problem". The filename above is an
   illustration, not a recorded finding.
3. State the exact blocked task or claim and its parameters, assumptions, and
   acceptance criteria. Separate evidence from conjecture and cite durable
   artifacts, proofs, reproductions, or source records.
4. Specify the evidence or capability needed to move past the obstruction and
   one concrete next experiment or proof obligation that can distinguish
   explanations. Link an existing execution protocol where applicable; this
   record does not replace it or impose a new campaign spending cap.
5. Add the index row and link the problem from the affected note, plan, audit
   finding, or PR. Update these pointers and the record together as work lands.

Use this directory for consequential or recurring blockers. Routine defects
and TODOs can stay in their existing trackers. Record a blocker even when its
cause is unknown; explicit unknowns are useful, invented explanations are not.

## Status and evidence

- `open`: the obstruction remains and its next action is recorded.
- `investigating`: work on that action is underway; link the task or run.
- `resolved`: linked evidence meets the stated exit criteria for this specific
  obstruction. State any remaining limitations.
- `superseded`: a linked successor now tracks this obstruction; this is not a
  claim that it was solved.

These statuses describe the problem record, not the scientific status of a
hypothesis or the completion of a campaign. Stopping work does not resolve
a problem or authorize abandoning the user's requirement.

Keep observations, calculations, hypotheses, resource/access/assistance gaps,
and proved limits distinct. Bound a negative result by its exact inputs and
search limits. A resource ceiling, missing implementation, or failed search
does not establish a mathematical obstruction or a security theorem.

Preserve historical evidence and append dated updates explaining corrections.
A resolution cites its evidence and the criteria met; a supersession cites its
successor and the reason for replacing the record. Follow the repository's
existing review, archival, and claim-promotion rules.
