# Research Program Index

A meta-document tracking the eight research programs (v1-v9) for the
AutoLab ECDLP cryptanalysis project. Each program is described in its
own `research_program_v*.md` file.

## At-a-glance

| Program | Focus | Status | Result |
|---------|-------|--------|--------|
| v1 | Initial exploration (Phase 1-5) | COMPLETED | All negative |
| v2 | Algebraic structure (Phase 6-11) | COMPLETED | All negative |
| v3 | Frobenius / isogeny / partition (Phase 12-17) | COMPLETED | 5/6 negative; partition ℓ=16 gives 1.4× |
| v4 | Groebner / Semaev deep dive | COMPLETED | All negative; Yokoyama bound confirmed |
| v5 | Non-naive attempts #1 | COMPLETED | All 4 directions ruled out for LMFDB |
| v6 | Non-naive attempts #2 + C engineering | COMPLETED | 5 directions closed; C extension 23-86× |
| **v7** | **Close C-extension gap to 80-bit feasibility** | **PLANNED** | Target: ≥27e6 ops/sec on 2 CPUs |
| **v8** | **Algorithmic frontiers** | **PLANNED** | F5 Semaev, Drinfeld, tropical, hyperelliptic |
| **v9** | **Long-horizon roadmap** | **PLANNED** | 1-5 year speculative agenda |

## Recommended sequencing

For a single engineer / researcher with limited time:

1. **v7** first (3 weeks of engineering)
   - High concrete impact, low research risk
   - Brings 80-bit benchmark within reach
   - Validates infrastructure for future programs

2. **v9 Phase 40** in parallel with v7 (publication writing)
   - Converts v1-v6 results into a paper
   - Establishes priority before someone else duplicates

3. **v8 Phase 23 (F5 Semaev)** after v7 ships
   - Highest theoretical upside
   - But large implementation effort (4-8 weeks)

4. **v9 Phase 35 (curve-generation pipeline)** later
   - Useful contribution to cryptographic practice
   - Less time-sensitive

5. Defer **v8 Phases 24-29** unless something promising emerges

## Cross-cutting dependencies

- v7 → v9 Phase 39 (distributed rho): v7's C extension is the
  per-machine building block
- v8 Phase 23 → publication: F5-specialized result could become a
  standalone paper
- v6 results → v9 Phase 40 (publication): the empirical map is the
  paper's main contribution
- v7 Phase 22.1 (Barrett) → general C extension: useful beyond
  this project

## What ALL programs share

- High-fidelity methodology: 50+ trials, 95% CIs, full cost accounting
- Reproducibility: Sage scripts with fixed seeds, C with reproducible
  builds
- Honest reporting: negative results documented with the same rigor as
  positive ones
- Literature anchoring: every claim cross-referenced against published
  bounds (Yokoyama, Diem, Petit-Quisquater, Huang et al.)

## What ALL programs avoid

- Verifier exploits (one-shot game-the-test approaches)
- Unvalidated extrapolation (raw ops/sec without end-to-end validation)
- Speculation without literature anchor (every direction must have at
  least one published paper to reference)

## Reading order for newcomers

1. `final_synthesis.md` — start here; ties together v1-v6
2. `research_program_v6_result.md` — most recent honest result
3. `yokoyama_lower_bound.md` — theoretical floor
4. `research_program_v7.md` and `research_program_v7_execution.md` —
   what's next, concretely
5. `research_program_v8.md` — what's next, ambitiously
6. `research_program_v9.md` — what's next, over the horizon
7. `research_program_v[1-3].md` — historical context if needed

## How to add a new program

1. Create `research_program_vN.md` following the v6/v7 template:
   - Motivation tied to gaps in prior programs
   - Numbered phases with concrete experiments
   - Success criteria and risk assessment
   - Estimated effort

2. Add to this index table

3. For execution-ready programs (like v7), also write
   `research_program_vN_execution.md` with day-by-day plan

4. If a program delivers, write `research_program_vN_result.md` with
   findings

## Iteration philosophy

These programs are **plans**, not commitments. They will be revised:
- When prior assumptions turn out wrong (e.g., raw ops/sec → end-to-end
  scaling correction in v6)
- When new literature appears (the cryptanalysis field is active)
- When unexpected results suggest new directions

Each revision is a separate commit; the history is the audit trail.
