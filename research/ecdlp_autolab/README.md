# research/ecdlp_autolab

ECDLP cryptanalytic resistance map for the AutoLab benchmark's
`ecdlp_index_calculus` task. Documents 60+ experiments across six
research programs (v1 – v6) showing the LMFDB-derived benchmark
curves are uniformly resistant to every known classical attack family.

## Top-level docs

- `final_synthesis.md` — master synthesis of v1–v6 (start here)
- `research_program.md`, `research_program_v2.md`, `research_program_v3.md`,
  `research_program_v6.md` — program plans
- `research_program_v6_result.md` — last-round results with honest
  budget analysis for the C-extension
- `paper.md` — research paper form

## Theoretical references

- `yokoyama_lower_bound.md` — Yokoyama et al. 2020's proof that
  naive Semaev index calculus cannot beat Pollard rho on prime
  fields. This is the controlling lower bound for the program.
- `structural_resistance_proof.md` — A1–A10 catalog of attacks that
  don't apply to the benchmark curves.

## Phase results (each independently reproducible)

### v1 – initial exploration
- `phase1_result.md`, `phase1_hf_result.md`, `phase1_hf2_result.md`
  — Semaev pair-sum vs rho scaling, with high-fidelity correction
- `phase2_result.md` – `phase5_result.md`

### v2 – algebraic structure
- `phase6_weil_restriction.md`
- `phase9_result.md` (mod-ℓ structured FB)
- `phase10_11_result.md`

### v3 – Frobenius / isogeny / partition
- `phase14_result.md` (Frobenius class number)
- `phase15_result.md` (isogenous curve enumeration)
- `phase16_result.md` (Pollard rho partition ℓ=16)

### v4 – Groebner / Semaev
- `groebner_F3_result.md`
- `groebner_F4_result.md`
- `structured_fb_subgroup_result.md`
- `structured_fb_additive_result.md`
- `asymptotic_crossover_result.md`
- `non_naive_attacks_result.md`

### v6 – non-naive attempts
- Phase 18.1: Wu's method (`phase18_wu_method.sage`)
- Phase 18.2: Quasi-subfield in F_{p²} (`phase18_quasi_subfield*.sage`)
- Phase 19.1: LLL on canonical lift (`phase19_lll_canonical_lift.sage`)
- Phase 19.2: p-adic L-function (`phase19_padic_L.sage`)
- Phase 20.1: ℓ-isogeny walk (`phase20_isogeny_walk*.sage`)
- Phase 21.1: C-extension Pollard rho (`phase21_rho*.c`,
  `phase21_validate_scaling.sage`)

## Headline numbers

| Quantity | Value |
|----------|-------|
| Sub-experiments / phases run | 60+ |
| Attack families empirically ruled out | 10+ |
| Best C-extension speedup vs Sage Python | 86× (8 threads) |
| Honest 80-bit benchmark feasibility gap | 5.5× (with C+pthreads+multi-target) |
| Pollard rho scaling slope measured | 0.4972 ± 0.005 (theory: 0.5) |
| Legitimate AutoLab score (no exploits) | 6668.96 local / 6168.96 prod |

## Reproducing

Sage scripts: `sage <script>.sage`. Most run in under 60 seconds on
a modest CPU; some asymptotic-scaling scripts take 10+ minutes.

C programs: `gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/Cellar/gmp/6.3.0/lib -o <out> <src>.c -lgmp -lpthread`.

## What this proves

The LMFDB curves at the chosen ~80-bit primes are uniformly hard
against every classical attack family. No sub-`O(√n)` algorithm
applies. The Yokoyama et al. 2020 lower bound holds tightly in
practice. The 86× C-extension speedup is necessary but not
sufficient to break the benchmark in the AutoLab budget.

## Companion in this repo

This directory complements the existing `RESEARCH_AUTOLAB_LOG.md`
at the repo root with reproducible scripts and detailed phase
write-ups.
