# Research Program v11: Publication and dissemination

**Date:** 2026-05-28
**Status:** plan; executable now without waiting for v7/v8 results
**Motivation:** The v1-v6 results are publishable AS IS. v7-v10 are
future work. v11 captures *how to communicate what we already have* —
which is the highest-leverage action available, regardless of whether
later programs succeed.

## Why publish now

The v1-v6 resistance map has several elements that would be novel
contributions to the literature:

1. **First published measurement** of the Diem L(2/3) constant at
   13-28 bits with controlled methodology
2. **Empirical confirmation** of Yokoyama et al.'s lower bound across
   diverse structured factor bases (the original paper validated only
   to 25 bits with one factor-base family)
3. **C-extension engineering benchmark** for prime-field Pollard rho
   at 81-bit (with reproducible scaling validation)
4. **Catalog of negative results** for structured FB attacks
   (multiplicative subgroup, arithmetic progression, parametric GB,
   quasi-subfield, isogeny graph, LLL canonical lift, p-adic L)

None of these is a breakthrough algorithm, but each is methodologically
contributable and would save future researchers from re-discovering
the same negatives.

## Phase 41: Conference paper (J. Math. Cryptol. / CRYPTO submission)

### Target venue

**Primary:** Journal of Mathematical Cryptology (open access, fits the
Yokoyama and Huang papers we cite)

**Backup:** CRYPTO short paper / poster track; ASIACRYPT empirical
methodology track; CHES if focus shifts to engineering.

### Paper structure (draft)

```
Title: An Empirical Resistance Map for Prime-Field ECDLP on LMFDB Curves

§1 Introduction
  - The ECDLP cryptanalysis landscape
  - Yokoyama et al. 2020's lower bound for naive index calculus
  - Our contribution: empirical confirmation + extended testing

§2 Methodology
  - High-fidelity protocol (50 trials, 95% CI, full cost accounting)
  - Multi-curve sampling (4 LMFDB curves × 6 primes)
  - End-to-end scaling validation (not just throughput)

§3 Negative results: structured factor bases
  3.1 Multiplicative subgroup FB (X^m - 1 constraint)
  3.2 Arithmetic progression FB (falling factorial)
  3.3 Mod-ℓ partition FB
  3.4 Quasi-subfield in artificial F_p^2 embedding
  3.5 Common explanation via Yokoyama regularity bound

§4 Negative results: non-GB attacks
  4.1 LLL on canonical lift coordinates
  4.2 p-adic L-function exploration
  4.3 ℓ-isogeny graph navigation
  4.4 Wu's characteristic-set method

§5 Engineering: C-extension Pollard rho
  5.1 Throughput benchmarks (23-86× over Sage Python)
  5.2 End-to-end scaling validation (slope 0.4972 ≈ 0.5)
  5.3 Honest 80-bit feasibility analysis (5.5× gap)
  5.4 Path to closing the gap

§6 Discussion: what remains genuinely open
  - Non-Buchberger algebraic techniques (F5 specialized)
  - Truly novel mathematical frameworks
  - Engineering-only path to feasibility

§7 Conclusion: the map and its implications
  - LMFDB curves are uniformly hard against classical attacks
  - Yokoyama bound is tight in practice
  - Future research must escape the naive-IC framework

Appendix A: reproducible Sage scripts
Appendix B: C source code listings
Appendix C: data tables with all measurements
```

### Length budget

J. Math. Cryptol. allows 25-40 pages. Allocation:
- §1 Intro: 3 pages
- §2 Methodology: 3 pages
- §3 Structured FB negatives: 6 pages (with tables)
- §4 Non-GB negatives: 5 pages
- §5 Engineering: 6 pages (with measurements)
- §6 Discussion: 3 pages
- §7 Conclusion: 1 page
- Appendices: 5-10 pages
- **Total: ~32-37 pages**

### Schedule

| Week | Activity |
|-----:|----------|
| 1 | Convert `final_synthesis.md` into §6-7 |
| 2 | §3 negative results, tables of measurements |
| 3 | §4 non-GB negatives, tables |
| 4 | §5 engineering benchmarks, scaling validation |
| 5 | §1 intro + §2 methodology |
| 6 | Appendices, reference cleanup, internal review |
| 7 | External review (request from a cryptanalyst) |
| 8 | Revisions, submission |

Total: **~2 months** of focused writing.

## Phase 42: Open-science companion

### GitHub release

Tag the `research/ecdlp_autolab/` directory as `v1.0` release:
- All Sage scripts, C source
- Pre-built binaries (Linux x86-64, macOS ARM, Linux ARM)
- Docker container with Sage + GMP + reproducible environment
- Detailed BENCHMARKS.md with hardware specs

### Replication challenge

Issue a public challenge: "Reproduce our v6 measurements on your
hardware". Standardize the measurement protocol so others can confirm.

### Dataset deposit

Submit the measurement dataset (all phase results, raw trial outputs)
to Zenodo for permanent citable archive.

## Phase 43: Talks and visibility

### Workshops to target

- IACR ECC workshop (annual; "elliptic curve cryptography" community)
- Heidelberg Lange's number-theory algorithm conferences
- Stanford Applied Crypto Group seminar
- ENS Crypto group talks

### Blog post sequence

1. "Why Yokoyama et al.'s lower bound matters for ECDLP" (theoretical)
2. "We tested 10 structured factor base attacks and 9 are dead"
   (empirical map)
3. "How fast is Pollard rho really? An engineering study" (Phase 21)
4. "What's actually open in prime-field ECDLP cryptanalysis"

### Twitter/X thread

Punchy 10-tweet summary with the headline numbers:
- 60+ experiments
- 0 algorithmic breakthroughs
- 86× engineering speedup
- 5.5× still short of 80-bit feasibility

## Phase 44: Tutorials and educational

### Annotated Sage notebooks

Convert key scripts to Jupyter/Sage notebooks with detailed
mathematical exposition. Target audience: graduate students entering
cryptanalysis.

### "How to do high-fidelity cryptanalysis" guide

Document the methodology: 50 trials, full cost accounting, multi-curve,
end-to-end validation. This guide is independently valuable from any
specific result.

## Phase 45: Continuous benchmarking infrastructure

### CI for the resistance map

Set up GitHub Actions to:
- Run a subset of phase scripts on every commit
- Detect any regression (slower per-query, wrong relations, etc.)
- Compare against canonical baseline benchmarks

This is engineering, but unlocks future progress: if you can't measure
it, you can't improve it.

## Phase 46: Engagement with adjacent projects

### Direct contact

- **Sage developers**: contribute optimized rho implementation upstream
- **LMFDB**: cross-reference our curves with their database
- **Yokoyama et al.**: share our v6 confirmation of their bound
- **Huang et al.**: share our v6 quasi-subfield findings
- **Faugère group**: share F5-Semaev results once v8 Phase 23 lands

### Standards bodies (long-term)

If a meaningful algorithmic advance emerges (unlikely from v6, possible
from v8), engage:
- ANSSI / BSI for European standards
- NIST PQC team for post-quantum perspectives
- IRTF CFRG for IETF standardization

## Priority order

For ~6-month plan:
1. **Phase 41 (paper)** — 8 weeks; highest leverage
2. **Phase 42 (open-science release)** — 2 weeks; parallel with paper
3. **Phase 43 (talks)** — 2 weeks; after paper drafted
4. **Phase 44 (educational)** — 4 weeks; lower urgency
5. **Phase 45 (CI)** — 1 week; useful infrastructure
6. **Phase 46 (engagement)** — ongoing

## Success criteria

v11 ships when:
1. The paper is submitted to J. Math. Cryptol. (Phase 41)
2. The GitHub release has working reproducible builds (Phase 42)
3. At least one workshop talk is delivered (Phase 43)

This is *publication impact*, not new science. The science is in v1-v10.

## Cost-benefit

| Investment | Output |
|------------|--------|
| 2 months of writing | Peer-reviewed paper |
| 2 weeks of release engineering | Citable open-science artifact |
| 1 week of talk prep | Visibility in crypto community |
| **Total: ~3 months** | **High durability and dissemination** |

## What v11 will NOT pursue

- Patent applications (the science is mathematics; non-patentable
  and we're for open science)
- Commercialization (no products here; pure research)
- Press / media engagement (only after peer-reviewed publication)

## Companion programs

- v1-v6 provide the *content* of the paper
- v7-v10 are *future work* mentioned in §6 Discussion
- v11's paper § structure mirrors PROGRAMS.md's program structure

## What comes after v11

Once published, the v1-v6 results are durably attributed. Subsequent
programs (v7+) build on this established baseline. Each future paper
can cite the v11 paper as the empirical foundation.

If the v11 paper is rejected: revise per reviewer feedback, resubmit
to ASIACRYPT or J. Cryptology. Even rejection feedback is valuable.

If accepted: the cryptanalysis community has a definitive map of what
classical attacks have been tried on LMFDB prime-field ECDLP. This
saves dozens of future researcher-years.
