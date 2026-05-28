# Research Program v9: Long-horizon roadmap

**Date:** 2026-05-28
**Status:** speculative; multi-year vision
**Motivation:** v7 attacks the constant-factor gap; v8 attacks
algorithmic frontiers. v9 asks: **what would the cryptanalytic
landscape look like in 5-10 years**, and what foundational
investigations should we start now to be ready?

This program is intentionally over-the-horizon. Each phase is a
multi-month or multi-year investigation that doesn't fit in any
single AutoLab round. It's a research agenda, not a sprint plan.

## Phase 31: Provable separation between EC and DH-style DLP

**Open question (classical):** Is ECDLP provably harder than the
discrete-log problem in `F_p^*`? IFP (factoring) and DLP in `F_p^*`
both have L(1/3) sub-exponential algorithms. ECDLP (for "good"
curves) has only O(√n) attacks. **Why?**

The standard answer is "subexponential index calculus doesn't apply
to elliptic curves" — but this is descriptive, not proof. A proof
would show: any classical algorithm requires Ω(n^{1/2}) operations.

**Concrete first step:** Formalize the gap in the *generic group
model*. Shoup 1997 gives Ω(√n) in the GGM, but the GGM is too
restrictive (treats group elements as opaque). We need a
"semi-algebraic" model that allows the curve's coordinate structure
but not other side info.

**Risk:** This is a major open problem in cryptography. Several
candidates exist (algebraic complexity, polynomial hierarchy
relativizations) but no progress in the last 20 years.

## Phase 32: Quantum-classical hybrid attacks

**Hypothesis:** Shor's algorithm requires a fault-tolerant quantum
computer with thousands of logical qubits. But a NISQ device with
50-100 qubits + classical post-processing might give *partial*
speedups (e.g., via quantum walks, Grover-amplified rho, or
quantum-classical Pollard variants).

**Concrete first step:** Survey the literature on:
- Grover-amplified rho (Brassard-Hoyer-Tapp 1998 — collision-finding
  speedup `n^{1/3}` instead of `n^{1/2}`)
- Quantum walk-based discrete log algorithms
- Adiabatic optimization for ECDLP

**Risk:** Most quantum speedups assume fault-tolerant QC. NISQ-era
results are mostly null. Worth surveying though.

## Phase 33: Probabilistic ECDLP at low-confidence

**Hypothesis:** Suppose we relax "find the discrete log" to "find
the discrete log with probability ε for some ε > 0". With ε = 0.01,
many "almost works" algorithms might apply. The information-theoretic
lower bound is still Ω(√n) but the practical constants might be lower.

**Concrete first step:** Define a precise probabilistic-ECDLP problem
and study which classical algorithms (lattice, sieve, MLP-classifier)
achieve non-trivial ε at sub-`O(√n)` cost.

**Risk:** "almost solve" might not be useful for cryptanalysis
in practice (one wants the right key, not 1% of the right key).

## Phase 34: Foundational complexity-theoretic separations

**Open question:** Is ECDLP NP-intermediate? Pseudorandom-permutation-
hard? `NP ∩ co-NP ∩ BPP`-hard? Different answers have different
implications.

**Concrete first step:** Survey the state of complexity-theoretic
ECDLP placement. Compare to FACTORING (similarly placed) and
GRAPH-ISOMORPHISM (NP-intermediate, polynomial-time on quasi-polynomial
sizes).

**Risk:** Pure theory; no immediate algorithmic payoff.

## Phase 35: Cryptanalysis-aware curve generation

**Hypothesis:** If we generate curves with specific algebraic
invariants known to resist attacks, we get provably hard ECDLP
instances. The current "use NIST curves or random Brainpool"
approach is folklore.

**Concrete first step:** Design a curve-generation pipeline that:
- Avoids all attack families in our resistance map
- Has provably hard ECDLP (assuming the GGM)
- Generates explicit certificates of attack resistance

**Risk:** No new mathematics; mostly engineering. But useful
contribution to cryptographic practice.

## Phase 36: Open conjectures with ECDLP corollaries

**Examples:**
- BSD conjecture: would give explicit rank determination → some MW
  attacks become applicable to specific curves
- Lang-Trotter conjecture: gives statistics on trace, but no
  algorithmic exploit known
- Sato-Tate conjecture: similar; statistical, not algorithmic

**Concrete first step:** For each major conjecture, identify the
*algorithmic* implication for ECDLP if proven (or violated for
specific curves).

**Risk:** Conjectures may take decades to resolve.

## Phase 37: Side-channel-resistant ECDLP analysis

**Hypothesis:** Some implementations of ECDLP-using protocols
(ECDSA, ECDH) leak partial information. A research program
focused on *what the leakage means* for the underlying ECDLP
hardness is distinct from the pure-cryptanalytic work above.

**Concrete first step:** Catalog the known side-channel attacks
(Brumley-Tuveri 2012, Aldaya et al. 2019, ROBOT 2017) and
characterize the leakage model. Then study whether any leakage
gives sub-`O(√n)` ECDLP for the abstract problem.

**Risk:** Most side-channel attacks are implementation bugs,
not ECDLP attacks per se. But the boundary is interesting.

## Phase 38: Cross-disciplinary methods (Wave 1)

**Examples:**
- **Random matrix theory** for L-function distribution and trace
  statistics → predicting weak points in curve families
- **Spectral analysis** of isogeny graphs → eigenvalue distribution,
  Cheeger constant, mixing time
- **Algebraic dynamics** (post-quantum framework) → could ECDLP
  be reformulated as a dynamical system fixed-point problem?
- **Information geometry** → curve manifold, geodesic between curves,
  attack-vulnerability metric

**Concrete first step:** Read seminal papers in each area; identify
points of contact with ECDLP. Write a position paper.

**Risk:** Highly speculative; most cross-disciplinary attempts fail.

## Phase 39: Distributed cryptanalysis at scale

**Concrete deliverable:** Build the infrastructure for a 100-machine
distributed Pollard rho computation. With 100 machines × 8 CPUs × 4 days
= 3200 CPU-days, we could solve much larger ECDLP instances than the
AutoLab budget allows.

**Concrete first step:** Take the `phase21_rho_v3.c` C-extension,
add a coordinator service (gRPC or HTTP) that distributes trails to
a central database, implements distinguished-point matching across
machines, and handles fault tolerance.

**Risk:** Engineering, not research. But unblocks larger-scale
experiments.

## Phase 40: Publication of resistance map

**Deliverable:** Convert the v1-v8 results into a formal research
paper for J. Math. Cryptol., CRYPTO, EUROCRYPT, or ASIACRYPT.

**Title (draft):** "An Empirical and Theoretical Map of Cryptanalytic
Resistance for Prime-Field ECDLP on LMFDB Curves"

**Sections:**
1. Background: ECDLP, Pollard rho, Semaev/Diem/Yokoyama
2. Methodology: 50 trials, full cost accounting, multi-curve
3. Results: every attack family closed (Tables 1-5)
4. Engineering: 86× C speedup, 5.5× gap analysis
5. Theoretical confirmation: Yokoyama lower bound matches empirics
6. Discussion: what remains genuinely open

**Risk:** Significant writing effort (3-6 months).

## Priority for the v9 horizon

Within 1 year:
- Phase 39 (distributed infrastructure) — engineering, high-leverage
- Phase 40 (publication) — establishes priority on resistance map
- Phase 35 (curve-generation pipeline) — useful contribution

Within 3 years:
- Phase 38 (cross-disciplinary) — broad survey, possible breakthroughs
- Phase 32 (quantum-classical hybrid) — important for post-quantum era
- Phase 31 (provable separation) — major open problem

Within 5+ years:
- Phase 34 (complexity-theoretic placement) — generational problem
- Phase 36 (conjecture-conditional results) — wait for math progress

## What v9 deliberately omits

- Specific implementation work — those are in v7 (engineering) and v8
  (algorithmic)
- Concrete experiments at the AutoLab scale — v6 is the empirical
  ground truth
- Speculation about specific yet-undiscovered attacks — we can't plan
  for what we don't know

## Companion structure

| Program | Focus | Time horizon |
|---------|-------|--------------|
| v1-v3 | Initial exploration | weeks |
| v4 | Groebner / Semaev deep dive | days |
| v5 | Non-naive attempt #1 | days |
| v6 | Non-naive attempt #2 + C engineering | days |
| v7 | Close the C-extension gap | weeks |
| v8 | Algorithmic frontiers | months |
| **v9** | **Long-horizon roadmap** | **years** |

## Connection to the rest of the crypto repo

This v9 plan complements existing `RESEARCH_*.md` documents at the
repo root by providing a research agenda rather than completed work.
Each phase here links forward to where actual investigations would
live (e.g., a Phase 39 implementation might live in
`distributed_rho/` directory).
