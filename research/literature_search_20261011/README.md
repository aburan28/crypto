# Literature search: 2025–26 papers relevant to ECDLP, index calculus and Gröbner solving

**Date.** 2026-10-11. **Type.** Literature search round (no experiment was run).
**Status of every statement.** *Proposed* unless marked otherwise. Summaries come from
abstract pages that were opened and checked; full papers were not read. No claim here is
derived, independently verified or measured, and no claim about a deployed curve is
made or withdrawn.

Figure: [`landscape.svg`](landscape.svg) (source [`plot_landscape.py`](plot_landscape.py),
data [`papers.csv`](papers.csv)). PDF: [`report.pdf`](report.pdf) (source [`report.typ`](report.typ)).

## Question and scope

Has anything posted on IACR ePrint or arXiv between 2025-01 and 2026-10 (plus important
late-2024 items) given an algorithmic improvement this repo could use for elliptic-curve
DLP? That covers index calculus and summation polynomials, Gröbner/XL/Crossbred
solving, sparse linear algebra, rho/BSGS/kangaroo, and indirect routes (isogeny, cover,
pairing, small-root oracles).

## Method

Three parallel searches, about 115 distinct queries in total, against:

- ePrint keyword search and recent listings;
- the arXiv API (cs.CR, cs.SC, math.AC, math.NT, math.AG);
- general web search.

Each listed paper's abstract page was opened. Companion versions were merged, for
example arXiv 2410.23211 with ePrint 2025/822. Papers already cited in the repo were
excluded: ePrint 2024/1985, 2024/924, 2022/975, 2020/1315 and 2019/313; arXiv
2610.06783, 2304.07820, 2208.00844, 2205.09888 and 2112.10506.

Each paper is classified by the repo component it touches and by its expected impact
in [`papers.csv`](papers.csv), N = 61. The classification is a judgement made from the
abstracts.

## Findings

**Headline (negative).** No paper in the window claims or plausibly delivers an
exponent-level improvement for ECDLP, on prime-field, binary, or Koblitz curves.

- **Searches that came up empty.** Searches for summation polynomials, point
  decomposition, Weil descent/restriction combined with ECDLP, the GHS attack,
  first/last fall degree, Koblitz attacks, and kangaroo found nothing new in
  2025–26. The newest hits are ePrint 2019/313 (summation polynomials) and ePrint
  2021/676 (Weil descent / GHS).
- **No trapdoor or classical break.** Nothing new was found for NIST P-curves,
  secp256k1, Brainpool or the sect curves.

### 1. Lower bounds that formalise the repo's "generic" verdict

- **Corrigan-Gibbs, Henzinger, Wu, *The Structured Generic-Group Model*, ePrint
  [2026/384](https://eprint.iacr.org/2026/384), EUROCRYPT 2026.**
  - Any DL algorithm that exploits structure on only a δ fraction of group elements
    needs time Ω(min(√q, 1/δ)).
  - The abstract states analogous bounds for algorithms that exploit elliptic-curve
    point structure.
  - This is the theorem-level version of `RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md`
    §0.1. EC index calculus can beat rho only if the decomposable set is dense,
    δ > q^(-1/2), and still admits an algebraic test that is cheaper than a group
    operation.
- **Dinur, Keller, Marmor**, arXiv [2505.00894](https://arxiv.org/abs/2505.00894) /
  ePrint 2025/783, STOC 2026.
  - Non-adaptive advice cannot beat online time √N unless S = Ω(√N).
  - This caps every precompute arm whose online queries do not depend on the
    challenge. It matches the earlier finding that multiplier closure is a
    precompute-only route.
- **Akshima et al.**, ePrint [2026/2253](https://eprint.iacr.org/2026/2253): the
  multi-instance generic baseline is √(k·q).
- **Hu**, ePrint [2024/1923](https://eprint.iacr.org/2024/1923), itself concludes that
  no subexponential prime-field ECDLP method is known.

### 2. Tools aimed at the repo's open quantity: the solving degree of Weil-descended S_{m+1}

None of these papers bounds the solving degree for summation systems. They supply
computable diagnostics:

- **Kudo–Yokoyama**, ePrint [2024/528](https://eprint.iacr.org/2024/528) (revised
  2026-09).
  - Bounds the maximal Gröbner-basis degree through the saturation exponent of the
    homogenising variable.
  - Shows that affine and homogenised computations stay in step up to the first
    degree fall, which needs no semi-regularity assumption.
  - Companion: arXiv [2410.23211](https://arxiv.org/abs/2410.23211) /
    ePrint 2025/822, on generalised cryptographic semi-regularity.
- **Ning–Ran–Samardjiska**, *Multi-homogeneous XL*, ePrint
  [2025/2060](https://eprint.iacr.org/2025/2060).
  - Gives regularity theory for multi-graded systems.
  - The descended S_{m+1} has one variable block of size ℓ per point. A multi-graded
    regularity could fall below the total-degree d_reg the repo currently measures.
- **Gorla–Trebiani**, arXiv [2609.14052](https://arxiv.org/abs/2609.14052).
  - Proves solving-degree bounds for a model in elementary symmetric functions over
    GF(2).
  - This is the closest template for proving bounds on symmetrised summation systems.
- **Diagnostics for where the system stops looking generic:**
  - Sakata–Takagi, arXiv [2505.10246](https://arxiv.org/abs/2505.10246): generic
    leading monomials computed without any reductions, so the first divergence marks
    where the structure takes over.
  - Gaggero–Gorla, ePrint [2024/1924](https://eprint.iacr.org/2024/1924):
    "algebraic randomness" test.
  - Bardet–Gilard, arXiv [2502.12721](https://arxiv.org/abs/2502.12721): an example
    of proving a Hilbert series for a structured, non-semi-regular ideal.

### 3. Constant-factor engineering worth a falsifier

| Paper | Mechanism | Repo stage |
|---|---|---|
| Demin, arXiv [2607.06372](https://arxiv.org/abs/2607.06372) (Groebner.jl) | SIMD-packed coefficients: one trace replays k instances at once | inherited-F4 replay |
| Horak–Semaev, arXiv [2603.24890](https://arxiv.org/abs/2603.24890) | inconsistency probability of sparse F2 systems; most decompositions have no solution | early abort when 1 ∈ I |
| Vidal–Delaplace–Ionica, ePrint [2026/360](https://eprint.iacr.org/2026/360); Baena et al. [2024/758](https://eprint.iacr.org/2024/758); Jaques et al. [2026/408](https://eprint.iacr.org/2026/408) | Crossbred preprocessing trimming; admissibility condition; XL "between" integer degrees | crossbred |
| Bouillaguet–Sauvage, ePrint [2026/1704](https://eprint.iacr.org/2026/1704); Liu et al. [2026/040](https://eprint.iacr.org/2026/040) | Boolean MQ hybrid; evaluation over Hamming balls | GF(2) exhaustive phase (after quadratisation) |
| Chou–Niederhagen, ePrint [2026/1787](https://eprint.iacr.org/2026/1787) | Wiedemann-XL: compute only the needed kernel coordinates (≈1.5×) | sparse linear algebra |
| Al Aswad–Pierrot–Thomé, ePrint [2026/560](https://eprint.iacr.org/2026/560) | order-k Galois automorphisms give about k² in TNFS linear algebra | IC linear algebra; probably already captured by the GGMP n² Frobenius saving: check |
| Dryło, ePrint [2026/1299](https://eprint.iacr.org/2026/1299) | a G-invariant compression equals f∘Φ with Φ: E→E/G | 2/4-torsion symmetrisation equals x-only S_m on E/G: structural, not new speed |
| Banegas–Sarkis–Smith, ePrint [2026/704](https://eprint.iacr.org/2026/704) | x-only binary isogenies at 5s M (about 4.2× over Vélu) | GHS weak-curve isogeny walks |
| Cenk et al., ePrint [2026/2102](https://eprint.iacr.org/2026/2102) | table-free BSGS with co-Z paired giant steps and windowed batch inversion | BSGS/rho arithmetic |

### 4. Prime-field small-root oracle (S5 at δ ≈ 1/5)

Two new tools give the S5 Coppersmith route a cheap feasibility check on paper:

- Ryan's **cuso**, ePrint [2024/1577](https://eprint.iacr.org/2024/1577) (EUROCRYPT 2025),
  which computes optimal shift sets.
- The LP formulation of automated-Coppersmith asymptotics, ePrint
  [2026/1027](https://eprint.iacr.org/2026/1027) (ASIACRYPT 2026).

Neither targets S5, and the repo's measured reach is 1/80.

### 5. Indirect routes: none open a new ECDLP path

The volcano, isogeny and module-action papers improve isogeny *recovery*:

- ePrint [2025/1243](https://eprint.iacr.org/2025/1243): ascending volcanoes.
- ePrint [2024/1556](https://eprint.iacr.org/2024/1556): module action.
- ePrint [2025/1706](https://eprint.iacr.org/2025/1706): Kani's lemma via Clifford
  algebra.
- ePrint [2026/1823](https://eprint.iacr.org/2026/1823): HyperSolver, supersingular
  only.

Trace and embedding degree are invariants of the isogeny class, so these walks cannot
reach a MOV-weak or anomalous curve. GHS cover weakness on binary curves is not such
an invariant, which is why ePrint 2026/704 matters only for binary-curve walks.

### 6. Unsupported claims (flagged)

- **Hu, ePrint [2025/015](https://eprint.iacr.org/2025/015).** The method ("IICA") is
  undefined. The only claim is about a 2× success probability, with no parameters.
- **Mahalanobis, arXiv [2607.09814](https://arxiv.org/abs/2607.09814).** Zero-minor
  search; the abstract states no complexity.
- **Lin et al., D-Wave ECDLP (ICN 2026, DOI 10.23919/ICN.2026.0017).** A 16-bit
  instance.
- **Ma–Liu–Jiao "Proper-Cover", arXiv [2607.09163](https://arxiv.org/abs/2607.09163).**
  Claims to beat F5, but not on cryptographic systems. Not reproduced.

## Proposed falsifiers, ranked by expected information per hour

1. **Saturation exponent and first degree fall vs measured solving degree** (2024/528)
   on descended S3/S4 at n = 11–23, using systems the repo already generates.
   - *Falsifier:* the bound is loose, or the systems fail the top-degree
     semi-regularity hypothesis. That would itself characterise the structure.
2. **Multi-graded Hilbert function** of descended S3 at n ≈ 15 with point blocks
   (2025/2060).
   - *Falsifier:* no multidegree (d₁, d₂) with d₁ + d₂ below the measured total
     solving degree.
3. **δ table for the structured-GGM bound** (2026/384). For each Koblitz and
   prime-field IC configuration, compute δ = |decomposable set| / r from existing
   records and place the measured cost against min(√r, 1/δ).
   - *Falsifier:* measured cost below the bound would mean an accounting bug.
4. **Early abort by refutation degree** (2603.24890). Log the degree at which F4
   proves 1 ∈ I, compared with the degree at which it finds solutions.
5. **Lane-batched inherited-F4 replay** (2607.06372). Run k = 4–8 target points
   through one step plan.
6. **E/G check for torsion symmetrisation** (2026/1299). Compare the degree and term
   count of the FGHR-symmetrised S3/S4 on E with plain x-only S3/S4 on E/⟨T⟩.
7. **S5 feasibility on paper** via the LP of 2026/1027, or cuso at p ≈ 2^32.
   - *Falsifier:* an asymptotic bound below p^(1/5) closes the Coppersmith route.

## Canonical graphs checked

These graphs were checked and left unchanged, because the search produced no new
measured value or verified edge:

- `docs/index-calculus-scoreboard.html`
- `docs/ic/progress-timeline.json`
- `docs/curves/`
- `docs/performance-gains/`
- `figures/`

## Unresolved checks and coverage caveats

- ePrint search matches titles and keywords, and some requests were rate-limited
  (HTTP 429). A relevant paper whose title avoids these terms may have been missed.
- Two papers were not opened: ePrint 2026/1198 and 2025/271.
- The claim attributed to arXiv 2508.14011 that the secp256k1 interval record is
  129 bits came from a search snippet and was not verified.
- Finite-field DL record lists were not checked.
